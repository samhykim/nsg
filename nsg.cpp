#include "gen_hashes.h"
#include "find_xz.h"
#include "nsg_util.h"  
#include "threadsafe_multimap.h"

#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <thread>
#include <algorithm>
#include <atomic>
#include <iterator> 
#include <cstring>   
#include <map>
#include <sstream>
#include <string>
#include <set>

using namespace std::chrono;

#define MIN_L 10
#define NUM_LOCKS 419
// Inverse load factor of hash table.  nominally = 1.4.
#define INV_LOAD_FACTOR 1.4

typedef ThreadSafeMultiMap<hash_t, read_t> s_map;

static vector<SpinLock> spinLocks(NUM_LOCKS);

static vector<hash_s> hashes;
static atomic_llong NSGcount;

static vector<Match> graph;
static vector<Match> reverse_graph;

static read_t indexReads;
static int readLength;
static read_t numReads;

// Memory map pointer for input fasta file
static char* mmapptr;

// Construct string graph by finding all greedy and nsg reads
// We process the prefixes of given length l for all the reads
// and find the reads whose suffix with length l matches each prefix
void ProcessPrefix(const s_map* sMap, const read_t rStart, const read_t rEnd,
                   const int l) {
  read_t suffixMatch;
  int overlap;
  for (int read_id = rStart + 1; read_id <= rEnd; read_id++) {
    const auto its = sMap->equal_range(hashes[read_id + indexReads * l].prefix);
    int k = 1;
    // loop through reads with matching suffixes
    for (auto it = its.first; it != its.second; ++it) {
      suffixMatch = (it->second);

      Match& left_match = reverse_graph[read_id];
      Match& right_match = graph[suffixMatch];

      spinLocks[suffixMatch % NUM_LOCKS].lock();

      if (!left_match.nsg_found) {
        if (!left_match.greedy_found) {
          if (!right_match.greedy_found) {
            right_match.SetGreedy(read_id, l);
            left_match.SetGreedy(suffixMatch, l);
          } else {
            if (!right_match.nsg_found) {
              right_match.SetNSG(read_id, l);
              left_match.SetGreedy(suffixMatch, l);
              NSGcount++;
            }
          }
        } else if (!right_match.greedy_found) {
          right_match.SetGreedy(read_id, l);
          left_match.SetNSG(suffixMatch, l);
        } else if (left_match.greedy_found) {
          overlap = readLength + l - left_match.greedy_weight;
          if (hashes[left_match.greedy_id + indexReads * overlap].prefix !=
              hashes[suffixMatch + indexReads * overlap].suffix) {
            if (right_match.greedy_id != read_id) {
              left_match.SetNSG(suffixMatch, l);
              right_match.SetNSG(read_id, l);
              NSGcount++;
            }
          }
        } else if (k == 2) {
          left_match.SetNSG(0, 0);
        }
      }
      k++;
      spinLocks[suffixMatch % NUM_LOCKS].unlock();
      if (k > 2) {
        break;
      }
    }
  }
}

// Add suffixes of length l to sMap
void AddSuffix(s_map* sMap, const read_t rStart, const read_t rEnd,
               const int l) {
  hash_t hashToInsert;
  for (read_t r = rStart; r < rEnd; r++) {
    // only need to add to the table if its nsg hasn't been found
    Match& right_match = graph[r];
    if (!right_match.nsg_found) {
      hashToInsert = hashes[r + indexReads * l].suffix;
      sMap->insert(hashToInsert, r);
    }
  }
}

// Clear and remove all elements in the suffix map
void ClearSuffixMap(s_map* sMap, const read_t bStart, const read_t bEnd) {
  sMap->clear_buckets(bStart, bEnd);
}

int main(int argc, const char* argv[]) {
  if (argc == 1) {
    cout << endl
         << "usage: NSG [-flags] [inputFilename] [outputFilename]" << endl
         << endl;
    cout << "valid flags are:" << endl;
    cout << "   t : run with t threads (default is max number of cores on your "
            "machine )" << endl;
    cout << "       must be followed by a number " << endl;
    cout << "   x : don't create output file" << endl
         << endl;
    cout << "   e : do read extension step" << endl
         << endl;

    cout << "example: './NSG -xs bact1.reads' will run in serial, disable "
            "output, and use reads from bact1.reads." << endl;

    return 0;
  }

  int num_threads = std::thread::hardware_concurrency();
  string DEFAULT_OUTPUT_FILE = "output.txt";
  const char* input_file;
  string output_file;
  bool output_flag = true;
  bool read_extend = false;
  int ip_index = 1, op_index;

  if (strncmp(argv[1], "-", 1) == 0) {
    // flags set
    string flags = argv[1];

    if (flags.find("x") != std::string::npos) {
      output_flag = false;
      cout << "output file disabled" << endl;
      ip_index = 2;
      op_index = 3;
    }
    if (flags.find("t") != std::string::npos) {
      num_threads = std::stoi(argv[2]);
      if (num_threads > std::thread::hardware_concurrency()) {
        std::cerr << num_threads << " threads not supported on this machine."
                  << endl;
        return -1;
      }
      ip_index = 3;
      op_index = 4;
    }
    if (flags.find("e") != std::string::npos) {
      read_extend = true;
      cout << "read extension step enabled" << endl;
    }
  } else {
    ip_index = 1;
    op_index = 2;
  }

  if (argc > op_index) {
    input_file = argv[ip_index];
    output_file = argv[op_index];
    cout << "Input file: " << input_file << endl;
    cout << "Output file: " << output_file << endl;
  } else {
    input_file = argv[ip_index];
    output_file = DEFAULT_OUTPUT_FILE;
    cout << "No output file specified. Writing to default file: "
         << DEFAULT_OUTPUT_FILE << endl;
  }

  std::ifstream readFile(input_file);
  if (!readFile.good()) {
    std::cerr << "Error opening read file '" << input_file << "'. Bailing out."
              << std::endl;
    return -1;
  }

  string readSeq;
  getline(readFile, readSeq);
  if (readSeq.find(">") != std::string::npos) {
    getline(readFile, readSeq);
  }
  readLength = (int)readSeq.length();

  // count number of reads in the file.
  // new lines will be skipped unless we stop it from happening:
  readFile.unsetf(std::ios_base::skipws);

  // count the newlines with an algorithm specialized for counting:
  indexReads = (read_t)std::count(std::istream_iterator<char>(readFile),
                                  std::istream_iterator<char>(),
                                  '\n') +
               1;  // need to add one for the line that we just read.
  readFile.close();

  cout << "Length of reads: " << readLength << endl;

  read_t num_buckets = indexReads * INV_LOAD_FACTOR;
  static s_map* p_suffix_map = new s_map(num_buckets);

  hashes.resize(indexReads * (readLength + 1));
  graph.resize(indexReads + 1);
  reverse_graph.resize(indexReads + 1);

  cout << "Number of threads = " << num_threads << std::endl;

  thread* t = new thread[num_threads - 1];
  read_t NSGcount = 0;

  clock_t e_suffix, e_prefix;
  clock_t s_suffix, s_prefix;
  double suffix_ms = 0, prefix_ms = 0, prefix_msClock = 0, suffix_msClock = 0;
  // Begin Recording Time for GenHashes
  clock_t start_cpu = clock();
  auto start0 = high_resolution_clock::now();
  auto startHash = high_resolution_clock::now();

  numReads = GenHashes(input_file, &hashes, indexReads, readLength);
  if (numReads < 0) {
    cout << "Problem reading input file.  No reads found." << endl;
    return -1;
  } else if (indexReads < numReads) {
    cout << "strange stuff: file line count was less than number of reads."
         << endl;
    return -1;
  }
  cout << "Number of reads: " << numReads << endl;

  auto endHash = high_resolution_clock::now();
  cout << "hash real time: "
       << duration_cast<duration<float>>(endHash - startHash).count() * 1000
       << "ms." << endl
       << endl;

  auto s_prefixClock = high_resolution_clock::now();
  auto e_prefixClock = high_resolution_clock::now();
  auto s_suffixClock = high_resolution_clock::now();
  auto e_suffixClock = high_resolution_clock::now();

  read_t readsPerThread = numReads / num_threads;
  read_t suffixPerThread = numReads / num_threads;
  read_t clearPerThread = numReads / num_threads;

  read_t suffixBound[num_threads + 1];
  read_t clearBound[num_threads + 1];
  read_t readBound[num_threads + 1];

  for (int i = 0; i < num_threads; i++) {
    suffixBound[i] = i * suffixPerThread;
    clearBound[i] = i * clearPerThread;
    readBound[i] = i * readsPerThread;
  }

  suffixBound[num_threads] = numReads;
  clearBound[num_threads] = num_buckets;
  readBound[num_threads] = numReads;

  NSGcount = 0;
  // Separate code block for ProcessPrefix
  // Ensures all threads finish before moving on the FindXZStructures
  //==========================================================================================
  // Begin iteration on overlaps l
  for (int l = readLength - 1; l >= MIN_L; l--) {
    //==========================================================================================
    // start clocks for tracking time spent recording/clearing suffixes in hash
    // table

    s_suffix = clock();
    s_suffixClock = high_resolution_clock::now();

    //==========================================================================================

    //==========================================================================================
    // Begin processing suffixes of length l

    if (l < readLength - 1) {
      for (int i = 0; i < num_threads - 1; i++) {
        t[i] = thread(ClearSuffixMap, p_suffix_map, clearBound[i],
                      clearBound[i + 1]);
      }
      ClearSuffixMap(p_suffix_map, clearBound[num_threads - 1],
                     clearBound[num_threads]);
      for (int i = 0; i < num_threads - 1; i++) {
        t[i].join();
      }
    }
    for (int i = 0; i < num_threads - 1; i++) {
      t[i] = thread(AddSuffix, p_suffix_map, suffixBound[i], suffixBound[i + 1],
                    l);
    }
    AddSuffix(p_suffix_map, suffixBound[num_threads - 1],
              suffixBound[num_threads], l);
    for (int i = 0; i < num_threads - 1; i++) {
      t[i].join();
    }

    // Done processing suffixes
    //==========================================================================================

    //==========================================================================================
    // stop clocks for tracking time spent recording/clearing suffixes in hash
    // table
    e_suffix = clock();
    e_suffixClock = high_resolution_clock::now();

    suffix_msClock +=
        duration_cast<duration<float>>(e_suffixClock - s_suffixClock).count();
    suffix_ms += (e_suffix - s_suffix) / double(CLOCKS_PER_SEC) * 1000;

    // start clocks for tracking time spent processing prefixes
    s_prefix = clock();
    s_prefixClock = high_resolution_clock::now();

    //==========================================================================================

    //==========================================================================================
    // Begin processing prefixes of length l

    // process the prefixes in parallel.  (Note: this code also works for
    // serial)
    for (int i = 0; i < num_threads - 1; i++) {
      t[i] = thread(ProcessPrefix, p_suffix_map, readBound[i], readBound[i + 1],
                    l);
    }
    ProcessPrefix(p_suffix_map, readBound[num_threads - 1],
                  readBound[num_threads], l);

    for (int i = 0; i < num_threads - 1; i++) {
      t[i].join();
    }

    // Done processing suffixes
    //==========================================================================================

    //==========================================================================================
    // stop clocks for tracking time spent processing prefixes
    e_prefix = clock();
    e_prefixClock = high_resolution_clock::now();

    prefix_msClock +=
        duration_cast<duration<float>>(e_prefixClock - s_prefixClock).count();
    prefix_ms += (e_prefix - s_prefix) / double(CLOCKS_PER_SEC) * 1000;
    //==========================================================================================

    if (NSGcount >= numReads) {  // all NSG's found.  Can exit loop.
      break;
    }

    //==========================================================================================
    // End of iteration on overlaps l
  }
  //==========================================================================================

  //==========================================================================================
  // On to Finding X and Z Structures Part of the Algorithm!
  //==========================================================================================
  NSGcount = 0;
  int minL = readLength;
  for (auto r = 0; r < numReads; ++r) {
    Match& right_match = graph[r];
    if (right_match.nsg_id != 0) {
      NSGcount++;
      minL = min(minL, right_match.nsg_id);
    }
  }
  cout << NSGcount << " of " << numReads << " NSG matches found." << endl;
  cout << "Min overlap: " << minL << endl
       << endl;
  cout << "Identifying X and Z Structures" << endl
       << endl;

  cout << "input file: " << input_file << endl;


  
  size_t filesize = GetFilesize(input_file);
  int fd = open(input_file, O_RDONLY);
  assert(fd != -1);
  // MEMORY MAPPED INPUT FILE. Create a memory map of the reads in order to quickly access 
  // individual nucleotides. Required for read extensions.
  mmapptr = (char*)mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, fd, 0);

  FindXZStructures(mmapptr, readLength, numReads, read_extend, &graph, &reverse_graph);

  // Cleanup memory map
  if (mmapptr != NULL) {
    munmap(mmapptr, filesize);
  }
  close(fd);

  //==========================================================================================

  auto end0 = high_resolution_clock::now();
  clock_t end_cpu = clock();

  cout << "suffix:   cpu  time: " << suffix_ms << "ms." << endl;
  cout << "          real time: " << suffix_msClock * 1000 << "ms." << endl
       << endl;
  cout << "prefix:   cpu  time: " << prefix_ms << "ms." << endl;
  cout << "          real time: " << prefix_msClock * 1000 << "ms." << endl
       << endl;

  cout << "combined  cpu  time: "
       << (end_cpu - start_cpu) / double(CLOCKS_PER_SEC) << "s." << endl;
  cout << "combined  real time: "
       << duration_cast<duration<float>>(end0 - start0).count() << "s." << endl
       << endl;

  cout << "Processing speed: "
       << (duration_cast<duration<float>>(end0 - start0).count()) /
              (numReads * readLength) * 1000000000 << "ns/char." << endl
       << endl;

  NSGcount = 0;
  for (auto r = 0; r < numReads; ++r) {
    Match& right_match = graph[r];
    if (right_match.nsg_id != 0) {
      NSGcount++;
      cout << "NSG remaining for: " << r << endl;
    }
  }
  cout << NSGcount << " of " << numReads << " NSG matches remaining." << endl;

  // Write to output file
  if (output_flag) {
    cout << "Writing output ...";
    ofstream myfile;
    myfile.open(output_file, ios::trunc);

    for (read_t read_id = 1; read_id <= numReads; read_id++) {
      const Match& match = graph[read_id];
      myfile << read_id << " " << match.greedy_id << " " << match.greedy_weight
             << " " << match.nsg_id << " " << match.nsg_weight << endl;
    }
    myfile.close();
    cout << " done!" << endl;
  }
  return 1;  // should have found all NSG-s
}
