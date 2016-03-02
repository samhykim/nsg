#ifndef __NSG__gen_hashes__
#define __NSG__gen_hashes__

#include "nsg_util.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <ctime>
#include <vector>
#include <thread>
#include <algorithm>
#include <stack>


#define A 13
#define B 2523345
#define P 1000000007

read_t GenHashes(const string fname, vector<hash_s> * hashes, read_t num_reads, int read_length) {
  hash_t hashA,hashB, fullHashA, fullHashB;
  hash_t aPow [read_length];
  hash_t bPow [read_length];
  
  aPow[0] = 1;
  bPow[0] = 1;
  for(int k = 1; k < read_length; k++) {
    aPow[k] = (aPow[k-1] * A) % P;
    bPow[k] = (bPow[k - 1] * B) % P;
  }
  
  std::ifstream readFile(fname);
  if(!readFile.good()){
    std::cerr << "Error opening read file '" << fname << "'. Bailing out." << std::endl;
    return -1;
  }
  
  string readSeq;
  // Start read index at 1
  read_t readID = 1;
  int k;
  while(getline( readFile, readSeq ).good()) {
    if (readSeq.find(">") != std::string::npos) {
      continue;
    }
      
    hashA=0;
    hashB=0;

    //roll up to get all suffix hashes
    for(k = read_length; k > 1; k--) { 
      hashA=((hashA+aPow[k-1]*readSeq[k-1])%P+P)%P; // so no hashes are negative
      hashB=((hashB+bPow[k-1]*readSeq[k-1])%P+P)%P; // so no hashes are negative
      (*hashes)[readID+num_reads*(read_length-k+1)].suffix=(hashA<<32) | hashB;
    }
    
    //save full hash of length L
    fullHashA = ((readSeq[0]+hashA)%P+P)%P; // sum_{i=0}^{L-1} a^i s_i mod P
    fullHashB = ((readSeq[0]+hashB)%P+P)%P; // sum_{i=0}^{L-1} b^i s_i mod P
    (*hashes)[readID+num_reads*(read_length)].suffix=(fullHashA<<32) | fullHashB;
    (*hashes)[readID+num_reads*(read_length)].prefix=(fullHashA<<32) | fullHashB;
    
    //roll down to get all prefix hashes
    hashA = ((A*(readSeq[0]+hashA-aPow[read_length-1]*readSeq[read_length-1]))%P+P)%P;
    hashB = ((B*(readSeq[0]+hashB-bPow[read_length-1]*readSeq[read_length-1]))%P+P)%P;
    (*hashes)[readID+num_reads*(read_length-1)].prefix=(hashA<<32) | hashB;;
    
    for(k = read_length - 2; k > 0; k--){
        hashA=((A*(hashA-aPow[read_length-1]*readSeq[k]))%P+P)%P;
        hashB=((B*(hashB-bPow[read_length-1]*readSeq[k]))%P+P)%P;
        (*hashes)[readID+num_reads*k].prefix=(hashA<<32) | hashB;
    }
    
    readID++;
    if (readID >= MAX_READS) {
      break;
    }
    
  }
  readFile.close();
  // Since we start index at 1
  // Return the number of reads
  return readID - 1;
}

#endif /* defined(__NSG__gen_hashes__) */
