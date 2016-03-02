#ifndef NSG_nsg_util_h
#define NSG_nsg_util_h

#include <sys/stat.h>
#include <sys/types.h>


using namespace std;
using read_t =  int ;
using hash_t = long long ; //64-bit hash.  32 bit is too short and will experience too many collisions

struct hash_s {
  hash_t prefix;
  hash_t suffix;
};

// Struct that contains greedy and nsg info
struct Match {
  void SetGreedy(read_t id, int weight) {
    greedy_id = id;
    greedy_weight = weight;
    e_greedy_weight = weight;
    greedy_found = true;
  }

  void SetNSG(read_t id, int weight) {
    nsg_id = id;
    nsg_weight = weight;
    e_nsg_weight = weight;
    nsg_found = true;
  }

  read_t greedy_id = 0;
  read_t nsg_id = 0;
  int greedy_weight = 0;
  // effective greedy weight - weight of edge if read is extended
  int e_greedy_weight = 0;
  int nsg_weight = 0;
  // effective nsg weight - weight of edge if read is extended
  int e_nsg_weight = 0;
  bool greedy_found = false;
  bool nsg_found = false;
  bool x_node = false;
};

// Struct for Edge in graph
struct Edge {
  Edge(read_t u, read_t v, int weight) : u(u), v(v), weight(weight) {}

  bool operator==(const Edge& other) { return u == other.u && v == other.v; }
  bool operator<(const Edge& other) { return weight < other.weight; }

  read_t u;
  read_t v;
  int weight;
};

size_t GetFilesize(const char* filename) {
  struct stat st;
  stat(filename, &st);
  return st.st_size;   
}

#endif
