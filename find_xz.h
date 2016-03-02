// Second part of the NSG algorithm: Transforming X and removing Z structures
#ifndef NSG_find_xz_h
#define NSG_find_xz_h

#include "nsg_util.h"

#include <map>
#include <vector>
#include <thread>
#include <algorithm>
#include <atomic>
#include <sstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <set>

namespace {
//==========================================================================================
// HELPER FUNCTIONS for FindXZStructures

// Count out-degrees of a given read node
int OutDegrees(Match m) {
  int out_degree = 0;
  if (m.greedy_id != 0) {
    out_degree++;
  }
  if (m.nsg_id != 0) {
    out_degree++;
  }
  return out_degree;
}

// Remove edge from the graph and its reverse_graph
void RemoveEdge(read_t u, read_t v, vector<Match>* graph,
                vector<Match>* reverse_graph) {
  if (v == (*graph)[u].greedy_id) {
    (*graph)[u].greedy_id = 0;
    (*graph)[u].greedy_weight = 0;
    (*graph)[u].e_greedy_weight = 0;
  } else {
    (*graph)[u].nsg_id = 0;
    (*graph)[u].nsg_weight = 0;
    (*graph)[u].e_nsg_weight = 0;
  }
  if (u == (*reverse_graph)[v].greedy_id) {
    (*reverse_graph)[v].greedy_id = 0;
    (*reverse_graph)[v].greedy_weight = 0;
    (*reverse_graph)[v].e_greedy_weight = 0;
  } else {
    (*reverse_graph)[v].nsg_id = 0;
    (*reverse_graph)[v].nsg_weight = 0;
    (*reverse_graph)[v].e_nsg_weight = 0;
  }
}

// Get the ith nucleotide from the read with id = read_id
char GetIthNucleotideFromRead(char* mmapptr, int read_length,
                              const read_t read_id, int i) {
  // We offset readLength by 1 to account for new line characters
  return mmapptr[(read_length + 1) * (read_id - 1) + i];
}

// Extend the read according to the Read Extension rule
void ReadExtend(bool right_oriented, char* mmapptr, int read_length,
                read_t read_id, vector<Match>* graph) {
  Match match = (*graph)[read_id];
  int extension = 0;
  for (int i = 0; i < read_length - match.e_greedy_weight; i++) {
    if (right_oriented) {
      if (GetIthNucleotideFromRead(mmapptr, read_length, match.greedy_id,
                                   i + match.e_greedy_weight) ==
          GetIthNucleotideFromRead(mmapptr, read_length, match.nsg_id,
                                   i + match.e_nsg_weight)) {
        extension++;
      } else {
        break;
      }
    } else {
      if (GetIthNucleotideFromRead(mmapptr, read_length, match.greedy_id,
                                   read_length - match.e_greedy_weight - i - 1) ==
          GetIthNucleotideFromRead(mmapptr, read_length, match.nsg_id,
                                   read_length - match.e_nsg_weight - i - 1)) {
        extension++;
      } else {
        break;
      }
    }
  }
  // Extend the effective overlap of the reads
  if ((*graph)[match.greedy_id].greedy_id == read_id) {
    (*graph)[match.greedy_id].e_greedy_weight += extension;
  } else if ((*graph)[match.greedy_id].nsg_id == read_id) {
    (*graph)[match.greedy_id].e_nsg_weight += extension;
  }
  if ((*graph)[match.nsg_id].greedy_id == read_id) {
    (*graph)[match.nsg_id].e_greedy_weight += extension;
  } else if ((*graph)[match.nsg_id].nsg_id == read_id) {
    (*graph)[match.nsg_id].e_nsg_weight += extension;
  }
}

// Find all Edges (u,v) such that u has outdegree 2 and
// v has indegree 2.
vector<Edge> FindSetU(bool read_extend, char* mmapptr, int read_length,
                      vector<Match>* graph, vector<Match>* reverse_graph) {
  vector<Edge> set_u;
  read_t read_id;
  for (read_id = 1; read_id <= graph->size(); read_id++) {
    Match node = (*graph)[read_id];
    if (!node.x_node) {
      if (OutDegrees(node) == 2) {
        if (read_extend) {
          ReadExtend(true, mmapptr, read_length, read_id, graph);
        }
        if (OutDegrees((*reverse_graph)[node.greedy_id]) == 2) {
          if (read_extend) {
            ReadExtend(false, mmapptr, read_length, node.greedy_id,
                       reverse_graph);
          }
          Edge e(read_id, node.greedy_id, node.e_greedy_weight);
          set_u.push_back(e);
        }

        if (OutDegrees((*reverse_graph)[node.nsg_id]) == 2) {
          if (read_extend) {
            ReadExtend(false, mmapptr, read_length, node.nsg_id, reverse_graph);
          }
          Edge e(read_id, node.nsg_id, node.e_nsg_weight);
          set_u.push_back(e);
        } 
      }
    }
  }
  return set_u;
}

}  // namespace

// Identify X and Z Structures for the second half of the NSG algorithm
void FindXZStructures(char* mmapptr, int read_length, int num_reads, bool read_extend,
                      vector<Match>* graph, vector<Match>* reverse_graph) {
  vector<Edge> set_u = FindSetU(read_extend, mmapptr, read_length, graph, reverse_graph);
  set<pair<int, int>> fixed_edges;
  int count = 0;
  while (set_u.size() > 0) {
    count++;
    cout << "size of set U: " << set_u.size() << endl;
    for (int i = 0; i < set_u.size(); i++) {
      read_t u2 = set_u[i].u;
      read_t v1 = set_u[i].v;
      int weight = set_u[i].weight;

      Match& u2_match = (*graph)[u2];
      Match& v1_reverse_match = (*reverse_graph)[v1];

      if (OutDegrees(u2_match) != 2 || OutDegrees(v1_reverse_match) != 2) {
        continue;
      }
      int v2_weight = 0;
      int u1_weight = 0;
      read_t v2, u1;
      if (v1 == u2_match.greedy_id) {
        v2_weight = u2_match.e_nsg_weight;
        weight = u2_match.e_greedy_weight;
        v2 = u2_match.nsg_id;
      } else {
        v2_weight = u2_match.e_greedy_weight;
        v2 = u2_match.greedy_id;
      }
      if (u2 == v1_reverse_match.greedy_id) {
        u1_weight = v1_reverse_match.e_nsg_weight;
        u1 = v1_reverse_match.nsg_id;
      } else {
        u1_weight = v1_reverse_match.e_greedy_weight;
        u1 = v1_reverse_match.greedy_id;
      }

      pair<int, int> u1v1_edge = make_pair(u1, v1);
      pair<int, int> u2v2_edge = make_pair(u2, v2);
      Match& u1_match = (*graph)[u1];
      Match& v2_reverse_match = (*reverse_graph)[v2];

      // Check if an X structure is identified
      if (u1_match.greedy_id == v2 || u1_match.nsg_id == v2) {
        int u1v2_weight = u1_match.e_greedy_weight;
        if (u1_match.nsg_id == v2) {
          u1v2_weight = u1_match.e_nsg_weight;
        }
        if (weight == u1_weight && u1_weight == v2_weight &&
            v2_weight == u1v2_weight) {
          cout << "X structure found among: " << endl;
          cout << "u1: " << u1 << " u2: " << u2 << " v1: " << v1
               << " v2: " << v2 << endl;

          num_reads++;
          Match x_node;
          x_node.SetGreedy(v1, u1_match.e_greedy_weight);
          x_node.SetNSG(v2, u1_match.e_greedy_weight);
          x_node.x_node = true;
          graph->push_back(x_node);
          u1_match.SetGreedy(num_reads, u1_match.e_greedy_weight);
          u1_match.SetNSG(0, 0);
          u2_match.SetGreedy(num_reads, u1_match.e_greedy_weight);
          u2_match.SetNSG(0, 0);
          // X node found but do not have same overlap
        } else {
          pair<int, int> u1v2_edge = make_pair(u1, v2);
          pair<int, int> u2v1_edge = make_pair(u2, v1);
          // Remove the cross edges in X structure
          if (weight <= v2_weight && weight <= u1_weight) {
            RemoveEdge(u1, v2, graph, reverse_graph);
            RemoveEdge(u2, v1, graph, reverse_graph);
            fixed_edges.insert(u1v1_edge);
            fixed_edges.insert(u2v2_edge);
          } else {
            RemoveEdge(u1, v1, graph, reverse_graph);
            RemoveEdge(u2, v2, graph, reverse_graph);
            fixed_edges.insert(u1v2_edge);
            fixed_edges.insert(u2v1_edge);
          }
        }
        // Check strictly less condition for edges and if
        // a Z structure can be found and removed
      } else if (weight != v2_weight || weight != u1_weight) {
        if (weight <= v2_weight && weight <= u1_weight ||
            fixed_edges.find(u1v1_edge) != fixed_edges.end() &&
                fixed_edges.find(u2v2_edge) != fixed_edges.end() ||
            fixed_edges.find(u1v1_edge) != fixed_edges.end() &&
                weight <= v2_weight ||
            fixed_edges.find(u2v2_edge) != fixed_edges.end() &&
                weight <= u1_weight) {
          RemoveEdge(u2, v1, graph, reverse_graph);
          fixed_edges.insert(u1v1_edge);
          fixed_edges.insert(u2v2_edge);
          // If Z structure is isolated (no zig-zag), remove middle edge
        } else if (OutDegrees(u1_match) == 1 && OutDegrees(v2_reverse_match) == 1) {
          RemoveEdge(u2, v1, graph, reverse_graph);
          fixed_edges.insert(u1v1_edge);
          fixed_edges.insert(u2v2_edge);
        } else if (OutDegrees(u1_match) == 1) {
          fixed_edges.insert(u1v1_edge);
        } else if (OutDegrees(v2_reverse_match) == 1) {
          fixed_edges.insert(u2v2_edge);
        }
        // There is a cyclical zigzag with all weights the same
      } else {
        RemoveEdge(u2, v1, graph, reverse_graph);
        fixed_edges.insert(u1v1_edge);
        fixed_edges.insert(u2v2_edge);
      }
    }
    set_u = FindSetU(read_extend, mmapptr, read_length, graph, reverse_graph);
  }
  cout << endl;
}

#endif
