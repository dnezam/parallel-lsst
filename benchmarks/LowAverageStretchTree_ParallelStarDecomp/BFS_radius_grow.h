#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {

template <class W>
struct BFS_rad {
  uintE* Parents;
  BFS_rad(uintE* _Parents) : Parents(_Parents){}
  inline bool update(uintE s, uintE d, W w) {
    if (Parents[d] == UINT_E_MAX) {
      Parents[d] = s;
      return 1;
    } else {
      return 0;
    }
  }
  inline bool updateAtomic(uintE s, uintE d, W w) {
    return (gbbs::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, s));
  }
  inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
};

template <class Graph>
inline std::pair<sequence<uintE>,std::vector<int>> BFS_grow_rad(Graph& G, uintE src, int offset) {
  using W = typename Graph::weight_type;
  /* Creates Parents array, initialized to all -1, except for src. */
  auto Parents = sequence<uintE>::from_function(G.n, [&](size_t i) { return UINT_E_MAX; });
  //std::cout << "Check Edges:\n";
  auto edge_list = to_edge_array<gbbs::empty>(G);
  // for(int i = 0; i < edge_list.E.size(); i++){
  //     std::cout << std::get<0>(edge_list.E[i]) << "_" << std::get<1>(edge_list.E[i]) << " ";
  // }
  // std::cout << "\n";


  Parents[src] = src;
  std::vector<int> dist(G.n,INT_MAX);
  dist[src] = 0;
  //std::cout << " offset " <<  offset;
  vertexSubset Frontier(G.n, src);
  size_t reachable = 0;
  int curr_dist = 1;
  //std::cout << " offset " <<  offset;
  while (!Frontier.isEmpty()) {
    //std::cout << Frontier.size() << "\n";
    reachable += Frontier.size();
    Frontier = edgeMap(G, Frontier, BFS_rad<W>(Parents.begin()), -1,
                       sparse_blocked | dense_parallel);
    // Also possible if we just loop over frontier
    //std::cout << "\nNOW n=" << G.n << "   ";

    // parallel_for(0, G.n, [&](size_t i) {
    //     if(Parents[i] != UINT_E_MAX && dist[i] == INT_MAX){
    //     dist[i] = curr_dist;
    //   }
    //             });

    for(uintE i = 0; i < G.n; i++){
      if(Parents[i] != UINT_E_MAX && dist[i] == INT_MAX){
        //std::cout << Parents[i] << " ";
        dist[i] = curr_dist;
      }
    }
    curr_dist++;
  }
  // std::cout << "Check Parents:\n";
  // for(int i = 0; i < Parents.size(); i++){
  //   std::cout << i << "_" << Parents[i] << " ";
  // }
  // std::cout << "\n";
  //std::cout << "Reachable: " << reachable << "\n";
  return std::make_pair(Parents,dist);
}

}  // namespace gbbs
