#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/edge_array.h"

#include "benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/BallGrowing/BallGrowing.h"

#include <random>
#include <tuple>

// TODO: find good value for c
#define c (0.05)

namespace gbbs
{
  namespace LowAvgStretch
  {

    template <class W>
    struct BFS_C
    {
      uintE *Parents;
      BFS_C(uintE *_Parents) : Parents(_Parents) {}
      inline bool update(uintE s, uintE d, W w)
      {
        if (Parents[d] == UINT_E_MAX)
        {
          Parents[d] = Parents[s];
          return 1;
        }
        else
        {
          return 0;
        }
      }
      inline bool updateAtomic(uintE s, uintE d, W w)
      {
        return (gbbs::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, Parents[s]));
      }
      inline bool cond(uintE d) { return (Parents[d] == UINT_E_MAX); }
    };

    template <class Graph>
    inline void BFS_clustering(Graph &G,
                               sequence<uintE> &srcs,
                               sequence<uintE> &cluster,
                               sequence<size_t> &cluster_index,
                               sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> &edges,
                               parlay::random_generator &gen,
                               double eps)
    {
      using W = typename Graph::weight_type;

      sequence<uintE> ball;
      sequence<uintE> ballShell;
      sequence<uintE> parents;

      sequence<size_t> radius;
      sequence<size_t> dists;
      sequence<size_t> r0s;

      debug(std::cout << "Entering BFS_clustering\n";);
      debug(assert(!srcs.empty()););
      debug(assert(cluster.size() == G.n););
      debug(assert(cluster_index.size() == G.n););
      debug(correctCluster(G, srcs, cluster););
      debug(correctClusterIndex(cluster_index, srcs););
      debug(noDoubleEdges(edges););

      BallGrowing(G, srcs, ball, ballShell, parents, cluster, cluster_index, edges, radius, gen, dists, r0s);
      if (srcs.empty())
      {
        debug(noDoubleEdges(edges););
        // std::cout << "We are finished\n";
        return;
      }

      debug(auto inBall = [&](const uintE &v)
      { return dists[v] <= r0s[cluster_index[cluster[v]]]; };
      auto inBallShell = [&](const uintE &v)
      { return (r0s[cluster_index[cluster[v]]] < dists[v] && dists[v] <= r0s[cluster_index[cluster[v]]] + 1); };

      for (uintE i = 0; i < G.n; i++) {
        if (inBall(i))
        {
          assert(contains(ball, i));
        }
        if (inBallShell(i))
        {
          assert(contains(ballShell, i));
        }
      });

      debug(assert(!srcs.empty()););
      debug(noDoubleEdges(edges););

      // debug(std::cout << "start_node[61]: " << start_nodes[61] << ", start_node[72]: " << start_nodes[72] << std::endl;);
      // debug(std::cout << "cluster[0]: " << cluster[0] << ", cluster[4]: " << cluster[4] << std::endl;);
      // debug(std::cout << "inBallShell(4): " << inBallShell(4) << ", contains(ballShell, 4): " << contains(ballShell, (uintE)4) << std::endl;);
      // debug(
      //     auto pars = BFS(G, 0);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable\n";
      //     } else {
      //       std::cout << "Not reachable\n";
      //     });

      // helper sequence (index to ballShell)
      sequence<uintE> vertices = sequence<uintE>::from_function(ballShell.size(), [&](uintE i)
                                                                { return i; });

      size_t *max_radius = parlay::max_element(radius);
      double beta = (c * log(ballShell.size()) / (eps * *max_radius));
      debug(std::cout << "beta: " << beta << std::endl;);

      // std::random_device rd;
      // std::mt19937 gen(rd());
      std::exponential_distribution<> d(beta);

      sequence<size_t> deltas = sequence<size_t>::from_function(ballShell.size(), [&](size_t i)
                                                                { auto r = gen[i]; return d(r); });

      // debug(std::cout << "ballshell size: " << ballShell.size() << std::endl;);
      // debug(for (uintE i = 0; i < deltas.size(); i++) {
      //   std::cout << "vertex " << ballShell[i] << " has delta: " << deltas[i] << " and cluster: " << cluster[ballShell[i]] << std::endl;
      // });

      size_t *max_delta = parlay::max_element(deltas);

      debug(std::cout << "max_delta: " << *max_delta << std::endl;);

      sequence<size_t> ads = parlay::map(deltas, [&](size_t i)
                                         { return (*max_delta) - i; });

      // for (uintE i = 0; i < deltas.size(); i++)
      // {
      //   std::cout << "vertex " << ballShell[i] << " has ad: " << ads[i] << std::endl;
      // }

      // implicitly denotes cluster (each node contains the node from which the BFS reaching it started)
      auto start_nodes =
          sequence<uintE>::from_function(G.n, [&](size_t i)
                                         { return UINT_E_MAX; });

      // Limit iterations
      max_delta = parlay::max_element(deltas);

      size_t iteration = 0;

      // Function to determine which startnodes to add to frontier in iteration
      auto to_add = [&](size_t i)
      { 
        if (ads[i] == iteration && start_nodes[ballShell[i]] == UINT_E_MAX){
          return ballShell[i];
        } else {return UINT_E_MAX;} };

      // auto new_starts = parlay::filter(ballShell, to_add);
      // new_starts = parlay::map(new_starts, [&](uintE i)
      //                          { return ballShell[i]; });

      auto new_starts = parlay::filter(sequence<uintE>::from_function(ballShell.size(), to_add), [&](uintE i)
                                       { return i != UINT_E_MAX; });

      // std::cout
      //     << "new_starts size: " << new_starts.size() << std::endl;
      // for (uintE i = 0; i < new_starts.size(); i++)
      // {
      //   std::cout << "vertex " << new_starts[i] << " added in iteration: " << iteration << std::endl;
      // }

      vertexSubset Frontier(G.n);
      add_to_vsubset(Frontier, new_starts.begin(), new_starts.size());
      parallel_for(0, new_starts.size(), [&](uintE i)
                   { start_nodes[new_starts[i]] = new_starts[i]; });

      // TODO: Fix this loop to actually work!!!
      while (iteration <= (*max_delta) || !Frontier.isEmpty())
      {
        // std::cout << "iteration\n";
        iteration++;
        if (!Frontier.isEmpty())
        {
          Frontier = edgeMap(G, Frontier, BFS_C<W>(start_nodes.begin()), -1,
                             sparse_blocked | dense_parallel);
        }
        new_starts = parlay::filter(sequence<uintE>::from_function(ballShell.size(), to_add), [&](uintE i)
                                    { return i != UINT_E_MAX; });

        // std::cout << "new_starts.size(): " << new_starts.size() << std::endl;
        if (new_starts.size() > 0)
        {
          // std::cout
          //     << "new_starts size: " << new_starts.size() << std::endl;
          // for (uintE i = 0; i < new_starts.size(); i++)
          // {
          //   std::cout << "vertex " << new_starts[i] << " added in iteration: " << iteration << std::endl;
          // }
          // std::cout << "calling add_to_vsubset\n";
          add_to_vsubset(Frontier, new_starts.begin(), new_starts.size());
          // std::cout << "called add_to_vsubset\n";
          parallel_for(0, new_starts.size(), [&](uintE i)
                       { start_nodes[new_starts[i]] = new_starts[i]; });
        }
      }

      // debug(std::cout << "start_nodes[4]: " << start_nodes[4] << ", UINT_E_MAX: " << UINT_E_MAX << std::endl;);

      // debug(
      //     std::cout << "iteration is now: " << iteration << std::endl; iteration = 62; if (to_add(4)) {
      //       std::cout << "should have added vertex 4\n";
      //     } else {
      //       std::cout << "to_add seems to be wrong!\n";
      //     });

      parallel_for(0, cluster.size(), [&](size_t i)
                   {
        if (start_nodes[i] != UINT_E_MAX)
        {
          cluster[i] = start_nodes[i];
        } });

      auto cut_edges = [&](const uintE &src, const uintE &ngh, const W &wgh)
      {
        if (start_nodes[src] != start_nodes[ngh])
        {
          return 1;
        }
        else
        {
          return 0;
        }
      };

      // debug(std::cout << "start_node[0]: " << start_nodes[0] << ", start_node[4]: " << start_nodes[4] << std::endl;);
      // debug(std::cout << "cluster[0]: " << cluster[0] << ", cluster[4]: " << cluster[4] << std::endl;);

      // debug(
      //     pars = BFS(G, 0);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 0\n";
      //     } else {
      //       std::cout << "Not reachable from 0\n";
      //     });

      // debug(
      //     pars = BFS(G, 25);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 25\n";
      //     } else {
      //       std::cout << "Not reachable from 25\n";
      //     });

      // debug(
      //     pars = BFS(G, 61);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 61\n";
      //     } else {
      //       std::cout << "Not reachable from 61\n";
      //     });

      // debug(
      //     pars = BFS(G, 69);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 69\n";
      //     } else {
      //       std::cout << "Not reachable from 69\n";
      //     });

      // debug(
      //     pars = BFS(G, 116);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 116\n";
      //     } else {
      //       std::cout << "Not reachable from 116\n";
      //     });

      // debug(
      //     pars = BFS(G, 119);
      //     if (isReachable(pars, 4)) {
      //       std::cout << "Is Reachable from 119\n";
      //     } else {
      //       std::cout << "Not reachable from 119\n";
      //     });

      filterEdges(G, cut_edges);

      auto unique_start_nodes = parlay::remove_duplicates_ordered(start_nodes, [&](uintE a, uintE b)
                                                                  { return a < b; });

      // Removes only the elements that are equal to UINT_E_MAX because it is literally flipped to the documentation...
      unique_start_nodes = parlay::remove_if(unique_start_nodes, [&](const uintE v)
                                             { return v != UINT_E_MAX; });

      auto new_edges = sequence<std::tuple<uintE, uintE, typename Graph::weight_type>>::from_function(unique_start_nodes.size(), [&](size_t i)
                                                                                                      { return std::make_tuple(unique_start_nodes[i], parents[unique_start_nodes[i]], gbbs::empty()); });

      debug(std::cout << "appending " << new_edges.size() << " edges\n";);
      edges.append(parlay::make_slice(new_edges));

      srcs.append(parlay::make_slice(unique_start_nodes));
      debug(auto unique_srcs = parlay::remove_duplicates_ordered(srcs, [&](uintE a, uintE b)
                                                                 { return a < b; });
            assert(unique_srcs.size() == srcs.size()););

      debug(assert(isDisjoint(G, srcs)););
      debug(correctCluster(G, srcs, cluster););
      debug(correctClusterIndex(cluster_index, srcs););
      debug(noDoubleEdges(edges););

      debug(std::cout << "exiting BFS_clustering\n";);

      return;
    }
  }
}