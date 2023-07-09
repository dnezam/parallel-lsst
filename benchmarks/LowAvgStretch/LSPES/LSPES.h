#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/contract.h"
// Dear all, please pay attention to this include path
// ALWAYS start from the project root(where the WORKSPACE file is defined)
#include "benchmarks/LowDiameterDecomposition/MPX13/LowDiameterDecomposition.h"
#include "benchmarks/Spanner/MPXV15/Spanner.h"
#include <iostream>

// LSPES - Low-Stretch Parallel Exponential Shift
namespace gbbs {
    namespace LowAvgStretch {
        using edge = gbbs::spanner::edge;

        template <
                template <class W> class vertex, class W,
                // This means W(weight) is not defined, thus unweighted graph
                typename std::enable_if<std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline gbbs::edge_array<gbbs::empty> LSPES(symmetric_graph<vertex, W>& G) {
            // TODO: Add parameter for beta
            // TODO: Figure out what the meaning of "permute" is in LDD
            // TODO: Look at the method signatures of the functions we use and actually
                // understand what the parameters are for (maybe changing the defaults could be useful)
            // TODO: Check whether the memory operations make sense (i.e. I don't copy stuff 
            // unnecessarily)

            
            using cluster_and_parent = gbbs::spanner::cluster_and_parent;
            using unweighted_edge = std::tuple<uintE, uintE, gbbs::empty>;
            using index = size_t;
            using indexed_edge = std::tuple<uintE, uintE, index>;  // start vertex, end vertex, index

            size_t original_n = G.n;

            sequence<unweighted_edge> original_edges = G.edges();
            auto initial_edges = sequence<indexed_edge>(G.m);
            auto current_edges = sequence<indexed_edge>(G.m);
            auto solution_ids = sequence<index>();

            size_t original_edges_size = original_edges.size();
            parallel_for(0, original_edges_size, [&](size_t i) {
                indexed_edge ie = std::make_tuple(std::get<0>(original_edges[i]),
                                                  std::get<1>(original_edges[i]),
                                                  i);
                initial_edges[i] = ie;
                current_edges[i] = ie;
            });

            // TODO: Parallelize the tasks
            while (G.n > 0) {
                // Determine clusters + implicit trees within clusters (taken from Spanner.h)
                auto cluster_and_parents = gbbs::spanner::LDD_parents(G, 0.2, true);

                // Determine edges within clusters (taken from Spanner.h)
                auto tree_edges_with_loops = parlay::delayed_seq<edge>(G.n, [&](size_t i) {
                    return std::make_pair(i, cluster_and_parents[i].parent);
                });
                auto tree_edges = parlay::filter(tree_edges_with_loops, [&](const edge& e) {
                    return e.first != e.second;
                });

                // Add tree edges to solution
                sequence<index> tree_edge_indices = parlay::map(tree_edges, [&](const edge& e) {
                    const auto& m = *std::lower_bound(current_edges.cbegin(), current_edges.cend(),
                                                    std::make_tuple(e.first, e.second, 0));
                    return std::get<2>(m);
                });
                solution_ids = parlay::append(solution_ids, tree_edge_indices);

                // Contract the graph
                auto clusters = parlay::map(cluster_and_parents, [](const cluster_and_parent& cp) {
                    return cp.cluster;
                });
                size_t num_clusters = contract::RelabelIds(clusters);
                
                auto c_out = contract::contract(G, clusters, num_clusters);

                // Prepare current_edges for next iteration
                current_edges = parlay::map(current_edges, [&](const indexed_edge& ie) {
                    uintE u = clusters[std::get<0>(ie)];
                    uintE v = clusters[std::get<1>(ie)];
                    if (u == v) {
                        return std::make_tuple(UINT_E_MAX, UINT_E_MAX, SIZE_MAX);
                    } else {
                        return std::make_tuple(u, v, std::get<2>(ie));
                    }
                });
                current_edges = parlay::filter(current_edges, [&](const indexed_edge& ie) {
                    return std::get<2>(ie) != SIZE_MAX;
                });
                current_edges = parlay::sort(current_edges);

                G = std::move(std::get<0>(c_out));
            }
            // Convert solution ids into an edge array (including the (empty) weight)
            sequence<unweighted_edge> solution_edges = parlay::map(solution_ids, [&](const index& id) -> unweighted_edge {
                return std::make_tuple(std::get<0>(initial_edges[id]), std::get<1>(initial_edges[id]), gbbs::empty());
            });

            return edge_array(std::move(solution_edges), original_n);
        }

        template <
                template <class W> class vertex, class W,
                // This means W(weight) is not undefined, thus weighted graph
                typename std::enable_if<!std::is_same<W, gbbs::empty>::value, int>::type = 0>
        inline parlay::sequence<edge> LSPES(symmetric_graph<vertex, W>& GA) {
            // Weighted Algo here

            return parlay::sequence<edge>();
        }

    }
}