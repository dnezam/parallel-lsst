#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/edge_array.h"

#include "benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved/BFS_clustering/BFS_clustering.h"

#include <tuple>

namespace gbbs
{
    namespace LowAvgStretch
    {
        template <class Graph>
        inline edge_array<typename Graph::weight_type> SSSP_LowStretchTree(Graph &G, uintE src)
        {
            sequence<uintE> srcs;
            sequence<uintE> clusters;
            sequence<size_t> cluster_index;
            sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> edges;
            parlay::random_generator gen;

            srcs.push_back(src);
            // // debug(assert(isConnected(G, src)););

            clusters = sequence<uintE>::from_function(G.n, [&](size_t i)
                                                      { return src; });
            // This assumes a connected graph and as such we start wiht one cluster
            cluster_index = sequence<size_t>::from_function(G.n, [&](size_t i)
                                                            { return src; });

            // debug(correctClusterIndex(cluster_index, srcs););
            // debug(correctCluster(G, srcs, clusters););
            int iteration = 0;
            while (!srcs.empty())
            {
                parlay::random_generator gen(iteration);
                BFS_clustering(G, srcs, clusters, cluster_index, edges, gen, 0.2);
                iteration++;
            }
            return edge_array<typename Graph::weight_type>(std::move(edges), G.n);
        }
    }
}