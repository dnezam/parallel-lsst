#pragma once

#include "gbbs/gbbs.h"
#include "gbbs/edge_array.h"

#include "benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition//BFS_clustering/BFS_clustering.h"

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
            debug(assert(isConnected(G, src)););

            clusters = sequence<uintE>::from_function(G.n, [&](size_t i)
                                                      { return src; });
            cluster_index = sequence<size_t>::from_function(G.n, [&](size_t i)
                                                           {
                if (i == src) {
                    return (size_t) 0;
                } 
                return SIZE_MAX; });

            debug(correctClusterIndex(cluster_index, srcs););
            debug(correctCluster(G, srcs, clusters););
            while (!srcs.empty())
            {
                BFS_clustering(G, srcs, clusters, cluster_index, edges, gen, 0.2);
            }
            return edge_array<typename Graph::weight_type>(std::move(edges), G.n);
        }
    }
}