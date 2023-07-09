#include "BFS_clustering.h"

namespace gbbs
{

    template <class Graph>
    double BFS_clustering_runner(Graph &G, commandLine P)
    {
        sequence<uintE> srcs;
        sequence<uintE> clusters;
        sequence<size_t> cluster_index;

        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);
        clusters.push_back(0);

        cluster_index.push_back(0);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        cluster_index.push_back(-1);
        uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
        srcs.push_back(src);

        sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> edges;

        parlay::random_generator gen(42);

        timer t;
        t.start();
        LowAvgStretch::BFS_clustering(G, srcs, clusters, cluster_index, edges, gen, 0.2);
        double tt = t.stop();

        std::cout << "### Running Time: " << tt << std::endl;
        std::cout << "### srcs size: " << srcs.size() << std::endl;

        std::cout << "cluster sources with the nodes they are connected to from other clusters\n";
        for (size_t i = 0; i < srcs.size(); i++)
        {
            std::cout << "vertex: " << srcs[i] << std::endl;
        }

        std::cout << "Edges:\n";
        for (int i = 0; i < edges.size(); i++) {
            std::cout << "{" << std::get<0>(edges[i]) << "," << std::get<1>(edges[i]) << "}\n";
        }
        return tt;
    }
}

generate_symmetric_once_main(gbbs::BFS_clustering_runner, false);