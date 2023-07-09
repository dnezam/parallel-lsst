#include <iostream>
#include "BallGrowing.h"

namespace gbbs
{
    namespace LowAvgStretch
    {

        template <class Graph>
        double BallGrowing_radius_runner(Graph &G, commandLine P)
        {
            sequence<size_t> radius;
            uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
            sequence<uintE> srcs;
            sequence<uintE> parents;
            sequence<size_t> dists;
            sequence<uintE> clusters;
            size_t max_rad;
            sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> edges;
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(10);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);

            srcs.push_back(src);
            srcs.push_back(10);
            timer t;
            t.start();
            LowAvgStretch::BFS_radius(G, srcs, radius, parents, clusters, dists, max_rad, edges);
            double tt = t.stop();

            std::cout << "### Running Time: " << tt << std::endl;
            std::cout << "### Radii: " << std::endl;
            for (int i = 0; i < radius.size(); i++)
            {
                std::cout << "Radius of src " << srcs[i] << " is " << radius[i] << std::endl;
            }

            std::cout << "Distances: \n";
            for (uintE i = 0; i < parents.size(); i++)
            {
                std::cout << "Vertex: " << i << ", distance: " << dists[i] << ", parent: " << parents[i] << std::endl;
            }

            return tt;
        }

        template <class Graph>
        double BallGrowing_runner(Graph &G, commandLine P)
        {

            uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
            sequence<uintE> srcs;
            sequence<uintE> ballShell;
            sequence<uintE> ball;
            sequence<uintE> parents;
            sequence<uintE> clusters;
            sequence<size_t> cluster_index;
            sequence<size_t> radii;
            size_t max_rad;

            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(10);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(0);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);
            clusters.push_back(10);

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
            cluster_index.push_back(1);
            cluster_index.push_back(-1);
            cluster_index.push_back(-1);
            cluster_index.push_back(-1);
            cluster_index.push_back(-1);
            cluster_index.push_back(-1);

            srcs.push_back(src);
            srcs.push_back(10);
            sequence<std::tuple<uintE, uintE, typename Graph::weight_type>> edges;
            parlay::random_generator gen(42);
            timer t;
            std::cout << "src node: " << src << std::endl;
            sequence<size_t> dists;
            sequence<size_t> r0s;
            t.start();
            LowAvgStretch::BallGrowing(G, srcs, ball, ballShell, parents, clusters, cluster_index, edges, radii, max_rad, gen, dists, r0s);
            double tt = t.stop();

            std::cout << "### Running Time: " << tt << std::endl;
            std::cout << "### Ball size: " << ball.size() << ", BallShell size: " << ballShell.size() << std::endl;

            std::cout << "Ball contains: \n";
            for (uintE i = 0; i < ball.size(); i++)
            {
                std::cout << "vertex: " << ball[i] << std::endl;
            }

            std::cout << "BallShell contains: \n";
            for (uintE i = 0; i < ballShell.size(); i++)
            {
                std::cout << "vertex: " << ballShell[i] << std::endl;
            }

            return tt;
        }
    }
}

// generate_symmetric_main(gbbs::LowAvgStretch::BallGrowing_radius_runner, false);
generate_symmetric_once_main(gbbs::LowAvgStretch::BallGrowing_runner, false);