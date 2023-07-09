#include "benchmarks/BFS/NonDeterministicBFS/BFS.h"
#include <random>

namespace gbbs {
    template <class Graph>
    auto LowAvgStretch_Algo_Spanner_Runner(Graph& G, commandLine P) {
        const int range_from  = 0;
        const int range_to    = G.n - 1;
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<gbbs::uintE> distr(range_from, range_to);
        const gbbs::uintE rand_id = distr(generator);
        const gbbs::uintE source{rand_id};
        auto seq = gbbs::BFS(G, source);
        std::vector<gbbs_io::Edge<gbbs::empty>> tree{};
        for (gbbs::uintE i = 0; i < seq.size(); i++) {
            if (i != seq[i]) {
                tree.push_back({i, seq[i]});
            }
        }
        // std::cout << seq.size() << std::endl;
        return tree;
    }
}

generate_custom_main(gbbs::LowAvgStretch_Algo_Spanner_Runner);