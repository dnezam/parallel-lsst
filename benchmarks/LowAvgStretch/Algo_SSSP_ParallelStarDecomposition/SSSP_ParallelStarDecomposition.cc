#include "SSSP_ParallelStarDecomposition.h"

namespace gbbs {

    template <class Graph>
    inline edge_array<typename Graph::weight_type> SSSP_ParallelStarDecomposition_runner(Graph& G, commandLine P) {
    // inline auto SSSP_ParallelStarDecomposition_runner(Graph& G, commandLine P) {

        uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
        edge_array<typename Graph::weight_type> edge_arr = LowAvgStretch::SSSP_LowStretchTree(G, src);

        // debug(std::cout << "Graph size: " << G.n << ", edges: " << edge_arr.E.size() << std::endl);
        // debug(for (std::size_t i = 0; i < edge_arr.E.size(); i++) {
        //     std::cout << "{" << std::get<0>(edge_arr.E[i]) << "," << std::get<1>(edge_arr.E[i]) << "}\n";
        // })
        return edge_arr;
    }

    template <class Graph>
    inline edge_array<typename Graph::weight_type> SSSP_ParallelStarDecomposition_test_runner(Graph& G) {
    // inline auto SSSP_ParallelStarDecomposition_runner(Graph& G, commandLine P) {

        uintE src = 0;
        edge_array<typename Graph::weight_type> edge_arr = LowAvgStretch::SSSP_LowStretchTree(G, src);

        // debug(std::cout << "Graph size: " << G.n << ", edges: " << edge_arr.E.size() << std::endl);
        // debug(for (std::size_t i = 0; i < edge_arr.E.size(); i++) {
        //     std::cout << "{" << std::get<0>(edge_arr.E[i]) << "," << std::get<1>(edge_arr.E[i]) << "}\n";
        // })
        return edge_arr;
    }
}

generate_custom_main(gbbs::SSSP_ParallelStarDecomposition_runner);
// generate_testing_main(gbbs::SSSP_ParallelStarDecomposition_test_runner);