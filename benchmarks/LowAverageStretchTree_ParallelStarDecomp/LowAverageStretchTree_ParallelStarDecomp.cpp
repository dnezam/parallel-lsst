#include "LowAverageStretchTree_ParallelStarDecomp.h"

// ./bazel-bin/benchmarks/LowAverageStretchTree_ParallelStarDecomp/LAS_ParaStar -s ./inputs/ErdosRenyi/ErdosRenyi_5000_2500582_BINgbbs::SSSP_ParallelStarDecomposition_runner_0
// ./bazel-bin/benchmarks/LowAverageStretchTree_ParallelStarDecomp/LAS_ParaStar -s ./inputs/ErdosRenyi/ErdosRenyi_20_39_BIN

namespace gbbs {
    
    using edge = std::tuple<uintE, uintE>;
    template <class Graph>
    gbbs::edge_array<gbbs::empty> LowAvgStretch_ParaStar(Graph& G, commandLine P) {
        int n = G.n;
        auto res = LowAvgStretch::ParaStarDecomp(G);
        
        auto edge_sequence =  parlay::map(res, [](auto edge) {
            return std::tuple_cat(edge, std::make_tuple(gbbs::empty()));
        });

        res = parlay::sort(res);

        // // std::cout << "Hand_in size: " << edge_sequence.size() << "\n";
        // // std::cout << "Check reverse edges\n";
        // for(int i = 0; i < res.size(); i++){
        //     if(std::find(begin(res),end(res),std::make_tuple(std::get<1>(res[i]),std::get<0>(res[i])))){
        //         //std::cout << std::get<0>(res[i]) << "_" << std::get<1>(res[i]) << " ";
        //     }
        // }
        //std::cout << "\n";

        return edge_array(std::move(edge_sequence), n);
    }

    template <class Graph>
    gbbs::edge_array<gbbs::empty> LowAvgStretch_ParaStar_Test(Graph& G) {
        int n = G.n;
        auto res = LowAvgStretch::ParaStarDecomp(G);
        
        auto edge_sequence =  parlay::map(res, [](auto edge) {
            return std::tuple_cat(edge, std::make_tuple(gbbs::empty()));
        });

        res = parlay::sort(res);

        // // std::cout << "Hand_in size: " << edge_sequence.size() << "\n";
        // // std::cout << "Check reverse edges\n";
        // for(int i = 0; i < res.size(); i++){
        //     if(std::find(begin(res),end(res),std::make_tuple(std::get<1>(res[i]),std::get<0>(res[i])))){
        //         //std::cout << std::get<0>(res[i]) << "_" << std::get<1>(res[i]) << " ";
        //     }
        // }
        //std::cout << "\n";

        return edge_array(std::move(edge_sequence), n);
    }
}

generate_custom_main(gbbs::LowAvgStretch_ParaStar);
//generate_symmetric_main(gbbs::LowAvgStretch_ParaStar, false);
//generate_testing_main(gbbs::LowAvgStretch_ParaStar_Test);