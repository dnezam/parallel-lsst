#include "LSPES.h"

namespace gbbs {

    template <class Graph>
    gbbs::edge_array<gbbs::empty> LowAvgStretch_LSPES_Runner(Graph& G, commandLine P) {
        return LowAvgStretch::LSPES(G);
    }

    template <class Graph>
    gbbs::edge_array<gbbs::empty> LowAvgStretch_LSPES_Test(Graph& G) {
        return LowAvgStretch::LSPES(G);
    }
    
}

generate_custom_main(gbbs::LowAvgStretch_LSPES_Runner);
// generate_testing_main(gbbs::LowAvgStretch_LSPES_Test);
