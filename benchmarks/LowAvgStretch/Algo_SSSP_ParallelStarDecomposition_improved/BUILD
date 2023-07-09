cc_library(
    name = "SSSP_ParallelStarDecompositionLib_improved",
    hdrs = ["SSSP_ParallelStarDecomposition.h"],
    deps = [
        "//gbbs",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved/BFS_clustering:BFS_clustering_improved",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved/BallGrowing:BallGrowing_improved",
    ],
)

cc_binary(
    name = "SSSP_ParallelStarDecomposition_improved",
    srcs = ["SSSP_ParallelStarDecomposition.cc"],
    deps = [":SSSP_ParallelStarDecompositionLib_improved"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)