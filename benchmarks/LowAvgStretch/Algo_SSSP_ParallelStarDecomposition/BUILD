cc_library(
    name = "SSSP_ParallelStarDecompositionLib",
    hdrs = ["SSSP_ParallelStarDecomposition.h"],
    deps = [
        "//gbbs",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/BFS_clustering:BFS_clustering",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/BallGrowing:BallGrowing",
    ],
)

cc_binary(
    name = "SSSP_ParallelStarDecomposition",
    srcs = ["SSSP_ParallelStarDecomposition.cc"],
    deps = [":SSSP_ParallelStarDecompositionLib"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)