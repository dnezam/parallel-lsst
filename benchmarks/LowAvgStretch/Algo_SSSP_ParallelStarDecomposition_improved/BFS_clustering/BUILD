cc_library(
    name = "BFS_clusteringLib_improved",
    hdrs = ["BFS_clustering.h"],
    deps = [
        "//gbbs",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved/BallGrowing:BallGrowing_improved",
        "//benchmarks/BFS/NonDeterministicBFS:BFS"
    ],
)

cc_binary(
    name = "BFS_clustering_improved",
    srcs = ["BFS_clustering.cc"],
    deps = [":BFS_clusteringLib_improved"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)