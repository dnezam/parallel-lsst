cc_library(
    name = "BFS_clusteringLib",
    hdrs = ["BFS_clustering.h"],
    deps = [
        "//gbbs",
        "//benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/BallGrowing:BallGrowing",
        "//benchmarks/BFS/NonDeterministicBFS:BFS"
    ],
)

cc_binary(
    name = "BFS_clustering",
    srcs = ["BFS_clustering.cc"],
    deps = [":BFS_clusteringLib"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)