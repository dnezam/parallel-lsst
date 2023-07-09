cc_library(
    name = "LowAvgStretch_ParallelSimple",
    hdrs = ["LowAverageStretchTree_ParallelStarDecomp.h", "BFS_radius_grow.h"],
    deps = [
        "//gbbs",
        "//benchmarks/BFS/NonDeterministicBFS:BFS",
        "//gbbs:contract",
    ],
)

cc_binary(
    name = "LAS_ParaStar",
    srcs = ["LowAverageStretchTree_ParallelStarDecomp.cpp"],
    deps = [":LowAvgStretch_ParallelSimple"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)