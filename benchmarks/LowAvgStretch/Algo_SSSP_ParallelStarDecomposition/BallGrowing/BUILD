cc_library(
    name = "BallGrowingLib",
    hdrs = ["BallGrowing.h"],
    deps = [
        "//gbbs",
        "//benchmarks/BFS/NonDeterministicBFS:BFS"
    ],
)

cc_binary(
    name = "BallGrowing",
    srcs = ["BallGrowing.cc"],
    deps = [":BallGrowingLib"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)