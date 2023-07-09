cc_library(
    name = "BallGrowingLib_improved",
    hdrs = ["BallGrowing.h"],
    deps = [
        "//gbbs",
        "//benchmarks/BFS/NonDeterministicBFS:BFS"
    ],
)

cc_binary(
    name = "BallGrowing_improved",
    srcs = ["BallGrowing.cc"],
    deps = [":BallGrowingLib_improved"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)