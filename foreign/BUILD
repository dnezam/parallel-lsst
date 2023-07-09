package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "foreign_lib",
    srcs = [
        "abraham_parallell.cpp",
        "abraham_seq.cpp",
        "alonkarp.cpp",
        "basics.cpp",
        "calcstretch.cpp",
        "correction.cpp",
        "lowstretchtree.cpp",
    ],
    hdrs = glob(["*.h"]),
    deps = [
        "//gbbs/helpers:parse_command_line",
        "@parlaylib//parlay:io",
        "@parlaylib//parlay:primitives",
        "@parlaylib//parlay:random",
        "//lemon",
        "//tools",
    ],
)

cc_binary(
    name = "test",
    srcs = ["par_main.cpp"],
    deps = [":foreign_lib"],
)

cc_binary(
    name = "seq_weighted_alon",
    srcs = ["seq_weighted_alon.cpp"],
    deps = [
        ":foreign_lib",
    ],
)

cc_binary(
    name = "seq_unweighted_alon",
    srcs = ["seq_unweighted_alon.cpp"],
    deps = [
        ":foreign_lib",
    ],
)

cc_binary(
    name = "seq_weighted_abraham",
    srcs = ["seq_weighted_abraham.cpp"],
    deps = [
        ":foreign_lib",
    ],
)

cc_binary(
    name = "par_weighted_abraham",
    srcs = ["par_weighted_abraham.cpp"],
    deps = [
        ":foreign_lib",
    ],
)

cc_binary(
    name = "seq_weighted_emel",
    srcs = ["seq_weighted_emel.cpp"],
    deps = [
        ":foreign_lib",
    ],
)

cc_binary(
    name = "seq_weighted_kruskal",
    srcs = ["seq_weighted_kruskal.cpp"],
    deps = [
        ":foreign_lib",
    ],
)
