package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "tools",
    srcs = ["common.cxx"],
    hdrs = glob(["*.hpp"]),
    deps = [
        "//gbbs",
        "//gbbs:graph_io",
        "//lemon",
        "//networkit",
    ],
)

cc_binary(
    name = "nx2gbbs",
    srcs = ["nx2gbbs.cxx"],
    deps = [":tools"],
)

cc_binary(
    name = "test_stretch_gbbs",
    srcs = ["test_stretch_gbbs.cxx"],
    deps = [":tools"],
)

cc_binary(
    name = "test_stretch_nx",
    srcs = ["test_stretch_nx.cxx"],
    deps = [":tools"],
)
