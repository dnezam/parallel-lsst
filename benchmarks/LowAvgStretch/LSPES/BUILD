cc_library(
    name = "LSPES",
    hdrs = ["LSPES.h"],
    deps = [
        "//gbbs",
        # I guess you are using this?
        "//benchmarks/LowDiameterDecomposition/MPX13:LowDiameterDecomposition",
        "//benchmarks/Spanner/MPXV15:Spanner",
        "//gbbs:contract",
    ],
)

cc_binary(
    name = "RunLSPES",
    srcs = ["LSPES.cpp"],
    deps = [":LSPES"],
)

package(
    # So that others can see your cc_library
    default_visibility = ["//visibility:public"],
)