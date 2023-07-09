# low-avg-stretch

## Files created/modified by us

* [low-avg-stretch/benchmarks/LowAvgStretch/LSPES](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAvgStretch/LSPES) (created) - *Contains implementation for our algorithm LSPES*
* [low-avg-stretch/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition) (created) - *Contains implementation for our Star Decomposition used in Benchmarking*
* [low-avg-stretch/benchmarks/LowAvgStretch/Algo_spanner](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAvgStretch/Algo_spanner) (created) - *Contains implementation for our random BFS*
* [low-avg-stretch/benchmarks/LowAverageStretchTree_ParallelStarDecomp](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAverageStretchTree_ParallelStarDecomp) (created) - *Contains implementation for our algorithm AKPW*
* [low-avg-stretch/gbbs/graph.h](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/blob/main/gbbs/graph.h) (modified) - *Fix bugs in GBBS library; see [GitHub Issue: Missing return value in graph.h (Undefined Behavior?)](https://github.com/ParAlg/gbbs/issues/77), [GitHub Issue: deletion_fn() is called too late, resulting in SEGFAULT](https://github.com/ParAlg/gbbs/issues/80)*
* [gbbs/benchmark.h](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/blob/main/gbbs/benchmark.h) (modified) - *Contains benchmarking infrastructure for algorithms in `/benchmarks/LowAvgStretch`*
* [low-avg-stretch/foreign/abraham_parallell.cpp](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/blob/main/foreign/abraham_parallell.cpp) (created) - *Contains implementation for our implementation of weighted Abraham (not evaluated)*
* [low-avg-stretch/foreign/seq_unweighted_alon.cpp](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/blob/main/foreign/seq_unweighted_alon.cpp) (created) - *Contains implementation of sequential Alon Karp*
* [low-avg-stretch/foreign/benchmark.h](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/blob/main/foreign/benchmark.h) (created) - *Contains benchmarking infrastructure for algorithms in `/foreign`*
* [low-avg-stretch/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved) (created) - *Contains improved implementation for our Star Decomposition resulting in approximately 2x speedup as discussed in the report.*

## Important Notes

Our AKPW algorithm has an unfortunate name in our codebase. This is due to some early confusion and was kept. As seen above it is located in [low-avg-stretch/benchmarks/LowAverageStretchTree_ParallelStarDecomp](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/benchmarks/LowAverageStretchTree_ParallelStarDecomp). Additionally in the benchmarking output it is named `StarDecomp`, while our Star Decomposition is named `SSSP_StarDecomp` there.

Also, the improved Star Decomposition is not included in the benchmarking script. If you want to use the script to benchmark it you would have to modify line 32 to the location of the improved version.

## Requirements

We used the following software:

* GCC version 8.1
* Python version 3.8
* Bazel version 1.15

Additionally we used Parlay in combination with GBBS where the source code is included in the repo.

This repo also contains a Dockerfile that sets up the execution environment and can be used with VS Code (with also [Remote Container Extension](https://code.visualstudio.com/docs/devcontainers/containers)).

Please note that if you are not using the Docker image, make sure that you install the correct bazel version or use the include executable `./bazel`.

## Building the code

All the code in this repository is to be built using Bazel. To do this use the commands listed below. Note that if you are using the provided executable, prepend a `./` before each command.

* ```sh
  bazel build //benchmarks/LowAvgStretch/Algo_spanner:RandomSpanner
    ```

* ```sh
  bazel build //benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition:SSSP_ParallelStarDecomposition
    ```

* ```sh
  bazel build //benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition_improved:SSSP_ParallelStarDecomposition_improved
    ```

* ```sh
  bazel build //benchmarks/LowAvgStretch/LSPES:RunLSPES
    ```

* ```sh
  bazel build //benchmarks/LowAverageStretchTree_ParallelStarDecomp:LAS_ParaStar
    ```

* ```sh
  bazel build --config=bench //foreign:seq_unweighted_alon
    ```

Please note that the following are also available, but were not benchmarked.

* ```sh
  bazel build --config=bench //foreign:par_weighted_abraham
    ```

* ```sh
  bazel build --config=bench //foreign:seq_weighted_alon
    ```

* ```sh
  bazel build --config=bench //foreign:seq_weighted_abraham
    ```

* ```sh
  bazel build --config=bench //foreign:seq_weighted_emel
    ```

* ```sh
  bazel build --config=bench //foreign:seq_weighted_kruskal
    ```

The following will have to be called to build the tools for calculating stretch:

```sh
bazel build --config=bench //tools:*
```

The executable binary will be placed in the directory `bazel-bin/`.

## Running code

### Inputs

For each graph we have two input formats: Binary(filename ends with `_BIN`), and GBBS(filename ends with `_GBBS`). This is due to the different input formats required by code from the directory `foreign` and GBBS. Any algorithm in the directory `foreign` takes the Binary format as input.

### Commands

With some of the algorithms we can run the stretch calculation seperately. As such they will need to specify an output directory to write the output graph into. Others will directly do the stretch calculation internally and print it to commandline. All algorithms will print the runtime to the commandline.

* ```sh
  ./bazel-bin/benchmarks/LowAvgStretch/Algo_spanner/RandomSpanner [-rounds x] <path_to_inputfile(GBBS)> <path_to_outdir>
    ```

* ```sh
  ./bazel-bin/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/SSSP_ParallelStarDecomposition [-rounds x] <path_to_inputfile(GBBS)> <path_to_outdir>
    ```

* ```sh
  ./bazel-bin/benchmarks/LowAvgStretch/LSPES/RunLSPES [-rounds x] <path_to_inputfile(GBBS)> <path_to_outdir>
    ```

* ```sh
  ./bazel-bin/benchmarks/LowAverageStretchTree_ParallelStarDecomp/LAS_ParaStar [-rounds x] <path_to_inputfile(GBBS)> <path_to_outdir>
    ```

* ```sh
  ./bazel-bin/foreign/seq_unweighted_alon [-rounds x] <path_to_inputfile(BIN)> <path_to_outdir>
    ```

For the not benchmarked algorithms:

* ```sh
  ./bazel-bin/foreign/par_weighted_abraham [-rounds x] <path_to_inputfile(BIN)>
    ```

* ```sh
  ./bazel-bin/foreign/seq_weighted_alon [-rounds x] <path_to_inputfile(BIN)>
    ```

* ```sh
  ./bazel-bin/foreign/seq_weighted_abraham [-rounds x] <path_to_inputfile(BIN)>
    ```

* ```sh
  ./bazel-bin/foreign/seq_weighted_emel [-rounds x] <path_to_inputfile(BIN)>
    ```

* ```sh
  ./bazel-bin/foreign/seq_weighted_kruskal [-rounds x] <path_to_inputfile(BIN)>
    ```

For the stretch calculation choose the one according to the input format of the algorithm used:

* ```sh
  ./bazel-bin/tools/test_stretch_gbbs <original_input_graph(GBBS)> <algorithm_output_dir>/
    ```

* ```sh
  ./bazel-bin/tools/test_stretch_nx <original_input_graph(BIN)> <algorithm_output_dir>/
    ```

Please note that the output directory should only hold trees produced from the same input graph. Also the script will calculate the stretch for all trees in the directory.

## Benchmarking inputs

All our benchmarking inputs are located in the branches `all_data` and `max_data` using `git-lsf`. In particular `max_data` holds our inputs for strong scaling in the directory [test_data](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/max_data/test_data), and `all_data` holds our inputs for weak scaling in the directory [weak_scaling_data](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/all_data/weak_scaling_data)

## Scripts

We have included several scripts in the directory [scripts/dphpc_benchmarking](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/scripts/dphpc_benchmarking). This includes a script for automated benchmarking called `run_benchmarks.py` and three scripts to produce plots from the output data called `plot_single_graph_multiple_algos.py`, `plot_snap2.py`, and `plot_weak_scaling-py`

## Results

We have included our benchmarking results in the directory [scripts/dphpc_benchmarking/benchmarks](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/scripts/dphpc_benchmarking/benchmarks) along with some plots in [scripts/dphpc_benchmarking/plot_output](https://gitlab.ethz.ch/dnezamabadi/low-avg-stretch/-/tree/main/scripts/dphpc_benchmarking/plot_output)
