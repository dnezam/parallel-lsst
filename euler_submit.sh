#!/bin/bash

#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=20960

export OMP_NUM_THREADS=24;
# ./bazel-bin/benchmarks/LowAvgStretch/Algo_spanner/RandomSpanner -rounds 10  -s inputs/ErdosRenyi_Bench/ErdosRenyi_5000_2500582_BIN
# bazel-bin/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/SSSP_ParallelStarDecomposition -s inputs/ErdosRenyi_Bench/ErdosRenyi_Bench_205000_5121926_BIN

./bazel-bin/foreign/seq_unweighted_alon -rounds 2 final_inputs/strong_scale_max/ErdosRenyi_10000000_100004838_BIN final_inputs/strong_scale_max/out/