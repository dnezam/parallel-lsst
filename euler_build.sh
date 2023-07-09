#!/bin/bash

#SBATCH -n 32
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=1024

./bazel build --config=openmp --compilation_mode=dbg //...
./bazel run @hedron_compile_commands//:refresh_all -- --config=openmp --compilation_mode=dbg