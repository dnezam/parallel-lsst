#!/bin/bash

#SBATCH -n 8
#SBATCH --time=0:05:00
#SBATCH --mem-per-cpu=1024

export OMP_NUM_THREADS=8

./bazel-bin/foreign/foreign_main