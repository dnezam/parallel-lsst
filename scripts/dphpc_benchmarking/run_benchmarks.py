import numpy as np
import pandas as pd
import scipy.stats as sp
import matplotlib.pyplot as plt
import subprocess
import time
import os
import glob
import logging

from math import ceil, floor, sqrt

# Keep LARGER THAN 7 (> 7)
MIN_ITERATIONS = 8
MAX_ITERATIONS = 50

PATH_TO_GRAPHS = 'test_data/' 
OUTPUT_PATH = 'out_dir/'

USE_MEDIAN = True

CONFIDENCE = 0.95
WITHIN = 0.05

AVAIL_PROCESSORS = 120

EVALUATION_PROCESSORS = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 120]
# EVALUATION_PROCESSORS = [1, 4, 8,16]
# list of algorithms with the binary path, and the build path, and name, and needsOpenMP(bool), testing method (1:Group 1, 2:Group2A, 3:Group2B)
PARA_ALGORITHMS = [
    ['benchmarks/LowAvgStretch/Algo_spanner/RandomSpanner', 'benchmarks/LowAvgStretch/Algo_spanner:RandomSpanner', 'RandBFS', False, 2],
    ['benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/SSSP_ParallelStarDecomposition', 'benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition:SSSP_ParallelStarDecomposition', 'SSSP_StarDecomp', False, 2],
    ['benchmarks/LowAvgStretch/LSPES/RunLSPES', 'benchmarks/LowAvgStretch/LSPES:RunLSPES', 'LSPES', False, 2],
    ['benchmarks/LowAverageStretchTree_ParallelStarDecomp/LAS_ParaStar', 'benchmarks/LowAverageStretchTree_ParallelStarDecomp:LAS_ParaStar', 'StarDecomp', True, 2],
#    ['foreign/par_weighted_abraham','foreign:par_weighted_abraham', 'ParAbrahamBartal', True, 1],
] 

SEQ_ALGORITHMS = [
 #   ['foreign/seq_weighted_alon', 'foreign:seq_weighted_alon', 'AlonKarp', True, 1],
    ['foreign/seq_unweighted_alon', 'foreign:seq_unweighted_alon', 'UnweightedAlonKarp', True, 3],
 #   ['foreign/seq_weighted_abraham', 'foreign:seq_weighted_abraham', 'AbrahamBartal', True, 1],
 #   ['foreign/seq_weighted_emel', 'foreign:seq_weighted_emel', 'EmelElkin', True, 1],
 #   ['foreign/seq_weighted_kruskal', 'foreign:seq_weighted_kruskal', 'Kruskal', True, 1],
]

GRAPH_TYPES = {
    'Bara' : [25000000],
    'ErdosDense' : [140000],
    'ErdosSparse' : [16000000],
    'Snap-LiveJournal' : [3997961],
    'Snap-road-CA' : [1957026],
    'Snap-citation' : [3764116],
    'Bad' : [1000000],
}

GRAPH_OUT_NAMES = {
    'Bara' : 'BarabasiAlbert',
    'ErdosDense' : 'DenseErdosRenyi',
    'ErdosSparse' : 'SparseErdosRenyi',
    'Snap-LiveJournal' : 'SNAP LiveJournal',
    'Snap-road-CA' : 'SNAP RoadCA',
    'Snap-citation' : 'SNAP Citation',
    'Bad' : 'BadBFS',
}

# TODO: Get correct location in output
def get_time(output_str):
    # logging.info(output_str)
    return output_str.replace('\n', ' ').split(' ')[-2]

def calculate_avg_stretch(alg, in_loc, in_graph_gbbs, in_graph_bin):
    if alg[4] == 2:
        run_cmd = f'numactl -i all ./bazel-bin/tools/test_stretch_gbbs {in_graph_gbbs} {in_loc}/{alg[2]}/'
    elif alg[4] == 3:
        run_cmd = f'numactl -i all ./bazel-bin/tools/test_stretch_nx {in_graph_bin} {in_loc}/{alg[2]}/'
    else:
        logging.info("Something went wrong")
    custom_env = os.environ.copy()
    custom_env["PARLAY_NUM_THREADS"] = str(AVAIL_PROCESSORS)
    custom_env["OMP_NUM_THREADS"] = str(AVAIL_PROCESSORS)
    p = subprocess.run(run_cmd.split(' '), capture_output=True)

    return p.stdout.decode("utf-8")

def get_avg_stretch(output_str):
    return output_str.replace('\n', ' ').split(' ')[-2]

def get_avg_stretch_group1(output_str):
    return output_str.replace('\n', ' ').split(' ')[-4]

def get_median_CI(value_arr):
    sorted_arr = np.sort(value_arr)
    z = sp.norm.ppf(1- ((1-CONFIDENCE)/2))
    n = len(sorted_arr)
    lower_index = floor((n-(z*sqrt(n)))/2)
    upper_index = ceil(1+(n+(z*sqrt(n)))/2)
    lower = sorted_arr[lower_index-1]
    upper = sorted_arr[upper_index-1]
    # logging.info(f'CI: [{lower}, {upper}]')
    return lower, upper

def get_mean_CI(value_arr):
    lower_perc = 100*((1-CONFIDENCE)/2)
    upper_perc = 100*(1-((1-CONFIDENCE)/2))
    lower = np.percentile(value_arr, lower_perc, method='median_unbiased')
    upper = np.percentile(value_arr, upper_perc, method='median_unbiased')
    # logging.info(f'CI: [{lower}, {upper}]')
    return lower, upper

def good_enough_CI(value_arr):
    if USE_MEDIAN:
        lower, upper = get_median_CI(value_arr)
        median = np.median(value_arr)
        lower_cutoff = (1-WITHIN)*median
        upper_cutoff = (1+WITHIN)*median
        # logging.info(f'median: {median}, lower_cutoff: {lower_cutoff}, upper_cutoff: {upper_cutoff}')
        if lower >= lower_cutoff and upper <= upper_cutoff:
            return True
        else:
            return False
    else:
        lower, upper = get_mean_CI(value_arr)
        mean = np.mean(value_arr)
        lower_cutoff = (1-WITHIN)*mean
        upper_cutoff = (1+WITHIN)*mean
        # logging.info(f'mean: {mean}, lower_cutoff: {lower_cutoff}, upper_cutoff: {upper_cutoff}')
        if lower >= lower_cutoff and upper <= upper_cutoff:
            return True
        else:
            return False


if __name__ == "__main__":

    output_dir = f"scripts/dphpc_benchmarking/bench_output_{time.strftime('%Y%m%d-%H%M%S')}"
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(filename=f'{output_dir}/output.log', level=logging.INFO)
    logging.info('Starting benchmark')

    p = subprocess.run('./bazel build --config=bench //tools:*'.split(' '), capture_output=True)
    if b'error' in p.stderr:
        logging.error(f'stderr: {p.stderr}')
    if b'Build completed successfully' in p.stderr:
        logging.info('Tools built success')

    for alg in PARA_ALGORITHMS:
        logging.info(f'running benchmarks of algorithm: {alg[0]}')
        build_cmd = f'./bazel build //{alg[1]}' 
        if alg[3]:
            build_cmd = f'./bazel build --config=bench //{alg[1]}' 
        # logging.info(build_cmd)
        
        alg_dir = f'./{output_dir}/{alg[2]}'
        os.makedirs(alg_dir, exist_ok=True)

        p = subprocess.run(build_cmd.split(' '), capture_output=True)
        if b'error' in p.stderr:
            logging.error(f'stderr: {p.stderr}')
        if b'Build completed successfully' in p.stderr:
            logging.info('Build success')

        for graph_type in GRAPH_TYPES.keys():
            for graph_size in GRAPH_TYPES[graph_type]:
                graph_dir = f'{PATH_TO_GRAPHS}{graph_type}/{graph_size}'
                gbbs_file = glob.glob(f'{graph_dir}/*_GBBS')
                bin_file = glob.glob(f'{graph_dir}/*_BIN')
                # Make dirs for outputs incase inexistant
                os.makedirs(f'{graph_dir}/{alg[2]}', exist_ok=True) 
                os.makedirs(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}', exist_ok=True)

                logging.info(f'Running benchmarks on Graphs: {graph_dir}')

                custom_env = os.environ.copy()

                for num_processors in EVALUATION_PROCESSORS:
                    time_arr = np.empty(0)
                    stretch_arr = np.empty(0)
                    logging.info(num_processors)
                    custom_env["PARLAY_NUM_THREADS"] = str(num_processors)
                    custom_env["OMP_NUM_THREADS"] = str(num_processors)
                    iteration = 0

                    # for iterations in range(MIN_ITERATIONS):
                    while iteration < MIN_ITERATIONS:
                        if alg[4] == 1:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]}'
                        if alg[4] == 2:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {gbbs_file[iteration % len(gbbs_file)]} {graph_dir}/{alg[2]}/'
                        if alg[4] == 3:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                        p = subprocess.run(run_cmd.split(' '), env=custom_env, capture_output=True)
                        output_str = p.stdout.decode("utf-8")
                        time_arr = np.append(time_arr, np.array(float(get_time(output_str))))
                        if alg[4] == 1:
                            stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch_group1(output_str))))
                        else:
                            stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch(calculate_avg_stretch(alg, graph_dir, gbbs_file[iteration % len(gbbs_file)], bin_file[iteration % len(bin_file)])))))
                        iteration += 1
                        
                    
                    while ((not good_enough_CI(time_arr)) or (not good_enough_CI(stretch_arr))) and iteration < MAX_ITERATIONS:
                        if alg[4] == 1:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]}'
                        if alg[4] == 2:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {gbbs_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                        if alg[4] == 3:
                            run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                        p = subprocess.run(run_cmd.split(' '), env=custom_env, capture_output=True)
                        output_str = p.stdout.decode("utf-8")
                        time_arr = np.append(time_arr, np.array(float(get_time(output_str))))
                        if alg[4] == 1:
                            stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch_group1(output_str))))
                        else:
                            stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch(calculate_avg_stretch(alg, graph_dir, gbbs_file[iteration % len(gbbs_file)], bin_file[iteration % len(bin_file)])))))
                        iteration += 1

                    logging.info(f'took {iteration} iterations')

                    df = pd.DataFrame({'Runtime':time_arr, 'AvgStretch':stretch_arr})
                    df.to_csv(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}/{num_processors}.csv')
                    plt.figure()
                    figure, axis = plt.subplots(1, 2)
                    axis[0].hist(time_arr)
                    axis[0].set_title(f'Runtime, {len(time_arr)} samples')
                    axis[1].hist(stretch_arr)
                    axis[1].set_title(f'Average Stretch, {len(stretch_arr)} samples')

                    plt.savefig(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}/{num_processors}.png')
                    plt.clf()

    for alg in SEQ_ALGORITHMS:
        logging.info(f'running benchmarks of algorithm: {alg[0]}')
        build_cmd = f'./bazel build //{alg[1]}' 
        if alg[3]:
            build_cmd = f'./bazel build --config=bench //{alg[1]}' 
        # logging.info(build_cmd)
        
        alg_dir = f'./{output_dir}/{alg[2]}'
        os.makedirs(alg_dir, exist_ok=True)

        p = subprocess.run(build_cmd.split(' '), capture_output=True)
        if b'error' in p.stderr:
            logging.error(f'stderr: {p.stderr}')
        if b'Build completed successfully' in p.stderr:
            logging.info('Build success')

        for graph_type in GRAPH_TYPES.keys():
            for graph_size in GRAPH_TYPES[graph_type]:
                graph_dir = f'{PATH_TO_GRAPHS}{graph_type}/{graph_size}'
                gbbs_file = glob.glob(f'{graph_dir}/*_GBBS')[0]
                bin_file = glob.glob(f'{graph_dir}/*_BIN')[0]
                # Make dirs for outputs incase inexistant
                os.makedirs(f'{graph_dir}/{alg[2]}', exist_ok=True) 
                os.makedirs(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}', exist_ok=True)

                logging.info(f'Running benchmarks on Graphs: {graph_dir}')

                custom_env = os.environ.copy()

                num_processors = 1
                time_arr = np.empty(0)
                stretch_arr = np.empty(0)
                logging.info(num_processors)
                custom_env["PARLAY_NUM_THREADS"] = str(num_processors)
                custom_env["OMP_NUM_THREADS"] = str(num_processors)

                iteration = 0

                while iteration < MIN_ITERATIONS:
                    if alg[4] == 1:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]}'
                    if alg[4] == 2:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {gbbs_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                    if alg[4] == 3:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                    p = subprocess.run(run_cmd.split(' '), env=custom_env, capture_output=True)
                    output_str = p.stdout.decode("utf-8")
                    time_arr = np.append(time_arr, np.array(float(get_time(output_str))))
                    if alg[4] == 1:
                        stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch_group1(output_str))))
                    else:
                        stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch(calculate_avg_stretch(alg, graph_dir, gbbs_file[iteration % len(gbbs_file)], bin_file[iteration % len(bin_file)])))))
                    iteration += 1
                    
                
                while ((not good_enough_CI(time_arr)) or (not good_enough_CI(stretch_arr))) and iteration < MAX_ITERATIONS:
                    if alg[4] == 1:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]}'
                    if alg[4] == 2:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {gbbs_file[iteration % len(gbbs_file)]} {graph_dir}/{alg[2]}/'
                    if alg[4] == 3:
                        run_cmd = f'numactl -i all ./bazel-bin/{alg[0]} -rounds 1 {bin_file[iteration % len(bin_file)]} {graph_dir}/{alg[2]}/'
                    p = subprocess.run(run_cmd.split(' '), env=custom_env, capture_output=True)
                    output_str = p.stdout.decode("utf-8")
                    time_arr = np.append(time_arr, np.array(float(get_time(output_str))))
                    if alg[4] == 1:
                        stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch_group1(output_str))))
                    else:
                        stretch_arr = np.append(stretch_arr, np.array(float(get_avg_stretch(calculate_avg_stretch(alg, graph_dir, gbbs_file[iteration % len(gbbs_file)], bin_file[iteration % len(bin_file)])))))
                    iteration += 1

                logging.info(f'took {iteration} iterations')

                df = pd.DataFrame({'Runtime':time_arr, 'AvgStretch':stretch_arr})
                df.to_csv(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}/{num_processors}.csv')
                plt.figure()
                figure, axis = plt.subplots(1, 2)
                axis[0].hist(time_arr)
                axis[0].set_title(f'Runtime, {len(time_arr)} samples')
                axis[1].hist(stretch_arr)
                axis[1].set_title(f'Average Stretch, {len(stretch_arr)} samples')

                plt.savefig(f'{alg_dir}/{GRAPH_OUT_NAMES[graph_type]}/{graph_size}/{num_processors}.png')
                plt.clf()
                    

               
