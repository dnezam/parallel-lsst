import os
import subprocess
import pandas as pd
graph_gen_algo = 'ErdosRenyi_Bench'
file_prefix = graph_gen_algo + '_'
base = 'inputs/'
msg_spanner = 'Is spanner: '
msg_timmer = '# time per iter: '
msg_stretch = '# avg stretch per iter: '

low_avg_stretch_bin = [
    'bazel-bin/benchmarks/LowAvgStretch/Algo_spanner/RandomSpanner',
    # 'bazel-bin/benchmarks/LowAvgStretch/Algo_SSSP_ParallelStarDecomposition/SSSP_ParallelStarDecomposition'
]

graph_algo_base = base + graph_gen_algo + '/'
statistics = pd.DataFrame(columns=['Algo', 'n', 'm', 'Time', 'Stretch'])
for graph_file in os.listdir(graph_algo_base):
    print(graph_file)
    full_path = graph_algo_base + graph_file
    info = graph_file.strip(file_prefix).split('_')
    n = int(info[0])
    m = int(info[1])
    for bin in low_avg_stretch_bin:
        out = subprocess.check_output([bin, '-rounds', '3', '-s', full_path])
        lines = out.split(b'\n')
        time = []
        stretch = []
        check_spanner = True
        for line in lines:
            line = line.decode()
            if line.find(msg_spanner) == 0:
                line = line.strip(msg_spanner)
                if line != '1':
                    breakpoint()
                    print('Not a spanner!')
                    check_spanner = False
                    break
        if check_spanner:
            for line in lines:
                line = line.decode()
                if line.find(msg_timmer) == 0:
                    time.append(float(line.strip(msg_timmer)))
                if line.find(msg_stretch) == 0:
                    stretch.append(float(line.strip(msg_stretch)))
        statistics.loc[len(statistics.index)] = [bin, n, m, sum(time)/len(time), sum(stretch)/len(stretch)]

statistics.to_csv('analysis/' + graph_gen_algo+'.csv')
    