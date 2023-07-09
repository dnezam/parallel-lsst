from itertools import cycle
import os
import subprocess
# from turtle import color
import numpy as np
from scipy import ndimage
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import signal
import sys
import glob
import pandas as pd
from math import ceil, floor, sqrt
import scipy.stats as sp

###################################################################################
#       If the GUI is very laggy and you get the following messages:
#   QSocketNotifier: Can only be used with threads started with QThread
#   qt.qpa.wayland: Wayland does not support QWindow::requestActivate()
#
#       Try running 'export XDG_SESSION_TYPE=x11' before running this script
####################################################################################

TITLE = '$What of size $Size x processors'
Y_LABEL = "$What vs #processors"
Y_LIMIT = False
LOWER_LIMIT_Y = -0.2 #Keep
UPPER_LIMIT_Y = 15
BASE_DIRECTORY = "scripts/dphpc_benchmarking/"
DIRECTORY = "benchmarks/weak_scaling"
USE_BOX = True

CONFIDENCE = 0.95
WITHIN = 0.05


par_algos = {
    'RandBFS' : '#b30000',
    # 'SSSP_StarDecomp' : "#7c1158",
    'LSPES' : "#4421af",
    'StarDecomp' : "#1a94ff",
    # 'ParAbrahamBartal' : "#00b7c7",
}

out_names = {
    'RandBFS' : 'RandBFS',
    'SSSP_StarDecomp': 'StarDecomp',
    'LSPES' : 'LSPES',
    'StarDecomp' : 'AKPW',
    'ParAbrahamBartal' : 'ParAbrahamBartal',
    'UnweightedAlonKarp' : 'SequentialAlonKarp'
}


# GRAPH_TYPE = 'BarabasiAlbert'
# RATIO = 12000000
# GRAPH_TYPE = 'DenseErdosRenyi'
# RATIO = 100000
GRAPH_TYPE = 'SparseErdosRenyi'
RATIO = 10000000



output_folder = f"{BASE_DIRECTORY}plot_output"
SAVE_AS_EPS = True
PLOT_STRETCH = False
# TIME_UNIT = 'hours'
# TIME_UNIT = 'minutes'
# TIME_UNIT = 'seconds'
TIME_UNIT = 'miliseconds'
# TIME_UNIT = 'microseconds'

if PLOT_STRETCH:
    Y_LABEL = Y_LABEL.replace('$What', 'Average Stretch')
else:
    Y_LABEL = Y_LABEL.replace('$What', f'Runtime [{TIME_UNIT}]')

TITLE = TITLE.replace('$Size', str(RATIO))
TITLE = TITLE.replace('$What', GRAPH_TYPE)


DIMENSIONS = [1, 2, 3, 4, 5, 6, 7, 8, 16, 32]#, 64, 128]
SSSP_DIMENSIONS = [1, 2, 3, 4, 5, 6, 7, 8]
MACHINE_SPEC = "32-Core AMD EPYC 7R13 3.6 GHz"
FONT_SIZE = 12
MARKERSIZE = 4
FONT_WEIGHT = 'normal'
NAME_OVERWRITE = f"{GRAPH_TYPE}_weak_scaling_Star"
NAME_OVERWRITE = ''

# The following was tested
if SAVE_AS_EPS:
    FONT_SIZE = 16
    MARKERSIZE = 5
    FONT_WEIGHT = 'bold'
    MACHINE_SPEC = "32-Core AMD EPYC 7R13 3.6 GHz"

def signal_handler(sig, frame):
    print('You pressed Ctrl+C, image not saved!')
    sys.exit(0)

def map_to_hours(x):
    return map_to_minutes(x)/60.0

def map_to_minutes(x):
    return map_to_seconds(x)/60.0

def map_to_seconds(x):
    return x/1000000.0

def map_to_miliseconds(x):
    return x/1000.0
signal.signal(signal.SIGINT, signal_handler)

def read_performance_data(algo, graph, num_processors):
    file_name = f"{BASE_DIRECTORY}{DIRECTORY}/{algo}/{graph}/{num_processors}.csv"
    df = pd.read_csv(file_name, index_col=0)
    # print(df)
    return df

def get_median_CI(value_arr):
    sorted_arr = np.sort(value_arr)
    z = sp.norm.ppf(1- ((1-CONFIDENCE)/2))
    n = len(sorted_arr)
    lower_index = floor((n-(z*sqrt(n)))/2)
    upper_index = ceil(1+(n+(z*sqrt(n)))/2)
    lower = sorted_arr[lower_index-1]
    upper = sorted_arr[upper_index-1]
    # print(f'CI: [{lower}, {upper}]')
    return lower, upper

####-----------------------------------------------------------------------------------------#####

title = "$\\bf{" + TITLE.replace(' ', "\ ") + "}$"# on " + MACHINE_SPEC
plt.rcParams["font.family"] = 'sans serif' # TODO fix better font

def my_legend(axis, lines):
    # Initial setting
    N = 256
    Nlines = len(lines)
    xmin, xmax = axis.get_xlim()
    ymin, ymax = axis.get_ylim()
    # the 'point of presence' matrix
    pop = np.zeros((Nlines, N, N), dtype=float)
    for l in range(Nlines):
        # get xy data and scale it to the NxN squares
        xy = lines[l].get_xydata()
        xy = (xy - [xmin,ymin]) / ([xmax-xmin, ymax-ymin]) * N
        xy = xy.astype(np.int32)
        # mask stuff outside plot        
        mask = (xy[:,0] >= 0) & (xy[:,0] < N) & (xy[:,1] >= 0) & (xy[:,1] < N)
        xy = xy[mask]
        # add to pop
        for p in xy:
            pop[l][tuple(p)] = 1.0
    # find whitespace, nice place for labels
    ws = 1.0 - (np.sum(pop, axis=0) > 0) * 1.0 
    # don't use the borders, modified for better placement with our data
    ws[:,:8]   = 0  # lower
    ws[:,N-8:] = 0  # upper
    ws[:48,:]   = 0  # left
    ws[N-48:,:] = 0   # right

    for l, key in zip(range(Nlines), par_algos):
        # positive weights for current line, negative weight for others....
        w = -0.3 * np.ones(Nlines, dtype=float)
        w[l] = 0.5
        # calculate a field         
        # blur the pop's
        for z in range(Nlines):
            pop[z] = ndimage.gaussian_filter(pop[z], sigma=N/8)
        p = ws + np.sum(w[:, np.newaxis, np.newaxis] * pop, axis=0)
        if False:
            plt.figure()
            plt.imshow(p, interpolation='nearest')
            plt.title(axis.lines[l][0].get_label())
            plt.show()

        pos = np.argmax(p)  # note, argmax flattens the array first 
        best_x, best_y =  (pos / N, pos % N) 
        x = xmin + (xmax-xmin) * best_x / N       
        y = ymin + (ymax-ymin) * best_y / N       
        # set already used label position
        best_x = int(best_x)
        ws[best_x-40:best_x+40, best_y-40:best_y+40] = 0
        # and add label
        an = axis.annotate(out_names[key], xy=(x, y), xycoords='data', ha="center", va="center", color=par_algos[key], weight=FONT_WEIGHT)
        an.draggable(True)


def create_plot():
    sns.set_theme()
    sns.set(rc={'figure.figsize':(8, 5)})
    sns.set(style="darkgrid")  

    fig, p = plt.subplots()

    colors = []

    for alg in par_algos:
        data = []
        CIS = []
        if alg == 'SSSP_StarDecomp':
            positions = [i for i in range(len(SSSP_DIMENSIONS))]
            for num_processors in SSSP_DIMENSIONS:
                df = read_performance_data(alg, f'{GRAPH_TYPE}', num_processors)
                if PLOT_STRETCH:
                    df = df.drop('Runtime', axis=1)
                else:
                    df = df.drop('AvgStretch', axis=1)
                # print(df)
                arr = df.to_numpy().flatten()
                # print(arr.shape)
                if not PLOT_STRETCH:
                    if TIME_UNIT == 'hours':
                        arr = np.vectorize(map_to_hours)(arr)
                    if TIME_UNIT == 'minutes':
                        arr = np.vectorize(map_to_minutes)(arr)
                    if TIME_UNIT == 'seconds':
                        arr = np.vectorize(map_to_seconds)(arr)
                    if TIME_UNIT == 'miliseconds':
                        arr = np.vectorize(map_to_miliseconds)(arr)
                # print(arr)
                data.append(arr)
                lower, upper = get_median_CI(arr)
                CIS.append([lower, upper])
            boxes = p.boxplot(data, notch=True, positions=positions, manage_ticks=False, conf_intervals=CIS, sym='')
            # print(len(boxes))
            colors.append(boxes['boxes'][0])

            for patch in boxes['boxes']:
                patch.set_color(par_algos[alg])
            for whisker in boxes['whiskers']:
                whisker.set_color(par_algos[alg])
            for cap in boxes['caps']:
                cap.set_color(par_algos[alg])
            continue

        positions = [i for i in range(len(DIMENSIONS))]
        for num_processors in DIMENSIONS:
            df = read_performance_data(alg, f'{GRAPH_TYPE}', num_processors)
            if PLOT_STRETCH:
                df = df.drop('Runtime', axis=1)
            else:
                df = df.drop('AvgStretch', axis=1)
            # print(df)
            arr = df.to_numpy().flatten()
            # print(arr.shape)
            if not PLOT_STRETCH:
                if TIME_UNIT == 'hours':
                    arr = np.vectorize(map_to_hours)(arr)
                if TIME_UNIT == 'minutes':
                    arr = np.vectorize(map_to_minutes)(arr)
                if TIME_UNIT == 'seconds':
                    arr = np.vectorize(map_to_seconds)(arr)
                if TIME_UNIT == 'miliseconds':
                    arr = np.vectorize(map_to_miliseconds)(arr)
            # print(arr)
            data.append(arr)
            lower, upper = get_median_CI(arr)
            CIS.append([lower, upper])


        boxes = p.boxplot(data, notch=True, positions=positions, manage_ticks=False, conf_intervals=CIS, sym='')
        # print(len(boxes))
        colors.append(boxes['boxes'][0])

        for patch in boxes['boxes']:
            patch.set_color(par_algos[alg])
        for whisker in boxes['whiskers']:
            whisker.set_color(par_algos[alg])
        for cap in boxes['caps']:
            cap.set_color(par_algos[alg])



    if Y_LIMIT:
        plt.ylim(LOWER_LIMIT_Y, UPPER_LIMIT_Y)
    # create xticks manually formatted
    positions = [i for i in range(len(DIMENSIONS))]
    plt.xticks(positions)
    xlabels = [f'{DIMENSIONS[x] if (x % 1 == 0) else ""}' for x in p.get_xticks()] # TODO: Fix names
    p.set_xticklabels(xlabels)
    # style labels & grid
    plt.ylabel('')
    p.yaxis.grid(True)
    bbox = p.get_yticklabels()[-1].get_window_extent()
    axis_aligned_coor,_ = p.transAxes.inverted().transform([bbox.x0, bbox.y0])
    axis_aligned_coor = axis_aligned_coor - 0.03 # Align similarly to example on the slides
    p.text(axis_aligned_coor,1.03, Y_LABEL, ha='left', transform=p.transAxes, fontsize=FONT_SIZE)
    plt.xlabel("", fontsize=FONT_SIZE)
    p.xaxis.grid(False)
    plt.title(title, x=axis_aligned_coor, y=1.07, ha='left', fontsize=FONT_SIZE)

    my_legend(p, colors)

    plot_name = [out_names[x] for x in par_algos.keys()]
    plt.tight_layout()
    plt.show(block=False)
    plot_name = "-".join(sorted(plot_name))
    if len(NAME_OVERWRITE) == 0:
        output_file = f'{output_folder}/weak-scaling-{GRAPH_TYPE}{"-stretch" if PLOT_STRETCH else ""}-{plot_name}.{"eps" if SAVE_AS_EPS else "png"}'
    else:
        output_file = f'{output_folder}/{NAME_OVERWRITE}.{"eps" if SAVE_AS_EPS else "png"}'
    input(f"Input any key to save image to {output_file}")
    if len(plot_name) > 200:
        plot_name = plot_name[:200] + "-rest"
    if not SAVE_AS_EPS:
        plt.savefig(output_file, dpi=500)
    else:
        plt.savefig(output_file, format='eps')
    print("Saved the figure.")
    plt.clf()

if __name__ == "__main__":
    create_plot()
