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

TITLE = "Other graphs"
Y_LABEL = "$What vs Graph type(size)"
Y_LIMIT = False
LOWER_LIMIT_Y = -0.2 #Keep
UPPER_LIMIT_Y = 6.5
BASE_DIRECTORY = "scripts/dphpc_benchmarking/"
DIRECTORY = "benchmarks/strong_scaling"
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

seq_algos = {
    # 'UnweightedAlonKarp' : '#ca6ab2',
}

algos = dict(**par_algos, **seq_algos)

out_names = {
    'RandBFS' : 'RandBFS',
    'SSSP_StarDecomp': 'StarDecomp',
    'LSPES' : 'LSPES',
    'StarDecomp' : 'AKPW',
    'ParAbrahamBartal' : 'ParAbrahamBartal',
    'UnweightedAlonKarp' : 'SequentialAlonKarp'
}

graphs = {
    'BadBFS' : ['1000000', -0.05, 2.2],
    'SNAP Citation' : ['3764116', -0.02, 0.9],
    'SNAP LiveJournal' : ['3997961', -0.05, 2.1],
    'SNAP RoadCA' : ['1957026', -0.045, 1.5],
}


# GRAPH_TYPE = 'BarabasiAlbert'
# RATIO = 12000000
# GRAPH_TYPE = 'DenseErdosRenyi'
# RATIO = 100000
# GRAPH_TYPE = 'SparseErdosRenyi'
# RATIO = 10000000



output_folder = f"{BASE_DIRECTORY}plot_output"
SAVE_AS_EPS = True
PLOT_STRETCH = True
# TIME_UNIT = 'hours'
# TIME_UNIT = 'minutes'
TIME_UNIT = 'seconds'
# TIME_UNIT = 'miliseconds'
# TIME_UNIT = 'microseconds'

if PLOT_STRETCH:
    Y_LABEL = Y_LABEL.replace('$What', 'Average Stretch')
    graphs = {
        'BadBFS' : ['1000000', -0.05, 500],
        'SNAP Citation' : ['3764116', -0.02, 15],
        'SNAP LiveJournal' : ['3997961', -0.05, 8],
        'SNAP RoadCA' : ['1957026', -0.045, 12.5],
    }
else:
    Y_LABEL = Y_LABEL.replace('$What', f'Runtime [{TIME_UNIT}]')




DIMENSIONS = [i for i in range(len(algos))]
MACHINE_SPEC = "32-Core AMD EPYC 7R13 3.6 GHz"
FONT_SIZE = 12
MARKERSIZE = 4
FONT_WEIGHT = 'normal'
# NAME_OVERWRITE = f"{GRAPH_TYPE}_weak_scaling_Star"
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

title = "$\\bf{" + TITLE.replace(' ', "\ ") + "}$ on " + MACHINE_SPEC
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

    for l, key in zip(range(Nlines), algos):
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
        if key == 'UnweightedAlonKarp':
            an = axis.annotate(out_names[key], xy=(x, y), xycoords='data', ha="center", va="center", color=seq_algos[key], weight=FONT_WEIGHT)
        else:
            an = axis.annotate(out_names[key], xy=(x, y), xycoords='data', ha="center", va="center", color=par_algos[key], weight=FONT_WEIGHT)
        an.draggable(True)


def create_plot():
    sns.set_theme()
    sns.set(rc={'figure.figsize':(8, 5)})
    sns.set(style="darkgrid")  

    # fig, p = plt.subplots(ncols=len(graphs))#, layout='constrained')
    print(len(graphs))
    fig = plt.figure()
    spec = fig.add_gridspec(ncols=len(graphs), left=0.075, right=0.95, wspace=0.5, bottom=0.125)
    p = []
    for i in range(len(graphs)):
        p.append(fig.add_subplot(spec[i]))

    colors = []

    position = 0
    xlabels = []
    for graph in graphs.keys():
        
        for alg in par_algos:

            df = read_performance_data(alg, f'{graph}/{graphs[graph][0]}', 32)
            if PLOT_STRETCH:
                df = df.drop('Runtime', axis=1)
            else:
                df = df.drop('AvgStretch', axis=1)
            arr = df.to_numpy().flatten()
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
            lower, upper = get_median_CI(arr)
            boxes = p[position].boxplot([arr], notch=True, positions=[0], manage_ticks=False, conf_intervals=[[lower, upper]], sym='')
            # print(len(boxes))
            colors.append(boxes['boxes'][0])

            for patch in boxes['boxes']:
                patch.set_color(par_algos[alg])
            for whisker in boxes['whiskers']:
                whisker.set_color(par_algos[alg])
            for cap in boxes['caps']:
                cap.set_color(par_algos[alg])

        for alg in seq_algos:

            df = read_performance_data(alg, f'{graph}/{graphs[graph][0]}', 1)
            if PLOT_STRETCH:
                df = df.drop('Runtime', axis=1)
            else:
                df = df.drop('AvgStretch', axis=1)
            arr = df.to_numpy().flatten()
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
            lower, upper = get_median_CI(arr)
            boxes = p[position].boxplot([arr], notch=True, positions=[0], manage_ticks=False, conf_intervals=[[lower, upper]], sym='')
            # print(len(boxes))
            colors.append(boxes['boxes'][0])

            for patch in boxes['boxes']:
                patch.set_color(seq_algos[alg])
            for whisker in boxes['whiskers']:
                whisker.set_color(seq_algos[alg])
            for cap in boxes['caps']:
                cap.set_color(seq_algos[alg])

        p[position].set_xticks([0])
        xlabels = ['' for x in p[position].get_xticks()]
        p[position].set_xticklabels(xlabels)
        p[position].set_xlabel(f'{graph}\n({graphs[graph][0]})')
        p[position].set_ylabel('')
        p[position].yaxis.grid(True)
        p[position].xaxis.grid(False)
        p[position].margins(0.15)
        p[position].set_ylim(graphs[graph][1], graphs[graph][2])
        position += 1


    # If len(graph)==4:
    p[0].text(-0.35,1.03, Y_LABEL, ha='left', transform=p[0].transAxes, fontsize=FONT_SIZE)
    plt.title(title, x=-4.85, y=1.07, ha='left', fontsize=FONT_SIZE)
    # If len(graph)==3:
    # p[0].text(-0.19,1.03, Y_LABEL, ha='left', transform=p[0].transAxes, fontsize=FONT_SIZE)
    # plt.title(title, x=-3.19, y=1.07, ha='left', fontsize=FONT_SIZE)

    my_legend(p[-1], colors)

    plot_name = [out_names[x] for x in algos.keys()]
    # plt.tight_layout()
    # fig.tight_layout(pad=0.25)
    plt.show(block=False)
    plot_name = "-".join(sorted(plot_name))
    if len(NAME_OVERWRITE) == 0:
        output_file = f'{output_folder}/others{"-stretch" if PLOT_STRETCH else ""}-{plot_name}.{"eps" if SAVE_AS_EPS else "png"}'
    else:
        output_file = f'{output_folder}/{NAME_OVERWRITE}.{"eps" if SAVE_AS_EPS else "png"}'
    input(f"Input any key to save image to {output_file}")
    if len(plot_name) > 200:
        plot_name = plot_name[:200] + "-rest"
    if not SAVE_AS_EPS:
        plt.savefig(output_file, dpi=1000)
    else:
        plt.savefig(output_file, format='eps')
    print("Saved the figure.")
    plt.clf()

if __name__ == "__main__":
    create_plot()
