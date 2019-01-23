from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import socket
import matplotlib
import numpy as np
import os
import argparse
machine = socket.gethostname()
if machine == "bsccv03":
    matplotlib.use('wxagg')
elif 'login' in machine:
    matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
try:
    # This might fail for older versions of matplotlib (e.g in life cluster)
    plt.style.use("ggplot")
except NameError:
    pass


def printHelp():
    """
        Create command line interface

        :returns: str -- Output filename ( if specified )
    """
    desc = "Program that prints the number of clusters of each threshold that"\
           "have been selected by the spwaning throughout an adaptive sampling simulation. "\
           "It must be run in the root folder. "
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-filename", type=str, default="", help="Output filename")
    parser.add_argument("-o", "--output", type=str, default="", help="Output folder")
    args = parser.parse_args()
    return args.filename, args.output


def main(filename, output_folder):
    print("FILENAME", filename)
    templateSummary = "%d/clustering/summary.txt"
    allFolders = os.listdir(".")
    numberOfEpochs = len([epoch for epoch in allFolders if epoch.isdigit() and
                          os.path.isfile(templateSummary % int(epoch))])

    clustering = np.loadtxt(templateSummary % (numberOfEpochs-1))
    clustersThres = list(set(clustering[:, 4]))
    clustersThres.sort(reverse=True)
    spawningPerThres = np.zeros((numberOfEpochs, len(clustersThres)))
    for i in range(numberOfEpochs):
        clustering = np.loadtxt(templateSummary % i)
        for j, threshold in enumerate(clustersThres):
            spawningPerThres[i, j] = clustering[clustering[:, 4] == threshold, 2].sum()
    line_objects = plt.plot(spawningPerThres)
    plt.legend(line_objects, tuple(["Cluster size %d" % x for x in clustersThres]),
               loc="best")
    plt.title("Processors spawned per epoch and cluster size")
    plt.xlabel("Epoch")
    plt.ylabel("Number of spawned processors")
    if output_folder and not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if filename != "":
        plt.savefig(os.path.join(output_folder, "%s_spawning.png" % filename))
    plt.show()

if __name__ == "__main__":
    filename, out_folder = printHelp()
    main(filename, out_folder)
