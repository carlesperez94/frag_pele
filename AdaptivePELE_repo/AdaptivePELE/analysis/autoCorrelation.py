from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import glob
import argparse
import numpy as np
from AdaptivePELE.freeEnergies import cluster
from AdaptivePELE.freeEnergies import utils
from AdaptivePELE.utilities import utilities
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def parse_arguments():
    """
        Create command-line interface
    """
    desc = "Calculate the autocorrelation function of a MSM discretization"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--lagtimes", type=int, help="Upper limit of the lagtimes to analyse")
    parser.add_argument("-dl", "--dlagtimes", type=int, default=20, help="Separation of the lagtimes to analyse")
    parser.add_argument("-n", "--n_clusters", type=int, help="Number of clusters")
    parser.add_argument("-o", default=None, help="Path of the folder where to store the plots")
    parser.add_argument("--savePlots", action="store_true", help="Save the plots to disk")
    parser.add_argument("--showPlots", action="store_true", help="Show the plots to screen")
    parser.add_argument("--dtrajs", type=str, help="Path to the folder with the discretized trajectories")
    parser.add_argument("--clusters", type=str, default=None, help="Path to the clustering file")
    parser.add_argument("--trajs", type=str, default=None, help="Path to the trajectories files")
    args = parser.parse_args()
    return args.clusters, args.lagtimes, args.o, args.savePlots, args.showPlots, args.dtrajs, args.trajs, args.n_clusters, args.dlagtimes


def __rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def __rmFiles(trajWildcard):
    allfiles = glob.glob(trajWildcard)
    for f in allfiles:
        __rm(f)


def __cleanupFiles(trajWildcard, cleanupClusterCenters=True):
    __rmFiles("clustering_object.pkl")
    # __rmFiles("MSM_object.pkl")
    # __rmFiles("discretized/traj_*")
    __rmFiles(trajWildcard)
    if cleanupClusterCenters:
        __rmFiles("discretized/clusterCenter*")


def calculateAutoCorrelation(lagtimes, dtrajs, nclusters, nLags):
    C = np.zeros((nclusters, nLags))
    Ci = np.zeros((nclusters, nLags))
    Cf = np.zeros((nclusters, nLags))
    autoCorr = np.zeros((nclusters, nLags))
    N = 0
    M = np.zeros(nLags)
    for trajectory in dtrajs:
        traj = np.loadtxt(trajectory, dtype=int)
        Nt = traj.size
        N += Nt
        for il, lagtime in enumerate(lagtimes):
            M[il] += Nt-lagtime
            for i in range(Nt-lagtime):
                autoCorr[traj[i], il] += (traj[i] == traj[i+lagtime])
                C[traj[i], il] += 1
                Ci[traj[i], il] += 1
                if i > lagtime:
                    Cf[traj[i], il] += 1
            for j in range(Nt-lagtime, Nt):
                C[traj[j], il] += 1
                Cf[traj[j], il] += 1

    mean = C/float(N)
    var = (N*C-(C**2))/float(N*(N-1))
    autoCorr += M*mean**2-(Ci+Cf)*mean
    autoCorr /= N
    autoCorr /= var
    return autoCorr


def create_plots(autoCorr, plots_path, save_plot, show_plot, nclusters, lagtimes, threshold=2, title=""):
    if threshold < 1:
        fig_filename = "autoCorr_thres_%s.png" % str(threshold).replace(".", "_")
        filtered = np.where(autoCorr[:, -1] > threshold)[0]
        if len(filtered) == 0:
            raise ValueError("The threshold specified is too strict, no states found above it")
    else:
        fig_filename = "autoCorr_no_thres.png"
        filtered = list(range(nclusters))
    plt.figure()
    axes = plt.plot(lagtimes, autoCorr.T[:, filtered])
    plt.xlabel("Lagtime")
    plt.title("Autocorrelation of membership function %s" % title)
    if len(filtered) < 20:
        for i, ax in zip(filtered, axes):
            ax.set_label("Cluster %d" % i)
        plt.legend()
    if save_plot:
        plt.savefig(os.path.join(plots_path, fig_filename))
    if show_plot:
        plt.show()


def main(lagtime, clusters_file, disctraj, trajs, n_clusters, plots_path, save_plot, show_plot, lagtime_resolution=20):
    lagtimes = list(range(1, lagtime, lagtime_resolution))
    n_lags = len(lagtimes)
    if disctraj is None:
        clusteringObject = cluster.Cluster(n_clusters, trajs, "traj*", alwaysCluster=False)
        if clusters_file is not None:
            # only assign
            utilities.makeFolder(clusteringObject.discretizedFolder)
            clusteringObject.clusterCentersFile = clusters_file
        clusteringObject.clusterTrajectories()
        disctraj = clusteringObject.discretizedFolder
        clusterCenters = clusteringObject.clusterCenters
    else:
        clusterCenters = np.loadtxt(clusters_file)
    if len(clusterCenters) != n_clusters:
        raise ValueError("Number of clusters specified in the -n parameter does not match the provided clusters")
    print("Calculating autocorrelation...")
    dtrajs = glob.glob(os.path.join(disctraj, "traj*"))
    dtrajs_loaded = [np.loadtxt(dtraj, dtype=int) for dtraj in dtrajs]

    autoCorr = utils.calculateAutoCorrelation(lagtimes, dtrajs_loaded, n_clusters, n_lags)
    np.save("autoCorr.npy", autoCorr)
    # __cleanupFiles(parameters.trajWildcard, False)

    utilities.write_PDB_clusters(np.vstack((clusterCenters.T, np.abs(autoCorr[:, -1]))).T, use_beta=True, title="cluster_autoCorr.pdb")
    print("Clusters over correlation time limit")
    correlation_limit = np.exp(-1)
    states2 = np.where(autoCorr[:, -1] > correlation_limit)[0]
    size2 = states2.size
    if len(states2):
        print(" ".join(map(str, states2)))
    print("Number of clusters:", size2, ", %.2f%% of the total" % (100*size2 / float(n_clusters)))
    print("Clusters with more than 0.1 autocorrelation")
    states1 = np.where(autoCorr[:, -1] > 0.1)[0]
    size1 = states1.size
    if len(states1):
        print(" ".join(map(str, states1)))
    print("Number of clusters:", size1, ", %.2f%% of the total" % (100*size1 / float(n_clusters)))
    if size2 > 0:
        print("Correlation time not achieved at lagtime %d" % lagtime)
    else:
        for i in range(len(lagtimes)):
            states = np.where(autoCorr[:, -i-1] > correlation_limit)[0]
            if len(states):
                string_states = ", ".join(map(str, states))
                print("Correlation time %d, for states: %s" % (lagtimes[-i], string_states))
                break

    if plots_path is None:
        plots_path = ""
    else:
        utilities.makeFolder(plots_path)
    create_plots(autoCorr, plots_path, save_plot, show_plot, n_clusters, lagtimes, threshold=2.0)

if __name__ == "__main__":
    clusters, lagtime_upper, output_plots_path, save_plots, show_plots, dtraj_path, trajs_path, num_clusters, lagtimes_res = parse_arguments()
    main(lagtime_upper, clusters, dtraj_path, trajs_path, num_clusters, output_plots_path, save_plots, show_plots, lagtime_resolution=lagtimes_res)
