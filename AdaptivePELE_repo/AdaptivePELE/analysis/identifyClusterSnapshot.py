from __future__ import print_function, unicode_literals
import os
import glob
import argparse
import numpy as np
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import atomset, RMSDCalculator


def parseArguments():
    """
        Parse the command-line options

        :returns: :py:class:`.Clustering`, int, int, int, str -- Clustering
            object, number of trajectory, number of snapshot, number of epoch,
            output path where to write the files
    """
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("epoch", type=int, help="Path to the epoch to search the snapshot")
    parser.add_argument("trajectory", type=int, help="Trajectory number")
    parser.add_argument("snapshot", type=int, help="Snapshot to select (in accepted steps)")
    parser.add_argument("resname", type=str, help="Name of the ligand in the pdb")
    parser.add_argument("clustering", type=str, help="Path to the clustering object")
    parser.add_argument("--top", type=str, default=None, help="Name of the pdb topology for loading non-pdb trajectories")
    args = parser.parse_args()
    return args.epoch, args.trajectory, args.snapshot, args.resname, args.clustering, args.top


def main(epoch_num, trajectory, snapshot_num, resname, clustering_object, topology):
    calc = RMSDCalculator.RMSDCalculator()
    clustering_object = utilities.readClusteringObject(clustering_object)
    n_clusters = np.loadtxt(os.path.join(str(epoch_num-1), "clustering", "summary.txt")).shape[0]
    if topology is not None:
        topology_contents = utilities.getTopologyFile(topology)
    else:
        topology_contents = None
    filename = glob.glob(os.path.join(str(epoch_num), "*traj*_%d.*" % trajectory))
    if not filename:
        raise ValueError("No file with the specified epoch and trajectory found")
    try:
        snapshots = utilities.getSnapshots(filename[0], topology=topology)[snapshot_num]
    except IndexError:
        raise IndexError("Snapshot number %d not found in trajectory %d for epoch %d, please check that the arguments provided are correct" % (snapshot_num, trajectory, epoch_num))
    pdb = atomset.PDB()
    pdb.initialise(snapshots, resname=resname)
    for i, cluster in enumerate(clustering_object[:n_clusters]):
        dist = calc.computeRMSD(pdb, cluster.pdb)
        if dist < cluster.threshold:
            print("Snapshot belongs to cluster", i)
            return
    print("Snapshot not assigned to any cluster! :(")

if __name__ == "__main__":
    epoch, traj_num, snapshot, ligand_resname, clustering, top = parseArguments()
    main(epoch, traj_num, snapshot, ligand_resname, clustering, top)
