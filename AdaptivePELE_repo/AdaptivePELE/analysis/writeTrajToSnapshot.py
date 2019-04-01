import os
import sys
import argparse
import glob
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import atomset


def parseArguments():
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("clusteringObject", type=str, help="Path to the clustering object")
    parser.add_argument("epoch", type=str, help="Path to the epoch to search the snapshot")
    parser.add_argument("trajectory", type=int, help="Trajectory number")
    parser.add_argument("snapshot", type=int, help="Snapshot to select (in accepted steps)")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    parser.add_argument("-top", type=str, default=None, help="Topology file for non-pdb trajectories")
    args = parser.parse_args()
    return args.clusteringObject, args.trajectory, args.snapshot, args.epoch, args.o, args.top


if __name__ == "__main__":
    clusteringObject, trajectory, snapshot, epoch, outputPath, topology = parseArguments()
    if outputPath is not None:
        outputPath = os.path.join(outputPath, "")
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    else:
        outputPath = ""

    if topology is not None:
        topology_contents = utilities.getTopologyFile(topology)
    else:
        topology_contents = None

    sys.stderr.write("Reading clustering object...\n")
    cl = utilities.readClusteringObject(clusteringObject)
    if cl.conformationNetwork is None:
        sys.exit("Clustering object loaded has no conformation network!!")
    conf = cl.conformationNetwork
    filename = glob.glob(epoch+"/*traj*_%d*" % trajectory)
    snapshots = utilities.getSnapshots(filename[0], topology)
    snapshots = snapshots[:snapshot+1]
    if not isinstance(snapshots[0], basestring):
        new_snapshots = []
        for snapshot in snapshots:
            PDB = atomset.PDB()
            PDB.initialise(snapshot, topology=topology_contents)
            new_snapshots.append(PDB.get_pdb_string())
        snapshots = new_snapshots

    procMapping = open(os.path.join(epoch, "processorMapping.txt")).read().rstrip().split(',')
    leaf = procMapping[trajectory-1]
    pathway = conf.createPathwayToCluster(int(leaf))
    cl.writePathwayTrajectory(pathway, outputPath+"pathway.pdb")
    with open(outputPath+"pathway.pdb", "a") as f:
        f.write("ENDMDL\n".join(snapshots))
