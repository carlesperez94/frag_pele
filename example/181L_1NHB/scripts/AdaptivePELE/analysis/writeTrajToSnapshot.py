import os
import sys
import argparse
import glob
from AdaptivePELE.utilities import utilities


def parseArguments():
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("clusteringObject", type=str, help="Path to the clustering object")
    parser.add_argument("trajectory", type=int, help="Trajectory number")
    parser.add_argument("snapshot", type=int, help="Snapshot to select (in accepted steps)")
    parser.add_argument("epoch", type=str, help="Path to the epoch to search the snapshot")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    args = parser.parse_args()
    return args.clusteringObject, args.trajectory, args.snapshot, args.epoch, args.o


if __name__ == "__main__":
    clusteringObject, trajectory, snapshot, epoch, outputPath = parseArguments()
    if outputPath is not None:
        outputPath = os.path.join(outputPath, "")
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    else:
        outputPath = ""
    sys.stderr.write("Reading clustering object...\n")
    cl = utilities.readClusteringObject(clusteringObject)
    if cl.conformationNetwork is None:
        sys.exit("Clustering object loaded has no conformation network!!")
    conf = cl.conformationNetwork
    filename = glob.glob(epoch+"/*traj*_%d.pdb" % trajectory)
    snapshots = utilities.getSnapshots(filename[0])
    snapshots = snapshots[:snapshot+1]
    procMapping  = open(os.path.join(epoch, "processorMapping.txt")).read().rstrip().split(',')
    leaf = procMapping[trajectory-1]
    pathway = conf.createPathwayToCluster(int(leaf))
    cl.writePathwayTrajectory(pathway, outputPath+"pathway.pdb")
    with open(outputPath+"pathway.pdb", "a") as f:
        f.write("ENDMDL\n".join(snapshots))
