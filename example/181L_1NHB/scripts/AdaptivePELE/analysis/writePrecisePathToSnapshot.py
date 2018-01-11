import os
import sys
import argparse
import glob
import itertools
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
    pathway = []
    # Strip out trailing backslash if present
    pathPrefix, epoch = os.path.split(epoch.rstrip("/"))
    sys.stderr.write("Creating pathway...\n")
    while epoch != "0":
        filename = glob.glob(os.path.join(pathPrefix,epoch,"*traj*_%d.pdb" % trajectory))
        snapshots = utilities.getSnapshots(filename[0])
        snapshots = snapshots[:snapshot+1]
        pathway.insert(0, snapshots)
        procMapping  = open(os.path.join(pathPrefix, epoch, "processorMapping.txt")).read().rstrip().split(':')
        epoch, trajectory, snapshot = map(int, procMapping[trajectory-1][1:-1].split(','))
        epoch = str(epoch)
    sys.stderr.write("Writing pathway...\n")
    with open(outputPath+"pathway.pdb", "a") as f:
        f.write("ENDMDL\n".join(itertools.chain.from_iterable(pathway)))
