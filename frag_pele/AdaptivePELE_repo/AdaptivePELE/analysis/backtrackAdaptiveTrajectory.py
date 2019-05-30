"""
    Recreate the trajectory fragments to the led to the discovery of a snapshot,
    specified by the tuple (epoch, trajectory, snapshot) and write as a pdb file
"""
from __future__ import print_function
import os
import sys
import argparse
import glob
import itertools
from AdaptivePELE.utilities import utilities


def parseArguments():
    """
        Parse the command-line options

        :returns: :py:class:`.Clustering`, int, int, int, str -- Clustering
            object, number of trajectory, number of snapshot, number of epoch,
            output path where to write the files
    """
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("epoch", type=str, help="Path to the epoch to search the snapshot")
    parser.add_argument("trajectory", type=int, help="Trajectory number")
    parser.add_argument("snapshot", type=int, help="Snapshot to select (in accepted steps)")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    parser.add_argument("--name", type=str, default="pathway.pdb", help="Name of the pdb to write the files")
    args = parser.parse_args()
    return args.trajectory, args.snapshot, args.epoch, args.o, args.name


def main(trajectory, snapshot, epoch, outputPath, out_filename):
    if outputPath is not None:
        outputPath = os.path.join(outputPath, "")
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    else:
        outputPath = ""
    if os.path.exists(outputPath+out_filename):
        # If the specified name exists, append a number to distinguish the files
        name, ext = os.path.splitext(out_filename)
        out_filename = "".join([name, "_%d", ext])
        i = 1
        while os.path.exists(outputPath+out_filename % i):
            i += 1
        out_filename %= i
    pathway = []
    # Strip out trailing backslash if present
    pathPrefix, epoch = os.path.split(epoch.rstrip("/"))
    sys.stderr.write("Creating pathway...\n")
    while True:
        filename = glob.glob(os.path.join(pathPrefix, epoch, "*traj*_%d.pdb" % trajectory))
        snapshots = utilities.getSnapshots(filename[0])
        snapshots = snapshots[:snapshot+1]
        pathway.insert(0, snapshots)
        if epoch == '0':
            # Once we get to epoch 0, we just need to append the trajectory
            # where the cluster was found and we can break out of the loop
            break
        procMapping = open(os.path.join(pathPrefix, epoch, "processorMapping.txt")).read().rstrip().split(':')
        epoch, trajectory, snapshot = map(int, procMapping[trajectory-1][1:-1].split(','))
        epoch = str(epoch)
    sys.stderr.write("Writing pathway...\n")
    with open(outputPath+out_filename, "a") as f:
        f.write("ENDMDL\n".join(itertools.chain.from_iterable(pathway)))


if __name__ == "__main__":
    trajectory, snapshot, epoch, outputPath, out_filename = parseArguments()
    main(trajectory, snapshot, epoch, outputPath, out_filename)
