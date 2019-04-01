from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.utilities import utilities
import argparse
import os


def parseArguments():
    desc = "Program that converts a trajectory from a non-pdb format into a pdb."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("trajectory", help="Trajectory files to convert")
    parser.add_argument("-o", type=str, default="", help="Output file")
    parser.add_argument("--dir", type=str, default="", help="Output dir")
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories")
    args = parser.parse_args()
    return args.trajectory, args.o, args.dir, args.top


def main(trajectory, output_file, output_path, topology):
    if output_path:
        utilities.makeFolder(output_path)
    if not output_file:
        output_file = "%s.pdb" % os.path.splitext(os.path.split(trajectory)[1])[0]
    utilities.convert_trajectory_to_pdb(trajectory, topology, output_file, output_path)


if __name__ == "__main__":
    traj_file, out_file, output_dir, top = parseArguments()
    main(traj_file, out_file, output_dir, top)
