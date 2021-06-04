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
import re
from ast import literal_eval
from AdaptivePELE.utilities import utilities
try:
    basestring
except NameError:
    basestring = str


def parseArguments():
    """
        Parse the command-line options

        :returns: str, str, str, str --  path to file to backtrack,
            path to the result files,
            output path where to write the files, name of the files
    """
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("file_to_backtrack", type=str, help="File of the selected_results folder that you want to"
                                                            "backtrack.")

    parser.add_argument("-p", type=str, default="growing_results", help="Path to the folder which contains the growing"
                                                                        "results of all growing steps.")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    parser.add_argument("--name", type=str, default="pathway.pdb", help="Name of the pdb to write the files")
    args = parser.parse_args()
    return args.file_to_backtrack, args.p, args.o, args.name


def main(file_to_backtrack, results_path, outputPath, out_filename):
    """

    :param file_to_backtrack: File of the selected_results folder that you want to backtrack.
    :type file_to_backtrack: str
    :param results_path: Path where the growing simulations are stored.
    :type results_path: str
    :param outputPath: Output folder path.
    :type outputPath: str
    :param out_filename: Output filename prefix.
    :type out_filename: str
    :return: None
    """
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
    sys.stderr.write("Creating pathway...\n")
    # Get information of the input structure
    trajectory, snapshot, growing_id = extract_info_from_selected_file(path_to_selected_file=file_to_backtrack)
    # Obtain trajectory from sampling folder
    sampling_file_from = "sampling_result_{}/trajectory_{}.pdb".format(growing_id, trajectory)
    sampling_models = utilities.getSnapshots(sampling_file_from)
    # Get all snapshots from the input snapshot and add it to pathway
    snapshots = sampling_models[1:snapshot + 1]
    pathway.insert(0, snapshots)
    growing_epochs = glob.glob(os.path.join(results_path, "{}_growing_output*".format(growing_id)))
    if not growing_epochs:
        raise ValueError("Trajectory %s not found!" % os.path.join(results_path))
    for n in range(len(growing_epochs)-1, -1, -1):
        # Reading mapping information
        procMapping = open(os.path.join(results_path, "{}_growing_output{}".format(growing_id, n),
                                        "processorMapping.txt")).read().split(":")
        # Extract the spawning structure of the previous simulation
        cluster_from = procMapping[trajectory-1].rstrip("\n")
        cluster_from_tup = literal_eval(cluster_from)
        # Update the trajectory and snapshot to be used
        trajectory = cluster_from_tup[1]
        snapshot = cluster_from_tup[2]
        # Extract the filename and the snapshot that has spawn
        filename = os.path.join(results_path, "{}_growing_output{}/trajectory_{}.pdb".format(growing_id, n, trajectory))
        print(n, filename, trajectory, snapshot, procMapping)
        snapshots = utilities.getSnapshots(filename)
        if n == 0:
            initial = 0
        else:
            initial = 1
        # Take all MODELS from the snapshot
        snapshots = snapshots[initial:snapshot+1]
        pathway.insert(0, snapshots)

    sys.stderr.write("Writing pathway...\n")
    with open(outputPath+out_filename, "a") as f:
        f.write("ENDMDL\n".join(itertools.chain.from_iterable(pathway)))


def extract_info_from_selected_file(path_to_selected_file):
    pattern_trajectory = re.compile('trajectory_[\d+]\.')
    pattern_snapshot = re.compile('\.[\d+]_')
    trajectory_and_others = re.findall(pattern_trajectory, path_to_selected_file)
    trajectory = re.findall('\d+', trajectory_and_others[0])
    snapshot_and_others = re.findall(pattern_snapshot, path_to_selected_file)
    snapshot = re.findall('\d+', snapshot_and_others[0])
    growing_id = path_to_selected_file.split("epochsampling_result_")[-1].split("_trajectory")[0]
    return int(trajectory[0]), int(snapshot[0]), growing_id


if __name__ == "__main__":
    traj, results_path, output_path, output_filename = parseArguments()
    main(traj, results_path, output_path, output_filename)
