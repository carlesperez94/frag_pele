from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import argparse
import numpy as np
import mdtraj as md
import multiprocessing as mp
from AdaptivePELE.utilities import utilities
from AdaptivePELE.analysis import correctRMSD
from AdaptivePELE.freeEnergies.extractCoords import getTopologyObject


def parseArguments():
    """
        Parse the command-line options

        :returns: object -- Object containing the options passed
    """
    desc = "Program that caculates the relative SASA of a ligand."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("resname", type=str, help="Ligand resname")
    parser.add_argument("--path", type=str, default=".", help="Path where the simulation is stored")
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories or path to Adaptive topology object")
    parser.add_argument("--out_name", type=str, default="fixedReport", help="Name of the modified report files (default is fixedReport)")
    parser.add_argument("--out_folder", type=str, default=None, help="Path where to store the report files (default is fixedReport)")
    parser.add_argument("-n", type=int, default=1, help="Number of processors to parallelize")
    parser.add_argument("--fmt_str", type=str, default="%.4f", help="Format of the output file (default is .4f which means all floats with 4 decimal points)")
    parser.add_argument("--new_report", action="store_true", help="Whether to create new report files instead of modifying existing ones")
    args = parser.parse_args()

    return args.resname, args.path, args.top, args.out_name, args.fmt_str, args.n, args.out_folder, args.new_report


def calculateSASA(trajectory, topology, res_name):
    """
        Calculate the SASA of a ligand in a trajectory

        :param trajectory: Name of the trajectory file
        :type trajectory: str
        :param topology: Topology of the trajectory (needed for non-pdb trajs)
        :type topology: str
        :param res_name: Ligand resname
        :type res_name: str
    """
    t = md.load(trajectory, top=topology)
    res_atoms = t.top.select("resname '%s'" % res_name)
    t2 = t.atom_slice(res_atoms)
    for atom in t2.top.atoms:
        # mdtraj complains if the ligand residue index is not 0 when
        # isolated
        atom.residue.index = 0
    res_index = t.top.atom(res_atoms[0]).residue.index
    sasa = md.shrake_rupley(t, mode="residue")
    sasa_empty = md.shrake_rupley(t2, mode="residue")
    return sasa[:, res_index]/sasa_empty[:, 0]


def process_file(traj, top_file, resname, report, outputFilename, format_out, new_report, epoch):
    sasa_values = calculateSASA(traj, top_file, resname)
    header = ""
    if not new_report:
        try:
            reportFilename = glob.glob(report)[0]
        except IndexError:
            raise IndexError("File %s not found" % report)

        with open(reportFilename) as f:
            header = f.readline().rstrip()
            if not header.startswith("#"):
                header = ""
            reportFile = utilities.loadtxtfile(f)

        fixedReport = correctRMSD.extendReportWithRmsd(reportFile, sasa_values)
    else:
        fixedReport = sasa_values

    with open(outputFilename, "w") as fw:
        if header:
            fw.write("%s\tSASA\n" % header)
        else:
            fw.write("# SASA\n")
        np.savetxt(fw, fixedReport, fmt=format_out, delimiter="\t")


def process_folder(epoch, folder, trajName, reportName, output_filename, top):
    if epoch is None:
        allTrajs = glob.glob(os.path.join(folder, trajName))
        full_reportName = os.path.join(folder, reportName)
    else:
        allTrajs = glob.glob(os.path.join(folder, epoch, trajName))
        full_reportName = os.path.join(folder, epoch, reportName)
        epoch = int(epoch)

    allFiles = []
    for traj in allTrajs:
        trajNum = utilities.getTrajNum(traj)
        if top is not None:
            top_file = top.getTopologyFile(epoch, trajNum)
        else:
            top_file = None
        report_file = full_reportName % trajNum
        allFiles.append((traj, report_file, top_file, epoch, output_filename % trajNum))
    return allFiles


def main(resname, folder, top, out_report_name, format_out, nProcessors, output_folder, new_report):
    """
        Calculate the relative SASA values of the ligand

        :param resname: Ligand resname
        :type resname: str
        :param folder: Path the simulation
        :type folder: str
        :param top: Path to the topology
        :type top: str
        :param out_report_name: Name of the output file
        :type out_report_name: str
        :param format_out: String with the format of the output
        :type format_out: str
        :param nProcessors: Number of processors to use
        :type nProcessors: int
        :param output_folder: Path where to store the new reports
        :type output_folder: str
        :param new_report: Whether to create new reports
        :type new_report: bool
    """
    # Constants
    if output_folder is not None:
        out_report_name = os.path.join(output_folder, out_report_name)
    outputFilename = "_".join([out_report_name, "%d"])
    trajName = "*traj*"
    reportName = "*report*_%d"
    if nProcessors is None:
        nProcessors = utilities.getCpuCount()
    nProcessors = max(1, nProcessors)
    print("Calculating SASA with %d processors" % nProcessors)
    pool = mp.Pool(nProcessors)
    epochs = utilities.get_epoch_folders(folder)
    if top is not None:
        top_obj = getTopologyObject(top)
    else:
        top_obj = None
    files = []
    if not epochs:
        # path does not contain an adaptive simulation, we'll try to retrieve
        # trajectories from the specified path
        files = process_folder(None, folder, trajName, reportName, os.path.join(folder, outputFilename), top_obj)
    for epoch in epochs:
        print("Epoch", epoch)
        files.extend(process_folder(epoch, folder, trajName, reportName, os.path.join(folder, epoch, outputFilename), top_obj))
    results = []
    for info in files:
        results.append(pool.apply_async(process_file, args=(info[0], info[2], resname, info[1], info[4], format_out, new_report, info[3])))
    for res in results:
        res.get()

if __name__ == "__main__":
    lig_name, path, topology_path, out_name, fmt_str, n_proc, out_folder, new_reports = parseArguments()
    main(lig_name, path, topology_path, out_name, fmt_str, n_proc, out_folder, new_reports)
