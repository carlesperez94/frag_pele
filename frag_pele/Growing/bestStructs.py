#!/usr/bin/python2.7

import os
import errno
import argparse
import pandas as pd
import glob
import re
import sys
import socket

"""

   Description: Parse all the reports found under 'path' and sort them all
   by the chosen criteria (Binding Energy as default) having into account the
   frequency our pele control file writes a structure through the -ofreq param
   (1 by default). To sort from higher to lower value use -f "max" otherwise
   will rank the structures from lower to higher criteria's values. The number
   of structures will be ranked is controlled by -i 'nstruct' (default 10).

   For any problem do not hesitate to contact us through the email address written below.

"""
machine = socket.getfqdn()
__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# DEFAULT VALUES
ORDER = "min"
CRITERIA = ["Binding", "Energy"]
OUTPUT = "Structure_{}.pdb"
N_STRUCTS = 1
FREQ = 1
REPORT = "report"
TRAJ = "trajectory"
ACCEPTED_STEPS = 'numberOfAcceptedPeleSteps'
if "bsccv" in machine:
    ACCEPTED_STEPS = 'AcceptedSteps'
DIR = os.path.abspath(os.getcwd())


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("filename", type=str, help="Criteria we want to rank and output the strutures for. Must be a clumn of the report. i.e: Binding Energy")
    parser.add_argument("crit", type=str, nargs='+', help="Criteria we want to rank and output the strutures for. Must be a clumn of the report. i.e: Binding Energy")
    parser.add_argument("--steps", "-as", type=str, help="Name of the accepted steps column in the report files. i.e: numberOfAcceptedPeleSteps", default=ACCEPTED_STEPS)
    parser.add_argument("--path", type=str, help="Path to Pele's results root folder i.e: path=/Pele/results/", default=DIR)
    parser.add_argument("--nst", "-n", type=int, help="Number of produced structures. i.e: 20" , default=N_STRUCTS)
    parser.add_argument("--sort", "-s", type=str, help="Look for minimum or maximum value --> Options: [min/max]. i.e: max", default=ORDER)
    parser.add_argument("--ofreq", "-f", type=int, help="Every how many steps the trajectory were outputted on PELE i.e: 4", default=FREQ)
    parser.add_argument("--out", "-o", type=str, help="Output Path. i.e: BindingEnergies_apo", default="".join(CRITERIA))
    parser.add_argument("--numfolders", "-nm", action="store_true", help="Not to parse non numerical folders")
    args = parser.parse_args()

    return args.filename, os.path.abspath(args.path), " ".join(args.crit), args.nst, args.sort, args.ofreq, args.out, args.steps, args.numfolders


def main(criteria, file_name, path=DIR, n_structs=10, sort_order="min", out_freq=FREQ, output="".join(CRITERIA), steps = ACCEPTED_STEPS, numfolders=False):
    """

      Description: Rank the traj found in the report files under path
      by the chosen criteria. Finally, output the best n_structs.

      Input:

         Path: Path to look for *report* files in all its subfolders.

         Criteria: Criteria to sort the structures.
         Needs to be the name of one of the Pele's report file column.
         (Default= "Binding Energy")

         n_structs: Numbers of structures to create.

         sort_order: "min" if sorting from lower to higher "max" from high to low.

         out_freq: "Output frequency of our Pele control file"

     Output:

        f_out: Name of the n outpu
    """

    reports = glob.glob(os.path.join(path, "*report_*"))
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")

    # Data Mining
    min_values, steps = parse_values(reports, n_structs, criteria, sort_order, steps)
    values = min_values[criteria].tolist()
    paths = min_values[DIR].tolist()
    epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in paths]
    file_ids = min_values.report.tolist()
    step_indexes = min_values[steps].tolist()
    files_out = ["epoch{}_trajectory_{}.{}_{}{}.pdb".format(epoch, report, int(step), criteria.replace(" ",""), value) \
       for epoch, step, report, value in zip(epochs, step_indexes, file_ids, values)]
    files_out_tmp_tuple = [("epoch{}_trajectory_{}.{}_{}{}.pdb".format(epoch, report, int(step),
                                                                                criteria.replace(" ",""),
                                                                                value), value) \
       for epoch, step, report, value in zip(epochs, step_indexes, file_ids, values)]
    files_out_best = sorted(files_out_tmp_tuple, key=lambda x: x[1])[0][0]
    for f_id, f_out, step, path, n in zip(file_ids, files_out, step_indexes, paths, range(0, n_structs)):

        # Read Trajetory from PELE's output
        f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.pdb".format(f_id)))
        if len(f_in) == 0:
            sys.exit("Trajectory {} not found. Be aware that PELE trajectories must contain the label \'trajectory\' in their file name to be detected".format("*trajectory*_{}".format(f_id)))
        f_in = f_in[0]
        with open(f_in, 'r') as input_file:
            file_content = input_file.read()
        trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' %int((step)/out_freq+1), file_content,re.DOTALL)

        traj = []
        with open(os.path.join(file_name, f_out), 'w') as f:
            traj.append("MODEL     %d" %int((step)/out_freq+1))
            try:
                traj.append(trajectory_selected.group(1))
            except AttributeError:
                raise AttributeError("Model not found. Check the -f option.")
            traj.append("ENDMDL\n")
            f.write("\n".join(traj))
        print("MODEL {} has been selected".format(f_out))
    return files_out_best, files_out


def parse_values(reports, n_structs, criteria, sort_order, steps):
    """

       Description: Parse the 'reports' and create a sorted array
       of size n_structs following the criteria chosen by the user.

    """

    INITIAL_DATA = [(DIR, []),
                    (REPORT, []),
                    (steps, []),
                    (criteria, [])
                    ]
    min_values = pd.DataFrame.from_dict(dict(INITIAL_DATA))
    for f in reports:
        report_number = os.path.basename(f).split("_")[-1]
        data = pd.read_csv(f, sep='    ', engine='python')
        # Skip first line not to get initial structure
        if not steps in list(data): steps="AcceptedSteps"
        selected_data = data.loc[1:, [steps, criteria]]
        if sort_order == "min":
                report_values =  selected_data.nsmallest(n_structs, criteria)
                report_values.insert(0, DIR, [f]*report_values[criteria].size)
                report_values.insert(1, REPORT, [report_number]*report_values[criteria].size)
                mixed_values = pd.concat([min_values, report_values])
                mixed_values.index = pd.RangeIndex(len(mixed_values.index))
                min_values = mixed_values.nsmallest(n_structs, criteria)

        else:
                report_values =  selected_data.nlargest(n_structs, criteria)
                report_values.insert(0, DIR, [f]*report_values[criteria].size)
                report_values.insert(1, REPORT, [report_number]*report_values[criteria].size)
                mixed_values = pd.concat([min_values, report_values])
                min_values = mixed_values.nlargest(n_structs, criteria)
    return min_values, steps


def filter_non_numerical_folders(reports, numfolders):
    """
    Filter non numerical folders among
    the folders to parse
    """
    if(numfolders):
        new_reports = [report for report in reports if(os.path.basename(os.path.dirname(report)).isdigit())]
        return new_reports
    else:
        return reports


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == "__main__":
    filename, path, crit, nst, sort, ofreq, out, steps, numfolders = parse_args()
    main(crit, filename, path,  nst, sort, ofreq, out, steps, numfolders)
