#!/usr/bin/python2.7

import os
import errno
import argparse
import pandas as pd
import glob
import re
import sys
from AdaptivePELE.analysis import splitTrajectory

"""

   Description: Parse all the reports found under 'path' and sort them all
   by the chosen criteria (Binding Energy as default) having into account the
   frequency our pele control file writes a structure through the -ofreq param
   (1 by default). To sort from higher to lower value use -f "max" otherwise
   will rank the structures from lower to higher criteria's values. The number
   of structures will be ranked is controlled by -i 'nstruct' (default 10).

   For any problem do not hesitate to contact us through the email address written below.

"""

__author__ = "Daniel Soler Viladrich"
__email__ = "daniel.soler@nostrumbiodiscovery.com"

# DEFAULT VALUES
ORDER = "min"
CRITERIA = ["Binding", "Energy"]
OUTPUT = "Structure_{}.pdb"
N_STRUCTS = 10
FREQ = 1
REPORT = "report"
TRAJ = "trajectory"
ACCEPTED_STEPS = 'numberOfAcceptedPeleSteps'
OUTPUT_FOLDER = 'BestStructs'
DIR = os.path.abspath(os.getcwd())
STEPS=3
HELP = "USE:\n\n- For xtc: python bestStructs.py 5 --top topology.pdb\n\n- For pdb:  python bestStructs.py 5"


def parse_args():

    parser = argparse.ArgumentParser(description=HELP, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("crit", type=str, nargs='+', help="Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: Binding Energy")
    parser.add_argument("--crit2", type=str, nargs='+', help="Second Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: SasaLig, i.e.2. 3", default="")
    parser.add_argument("--path", type=str, help="Path to Pele's results root folder i.e: path=/Pele/results/", default=DIR)
    parser.add_argument("--nst", "-n", type=int, help="Number of produced structures. i.e: 20" , default=N_STRUCTS)
    parser.add_argument("--sort", "-s", type=str, help="Look for minimum or maximum value --> Options: [min/max]. i.e: max", default=ORDER)
    parser.add_argument("--ofreq", "-f", type=int, help="Every how many steps the trajectory were outputted on PELE i.e: 4", default=FREQ)
    parser.add_argument("--out", "-o", type=str, help="Output Path. i.e: BindingEnergies_apo", default=OUTPUT_FOLDER)
    parser.add_argument("--top", "-t", type=str, help="PDB topology file. Mandatory when using xtf", default=None)
    parser.add_argument("--numfolders", "-nm", action="store_true", help="Not to parse non numerical folders")
    args = parser.parse_args()

    return os.path.abspath(args.path), " ".join(args.crit), " ".join(args.crit2), args.nst, args.sort, args.ofreq, args.out, args.numfolders, args.top


def main(criteria, path=DIR, n_structs=10, sort_order="min", out_freq=FREQ, output=OUTPUT_FOLDER, numfolders=False, criteria2=None, topology=None):
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
    #Get reports files
    reports = glob.glob(os.path.join(path, "*/*report*"))
    reports = glob.glob(os.path.join(path, "*report*")) if not reports else reports
    reports = filter_non_numerical_folders(reports, numfolders)
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")

    #Get metrics of reports
    if criteria.isdigit():
        steps, criteria = get_column_names(reports, STEPS, criteria)
    else:
        steps = get_column_names(reports, STEPS, criteria)

    #Data Mining
    if not criteria2:
        min_values = parse_values(reports, n_structs, criteria, sort_order, steps)
    else:
        min_values = parse_values(reports, n_structs, criteria, criteria2, sort_order, steps)

    #Process data to be outputted
    values = min_values[criteria].tolist()
    paths = min_values[DIR].tolist()
    epochs = [os.path.basename(os.path.normpath(os.path.dirname(Path))) for Path in paths]
    file_ids = min_values.report.tolist()
    step_indexes = min_values[steps].tolist()
    files_out = ["epoch{}_trajectory_{}.{}_{}{}.pdb".format(epoch, report, int(step), criteria.replace(" ",""), value) \
       for epoch, step, report, value in zip(epochs, step_indexes, file_ids, values)]

    #Read traj and output sanpshot  
    for f_id, f_out, step, path in zip(file_ids, files_out, step_indexes, paths):
        if not topology:
            extract_snapshot_from_pdb(path, f_id, output, topology, step, out_freq, f_out)
        else:
            extract_snapshot_from_xtc(path, f_id, output, topology, step, out_freq, f_out)

    #Return data
    return files_out, epochs, file_ids, step_indexes


def extract_snapshot_from_pdb(path, f_id, output, topology, step, out_freq, f_out):
    f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.pdb".format(f_id)))
    if not f_in:
        f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.*".format(f_id)))
    if len(f_in) == 0:
        sys.exit("Trajectory {} not found. Be aware that PELE trajectories must contain the label \'trajectory\' in their file name to be detected".format("*trajectory*_{}".format(f_id)))
    f_in = f_in[0]
    with open(f_in, 'r') as input_file:
        file_content = input_file.read()
    trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' % int((step) / out_freq+1), file_content,re.DOTALL)

    # Output Snapshot
    try:
        mkdir_p(output)
    except OSError:
        pass

    traj = []
    with open(os.path.join(output, f_out), 'w') as f:
        traj.append("MODEL     %d" % int((step)/out_freq+1))
        try:
            traj.append(trajectory_selected.group(1))
        except AttributeError:
            raise AttributeError("Model not found. Check the -f option.")
        traj.append("ENDMDL\n")
        f.write("\n".join(traj))
    print("Model {} selected".format(os.path.join(output, f_out)))

def extract_snapshot_from_xtc(path, f_id, output, topology, step, out_freq, f_out):
    f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.xtc".format(f_id)))
    if not f_in:
        f_in = glob.glob(os.path.join(os.path.dirname(path), "*trajectory*_{}.*".format(f_id)))
    if len(f_in) == 0: 
        sys.exit("Trajectory {} not found. Be aware that PELE trajectories must contain the label \'trajectory\' in their file name to be detected".format("*trajectory*_{}".format(f_id)))
    splitTrajectory.main(output, [f_in[0], ], topology, [(step)/out_freq+1, ], template= f_out)
    print("Model {} selected".format(f_out))
    
 


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
    min_values = pd.DataFrame.from_items(INITIAL_DATA)
    for file in reports:
        report_number = os.path.basename(file).split("_")[-1]
        data = pd.read_csv(file, sep='    ', engine='python')
        selected_data = data.loc[:, [steps, criteria]]
        if sort_order == "min":
                report_values = selected_data.nsmallest(n_structs, criteria)
                report_values.insert(0, DIR, [file]*report_values[criteria].size)
                report_values.insert(1, REPORT, [report_number]*report_values[criteria].size)
                mixed_values = pd.concat([min_values, report_values])
                min_values = mixed_values.nsmallest(n_structs, criteria)

        else:
                report_values = selected_data.nlargest(n_structs, criteria)
                report_values.insert(0, DIR, [file]*report_values[criteria].size)
                report_values.insert(1, REPORT, [report_number]*report_values[criteria].size)
                mixed_values = pd.concat([min_values, report_values])
                min_values = mixed_values.nlargest(n_structs, criteria)
    return min_values


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


def get_column_names(reports, steps, criteria):
    data = pd.read_csv(reports[0], sep='    ', engine='python')
    data = list(data)
    if criteria.isdigit():
        return data[int(steps)-1], data[int(criteria)-1]
    else:
        return data[int(steps)-1]


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == "__main__":
    path, criteria, criteria2, interval, sort_order, out_freq, output, numfolders, top = parse_args()
    main(criteria, path, interval, sort_order, out_freq, output, numfolders, criteria2, top)
