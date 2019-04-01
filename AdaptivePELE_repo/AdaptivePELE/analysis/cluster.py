import matplotlib
import mdtraj
import time
matplotlib.use('TkAgg')
import signal
import matplotlib.pyplot as plt
from matplotlib.widgets  import RectangleSelector
import numpy as np
from multiprocessing import Pool
import os
import hdbscan
import sklearn.metrics as mt
from tqdm import tqdm, trange
import prody
import errno
import argparse
import pandas as pd
import glob
import re
import sys
import warnings
import AdaptivePELE.analysis.splitTrajectory as st
import AdaptivePELE.analysis.backtrackAdaptiveTrajectory as bk

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
FREQ = 1
REPORT = "report"
TRAJ = "trajectory"
ACCEPTED_STEPS = 'numberOfAcceptedPeleSteps'
OUTPUT_FOLDER = 'BestStructs'
DIR = os.path.abspath(os.getcwd())
STEPS = 3


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("crit1", type=int, help="First Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: Binding Energy")
    parser.add_argument("crit2", type=int, help="Second Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: Binding Energy")
    parser.add_argument("ad_steps", type=int, help="Adaptive Steps i.e: self.ad_steps")
    parser.add_argument("--path", type=str, help="Path to Pele's results root folder i.e: path=/Pele/results/", default=DIR)
    parser.add_argument("--ofreq", "-f", type=int, help="Every how many steps the trajectory were outputted on PELE i.e: self.ad_steps", default=FREQ)
    parser.add_argument("--out", "-o", type=str, help="Output Path. i.e: BindingEnergies_apo", default=OUTPUT_FOLDER)
    parser.add_argument("--numfolders", "-nm", action="store_true", help="Not to parse non numerical folders")
    parser.add_argument("--top", "-t", type=str, help="Topology file for xtc", default=None)
    parser.add_argument("--first", action="store_true", help="Skip first line")
    parser.add_argument("--zcol", type=int, help="Thirs Criteria we want to rank and output the strutures for. Must be a column of the report. i.e: SASA", default=2)
    parser.add_argument("--resname", type=str, help="Resname of the ligand. Resquested for clusterization", default="LIG")
    parser.add_argument("--percentage", type=int, help="Percentage of snapshots taken based on the first criteria to clusterize", default=30)
    parser.add_argument("--thresh", type=float, help="Treshold of the first criteria below which it will be clusterize. i.e --thresh -60 (binding energy)", default=None)
    parser.add_argument("--cpus", type=int, help="Number of workers to use", default=1)
    args = parser.parse_args()

    return args.crit1, args.crit2, args.zcol, args.ad_steps, os.path.abspath(args.path), args.ofreq, args.out, args.numfolders, args.top, args.first, args.resname, args.percentage, args.thresh, args.cpus


def is_adaptive():
    folders = glob.glob("{}/*/".format(DIR))
    folders_numerical = [os.path.basename(os.path.normpath(folder)) for folder in folders]
    if len(folders_numerical) > 1:
        return True
    else:
        return False


def extract_ligand_coords(info):
    path, snapshot = info
    traj = mdtraj.load_frame(path, snapshot, top="topology.pdb")
    atoms_info = traj.topology.to_dataframe()[0]
    condition = atoms_info['resName'] == resname
    atom_numbers_ligand = atoms_info[condition].serial.tolist()
    coords = []
    for atom_num in atom_numbers_ligand:
        try:
            coords.extend(traj.xyz[0, atom_num-1].tolist())
        except IndexError:
            continue
    return coords



def cluster(data, resName, out_freq=1, cpus=1):
    print("Clusterizing")

    # Extract data
    epoch = data[DIR].tolist()
    trajectory = data[REPORT].tolist()
    snapshot = data.iloc[:, 4].tolist()
    snapshots = ["{}/*trajectory_{}.*".format(os.path.basename(os.path.dirname(e)), traj) for e, traj in zip(epoch, trajectory)]

    # Get Files
    paths = []
    for s in snapshots:
        if glob.glob(s):
            paths.extend(glob.glob(s))
        else:
            paths.extend(glob.glob(os.path.basename(s)))

    # Extract atom coordinates from files
    # Could be parallelize in a future
    t0 = time.time()
    print("Using {} cpus".format(cpus))
    p = Pool(processes=cpus)
    #Python2.7 compatibility
    pool_input = [[pt, s] for pt, s in zip(paths, snapshot)]
    all_coords = p.map(extract_ligand_coords, pool_input)
    p.close()
    t1 = time.time()
    print("Time extract atom coords")
    print(t1-t0)

    # Clusterize
    cluster_with_dbscan(paths, snapshot, all_coords, topology=topology)


def main(criteria1, criteria2, criteria3, ad_steps, path=DIR, out_freq=FREQ, output=OUTPUT_FOLDER, numfolders=False, topology=None, skip_first=False,
    resname="LIG", percentage=30, threshold=None, cpus=1):
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
    # Check whether is adaptive simulation or not
    adaptive = is_adaptive()

    # Find reports
    reports = find_reports(path, numfolders)

    # Retrieve Column Names from report
    steps, crit1_name, crit2_name, crit3_name = get_column_names(reports, STEPS, criteria1, criteria2, criteria3)

    # Retrieve Data from reports
    data = parse_values(reports, criteria1, criteria2, steps, crit1_name, crit2_name, skip_first)

    # Filter values
    if threshold:
        filter_data = Filter(data, percentage, crit1_name, thresh=threshold)
    else:
        filter_data = Filter(data, percentage, crit1_name)
        
    # cluster
    cluster(filter_data, resname, out_freq, cpus)


def Filter(values, percentage, column_name, thresh=None):
    if thresh:
        condition = (values[column_name] < thresh)
        data_filtered = values[condition]
    else:
        n_values_to_take = int(values.shape[0]*percentage/100)
        data_filtered = values.nlargest(n_values_to_take, column_name)
    print("Filtered simulation data {}".format(data_filtered.shape))
    return data_filtered


def find_reports(path, numfolders):
    reports = glob.glob(os.path.join(path, "*/*report_*"))
    reports = glob.glob(os.path.join(path, "*report_*")) if not reports else reports
    reports = filter_non_numerical_folders(reports, numfolders)
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")
    return reports


def parse_values(reports, criteria1, criteria2, steps, crit1_name, crit2_name, first=False, cpus=1):
    """

       Description: Parse the 'reports' and create a sorted array
       of size n_structs following the criteria chosen by the user.

    """
    
    info_reports = [ retrieve_report_data(report) for report in reports]
    data = pd.concat(info_reports)
    data.drop_duplicates(subset=[crit1_name, crit2_name], inplace=True)
    print("Simulation data {}".format(data.shape))
    return data


def retrieve_report_data(report):
    # Get report
    report_number = os.path.basename(report).split("_")[-1]
    # Read data
    try:
        data = pd.read_csv(report, sep='    ', engine='python')
    except pd.errors.EmptyDataError:
        warnings.warn("Report {} corrupted".format(report), UserWarning)
    # Skip first line if asked
    # if first and os.path.basename(os.path.dirname(report)):
        # data = data.iloc[1:]
    # Insert path and filename
    data.insert(0, DIR, [report]*data.shape[0])
    data.insert(1, REPORT, [report_number]*data.shape[0])
    return data


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


def get_column_names(reports, steps, criteria1, criteria2, criteria3):
    data = pd.read_csv(reports[0], sep='    ', engine='python')
    data = list(data)

    return data[int(steps)-1], data[criteria1-1], data[criteria2-1], data[criteria3-1]


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def cluster_with_dbscan(paths, snapshots, all_coordinates, out_freq=1, topology=None):

    """
    Use high performance computing hdbscan
    to do an all-atom cluster of the chosen
    plot structures
    """

    n_samples = len(snapshots)

    # Clusterize
    labels = []
    results = []
    t0 = time.time()
    try:
        db = hdbscan.HDBSCAN(min_samples=int(n_samples*0.10)+1).fit(all_coordinates)
    except ValueError:
        raise ValueError("Ligand not found check the option --resname. i.e python interactive.py 5 6 7 --resname LIG")
    result = db.labels_
    labels.append(len(set(result)))
    results.append(result)
    t1 = time.time()
    print("time clustering")
    print(t1-t0)

    # Get Best Result
    t0 = time.time()
    mx_idx = np.argmax(np.array(labels))
    final_result = results[mx_idx]
    try:
        silhouette_samples = mt.silhouette_samples(all_coordinates, final_result)
    except ValueError:
        raise ValueError("Clustering failed. Structures do not follow any pattern or they are not enough")
    max_clust = {label: [path, snap, sil] for (path, snap, label, sil) in zip(paths, snapshots, final_result, silhouette_samples)}

    # Get representative
    for path, snapshot, label, sil in zip(paths, snapshots, final_result, silhouette_samples):
        if sil > max_clust[label][2]:
            max_clust[label] = [path, snapshot, sil]

    # Get Structures
    for i, (label, info) in enumerate(max_clust.items()):
        # if label == -1: continue
        output = "Clusters"
        if not os.path.exists(output):
            os.mkdir(output)
        f_out = "cluster_{}.pdb".format(label+1)
        f_in, snapshot, _ = info

        #XTC
        if topology:
            found = st.main(output, [f_in, ], topology, [snapshot, ], template=f_out)
            if found:
                print("MODEL {} has been selected as {}".format(f_in, f_out))
            else:
                print("MODEL {} not found. Check -f option".format(f_in))
        #PDB
        else:
            traj = []
            model = (snapshot)/out_freq+1

            with open(f_in, 'r') as input_file:
                file_content = input_file.read()
                trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' %int(model), file_content, re.DOTALL)
            with open(os.path.join(output, f_out),'w') as f:
                traj.append("MODEL     %d" % int(model))
                try:
                    traj.append(trajectory_selected.group(1))
                except AttributeError:
                    raise AttributeError("Model {} not found. Check the -f option.".format(f_in))
                traj.append("ENDMDL\n")
                f.write("\n".join(traj))
            print("MODEL {} has been selected as {}".format(f_in, f_out))
    t1 = time.time()
    print("Time post processing")
    print(t1-t0)


if __name__ == "__main__":
    criteria1, criteria2, criteria3, ad_steps, path, out_freq, output, numfolders, topology, skip_first, resname, percentage, thresh, cpus = parse_args()
    main(criteria1, criteria2, criteria3, ad_steps, path, out_freq, output, numfolders, topology, skip_first, resname, percentage, thresh, cpus)
