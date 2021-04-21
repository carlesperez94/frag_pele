import sys
import os
import re
import numpy as np
import joblib
import mdtraj
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-pdbr", "--pdb_reference",required=True,
                                help="""Reference pdb file.""")
    required_named.add_argument("-p", "--fol_path", required=True,
                                help="""Path to folder to be analyzed. """)

    parser.add_argument("-r", "--resname", default="GRW",
                        help="Resame of the ligand.")
    parser.add_argument("-if", "--index_relation_file", default=None,
                        help="File containting atom-pairing indexing to"
                             " correct the reference structure.")
    parser.add_argument("-rp", "--report_pref", default='report_',
                        help="Prefix string to find report files")
    parser.add_argument("-np", "--processors", type=int, default=4,
                        help="Number of processors to paralelize the analysis")


    # Plot arguments
    parser.add_argument("-y", "--ycol", default="RMSD",
                        help="Column name of the report file that you want to use as X axis in the plot.")
    parser.add_argument("-op", "--plot", default=False, action="store_true",
                        help="Flag to plot the results.")

    args = parser.parse_args()

    return args.pdb_reference, args.fol_path, args.resname, \
           args.index_relation_file, args.ycol, args.plot, \
           args.report_pref, args.processors


def check_common_residues(pdb_target_prody, pdb_reference_prody):
    list_of_types_and_resnums_target= []
    list_of_types_and_resnums_reference = []
    for resnum, resname in zip(pdb_target_prody.getResnums(), pdb_target_prody.getResnames()):
        list_of_types_and_resnums_target.append((resnum, resname))
    for resnum, resname in zip(pdb_reference_prody.getResnums(), pdb_reference_prody.getResnames()):
        list_of_types_and_resnums_reference.append((resnum, resname))
    common_residues = sorted(list(set(list_of_types_and_resnums_target).intersection(list_of_types_and_resnums_reference)))
    only_resnums = []
    for residue in common_residues:
        resnum = residue[0]
        only_resnums.append(resnum)
    return only_resnums


def compute_rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)


def ligands_rmsd_calculator(pdb_target, pdb_reference, resname="GRW", index_relation=None):
    """

    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param write2report: if True export results in a file
    :param write2pdb: pdb file with the result of the superposition between "pdb_target" and "pdb_reference"
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD of the ligand
    """
    # Load the data to mdtraj
    target = mdtraj.load(pdb_target)
    reference = mdtraj.load(pdb_reference)
    # Get the indexes of heavy atoms of the ligands
    ligand = reference.topology.select('resname "{}" and symbol != H'.format(resname))
    ligand_target = target.topology.select('resname "{}" and symbol != H'.format(resname))
    print("INDEXES REFERENCE")
    print(ligand)
    print("INDEXES TARGET")
    print(ligand_target)
    ligand_traj_ref = reference.atom_slice(ligand)
    ligand_traj_target = target.atom_slice(ligand_target)
    ref_coords = ligand_traj_ref.xyz
    target_coords = ligand_traj_target.xyz
    if index_relation:
        with open(index_relation, "r") as in_f:
            indexes = in_f.readlines()
        list_of_pairs = []
        indexes_target = []
        indexes_ref = []
        # Read
        for idx_pair in indexes:
            idx_ref, idx_tar = idx_pair.split()
            list_of_pairs.append([int(idx_ref), int(idx_tar)])
            indexes_target.append(int(idx_tar))
            indexes_ref.append(int(idx_ref))
        # Get ordering for indexes from the target
        index_target_order = {}
        for n, index in enumerate(sorted(indexes_target)):
            index_target_order[index] = n
        # Get ordering for indexes from the reference
        index_ref_order = {}
        for n, index in enumerate(sorted(indexes_ref)):
            index_ref_order[index] = n
        new_coords = {}
        print("NEW LIGAND TARGET INDEXING RELATION:")
        for pair in sorted(list_of_pairs):
            idx_ref, idx_tar = pair
            print(index_ref_order[idx_ref], index_target_order[idx_tar])
            print(ligand_traj_ref.topology.atom(index_ref_order[idx_ref]), 
                  ligand_traj_target.topology.atom(index_target_order[idx_tar]))
            new_coords[index_target_order[idx_tar]] = ref_coords[0][index_ref_order[idx_ref]]
        ref_coords_new = []
        for key in sorted(new_coords):
            coords = new_coords[key]
            ref_coords_new.append(coords)
        ref_coords = []
        ref_coords.append(ref_coords_new)
            
    rmsd_list = []
    for snap in target_coords:
        rmsd = compute_rmsd(ref_coords[0], snap) * 10 # from nm to A
        rmsd_list.append(rmsd)
    return rmsd_list


def repair_pdbs(pdb_file):
    with open(pdb_file) as pdb:
        pdb_content = pdb.readlines()
    pdb_new = []
    for line in pdb_content:
        atom_name = line[12:16]
        if line.startswith("HETATM") and atom_name.strip().startswith("G"):
            line_new = list(line)
            size = len(atom_name.strip())
            white_spaces = 3-size
            line_new[12:16] = line[77:78] + atom_name.strip() + " "*white_spaces
            modified = "".join(line_new)
            pdb_new.append(modified)
        else:
            pdb_new.append(line)
    pdb_modified = "".join(pdb_new)
    with open(pdb_file, "w") as write_output:
        write_output.write(pdb_modified)


def sort_atoms_by_name(pdb_file):
    with open(pdb_file) as pdb:
        pdb_content = pdb.readlines()
    pdb_content_sorted = (sorted(pdb_content, key=lambda atom_name: atom_name[12:16]))
    new_pdb_content = "".join(pdb_content_sorted)
    with open(pdb_file, "w") as write_output:
        write_output.write(new_pdb_content)

def read_reports(path_to_folder, report_pref):
    """
        This function merge the content of different report for PELE simulations in a single file pandas Data Frame.
    """
    data = []
    report_list = glob.glob('{}[0-9]*'.format(os.path.join(path_to_folder, report_pref)))
    for report in report_list:
        tmp_data = pd.read_csv(report, sep='    ', engine='python')
        tmp_data = tmp_data.iloc[1:]  # We must discard the first row
        processor = re.findall('\d+$'.format(path_to_folder), report)
        tmp_data['Processor'] = processor[0]
        data.append(tmp_data)
    result = pd.concat(data)
    return result

def compute_trajectory_rmsd(traj, fol_path, pdb_reference, resname, 
                            index_relation, report_pref):
    # Compute RMSD
    result = ligands_rmsd_calculator(traj, pdb_reference, resname=resname,
                                     index_relation=index_relation)
    processor = re.findall(r'\d+', traj.split("/")[-1])[0]
    # Store data in reports
    report = "{}".format(os.path.join(fol_path,
                         "{}".format(report_pref + processor)))
    print(report, traj)
    try:
        new_lines = []
        with open(report) as rep:
            rep_lines = rep.readlines()
            rep_lines = [x.strip("\n") for x in rep_lines]
        for ind, line in enumerate(rep_lines):
            new_content = list(line.split("    "))
 
            if new_content[-1] == '':
                new_content = new_content[:-1]
 
            if ind == 0:
                new_content.append('rmsd')
            else:
                value = "{:.3f}".format(result[ind-1])
                new_content.append(value)
            new_line = "    ".join(new_content)
            new_lines.append(new_line)
        new_report = "\n".join(new_lines)
        new_rep_name = report.split("/")
        new_rep_name[-1] = "rmsd_" + new_rep_name[-1]
        new_rep_name = "/".join(new_rep_name)
        with open(new_rep_name, "w") as out:
            out.write(new_report)
    except Exception as e:
        print(f"{e} error in {report}!")

def compute_rmsd_in_serie(pdb_reference, fol_path, resname="GRW", 
                          plot=False, y=None, index_relation=None,
                          report_pref='report_', processors=4):
    equilibration_files = sorted(glob.glob("{}/trajectory_[0-9]*.pdb".format(fol_path)))
    # Compute RMSD
    joblib.Parallel(n_jobs=processors)(joblib.delayed(compute_trajectory_rmsd)
                   (traj, fol_path, pdb_reference, resname, index_relation, report_pref) 
                   for traj in equilibration_files)
    # Store a CSV
    df = read_reports(fol_path, "rmsd_{}".format(report_pref))
    # Overwriting the initial report
    df.to_csv(os.path.join(fol_path, 'rmsd_summary.csv'))
    if plot:
        filename = os.path.join(fol_path, "{}_rmsd.png".format(y))
        plot_results(df, xcolumn="rmsd", ycolumn=y, outfile=filename)


def check_if_multimodels_and_split(pdb_file):
    with open(pdb_file) as pdb:
        content = pdb.readlines()
    pdbs = []
    lines_in_model = []
    for line in content:
        lines_in_model.append(line)
        if "ENDMDL" in line:
            pdbs.append("".join(lines_in_model))
            lines_in_model = []
    return pdbs


def get_report_file_from_trajectory(trajectory_path, report_prefix):
    regex = re.compile(r'\d+')
    sufix = regex.findall(trajectory_path)[-1]
    prefix = trajectory_path.split("/")[0:-1]
    new_prefix = "/".join(prefix)
    path_completed = os.path.join(new_prefix, report_prefix+sufix)
    return path_completed


def plot_results(dataframe, xcolumn, ycolumn, outfile):
    dataframe.plot(xcolumn, ycolumn, kind="scatter")
    plt.xlabel(xcolumn)
    plt.ylabel(ycolumn)
    plt.savefig(outfile)
    plt.show()


def get_snapshot(trajectory_file, model, output_file):
    pdb = check_if_multimodels_and_split(trajectory_file)[model]
    with open(output_file,"w") as pdb2write:
        pdb2write.write(pdb)


def get_min_row(data, column_number):
    minimum = data.loc[data[column_number] == data[column_number].min()]
    return minimum


def get_pdb_by_min_value(report_file, column_number):
    data = pd.read_csv(report_file, sep="\t", header=None)
    minium_row = get_min_row(data, column_number)
    return minium_row[0].values[0], minium_row[1].values[0]


if __name__ == '__main__':
    pdb_reference, fol_path, resname, indx_rel, ycol, plot, \
    report_pref, processors = parse_arguments()
    compute_rmsd_in_serie(pdb_reference, fol_path, resname=resname, 
                          plot=plot, y=ycol, index_relation=indx_rel, 
                          report_pref=report_pref, processors=processors)

