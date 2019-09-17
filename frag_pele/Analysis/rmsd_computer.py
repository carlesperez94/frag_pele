import sys
import os
import re
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
    required_named.add_argument("-p", "--path", required=True,
                                help="""Path to folder to be analyzed. """)

    parser.add_argument("-r", "--resname", default="GRW",
                        help="Resame of the ligand.")
    parser.add_argument("-c", "--csv", default="summary/summary_report_1.csv",
                        help="Relative path to csv file")


    # Plot arguments
    parser.add_argument("-y", "--ycol", default="RMSD",
                        help="Column name of the report file that you want to use as X axis in the plot.")
    parser.add_argument("-op", "--plot", default=False, action="store_true",
                        help="Flag to plot the results.")

    args = parser.parse_args()

    return args.pdb_reference, args.path, args.resname, args.csv, args.ycol, args.plot


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


def ligands_rmsd_calculator(pdb_target, pdb_reference, resname="GRW"):
    """

    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param write2report: if True export results in a file
    :param write2pdb: pdb file with the result of the superposition between "pdb_target" and "pdb_reference"
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD of the ligand
    """
    # Reparation of PDB files if it is needed
    repair_pdbs(pdb_target)
    # Load the data to mdtraj
    target = mdtraj.load(pdb_target)
    reference = mdtraj.load(pdb_reference)
    # Get the indexes of heavy atoms of the ligands
    ligand = reference.topology.select('resname "{}"'.format(resname))
    rmsd = mdtraj.rmsd(target, reference=reference, atom_indices=ligand)
    return rmsd


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


def compute_rmsd_in_serie(pdb_reference, path, resname="GRW", pattern_to_csv="summary/summary_report_1.csv", plot=False,
                          y=None):
    equilibration_files = sorted(glob.glob("{}/trajectory_[0-9]*.pdb".format(path)))
    csv = os.path.join(path, pattern_to_csv)
    df = pd.read_csv(csv)
    for file in equilibration_files:
        result = ligands_rmsd_calculator(file, pdb_reference, resname=resname)
        processor = re.findall(r'\d+', file.split("/")[-1])[0]
        try:
            indexes = search_indexes_with_certain_value_in_column(df, int(processor), "Processor")
            for rmsd, index in zip(result, indexes):
                add_data_to_dataframe(dataframe=df, value_to_add=rmsd, row_to_add=index, column_name="rmsd")
        except Exception as e:
            print(e)
    # Overwriting the initial report
    df.to_csv(csv)
    if plot:
        filename = os.path.join(path, "{}_rmsd.png".format(y))
        plot_results(df, xcolumn="rmsd", ycolumn=y, outfile=filename)


def add_data_to_dataframe(dataframe, value_to_add, row_to_add, column_name="new_column"):
    dataframe.loc[row_to_add, column_name] = value_to_add


def search_indexes_with_certain_value_in_column(dataframe, value_to_seach, column_to_seach):
    idx = dataframe.index[dataframe[column_to_seach] == value_to_seach]
    return idx


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
    pdb_reference, path, resname, csv, ycol, plot = parse_arguments()
    compute_rmsd_in_serie(pdb_reference, path, pattern_to_csv=csv, resname=resname, plot=plot, y=ycol)

