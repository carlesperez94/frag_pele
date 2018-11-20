import sys
import os
import re
import ntpath
import prody
import glob
import argparse
import pandas as pn
import matplotlib.pyplot as plt
from analyser import get_min_value


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

    parser.add_argument("-pdbt", "--pdb_target",
                                help="""Target pdb file. """)

    parser.add_argument("-wp", "--w2pdb", default=False,
                        help="If true export a pdb file with the alignment result.")
    parser.add_argument("-ch", "--chain", default="L",
                        help="Name of the ligand's chain")
    parser.add_argument("-o", "--out", default="report_rmsd.txt",
                        help="Output file name.")
    parser.add_argument("-c", "--column_to_add", default="Binding Energy",
                        help="Column name of the report file that you want to add to the report")
    # Plot arguments
    parser.add_argument("-x", "--xcol", default="RMSD",
                        help="Column name of the report file that you want to use as X axis in the plot.")
    parser.add_argument("-op", "--out_plot", default="plot.png",
                        help="Filename of the output graph.")

    args = parser.parse_args()

    return args.pdb_reference, args.path, args.pdb_target, args.w2pdb, args.chain, args.out, args.column_to_add, \
           args.xcol, args.out_plot


def superimpose_backbones(pdb_target, pdb_reference):
    # Loading PDB files
    target = prody.parsePDB(pdb_target)
    reference = prody.parsePDB(pdb_reference)
    # Selection of the backbone
    target_backbone = target.select("backbone")
    reference_backbone = reference.select("backbone")
    # Superimpose the target to the reference backbone
    transformation = prody.calcTransformation(target_backbone, reference_backbone)
    prody.applyTransformation(transformation, target)
    return target, reference


def ligands_rmsd_calculator(pdb_target, pdb_reference, column_to_add=None, write2pdb=False, ligand_chain="L", ):
    """

    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param write2report: if True export results in a file
    :param write2pdb: pdb file with the result of the superposition between "pdb_target" and "pdb_reference"
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD of the ligand
    """
    repair_pdbs(pdb_target)
    target, reference = superimpose_backbones(pdb_target, pdb_reference)
    target_ligand = target.select("chain {} and heavy".format(ligand_chain))
    reference_ligand = reference.select("chain {} and heavy".format(ligand_chain))

    # Creating temporary files to perform RMSD computations in futher steps
    prody.writePDB("structure_superposed_target_tmp.pdb", target_ligand)
    prody.writePDB("structure_superposed_reference_tmp.pdb", reference_ligand)
    sort_atoms_by_name("structure_superposed_reference_tmp.pdb")
    results = []
    for i, pdb in enumerate(check_if_multimodels_and_split("structure_superposed_target_tmp.pdb")):
        with open("structure_superposed_target_tmp.pdb", "w") as write_pdb:
            write_pdb.write(pdb)
        sort_atoms_by_name("structure_superposed_target_tmp.pdb")
        target_sorted = prody.parsePDB("structure_superposed_target_tmp.pdb")
        reference_sorted = prody.parsePDB("structure_superposed_reference_tmp.pdb")
        RMSD = prody.calcRMSD(reference_sorted, target_sorted)
        print(RMSD)
        if column_to_add:
            column = get_report_column_from_trajectory(pdb_target, "report_", column_to_add, i)
            result_sumary = "{}\t{}\t{:2.3f}\t{}".format(pdb_target, i, RMSD, column)
        else:
            result_sumary = "{}\t{}\t{:2.3f}\t ".format(pdb_target, i, RMSD)
        results.append(result_sumary)

    if write2pdb:
        prody.writePDB("structure_superposed.pdb", target+reference)
    # Deleting temporary files
    os.remove("structure_superposed_target_tmp.pdb")
    os.remove("structure_superposed_reference_tmp.pdb")
    return results


def sidechains_rmsd_calculator(pdb_target, pdb_reference, radii, path, write2report=True, ligand_chain="L"):
    """
    :param pdb_target: problem pdb file
    :param pdb_reference: reference pdb file
    :param radii: area that we want to select around the ligand
    :param path: output path
    :param write2report: if true extract a report file
    :param ligand_chain: name of the chain of the ligand
    :return: superpose the backbone of the pdb_target to the pdb_reference and computes the RMSD for each side
    chain in the selection area
    """
    target, reference = superimpose_backbones(pdb_target, pdb_reference)
    selected_area_target = reference.select("protein and (within {} of chain {})".format(radii, ligand_chain))
    unique_residues_target = set(selected_area_target.getResnums())
    list_of_results = []
    for residue_target in unique_residues_target:
        res_selected_target = target.select("protein and resnum {}".format(residue_target))
        res_selected_reference = reference.select("protein and resnum {}".format(residue_target))
        target_CA = target.select("protein and resnum {} and name CA".format(residue_target))
        reference_CA = reference.select("protein and resnum {} and name CA".format(residue_target))
        RMSD = prody.calcRMSD(res_selected_reference, res_selected_target)
        distance_bet_CA = prody.calcRMSD(reference_CA, target_CA)
        residue_information = (residue_target, res_selected_target.getResnames()[0], RMSD, distance_bet_CA)
        list_of_results.append(residue_information)
    list_of_results = sorted(list_of_results)
    pdb_reference_name = os.path.splitext(ntpath.basename(pdb_reference))[0]
    pdb_target_name = os.path.splitext(ntpath.basename(pdb_target))[0]
    if write2report:
        filename = "{}_to_{}_sidech.csv".format(pdb_target_name, pdb_reference_name)
        filename = os.path.join(path, filename)
        with open(filename, "w") as report:
            for result in list_of_results:
                report.write("{:4d}\t{}\t{:5.3f}\t{:5.3f}\t{:5.3f}\n".format(result[0], result[1], float(result[2]),
                             float(result[3]), (float(result[2]) - float(result[3]))))


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


def compute_rmsd_in_serie(pdb_reference, path, column_to_add, ligand_chain="L", output_file="report_rmsd.txt"):
    equilibration_files = glob.glob("{}/equilibration*/trajectory*.pdb".format(path))
    if os.path.exists(os.path.join(path, output_file)):
        os.remove(os.path.join(path, output_file))
    for file in equilibration_files:
        result = ligands_rmsd_calculator(file, pdb_reference, column_to_add, write2pdb=False, ligand_chain=ligand_chain)
        for line in result:
            if os.path.exists(os.path.join(path, output_file)):
                with open(os.path.join(path, output_file), "a") as writing_file:
                    writing_file.write("{}\n".format(line))
            else:
                with open(os.path.join(path, output_file), "w") as writing_file:
                    writing_file.write("{}\n".format(line))


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


def get_report_column_from_trajectory(trajectory_path, report_prefix, column, model):
    report_path = get_report_file_from_trajectory(trajectory_path, report_prefix)
    report_data = pn.read_csv(report_path, sep="    ")
    return report_data[column][model]


def plot_results(report_file, x, y, xcolumn, ycolumn, output_name):
    data = pn.read_csv(report_file, sep="\t")
    data.plot(x, y, kind="scatter")
    plt.xlabel(xcolumn)
    plt.ylabel(ycolumn)
    print(output_name)
    plt.savefig(output_name)
    plt.show()


def get_snapshot(trajectory_file, model, output_file):
    pdb = check_if_multimodels_and_split(trajectory_file)[model]
    with open(output_file,"w") as pdb2write:
        pdb2write.write(pdb)


def get_min_row(data, column_number):
    minimum = data.loc[data[column_number] == data[column_number].min()]
    return minimum


def get_pdb_by_min_value(report_file, column_number):
    data = pn.read_csv(report_file, sep="\t", header=None)
    minium_row = get_min_row(data, column_number)
    return minium_row[0].values[0], minium_row[1].values[0]


if __name__ == '__main__':
    pdb_reference, path, pdb_target, w2pdb, chain, output_file, column_to_add, xcol, out_plot = parse_arguments()
    compute_rmsd_in_serie(pdb_reference, path, column_to_add, ligand_chain="L", output_file=output_file)
    get_snapshot(get_pdb_by_min_value(os.path.join(path, output_file), 2)[0],
                 get_pdb_by_min_value(os.path.join(path, output_file), 2)[1],
                 os.path.join(path, "minimum_rmsd.pdb"))
    get_snapshot(get_pdb_by_min_value(os.path.join(path, output_file), 3)[0],
                 get_pdb_by_min_value(os.path.join(path, output_file), 3)[1],
                 os.path.join(path, "minimum_{}.pdb".format(column_to_add.replace(" ", ""))))
    plot_results(os.path.join(path, output_file), 2, 3, xcol, column_to_add, os.path.join(path, out_plot))

