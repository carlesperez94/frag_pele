import matplotlib.pyplot as plt
import numpy as np
import argparse


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Script to perform plots. You can choose between different
    type of plots (depending on your input file): 'boxplot_single' if your input is the result of the analysis
    between two structures, for which you have computed the RMSD and the CA distance between them; 'boxplot_differences'
    if you want to show the comparison of RMSD and CA distance in two different files; and 'boxplot_atom_distances'
    if you have calculated atom-atom distances in three structures in order to show the differences.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-t", "--type", required=True, choices=['boxplot_single','boxplot_differences',
                                                                        'boxplot_atom_distances'],
                                help="""Plot type that you want to do. Choose between: boxplot_differences, 
                                boxplot_single. """)
    required_named.add_argument("-i", "--input", required=True,
                                help="""Path of the input file.""")

    parser.add_argument("-i2", "--input_2", default=False,
                        help=""""Path of the input file 2, only set if you want to do comparisons.""")

    parser.add_argument("-i3", "--input_3", default=False,
                        help=""""Path of the input file 2, only set if you want to do comparisons.""")

    args = parser.parse_args()

    return args.type, args.input, args.input_2, args.input_3


def plot_aminoacid_single(report_file1):
    """
    Given a single report file the function create a boxplot with two bars for n groups.
    :param report_file1: report file path
    :return: boxplot
    """
    # Reading report file 1
    with open(report_file1) as results:
        lines_file1 = results.readlines()
    n_groups_1 = len(lines_file1)
    residue_and_type_1 = []
    total_rmsd_1 = []
    distances_cas_1 = []
    for line in lines_file1:
        resid, type, rmsd_total, distance_ca, diff = line.split("\t")
        residue_and_type_1.append("{}{}".format(resid, type))
        total_rmsd_1.append(float(rmsd_total))
        distances_cas_1.append(float(distance_ca))

    fig, ax = plt.subplots()

    index = np.arange(n_groups_1)
    bar_width = 0.35

    opacity = 0.5

    rects1 = ax.bar(index, total_rmsd_1, bar_width,
                    alpha=opacity, color='r',
                    label='RMSD')

    rects2 = ax.bar(index + bar_width, distances_cas_1, bar_width,
                    alpha=opacity, color='b',
                    label=r'C$\alpha$ distance')

    ax.set_xlabel('Residues', fontsize=12)
    ax.set_ylabel(r'$\AA$',fontsize=12)
    ax.set_title('Comparison of sidechains between Initial and Final structures', fontsize=14)
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(residue_and_type_1, fontsize=10)
    ax.legend()

    fig.tight_layout()
    plt.show()


def plot_aminoacid_differences(report_file1, report_file2):
    """
    Given two report files the function create a boxplot with four bars for n groups (two bars for each file in order
    to compare them.).
    :param report_file1: report file path
    :param report_file2: report file path. This report requires the same labels and metrics than in the first one.
    :return: boxplot
    """
    # Reading report file 1
    with open(report_file1) as results:
        lines_file1 = results.readlines()
    n_groups_1 = len(lines_file1)
    residue_and_type_1 = []
    total_rmsd_1 = []
    distances_cas_1 = []
    for line in lines_file1:
        resid, type, rmsd_total, distance_ca, diff = line.split("\t")
        residue_and_type_1.append("{}{}".format(resid, type))
        total_rmsd_1.append(float(rmsd_total))
        distances_cas_1.append(float(distance_ca))

    # Reading report file 2
    with open(report_file2) as results:
        lines_file2 = results.readlines()
    n_groups_2 = len(lines_file2)
    residue_and_type_2 = []

    total_rmsd_2 = []
    distances_cas_2 = []

    for line in lines_file2:
        resid, type, rmsd_total, distance_ca, diff = line.split("\t")
        residue_and_type_2.append("{}{}".format(resid, type))
        total_rmsd_2.append(float(rmsd_total))
        distances_cas_2.append(float(distance_ca))

    # Simple check to ensure equal groups
    if not n_groups_1 == n_groups_2:
        exit("Different group lengths")
    if not residue_and_type_1 == residue_and_type_2:
        exit("The residues differ from the file 1 and file 2")

    fig, ax = plt.subplots()

    index = np.arange(n_groups_1)
    bar_width = 0.20

    opacity = 0.7

    rects1 = ax.bar(index, total_rmsd_1, bar_width,
                    alpha=opacity, color='b',
                    label='total RMSD initial')

    rects2 = ax.bar(index + bar_width, distances_cas_1, bar_width,
                    alpha=opacity, color='g',
                    label=r'C$\alpha$ distance initial')
    rects3 = ax.bar(index + bar_width*2, total_rmsd_2, bar_width,
                    alpha=opacity, color='r',
                    label='total RMSD final')

    rects4 = ax.bar(index + bar_width*3, distances_cas_2, bar_width,
                    alpha=opacity, color='y',
                    label=r'C$\alpha$ distance final')

    ax.set_xlabel('Residues', fontsize=12)
    ax.set_ylabel(r'$\AA$',fontsize=12)
    ax.set_title(r'RMSDs and C$\alpha$ distances for amino acid', fontsize=14)
    ax.set_xticks(index + bar_width / 4)
    ax.set_xticklabels(residue_and_type_1, fontsize=10)
    ax.legend()

    fig.tight_layout()
    plt.show()


def plot_atom_distances(input_initial, input_final, input_cristal):
    """
    Given three report files, the function create a boxplot with three bars for n groups (one bar for each file in order
    to compare them.).
    :param input_initial: report file (with distances) path of the initial structure
    :param input_final: report file (with distances) path of the final structure (model)
    :param input_cristal: report file (with distances) path of the crystal
    :return: boxplot
    """
    # Reading report file 1
    with open(input_initial) as results:
        lines_file1 = results.readlines()
    n_groups_initial = len(lines_file1)
    atoms_initial = []
    distances_initial = []
    for line in lines_file1:
        resid, atoms1, atoms2, distance = line.split()
        atoms_initial.append("{}-{}{}".format(atoms1, resid, atoms2))
        distances_initial.append(float(distance))

    # Reading report file 2
    with open(input_final) as results:
        lines_file2 = results.readlines()

    n_groups_final = len(lines_file2)
    atoms_final = []
    distances_final = []
    for line in lines_file2:
        resid, atoms1, atoms2, distance = line.split()
        atoms_final.append("{}-{}{}".format(atoms1, resid, atoms2))
        distances_final.append(float(distance))

    # Reading report file 3
    with open(input_cristal) as results:
        lines_file3 = results.readlines()

    n_groups_cristal = len(lines_file3)
    atoms_cristal = []
    distances_cristal = []
    for line in lines_file3:
        resid, atoms1, atoms2, distance = line.split()
        atoms_cristal.append("{}-{}{}".format(atoms1, resid, atoms2))
        distances_cristal.append(float(distance))
    # Simple check to ensure equal groups

    if not n_groups_initial == n_groups_final == n_groups_cristal:
        print(n_groups_initial, n_groups_final, n_groups_cristal)
        exit("Different group lengths")
    if not atoms_initial == atoms_final == atoms_cristal:
        exit("The atoms differ between files")

    fig, ax = plt.subplots()

    index = np.arange(n_groups_initial)
    bar_width = 0.20

    opacity = 0.5

    rects1 = ax.bar(index, distances_initial, bar_width,
                    alpha=opacity, color='r',
                    label='Initial')

    rects2 = ax.bar(index + bar_width, distances_final, bar_width,
                    alpha=opacity, color='b',
                    label='Model')
    rects3 = ax.bar(index + bar_width*2, distances_cristal, bar_width,
                    alpha=opacity, color='g',
                    label='Crystal')

    ax.set_xlabel('Interaction', fontsize=12)
    ax.set_ylabel(r'Distance ($\AA$)', fontsize=12)
    ax.set_title('Comparison of distances in EGFR', fontsize=14)
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(atoms_initial, fontsize=10)
    ax.legend()

    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    type, input, input_file2, input_file3= parse_arguments()
    if type == "boxplot_single":
        plot_aminoacid_single(input)
    if type == "boxplot_differences":
        plot_aminoacid_differences(input, input_file2)
    if type == "boxplot_atom_distances":
        plot_atom_distances(input, input_file2, input_file3)


