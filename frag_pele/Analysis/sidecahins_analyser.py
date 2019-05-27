import sys
import os
import numpy as np
import prody
import argparse


# Local import
from rmsd_computer import superimpose_backbones
from interaction_detector import pdb2prody
from interaction_detector import select_atom_given_name_type_and_num


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-t", "--type", required=True, choices=['sidechains', 'atom_distances'],
                                help="""Computation type that you want to do. Choose between: sidechains, 
                                atom_distances. """)

    required_named.add_argument("-pdbt", "--pdb_target", required=True,
                                help="""Target pdb file. """)


    parser.add_argument("-pdbr", "--pdb_reference",
                        help="""Reference pdb file.""")

    parser.add_argument("-a", "--area", default=False,
                        help=""""The area (volume) that will be used to select which amino acids are close to the 
                        ligand. THE USAGE OF THIS FLAG HAS NOT BEEN DEBUG SO IS NOT RECOMMENDED TO ACTIVATE IT.""")

    parser.add_argument("-ch", "--chain", default="L",
                        help="Name of the ligand's chain")
    parser.add_argument("-o", "--out", default="report_rmsd.txt",
                        help="Output file name.")
    parser.add_argument("-rf", "--res_file", default=False,
                        help="""Filename that contains which amino acids do you want to select to do the rmsd
                        computations.""")


    args = parser.parse_args()

    return args.type, args.pdb_reference, args.pdb_target, args.area, args.out, args.chain, args.res_file


def read_selecteds_from_file(filepath):
    residues_id_list = []
    with open(filepath, "r") as in_file:
        for aminoacid in in_file.readlines():
            residues_id_list.append(aminoacid)
    return [a.strip("\n") for a in residues_id_list]


def sidechains_rmsd_calculator(pdb_target, pdb_reference, res_file=False, area=False, write2report=False, ligand_chain="L"):
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
    if area:
        print("Selection set of {} Amstrongs".format(area))
        selected_area_target = reference.select("protein and (within {} of chain {})".format(area, ligand_chain))
        unique_residues_target = sorted(set(selected_area_target.getResnums()))
    elif res_file:
        aminoacids_list = read_selecteds_from_file(res_file)
        print("Searching the following amino acids: {}".format(aminoacids_list))
        selected_area_target = reference.select("resnum {}".format(' '.join(aminoacids_list)))
        unique_residues_target = sorted(set(selected_area_target.getResnums()))
    else:
        print("Please, set an input file or a radii to determine which amino acids will be used to compute the RMSD.")
    list_of_results = []
    for residue_target in unique_residues_target:
        res_selected_target = target.select("protein and resnum {} and heavy".format(residue_target))
        res_selected_reference = reference.select("protein and resnum {} and heavy".format(residue_target))
        target_CA = target.select("protein and resnum {} and name CA".format(residue_target))
        reference_CA = reference.select("protein and resnum {} and name CA".format(residue_target))
        try:
            RMSD = prody.calcRMSD(res_selected_reference, res_selected_target)
            distance_bet_CA = prody.calcRMSD(reference_CA, target_CA)
        except:
            print("ERROR because different number of atoms in residue {}".format(residue_target))
            print("ATOMS of the TARGET: {}".format(res_selected_target.getNames()))
            print("ATOMS of the REFERENCE: {}".format(res_selected_reference.getNames()))
        residue_information = (residue_target, res_selected_target.getResnames()[0], RMSD, distance_bet_CA)
        list_of_results.append(residue_information)
        print(residue_information)

    if write2report:
        filename = write2report
        with open(filename, "w") as report:
            for result in list_of_results:
                report.write("{:4d}\t{}\t{:5.3f}\t{:5.3f}\t{:5.3f}\n".format(result[0], result[1], float(result[2]),
                             float(result[3]), (float(result[2]) - float(result[3]))))


def compute_atom_distances(pdb_target, res_file, output_report, chain="L"):
    """
    This function calculate atom-atom distances for ligand and residue atoms. The residue number and atom names
    (for both, ligand and residue) must be specified in a file ('res_file').
    :param pdb_target: input PDB file path
    :param res_file: file with instructions. This file must have n rows with three format: RESNUM LATOMNAME/SLIGAND
    ATOMNAMES/SRESIDUE. If you want to calculate a distance using a center of mass write [ATOM1,ATOM2,ATOMN]

    :param output_report:
    :param chain:
    :return:
    """
    # Load PDB files
    target = pdb2prody(pdb_target)
    ligand = target.select("chain {}".format(chain))
    print(ligand.getNames())
    # Reading instructions from file
    list_of_instructions = read_selecteds_from_file(res_file)
    # Select the input atoms
    report = []
    for line in list_of_instructions:
        resnum, atom_name_ref, atom_name_tar = line.split()
        # If the user wants to select more than one atom he has to put them in a string with this format: [atom1,atomN...]
        if "[" in atom_name_ref or "]" in atom_name_ref:
            print("Multiple atom selection for the ligand")
            atom_string_with_comas = atom_name_ref.strip("[").strip("]")
            atom_list = atom_string_with_comas.split(",")
            atom_string = ' '.join(atom_list)
            atom_ref_selected = ligand.select("name {}".format(atom_string))
            print("Selected atoms: {}".format(atom_ref_selected.getNames()))
        else:
            print("Single atom selection for the ligand")
            atom_ref_selected = ligand.select("name {}".format(atom_name_ref))
            print("Selected atom: {}".format(atom_ref_selected.getNames()))

        if "[" in atom_name_tar or "]" in atom_name_tar:
            print("Multiple atom selection for the system")
            atom_string_with_comas = atom_name_tar.strip("[").strip("]")
            atom_list = atom_string_with_comas.split(",")
            atom_string = ' '.join(atom_list)
            atom_tar_selected = select_atom_given_name_type_and_num(target, resnum, atom_string)
            print("Selected atoms: {}".format(atom_tar_selected.getNames()))
        else:
            print("Single atom selection for the system")
            atom_tar_selected = select_atom_given_name_type_and_num(target, resnum, atom_name_tar)
            print("Selected atom: {}".format(atom_tar_selected.getNames()))
        try:
            number_of_selected_atoms_ref = len(atom_ref_selected.getNames())
        except AttributeError:
            exit("None atoms where selected. Please, check if the selected atoms exists in the ligand in {}".format(pdb_target))
        try:
            number_of_selected_atoms_tar = len(atom_tar_selected.getNames())
        except AttributeError:
            exit("None atoms where selected. Please, check if the selected atoms exists in the residue {} in {}".format(resnum, pdb_target))
        # Now there are four possibilities: len 1 in target and ref, len > 1 in one of both and len >1 in both.
        # If the len is more than 1 we will use the center of mass as a point to compute the distance.
        if number_of_selected_atoms_ref <= 1 and number_of_selected_atoms_tar <= 1:
            distance = prody.calcDistance(atom_tar_selected, atom_ref_selected)
        elif number_of_selected_atoms_ref <= 1 and number_of_selected_atoms_tar > 1:
            center_tar = prody.calcCenter(atom_tar_selected)
            atom_coords = atom_ref_selected.getCoords()[0]
            distance = np.linalg.norm(center_tar - atom_coords)
        elif number_of_selected_atoms_ref > 1 and number_of_selected_atoms_tar <= 1:
            center_ref = prody.calcCenter(atom_ref_selected)
            atom_coords = atom_tar_selected.getCoords()[0]
            distance = np.linalg.norm(atom_coords - center_ref)
        else:
            center_tar = prody.calcCenter(atom_tar_selected)
            center_ref = prody.calcCenter(atom_ref_selected)
            distance = np.linalg.norm(center_tar-center_ref)
        report_line = "{:4} {:10} {:10} {:6.3f}\n".format(resnum, ''.join(atom_ref_selected.getNames()),
                                                          ''.join(atom_tar_selected.getNames()), distance)
        report.append(report_line)

    report_final = ''.join(report)

    if output_report:
        with open(output_report, "w") as report_file:
            report_file.write(report_final)
    print(report_final)
    return report_final


if __name__ == '__main__':
    type, pdb_reference, pdb_target, area, write2report, chain, res_file = parse_arguments()
    if type == "sidechains":
        sidechains_rmsd_calculator(pdb_target, pdb_reference, res_file, area, write2report, chain)
    if type == "atom_distances":
        compute_atom_distances(pdb_target, res_file, write2report, chain)