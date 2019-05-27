import sys
import os
import glob
import pandas as pd
import argparse
import prody


def parse_arguments():
    """
            Parse user arguments

            Output: list with all the user arguments
        """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-tpdb", "--tar_pdb", required=True,
                                help="""Target PDB file.""")
    required_named.add_argument("-rpdb", "--ref_pdb", required=True,
                                help="""Reference PDB file.""")

    parser.add_argument("-ch", "--chain", default="L",
                        help="Ligand chain.")
    parser.add_argument("-d", "--dist", default=4,
                        help="Interaction distance threshold.")
    parser.add_argument("-w2r", "--write2report", default=False, action='store_true',
                        help="Write a report file with the interactions.")
    parser.add_argument("-ah", "--addh", default=False, action='store_true',
                        help="If set it takes into account ligand H's in the interaction.")
    parser.add_argument("-o", "--outfile", default=False,
                        help="If set it takes into account ligand H's in the interaction.")

    args = parser.parse_args()

    return args.tar_pdb, args.ref_pdb, args.chain, args.dist, args.write2report, args.addh, args.outfile


def pdb2prody(pdb_file):
    pdb_in_prody = prody.parsePDB(pdb_file)
    return pdb_in_prody


def get_atoms_of_chain(prody_object, ligand_chain, add_H):
    if not add_H:
        ligand = prody_object.select("heavy and chain {}".format(ligand_chain))
    else:
        ligand = prody_object.select("chain {}".format(ligand_chain))
    return ligand.getNames()


def get_elements_of_chain(prody_object, ligand_chain, add_H):
    if not add_H:
        ligand = prody_object.select("heavy and chain {}".format(ligand_chain))
    else:
        ligand = prody_object.select("chain {}".format(ligand_chain))
    return ligand.getElements()


def get_resname_of_chain(prody_object, ligand_chain):
    ligand = prody_object.select("chain {}".format(ligand_chain))
    return ligand.getResnames()[0]


def get_resnum_of_chain(prody_object, ligand_chain):
    ligand = prody_object.select("chain {}".format(ligand_chain))
    return ligand.getResid()[0]


def get_protein_atoms_at_interaction_distance(atom, pdb_prody, distance=4):
    atoms_close = pdb_prody.select("(protein or water) and within {} of name {}".format(distance, atom))
    if atoms_close == None:
        print("interaction not found")
        return [""]
    else:
        return atoms_close.getNames()


def get_protein_residues_at_interaction_distance(atom, pdb_prody, distance=4):
    atoms_close = pdb_prody.select("(protein or water) and within {} of name {}".format(distance, atom))
    residues = []
    if atoms_close == None:
        print("interaction not found")
        return residues
    else:
        resnums_close = atoms_close.getResnums()
        resnames_close = atoms_close.getResnames()
        for num, resname in zip(resnums_close, resnames_close):
            residues.append("{}{}".format(num, resname))
        return residues


def get_protein_residues_and_atoms_at_interaction_distance(atom, pdb_prody, distance=4):
    atoms_list = get_protein_atoms_at_interaction_distance(atom, pdb_prody, distance)
    residues_list = get_protein_residues_at_interaction_distance(atom, pdb_prody, distance)
    result = []
    for atom_at_dist, residue in zip(atoms_list, residues_list):
        result.append((atom, atom_at_dist, residue))
    return result


def write_report(atom_map, pdb_file):
    lines = []
    for contact in atom_map:
        line = "{:4} {:2} {:4} {}\n".format(contact[0], contact[1], contact[2], contact[3])
        lines.append(line)
    with open("liginteraction_report_{}.txt".format(os.path.splitext(pdb_file)[0]), "w") as report_file:
        report_file.writelines(lines)


def get_interaction_map(pdb_input, ligand_chain, distance, write2report, addh=False):
    pdb = pdb2prody(pdb_input)
    ligand_atoms = get_atoms_of_chain(pdb, ligand_chain, addh)
    elements = get_elements_of_chain(pdb, ligand_chain, addh)
    atom_map = []
    for element, atom in zip(elements, ligand_atoms):
        if not addh:
            if element == "H":
                pass
            else:
                for line in get_protein_residues_and_atoms_at_interaction_distance(atom, pdb, distance):
                    ligatomname, aaatomname, aatype = line
                    print("{:4} {:2} {:4} {}".format(ligatomname, element, aaatomname, aatype))
                    atom_map.append((ligatomname, element, aaatomname, aatype))
        else:
            for line in get_protein_residues_and_atoms_at_interaction_distance(atom, pdb, distance):
                ligatomname, aaatomname, aatype = line
                print("{:4} {:2} {:4} {}".format(ligatomname, element, aaatomname, aatype))
                atom_map.append((ligatomname, element, aaatomname, aatype))
    if write2report:
        write_report(atom_map, pdb_input)
    return atom_map


def get_aminoacids_at_interaction_dist(list_of_interacting_atoms):
    interacting_residues_list = []
    for interacting_atom in list_of_interacting_atoms:
        ligatom = interacting_atom[0]
        ligelement = interacting_atom[1]
        residue = interacting_atom[3]
        interacting_residues_list.append((ligatom, ligelement, residue))
    interacting_residues_list_unique = set(interacting_residues_list)
    return interacting_residues_list_unique


def add_origen_of_differents(list_of_target, list_of_differents):
    new_list = []
    for interaction in list_of_differents:
        if interaction in list_of_target:
            interaction = list(interaction)
            interaction.append("T")
            new_list.append(interaction)
        else:
            interaction = list(interaction)
            interaction.append("R")
            new_list.append(interaction)
    return new_list


def write2report_of_differences(differences_list, output_name):
    output_lines = []
    for interaction in sorted(differences_list):
        line = "{:4} {:2} {:7} {}\n".format(interaction[0], interaction[1], interaction[2], interaction[3])
        output_lines.append(line)
    with open(output_name, "w") as out_file:
        for line in output_lines:
            out_file.writelines(line)


def select_atom_given_name_type_and_num(prody_object, resnum, atom_name):
    selection = prody_object.select("resnum {} and name {}".format(resnum, atom_name))
    return selection


def read_atoms_to_get_distance(filename):
    atoms_to_compute = []
    with open(filename, "r") as instructions:
        lines = instructions.readlines()
        for line in lines:
            atom_name_ligand = line[0]
            atom_name_aa = line[1]
            resnum = line[2]
            resname = line[3]
            atoms_to_compute.append((atom_name_ligand, atom_name_aa, resnum, resname))
    return atoms_to_compute


def compute_distances_given_instruction(pdb_object, instruction, chain = "L"):
    ligand_resname = get_resname_of_chain(pdb_object, chain)
    ligand_resnum = get_resnum_of_chain(pdb_object, chain)
    atom_name_ligand, atom_name_aa, resnum_aa, resname_aa = instruction
    atom_of_ligand = select_atom_given_name_type_and_num(pdb_object, ligand_resnum, atom_name_ligand)
    atom_of_protein = select_atom_given_name_type_and_num(pdb_object, ligand_resnum, atom_name_ligand)
    distance = prody.calcDistance(atom_of_ligand, atom_of_protein)
    return atom_name_ligand, atom_name_aa, resnum_aa, resname_aa, distance


def main(pdb_target, pdb_reference, ligand_chain, distance, write2report, distances_file = False, addh=False, output_name=False):
    target_interactions = set(get_interaction_map(pdb_target, ligand_chain, distance, write2report, addh))
    reference_intertactions = set(get_interaction_map(pdb_reference, ligand_chain, distance, write2report, addh))
    differences = list(target_interactions - reference_intertactions)
    amount_t_interactions = len(target_interactions)
    amount_r_interactions = len(reference_intertactions)
    print("DIFFERENT ATOMS CONTACTS")
    differences = sorted(differences)
    differences = add_origen_of_differents(target_interactions, differences)
    differences_in_num = len(differences)
    for diff in differences:
        print("{:4} {:2} {:4} {:7} {}".format(diff[0], diff[1], diff[2], diff[3], diff[4]))
    print("Initial atoms at interaction distance: {}\nFinal atoms at interaction distance: {}\nDifferences: {}".format(
                                                                                  amount_r_interactions,
                                                                                  amount_t_interactions,
                                                                                  differences_in_num))
    target_aa_close = get_aminoacids_at_interaction_dist(target_interactions)
    ref_aa_close = get_aminoacids_at_interaction_dist(reference_intertactions)
    differences_aa = list(set(target_aa_close) - set(ref_aa_close))
    differences_aa = add_origen_of_differents(target_aa_close, differences_aa)
    print("DIFFERENT AA CONTACTS")
    for interaction in sorted(differences_aa):
        print("{:4} {:2} {:7} {}".format(interaction[0], interaction[1], interaction[2], interaction[3]))
    amount_aa_t_interactions = len(target_aa_close)
    amount_aa_r_interactions = len(ref_aa_close)
    differences_aa_in_num = len(differences_aa)
    print("Initial amino acids at interaction distance: {}\nFinal amino acids at interaction distance: {}\nDifferences: {}".format(
        amount_aa_r_interactions,
        amount_aa_t_interactions,
        differences_aa_in_num))
    if distances_file:
        for instruction in read_atoms_to_get_distance(distances_file):
            compute_distances_given_instruction(pdb_target, instruction, chain)
    if output_name:
        write2report_of_differences(differences_aa, output_name)


if __name__ == '__main__':
    pdb_target, pdb_reference, chain, dist, write2report, addh, outfile = parse_arguments()
    main(pdb_target, pdb_reference, chain, dist, write2report, addh, outfile)
