import sys
import re
import os
import logging
# Local imports 
import frag_pele.constants as c
import frag_pele.Helpers.checker as ch
# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def read_instructions_from_file(file):
    """
    It reads an "instruction file". This file contains the information of all the growing's that the user wants to
    perform, each one separated by newlines. Each instruction must have at least 3 columns, separated by tabulations:
    (1) name of the fragment's PDB file,
    (2) atom name of the core (if the user wants to add the H it must be separated by "-", example: C4-H2),
    (3) atom name of the fragment.
    If the line has more than 3 columns, it will be considered as successive growing's simulations, following the same
    structure previously mentioned but several times.
    :param file: filename of the instructions file.
    :return: list with all instructions processed.
    """
    list_of_instructions = []
    with open(file) as sf:
        instructions = sf.readlines()
        for line in instructions:
            if 0 < len(line.split()) <= 3:
                # Read from file the required information
                try:
                    fragment_pdb = line.split()[0]
                    core_atom = line.split()[1]
                    fragment_atom = line.split()[2]
                    ID = "{}{}{}".format(os.path.splitext(fragment_pdb)[0], core_atom, fragment_atom)
                    task = (fragment_pdb, core_atom, fragment_atom, ID)
                    list_of_instructions.append(task)
                except IndexError:
                    logger.critical("Check that the serie file {} contains: 'PDB_fragment_file'\t'PDB_core_atom_name'\t'PDB_fragment_atom_name' ".format(file))
            elif len(line.split()) > 3:
                growing_counter = len(line.split()) / 3
                successive_tasks = []
                for i in range(int(growing_counter)):
                    try:
                        fragment_pdb = line.split()[i * 3]
                        core_atom = line.split()[(i * 3) + 1]
                        if "*" in core_atom:
                            # Get fragment number and remove part of string
                            fragment_number = re.findall(r'[*]\d[*]', core_atom)[0].strip("*")
                            core_atom = core_atom.replace("*{}*".format(fragment_number), "")
                        else:
                            fragment_number = None
                        fragment_atom = line.split()[(i * 3) + 2]
                        ID = "{}{}{}".format(os.path.splitext(fragment_pdb)[0], core_atom, fragment_atom)
                        task = (fragment_pdb, core_atom, fragment_atom, ID, fragment_number)
                        successive_tasks.append(task)
                    except IndexError:
                        logger.critical("Check that the serie file {} contains: 'PDB_fragment_file'\t'PDB_core_atom_name'\t'PDB_fragment_atom_name' ".format(file))
                list_of_instructions.append(successive_tasks)
    return list_of_instructions


def get_pdb_fragments_and_atoms_from_instructions(list_of_instructions):
    """
    Given a list with the instructions processed it returns a list with containing the following elements:
    [fragment_pdb_file, pdb_atom_name_of_the_core (and H if required), pdb_atom_name_of_the_fragment (and H if required)]
    :param list_of_instructions: list with the instructions read from the instructions file. list
    :return: list with all fragments and atoms that will be used in the future. list
    """
    fragments_pdbs_and_atoms = []
    for instructions in list_of_instructions:
        if type(instructions) == list:
            for single_instruction in instructions:
                fragment_pdb, atom_core, atom_fragment = single_instruction[0], single_instruction[1], single_instruction[2]
                fragments_pdbs_and_atoms.append((fragment_pdb, atom_core, atom_fragment))
        else:
            fragment_pdb, atom_core, atom_fragment = instructions[0], instructions[1], instructions[2]
            fragments_pdbs_and_atoms.append((fragment_pdb, atom_core, atom_fragment))

    return fragments_pdbs_and_atoms

    
def check_instructions(instructions, complex_pdb, c_chain = "L", f_chain="L"):
    """
    It checks if the selected atoms exists in their correspondent PDB file and also checks if there are repeated
    PDB-atom-names in the PDB file.
    :param list_of_instructions: list with the instructions read from the instructions file. list
    :param complex: PDB file with the complex that contains the core
    :return: if something is wrong it raises an exception.
    """
    fragments_and_atoms = get_pdb_fragments_and_atoms_from_instructions([instructions])
    for fragment, atom_core, atom_fr in fragments_and_atoms:
        atoms_if_bond = extract_hydrogens_from_instructions([fragment, atom_core, atom_fr])
        if atoms_if_bond:
            ch.check_if_atom_exists_in_ligand(complex_pdb, atoms_if_bond[0], c_chain)
            ch.check_if_atom_exists_in_ligand(complex_pdb, atoms_if_bond[1], c_chain)
            ch.check_if_atom_exists_in_ligand(fragment, atoms_if_bond[2], f_chain)
            ch.check_if_atom_exists_in_ligand(fragment, atoms_if_bond[3], f_chain)
        else:
            ch.check_if_atom_exists_in_ligand(fragment, atom_fr, f_chain)
            ch.check_if_atom_exists_in_ligand(complex_pdb, atom_core, c_chain)
        with open(fragment) as content:
            fr_content = content.readlines()
        ch.check_duplicated_pdbatomnames(fr_content)


def extract_hydrogens_from_instructions(instruction):
    """
    If the core or the fragment atom contains a "-" means that the user is selecting an specific H to be bonded with
    the heavy atom. For this reason, this detects if this option has been selected by the user and extract the PDB-atom-name
    of each element.
    :param instruction: list that follow this structure: [fragment_pdb_file, pdb_atom_name_of_the_core
    (and H if required), pdb_atom_name_of_the_fragment (and H if required)]
    :return: if the "-" is found in the pdb_atom_name_of_the_core or pdb_atom_name_of_the_fragment it returns a list with
    PDB-atom-names split following this order: heavy_atom_core, hydrogen_core, heavy_atom_fragment, hydrogen_fragment.
    Otherwise it returns False.
    """
    if "-" in instruction[1] or "-" in instruction[2]:
        try:
            heavy_core = instruction[1].split("-")[0]
            hydrogen_core = instruction[1].split("-")[1]
            heavy_fragment = instruction[2].split("-")[0]
            hydrogen_fragment = instruction[2].split("-")[1]
            return heavy_core, hydrogen_core, heavy_fragment, hydrogen_fragment
        except IndexError:
            raise IndexError("To use steriochemistry the hydrogen from the core AND fragmet must be specified via controlfile")
    else:
        return False

