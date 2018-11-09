import sys
import os
import logging
import constants as c
import Helpers.checker as ch
# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def read_instructions_from_file(file):
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
                        fragment_atom = line.split()[(i * 3) + 2]
                        ID = "{}{}{}".format(os.path.splitext(fragment_pdb)[0], core_atom, fragment_atom)
                        task = (fragment_pdb, core_atom, fragment_atom, ID)
                        successive_tasks.append(task)
                    except IndexError:
                        logger.critical("Check that the serie file {} contains: 'PDB_fragment_file'\t'PDB_core_atom_name'\t'PDB_fragment_atom_name' ".format(file))
                list_of_instructions.append(successive_tasks)
    return list_of_instructions


def get_pdb_fragments_and_atoms_from_instructions(list_of_instructions):
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

    
def check_instructions(list_of_instructions, complex, c_chain = "L", f_chain="L"):
    fragments_and_atoms = get_pdb_fragments_and_atoms_from_instructions(list_of_instructions)
    complex = os.path.join(c.PRE_WORKING_DIR, complex)
    for fragment, atom_core, atom_fr in fragments_and_atoms:
        fragment = os.path.join(c.PRE_WORKING_DIR, fragment)
        ch.check_if_atom_exists_in_ligand(fragment, atom_fr, f_chain)
        ch.check_if_atom_exists_in_ligand(complex, atom_core, c_chain)
        with open(fragment) as content:
            fr_content = content.readlines()
        ch.check_duplicated_pdbatomnames(fr_content)



