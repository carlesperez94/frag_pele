import sys
import logging
import prody


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_duplicated_pdbatomnames(pdb):
    lines = pdb.readlines()
    pdb_atom_names_list = []
    for line in lines:
        if line.startswith("HETATM"):
            pdb_atom_name = line.split()[2]
            pdb_atom_names_list.append(pdb_atom_name)
    set_to_check = set(pdb_atom_names_list)
    list_to_check = sorted(list(set_to_check))
    sorted_list_names = sorted(pdb_atom_names_list)
    if list_to_check != sorted_list_names:
        sys.exit("REPEATED PDB ATOM NAMES IN PDB FILES!!")


def check_and_fix_pdbatomnames(pdb_file):
    with open(pdb_file) as pdb:
        check_duplicated_pdbatomnames(pdb)
