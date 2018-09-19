import sys
import logging
import glob
import prody


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_duplicated_pdbatomnames(pdb_content):
    pdb_atom_names_list = []
    for line in pdb_content:
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
        content = pdb.readlines()
        check_duplicated_pdbatomnames(content)
        for i, line in enumerate(content):
            if line.startswith("HETATM"):
                atom_name = line[12:15]
                if atom_name.strip().startswith("G"):
                    new_atom_name = line[77:78]+atom_name.strip()
                    line_to_list = list(line)
                    line_to_list[12:15] = new_atom_name
                    line_to_list = "".join(line_to_list)
                    content[i] = line_to_list
        check_duplicated_pdbatomnames(content)
        new_pdb = "".join(content)
    with open(pdb_file, "w") as writepdb:
        writepdb.write("{}".format(new_pdb))

