import prody
import pybel
import logging
import sys

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

def pdb_parser_ligand(pdb_file, ligand_chain="L"):
    """
    :param pdb_file: input PDB file of the complex that we want to get their ligand
    :param ligand_chain: chain where the ligand is placed
    :return: PRODY object with the atoms of the ligand of the input PDB
    """
    pdb = prody.parsePDB(pdb_file)
    ligand = pdb.select("chain {}".format(ligand_chain))
    if not ligand:
        logger.critical("Wrong chain selected!")
    elif ligand.ishetero:
        return ligand
    else:
       logger.critical("The selected chain does not contain heteroatoms!")


def selection_to_pdb(selection):
    """
    :param selection: prody selection
    :return: writes a PDB file containing the selection, named using the residue name of the first atom of the selection.
    """
    selection_pdb = prody.writePDB(selection.getResnames()[0], selection)
    return selection_pdb


def pdb_to_smile(pdb_file):
    """
    :param pdb_file: input pdb file
    :return: string with the translation of the pdb content to smile format
    """
    pdb = pybel.readfile("pdb", pdb_file).__next__()
    smile_and_name = pdb.write("smi")
    logger.info("SMILE of the pdb: {}".format(smile_and_name))
    smile = smile_and_name.split("\t")[0]
    return smile


