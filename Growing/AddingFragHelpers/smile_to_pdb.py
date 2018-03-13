import pybel
import logging
import sys

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def molecule_to_h_and_3d(smi_object):
    """
    Add H and 3D parameters to a pybel molecule object.
    :param smi_object: pybel object from smile.
    :return: pybel object with hydrogens and 3D parameters.
    """
    smi_object.addh()
    smi_object.make3D()  # Forcefield used by default: mmff94
    return smi_object


def molecule_to_pdb(molecule, filename):
    """
    :param molecule: pybel molecule that we want to transform into PDB file.
    :param filename: filename of the PDB.
    :return: PDB file of the molecule
    """
    output = pybel.Outputfile("pdb", filename)
    output.write(molecule)

