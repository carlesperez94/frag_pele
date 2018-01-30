import sys
import re
import os
import string
from shutil import copyfile
import logging


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def template_reader(template_name, path_to_template="DataLocal/Templates/OPLS2005/HeteroAtoms/"):
    """
    This function reads the content of a PELE's template and return it
    """

    with open(os.path.join(path_to_template, template_name), "r") as template_file:
        template_content = template_file.read()
        if not template_content:
            logger.critical("Template file {} is empty!".format(template_name))

    return template_content


def section_selector(template, pattern_1, pattern_2):
    """
    From a template string, this function return a section between two patterns.
    :param template: input template (string).
    :param pattern_1: pattern which sets the begining of the section.
    :param pattern_2: pattern which sets the end of the section.
    :return: string with the content of the section.
    """

    section_selected = re.search("{}\n(.*?){}".format(pattern_1, pattern_2), template, re.DOTALL)

    return section_selected.group(1)


def atoms_selector(template):
    """
    Given a template, it returns a dictionary with the atoms found.
    :param template: input template (string)
    :return: dictionary {"PDB atom name":"index"}
    """
    ROW_PATTERN = "\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)"

    # Select the section of the templates where we have the atoms defined
    atoms_section = section_selector(template, "\*", "NBON")
    # Obtain all rows (list of lists)
    rows = re.findall(ROW_PATTERN, atoms_section)
    # Get the atom name from all rows and insert it in a dictionary
    atoms = {}
    for row in rows:
        atoms[row[4]] = row[0]
    return atoms


# Temporary function
def new_atoms_detector(initial_atoms, final_atoms):
    """
    :param initial_atoms: initial dictionary with {"PDB atom names" : "index"}
    :param final_atoms: final dictionary with {"PDB atom names" : "index"}
    :return: dictionary with {"PDB atom names" : "index"} of the new atoms found.
    """
    # Find the differences between keys in both dictionaries
    differences = final_atoms.keys()-initial_atoms.keys()

    # Now, transform the object to dictionary in order to find the keys
    diff_dictionary = {}
    for atom_name in differences:
        diff_dictionary[atom_name] = final_atoms[atom_name]
    return diff_dictionary


def get_atom_properties(atoms_dictionary, template):
    """
    :param atoms_dictionary: dictionary with {"PDB atom names" : "index"}
    :param template: input template (string)
    :return: dictionary {"index" : ("vdw", "charge")}
    """
    # Definition of the pattern correspondent to NBON section
    NBON_PATTERN = "\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    # Selection of NBON section of the template
    nbon_section = section_selector(template, "NBON", "BOND")
    # Obtain all rows (list of lists)
    rows = re.findall(NBON_PATTERN, nbon_section)
    # Obtain a list with the indexes of atoms
    atom_indexes = []
    for atom_name, index in atoms_dictionary.items():
        atom_indexes.append(index)
    # Get the properties for each index
    properties = {}
    for row in rows:
        index = row[0]
        if index in atom_indexes:
            properties[index] = (row[1], row[3])
        else:
            pass

    return properties


initial_template = template_reader("mbez")
final_template = template_reader("pyjz")

atoms_selected_1 = atoms_selector(initial_template)
atoms_selected_2 = atoms_selector(final_template)

new_atoms = new_atoms_detector(atoms_selected_1, atoms_selected_2)
properties = get_atom_properties(new_atoms, final_template)

print(properties)

