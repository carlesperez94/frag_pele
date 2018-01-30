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
    :return:
    """
    ROW_PATTERN = "\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)"

    # Select the section of the templates where we have the atoms defined
    atoms_section = section_selector(template, "\*", "NBON")
    # Obtain all rows (list of lists)
    rows = re.findall(ROW_PATTERN, atoms_section)
    # Get the atom name from all rows and insert it in a dictionary
    atoms = {}
    for row in rows:
        atoms[row[0]] = row[4]
    logger.info("The following atoms have been found in the initial template:\n{}".format(atoms))
    return atoms


initial_template = template_reader("mbez")
atoms_selected = atoms_selector(initial_template)
print(atoms_selected)
