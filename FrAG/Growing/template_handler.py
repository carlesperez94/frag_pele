import sys
import re
import os
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Definition of reggex patterns
ATOM_PATTERN = "\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)"
NBON_PATTERN = "\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
BOND_PATTERN = "\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
WRITE_NBON_PATTERN = " {:5d}   {:3.4f}   {:3.4f}  {: 3.6f}   {:3.4f}   {:3.4f}   {:3.9f}  {: 3.9f}\n"
WRITE_BOND_PATTERN = " {:5d} {:5d}   {:5.3f} {: 2.3f}\n"


# Class definitions

class Atom:
    def __init__(self, index, PDB_atom_name, atom_type, sigma, epsilon, charge, rad_sgb,
                 rad_type, sgb_gamma, sgb_type, is_fragment=False, is_linker=False):
        self.index = index
        self.PDB_atom_name = PDB_atom_name
        self.atom_type = atom_type
        self.sigma = sigma
        self.epsilon = epsilon
        self.charge = charge
        self.rad_sgb = rad_sgb
        self.rad_type = rad_type
        self.sgb_gamma = sgb_gamma
        self.sgb_type = sgb_type
        self.is_fragment = is_fragment
        self.is_linker = is_linker

    def reduce_property(self, property, step, total_steps):
        if property == "sigma":
            self.sigma = float(self.sigma) * (step/(total_steps+1))
            return self.sigma
        elif property == "epsilon":
            self.epsilon = float(self.epsilon) * (step/(total_steps+1))
            return self.epsilon
        elif property == "charge":
            self.charge = float(self.charge) * (step/(total_steps+1))
            return self.charge
        elif property == "rad_sgb":
            self.rad_sgb = float(self.rad_sgb) * (step/(total_steps+1))
            return self.rad_sgb
        elif property == "rad_type":
            self.rad_type = float(self.rad_type) * (step/(total_steps+1))
            return self.rad_type
        elif property == "sgb_gamma":
            self.sgb_gamma = float(self.sgb_gamma) * (step/(total_steps+1))
            return self.sgb_gamma
        elif property == "sgb_type":
            self.sgb_type = float(self.sgb_type) * (step/(total_steps+1))
            return self.sgb_type
        else:
            raise ValueError(
                "Unknown property selected. Please, choose between the following properties: sigma, epsilon, charge, rad_sgb, rad_type, sgb_gamma, sgb_type")

    def set_fragment(self):
        self.is_fragment = True

    def set_linker(self):
        self.is_linker = True

        
class Template:
    def __init__(self, topology, nbon, bond, angles, is_initial=False):
        self.topology = topology
        self.nbon = nbon
        self.bond = bond
        self.angles = angles
        self.is_initial = is_initial

    def write_template(self, filename, path):
        content = []
        content.append("* LIGAND DATABASE FILE (OPLS2005)\n")
        content.append(self.topology[:-4])
        content.append(self.nbon[:-4])
        content.append(self.bond[:-4])
        content.append(self.angles)
        with open(os.path.join(path, filename), "w") as template_to_write:
            template_to_write.write("".join(content))


def section_selector(template, pattern_1, pattern_2):
    """
    From a template string, this function return a section between two patterns.
    :param template: input template (string).
    :param pattern_1: pattern which sets the begining of the section.
    :param pattern_2: pattern which sets the end of the section.
    :return: string with the content of the section.
    """

    section_selected = re.search("{}\n(.*?){}".format(pattern_1, pattern_2), template,  re.DOTALL)

    return section_selected.group(0)


def read_opls2005(file, path):
    """
    This function reads the content of a PELE's template and returns a Template object.
    :param file: name of the template file. string
    :param path: path to the file. string
    :return: Template object built with the content of the template.
    """

    with open(os.path.join(path, file), "r") as template_file:
        template_content = template_file.read()
        if not template_content:
            logger.critical("Template file {} is empty!".format(file))
    topology = section_selector(template_content, "\*", "NBON")
    nbon = section_selector(template_content, "NBON", "BOND")
    bond = section_selector(template_content, "BOND", "THET")
    angles = section_selector(template_content, "THET", "END")

    return Template(topology, nbon, bond, angles)


def pdb_atom_name_to_index(topology, pdb_atom_name):
    """
    From a topology (str), given a PDB atom name (str), it returns the correspondent index for the atom.
    :param topology: topology part of a template. string
    :param pdb_atom_name: name of an atom. string (len 4)
    :return: index. string object
    """
    rows = re.findall(ATOM_PATTERN, topology)
    for row in rows:
        if row[4] == pdb_atom_name:
            index = row[0]
            return index


def index_to_property(index, nbon, property):
    """
    Returns the selected property of an atom.
    :param index: index of the selected atom. string
    :param nbon: nbon part of a template. string
    :param property: choose between the following properties: sigma, epsilon, charge, rad_sgb, rad_type, sgb_gamma, sgb_type
    :return: depending on the selected propertie, returns an string with its content. string
    """
    rows = re.findall(NBON_PATTERN, nbon)
    for row in rows:
        if row[0] == index:
            if property == "sigma":
                sigma = row[1]
                return sigma
            elif property == "epsilon":
                epsilon = row[2]
                return epsilon
            elif property == "charge":
                charge = row[3]
                return charge
            elif property == "rad_sgb":
                rad_sgb = row[4]
                return rad_sgb
            elif property == "rad_type":
                rad_type = row[5]
                return rad_type
            elif property == "sgb_gamma":
                sgb_gamma = row[6]
                return sgb_gamma
            elif property == "sgb_type":
                sgb_type = row[7]
                return sgb_type
            else:
                raise ValueError("Unknown property selected. Please, choose between the following properties: sigma, epsilon, charge, rad_sgb, rad_type, sgb_gamma, sgb_type")


def template_to_atoms(template):
    rows_topology = re.findall(ATOM_PATTERN, template.topology)
    rows_nbon = re.findall(NBON_PATTERN, template.nbon)
    atom_list = []
    for row_top, row_nbon in zip(rows_topology, rows_nbon):
        index = row_top[0]
        atom_type = row_top[3]
        pdb_atom_name = row_top[4]
        sigma = row_nbon[1]
        epsilon = row_nbon[2]
        charge = row_nbon[3]
        rad_sgb = row_nbon[4]
        rad_type = row_nbon[5]
        sgb_gamma = row_nbon[6]
        sgb_type = row_nbon[7]
        atom_list.append(Atom(index, pdb_atom_name, atom_type, sigma, epsilon, charge, rad_sgb, rad_type, sgb_gamma, sgb_type))
    return atom_list

template = read_opls2005("ligz", "/home/carlespl/project/growing/grow/4CC5_frag8")
for atom in template_to_atoms(template):
    print(atom.reduce_property("sigma", 1, 10))
    print(atom.sigma)