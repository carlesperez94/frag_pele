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
    # Select the section of the templates where we have the atoms defined
    atoms_section = section_selector(template, "\*", "NBON")
    # Obtain all rows (list of lists)
    rows = re.findall(ATOM_PATTERN, atoms_section)
    # Get the atom name from all rows and insert it in a dictionary
    atoms = {}
    for row in rows:
        atom_name = row[4]
        atom_name_curated = atom_name.replace("_", "")
        atoms[atom_name_curated] = row[0]
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


# We will need to change this function in further updates
def transform_properties(original_atom, final_atom, initial_atm_dictionary, final_atm_dictionary,
                         initial_prp_dictionary, final_prp_dictionary):
    """
    :param original_atom: atom that we want to transform into another.
    :param final_atom: atom that we want to finally get.
    :param initial_atm_dictionary: initial dictionary {"PDB atom name" : "index"} which contain the atom.
    :param final_atm_dictionary: final dictionary {"PDB atom name" : "index"} which contain the atom.
    :param initial_prp_dictionary: initial dictionary {"index" : ("vdw", "charge")} which contain the properties
    of the template.
    :param final_prp_dictionary: final dictionary {"index" : ("vdw", "charge")} which contain the properties of
    the template.
    :return: dictionary with the properties modified.
    """
    # Collect indexes from original and final dictionaries
    index_original = initial_atm_dictionary[original_atom]
    index_final = final_atm_dictionary[final_atom]
    # Use this indexes to transform the properties of the final dictionary to the original ones
    initial_properties = initial_prp_dictionary[index_original]
    final_prp_dictionary[index_final] = initial_properties

    return final_prp_dictionary


def get_bonds(template):
    """
    :param template: template (string) that we want to extract bonding information
    :return: dictionary {("index_1", "index_2"): "bond length" }
    """
    # Selecting BOND information
    bonds_section = section_selector(template, "BOND", "THET")
    # Find data (list of lists)
    rows = re.findall(BOND_PATTERN, bonds_section)
    # Creating a dictionary having as key a tuple of indexes (atoms in bond) and the bond length
    bonds = {}
    for row in rows:
        bonds[(row[0], row[1])] = row[3]

    return bonds


def get_specific_bonds(atoms_dictionary, bonds_dictionary):
    """
    :param atoms_dictionary: dictionary with {"PDB atom names" : "index"} of the atoms that we want to get
    their bond length.
    :param bonds_dictionary: dictionary {("index_1", "index_2"): "bond length" } to obtain the information.
    :return: dictionary {("index_1", "index_2"): "bond length" }
    """
    # Get indexes of the dictionary with all the atoms
    atom_indexes = []
    for atom_name, index in atoms_dictionary.items():
        atom_indexes.append(index)
    # Get indexes of the dictionary with bonding data
    bonded_indexes = []
    for bond_indexes, length in bonds_dictionary.items():
        bonded_indexes.append(bond_indexes)
    # Create a dictionary where we are going to select the bonds for the atoms of the atom_dictionary
    selected_bonds_dictionary = {}
    for index in atom_indexes:
        for bond in bonded_indexes:
            # If we want to obtain only bonds that correspond to our atoms dictionary indexes
            # we have to apply this criteria (bond[1] is the atom that "recives" the bond)
            if bond[1] == index:
                selected_bonds_dictionary[bond] = bonds_dictionary[bond]
            else:
                pass
    return selected_bonds_dictionary


# We will need to change this function in further updates
def transform_bond_length(original_atom, initial_atm_dictionary, initial_bnd_dictionary, final_bnd_dictionary):
    """
    :param original_atom: name of the atom of the original template whose bond length will be used to be replaced
    for the final bond
    :param initial_atm_dictionary: dictionary with {"PDB atom names" : "index"} of the atoms of the original template
    :param initial_bnd_dictionary: dictionary {("index_1", "index_2"): "bond length" } of the initial template
    :param final_bnd_dictionary: dictionary {("index_1", "index_2"): "bond length" } of the final template
    :return: final_bnd_dictionary will be modified, getting the bond length of the original bonded atom in order to
    replace it in the final dictionary.
    """
    # Collect indexes from original and final dictionaries
    index_original = initial_atm_dictionary[original_atom]

    # Use this indexes to transform the bonds of the final dictionary to the original ones
    # We are only using original atom and initial dictionaries because the keys of the dictionary
    # (bonds) are repeated in both, initial and final, so we can replace the value of the correspondent key.
    initial_bonds = initial_bnd_dictionary.keys()
    for bonds in initial_bonds:
        if bonds[1] == index_original:
            final_bnd_dictionary[bonds] = initial_bnd_dictionary[bonds]
    return final_bnd_dictionary


def modify_properties(properties_dict, new_atoms_properties_dict, step):
    """
    :param properties_dict: dictionary {"index" : ("vdw", "charge")} that we want to modify.
    :param new_atoms_properties_dict: dictionary {"index" : ("vdw", "charge")} that we are going to use
    to know which atoms modify.
    :param step: integer with the number of the step of growing that we are in.
    :return: modification of properties_dict adding VDW and Charge to the atoms that we want to grow.
    """
    for index in properties_dict.keys():
        # In fact, we are only using the index of the new_atoms_properties_dict to get the changing atoms
        for new_index in new_atoms_properties_dict.keys():
            if index == new_index:
                if float(properties_dict[index][1]) != 0:
                    # We are adding value/steps to the current value of VDW and charge
                    properties_dict[index] = (float(properties_dict[index][0]) +
                                             (float(properties_dict[index][0])*step),
                                              float(properties_dict[index][1]) +
                                             (float(properties_dict[index][1])*step))
                else:
                    # We expect always to find positives or negatives values, otherwise we put this warning
                    # just in case...
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign")
            else:
                pass
    return properties_dict


def write_nbon_section(template, properties_dict):
    nbon_section = section_selector(template, "NBON", "BOND")
    rows = re.findall(NBON_PATTERN, nbon_section)
    section_modified = []
    sorted_keys_properties = sorted([int(x) for x in properties_dict.keys()])
    for key in sorted_keys_properties:
        for row in rows:
            if int(row[0]) == key:
                section_modified.append(WRITE_NBON_PATTERN.format(int(row[0]), float(properties_dict[str(key)][0]),
                                                                  float(row[2]), float(properties_dict[str(key)][1]),
                                                                  float(row[4]), float(row[5]),
                                                                  float(row[6]), float(row[7])))
            else:
                pass
    section_modified = "".join(section_modified)
    return section_modified


def modify_bonds(bonds_dictionary, new_atoms_properties_dict, step):
    """
    :param bonds_dictionary: dictionary {("index1", "index2") : "bond length"} that we want to modify.
    :param new_atoms_properties_dict: dictionary {"index" : ("vdw", "charge")} that we are going to use
    to know which bonds modify.
    :param step: integer with the number of the step of growing that we are in.
    :return: modification of bonds_dictionary adding length to the initial distance.
    """
    for index1, index2 in bonds_dictionary.keys():
        # We are only using the index of the new_atoms_properties_dict to get the changing atoms
        for new_index in new_atoms_properties_dict.keys():
            if index2 == new_index:
                # We are adding value/steps to the current value of length
                bonds_dictionary[(index1, index2)] = (float(bonds_dictionary[(index1, index2)]) +
                                                     (float(bonds_dictionary[(index1, index2)])*step))
            else:
                pass
    return bonds_dictionary


def write_bond_section(template, bonds_dict):
    bond_section = section_selector(template, "BOND", "THET")
    rows = re.findall(BOND_PATTERN, bond_section)
    section_modified = []
    # Sort keys (tuple) of the dictionary
    list_of_keys = list(bonds_dict.keys())
    list_of_keys.sort(key=lambda list: (int(list[0]), int(list[1])))
    for key in list_of_keys:
        for row in rows:
            if (row[0], row[1]) == key:
                section_modified.append(WRITE_BOND_PATTERN.format(int(row[0]), int(row[1]),
                                                                  float(row[2]), float(bonds_dict[tuple(key)])))
            else:
                pass
    section_modified = "".join(section_modified)
    return section_modified


def set_properties(properties_dict, new_atoms_properties_dict, steps):
    """
    :param properties_dict: dictionary {"index" : ("vdw", "charge")} that we want to set.
    :param new_atoms_properties_dict: dictionary {"index" : ("vdw", "charge")} that we are going to use
    to know which atoms modify.
    :param steps: integer with the number of steps that we are going to do to grow the template.
    :return: setting properties_dict VDW and Charge to the atoms that we want to grow.
    """
    for index in properties_dict.keys():
        # In fact, we are only using the index of the new_atoms_properties_dict to get the changing atoms
        for new_index in new_atoms_properties_dict.keys():
            if index == new_index:
                if float(properties_dict[index][1]) != 0:
                    # We are setting value/steps to the current data of VDW and charge
                    properties_dict[index] = ((float(properties_dict[index][0]) / steps),
                                             (float(properties_dict[index][1]) / steps))
                else:
                    # We expect always to find positives or negatives values, otherwise we put this warning
                    # just in case...
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign")
            else:
                pass
    return properties_dict


def set_bonds(bonds_dictionary, new_atoms_properties_dict, steps):
    """
    :param bonds_dictionary: dictionary {("index1", "index2") : "bond length"} that we want to set.
    :param new_atoms_properties_dict: dictionary {"index" : ("vdw", "charge")} that we are going to use
    to know which bonds modify.
    :param steps: integer with the number of steps that we are going to do to grow the template.
    :return: modification of bonds_dictionary setting the initial length.
    """
    for index1, index2 in bonds_dictionary.keys():
        # We are only using the index of the new_atoms_properties_dict to get the changing atoms
        for new_index in new_atoms_properties_dict.keys():
            if index2 == new_index:
                # We are adding value/steps to the current value of length
                bonds_dictionary[(index1, index2)] = (float(bonds_dictionary[(index1, index2)]) / steps)
            else:
                pass
    return bonds_dictionary


def write_template(reference_template, output_filename, nbon_content, bond_content,
                   output_path="DataLocal/Templates/OPLS2005/HeteroAtoms/"):
    """
    :param reference_template: string containing the whole template that we want to replace.
    :param output_filename: name of the file of the output template.
    :param nbon_content: string which contain the nbon section that we will replace to the reference template.
    :param bond_content: string which contain the bond section that we will replace to the reference template.
    """
    content_list = []
    atoms_section = section_selector(reference_template, "\*", "NBON")
    angles_section = section_selector(reference_template, "THET", "END")
    content_list.append("* LIGAND DATABASE FILE (OPLS2005)\n")
    content_list.append("*\n")
    content_list.append(atoms_section)
    content_list.append("NBON\n")
    content_list.append(nbon_content)
    content_list.append("BOND\n")
    content_list.append(bond_content)
    content_list.append("THET\n")
    content_list.append(angles_section)
    content_list.append("END")
    with open(os.path.join(output_path, output_filename), "w") as template_to_write:
        template_to_write.write("".join(content_list))


# These last two functions are the main algorithm to modify templates (we will modify them in further updates
# in order to avoid final templates)
def generate_starting_template(initial_template_file, final_template_file, original_atom_to_mod, final_atom_to_mod,
                               output_template_filename, path, steps=10):
    """
    :param initial_template_file: template file of the initial ligand.
    :param final_template_file: template file of the ligand with the fragment that we want to add.
    :param original_atom_to_mod: PDB atom name of the atom that we want to transform into another (in initial template).
    :param final_atom_to_mod: PDB atom name of the atom that will be transformed (in final template).
    :param output_template_filename: name of the output template.
    :param steps: number of growing steps.
    :return: template modified that will be used as starting point to do the growing process.
    """
    # Reading initial and final templates and convert them in strings
    initial_template = template_reader(initial_template_file, path)
    final_template = template_reader(final_template_file, path)
    # Select the atoms for this templates and convert them into dictionaries objects
    atoms_selected_initial = atoms_selector(initial_template)
    atoms_selected_final = atoms_selector(final_template)
    # Use this dictionaries in order to find differences in atoms to determine which are new ones
    new_atoms = new_atoms_detector(atoms_selected_initial, atoms_selected_final)
    # Get also the properties (VDW and Charge) of all the dictionaries
    properties_initial = get_atom_properties(atoms_selected_initial, initial_template)
    properties_final = get_atom_properties(atoms_selected_final, final_template)
    new_atoms_properties = get_atom_properties(new_atoms, final_template)
    # Get the bonding data for initial and final dictionaries
    bonds_initial = get_bonds(initial_template)
    bonds_final = get_bonds(final_template)
    # We want to generate a starting template, so we will set the properties correspondent to the first step of growing
    set_properties(properties_final, new_atoms_properties, steps)
    # As we are transforming a H to a heavy atom (C,O,S...), we want to keep the H properties in the final template so
    # given their PDB names we will transform these properties of the heavy atom into the H ones
    transform_properties(original_atom_to_mod, final_atom_to_mod, atoms_selected_initial, atoms_selected_final,
                         properties_initial, properties_final)
    # Now, we will repeat the same process for bonding data
    set_bonds(bonds_final, new_atoms_properties, steps)
    transform_bond_length(original_atom_to_mod, atoms_selected_initial, bonds_initial, bonds_final)
    # Once we have all data in place, we will replace the current content of the final template for the starting
    # values needed to grow
    nbon_section = write_nbon_section(final_template, properties_final)
    bond_section = write_bond_section(final_template, bonds_final)
    # Finally, join everything and write a file with the output template
    write_template(final_template, output_template_filename, nbon_section, bond_section, path)


# We could save a lot of steps and useless variables in this function if we use variables of the previous one,
# but I do not want to take the risk creating global variables
def grow_parameters_in_template(starting_template_file, initial_template_file, final_template_file,
                                original_atom_to_mod, final_atom_to_mod, output_template_filename,
                                path, step):
    """
    :param starting_template_file: template file resultant of generate_starting_template() which contain the properties
    of the atoms and bonds modified to start the growing.
    :param initial_template_file: template file of the ligand with the fragment that we want to add.
    :param final_template_file: template file of the ligand with the fragment that we want to add.
    :param original_atom_to_mod: PDB atom name of the atom that we want to transform into another (in initial template).
    :param final_atom_to_mod: PDB atom name of the atom that will be transformed (in final template).
    :param output_template_filename: name of the output template.
    :param step: current step of the growing loop.
    :return: template with the parameters increased linearly.
    """
    # Reading initial and final templates and convert them in strings
    starting_template = template_reader(starting_template_file, path)
    initial_template = template_reader(initial_template_file, path)
    final_template = template_reader(final_template_file, path)
    # Select the atoms for this templates and convert them into dictionaries objects
    atoms_selected_initial = atoms_selector(initial_template)
    atoms_selected_starting = atoms_selector(final_template)
    # Use this dictionaries in order to find differences in atoms to determine which are new ones
    new_atoms = new_atoms_detector(atoms_selected_initial, atoms_selected_starting)
    # Get also the properties (VDW and Charge) of the dictionaries
    properties_initial = get_atom_properties(atoms_selected_initial, initial_template)
    new_atoms_properties = get_atom_properties(new_atoms, final_template)
    properties_starting = get_atom_properties(atoms_selected_starting, starting_template)
    # Get the bonding data for initial and starting dictionaries
    bonds_initial = get_bonds(initial_template)
    bonds_starting = get_bonds(starting_template)
    # We want to add values to a starting template
    modify_properties(properties_starting, new_atoms_properties, step)
    transform_properties(original_atom_to_mod, final_atom_to_mod, atoms_selected_initial, atoms_selected_starting,
                         properties_initial, properties_starting)
    # Now, we will repeat the same process for bonding data
    modify_bonds(bonds_starting, new_atoms_properties, step)
    transform_bond_length(original_atom_to_mod, atoms_selected_initial, bonds_initial, bonds_starting)
    # Once we have all data in place, we will replace the current content of the final template for the starting
    # values needed to grow
    nbon_section = write_nbon_section(final_template, properties_starting)
    bond_section = write_bond_section(final_template, bonds_starting)
    # Finally, join everything and write a file with the output template
    write_template(final_template, output_template_filename, nbon_section, bond_section, path)
