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

    section_selected = re.search("{}\n(.*?){}".format(pattern_1, pattern_2), template,  re.DOTALL)

    return section_selected.group(0)


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
def transform_bond_length(original_atom, heavy_atom, atom_to_transform, initial_atm_dictionary, final_atm_dictionary,
                          initial_bnd_dictionary, final_bnd_dictionary, step, total_steps, modification = False):
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
    index_heavyatom_final = final_atm_dictionary[heavy_atom]
    index_heavyatom_initial = initial_atm_dictionary[heavy_atom]
    index_atom_to_transform = final_atm_dictionary[atom_to_transform]

    # Use this indexes to transform the bonds of the final dictionary to the original ones
    # We are only using original atom and initial dictionaries because the keys of the dictionary
    # (bonds) are repeated in both, initial and final, so we can replace the value of the correspondent key.
    initial_bonds = initial_bnd_dictionary.keys()
    final_bonds = final_bnd_dictionary.keys()
    step = step+1
    for bonds in initial_bonds:
        if bonds[1] == index_original and bonds[0] == index_heavyatom_initial:
            # Check this condition because in the final dictionary is not necessary to have the same key in bonds
            original_bond_length = initial_bnd_dictionary[bonds]
            for final_bond in final_bonds:
                if final_bond[0] == index_heavyatom_final and final_bond[1] == index_atom_to_transform:
                    # Add linearly bond length to the transformed atom
                    print("STEP : {}".format(step))
                    print("TOTAL STEPS : {}".format(total_steps))
                    if modification:
                        final_bond_length = float(final_bnd_dictionary[final_bond])
                    else:
                        final_bond_length = float(final_bnd_dictionary[final_bond]) * (total_steps + 1)
                    print("FINAL BOND LENGTH : {}".format(final_bond_length))
                    difference_of_len = abs(final_bond_length - float(original_bond_length))
                    print("DIFFERENCE OF LENGTH: {}".format(difference_of_len))
                    increase_of_len = difference_of_len / (total_steps+1)
                    print("INCREASE OF LENGTH: {}".format(increase_of_len))
                    original_bond_length = float(original_bond_length) + (step*increase_of_len)
                    print("STEP BOND LENGTH : {}".format(original_bond_length))
                    print(final_bnd_dictionary[final_bond])
                    final_bnd_dictionary[final_bond] = original_bond_length
                    print(final_bnd_dictionary[final_bond])

    return final_bnd_dictionary


# We will need to change this function in further updates
def transform_bond_length_grow(original_atom, heavy_atom, atom_to_transform, initial_atm_dictionary, final_atm_dictionary,
                               initial_bnd_dictionary, starting_bnd_dictionary, final_bnd_dictionary, step, total_steps,
                               modification = False):
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
    index_heavyatom_final = final_atm_dictionary[heavy_atom]
    index_heavyatom_initial = initial_atm_dictionary[heavy_atom]
    index_atom_to_transform = final_atm_dictionary[atom_to_transform]

    # Use this indexes to transform the bonds of the final dictionary to the original ones
    # We are only using original atom and initial dictionaries because the keys of the dictionary
    # (bonds) are repeated in both, initial and final, so we can replace the value of the correspondent key.
    initial_bonds = initial_bnd_dictionary.keys()
    final_bonds = final_bnd_dictionary.keys()
    step = step+1
    for bonds in initial_bonds:
        if bonds[1] == index_original and bonds[0] == index_heavyatom_initial:
            # Check this condition because in the final dictionary is not necessary to have the same key in bonds
            original_bond_length = initial_bnd_dictionary[bonds]
            print(original_bond_length)
            for final_bond in final_bonds:
                if final_bond[0] == index_heavyatom_final and final_bond[1] == index_atom_to_transform:
                    # Add linearly bond length to the transformed atom
                    print("STEP : {}".format(step))
                    print("TOTAL STEPS : {}".format(total_steps))
                    if modification:
                        final_bond_length = float(final_bnd_dictionary[final_bond])
                    else:
                        final_bond_length = float(final_bnd_dictionary[final_bond]) * (total_steps + 1)
                        print("FINAL BOND LENGTH : {}".format(final_bond_length))
                    difference_of_len = abs(final_bond_length - float(original_bond_length))
                    print("DIFFERENCE OF LENGTH: {}".format(difference_of_len))
                    increase_of_len = difference_of_len / (total_steps+1)
                    print("INCREASE OF LENGTH: {}".format(increase_of_len))
                    original_bond_length = float(original_bond_length) + (step*increase_of_len)
                    print("STEP BOND LENGTH : {}".format(original_bond_length))
                    print(final_bnd_dictionary[final_bond])
                    starting_bnd_dictionary[final_bond] = original_bond_length
                    print(starting_bnd_dictionary[final_bond])

    return starting_bnd_dictionary


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
    """
    Given a template and a dictionary of NBON properties the function rewrites the NBON section with the data in the
    dictionary.
    :param template: OPLS2005 template of a ligand. string
    :param properties_dict: dictionary {"index" : ("vdw", "charge")}
    :return: NBON section modified. string
    """
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
    """
    Given a template and a dictionary of BOND information the function rewrites the BOND section with the data in the
    dictionary.
    :param template: OPLS2005 template of a ligand. string
    :param bonds_dict: dictionary {("index1", "index2") : "bond length"}
    :return: BOND section modified. string
    """
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
                try:
                    properties_dict[index] = ((float(properties_dict[index][0]) / (steps + 1)),
                                              (float(properties_dict[index][1]) / (steps + 1)))
                except ZeroDivisionError:
                # We expect always to find positives or negatives values, otherwise we put this warning
                # just in case...
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign...")
                    continue
            else:
                pass
    return properties_dict


def set_pretemplate_properties(properties_dict, new_atoms_properties_dict, selected_atom_properties_dict, at_index, steps):
    """

    :param properties_dict:
    :param new_atoms_properties_dict:
    :param selected_atom_properties_dict:
    :param at_index:
    :param steps:
    :return:
    """
    for index in properties_dict.keys():
        # In fact, we are only using the index of the new_atoms_properties_dict to get the changing atoms
        for new_index in new_atoms_properties_dict.keys():
            if index == new_index:
                try:
                    # We are setting value/steps to the current data of VDW and charge
                    properties_dict[index] = ((float(properties_dict[index][0]) / (steps+1)),
                                             ((float(selected_atom_properties_dict[at_index][1]) / (steps+1))
                                               /
                                               len(new_atoms_properties_dict)))
                except ZeroDivisionError:
                # We expect always to find positives or negatives values, otherwise we put this warning
                # just in case...
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign...")
                continue
            else:
                pass
    return properties_dict


def set_properties_initial(properties_dict, new_atoms_properties_dict, steps):
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
                try:
                # We are setting value/steps to the current data of VDW and charge
                    properties_dict[index] = ((float(properties_dict[index][0]) / (steps)),
                                             (float(properties_dict[index][1]) / (steps)))
                except ZeroDivisionError:
                # We expect always to find positives or negatives values, otherwise we put this warning
                # just in case...
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign")
                    continue
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
                try:
                # We are adding value/steps to the current value of length
                    bonds_dictionary[(index1, index2)] = (float(bonds_dictionary[(index1, index2)]) / (steps+1))
                except ZeroDivisionError:
                    logger.warning("Charges of the atom are 0, we can not add charge if we do not know the sign")
                    continue
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
    content_list.append(atoms_section)
    content_list.append("\n")
    content_list.append(nbon_content)
    content_list.append("BOND\n")
    content_list.append(bond_content)
    content_list.append(angles_section)
    with open(os.path.join(output_path, output_filename), "w") as template_to_write:
        template_to_write.write("".join(content_list))


def get_specific_atom_properties(atom_pdb_names, atom_dictionary, properties_dictionary):
    """
    :param atom_pdb_names: list of strings with the atom names that we want to select.
    :param atom_dictionary: dictionary { "PDB atom name" : "index"} ("VDW", "CHARGE") } that we will use as reference to
    find the correct index of given atoms.
    :param properties_dictionary: dictionary {"index" : ("vdw", "charge")} that we are going to use to extract properties
    from the index obtained.
    :return: dictionary {"index":("VDW", "CHARGE")} of selected atoms.
    """
    dictionary = {}
    for name in atom_pdb_names:
        index = atom_dictionary[name]
        dictionary[index] = properties_dictionary[index]
    return dictionary


# These last three functions are the main algorithm to modify templates (we will modify them in further updates
# in order to avoid final templates)


def create_initial_template(initial_template, final_template, original_atom_to_mod, heavy_atom,  atom_to_transform,
                            output_template_filename="grwz_0", path="DataLocal/Templates/OPLS2005/HeteroAtoms/",
                            steps=10):
    """
    This function creates a pregrowing template. This template will be used to set the initial properties in a template
    in order to perform the growing. It will reduce the vdw to vdw/steps+1, the bond length to b.length/steps+1 and the
    charge to (charge_of_original_atom/steps+1)/number_of_new_atoms_ for all new atoms detected in the final template.
    :param initial_template: template file of the core ligand. filename string
    :param final_template: template file of the ligand with the fragment that we want to add. filename string
    :param original_atom_to_mod: PDB atom names of the atom of the initial template that will be used to extract the
    charge. list of strings
    :param heavy_atom: PDB atom name of the heavy atom bonded to the original_atom_to_mod. It will be use to replace
    the bond length in the final template for the original length. string
    :param output_template_filename: name of the output file with the template information. string
    :param steps: number of total steps to do the growing. int (10 by default)
    :return:
    """
    atoms_in_templates = []
    templates = []
    for template in [initial_template, final_template]:
        template_content = template_reader(template, path)
        templates.append(template_content)
        selected_atoms = atoms_selector(template_content)
        atoms_in_templates.append(selected_atoms)
    initial_atoms, final_atoms = atoms_in_templates
    initial_template_content, final_template_content = templates

    new_atoms = new_atoms_detector(initial_atoms, final_atoms)

    nbon_properties_new = get_atom_properties(new_atoms, final_template_content)
    nbon_properties_final = get_atom_properties(final_atoms, final_template_content)
    nbon_properties_initial = get_atom_properties(initial_atoms, initial_template_content)

    original_atom_prop = get_specific_atom_properties(original_atom_to_mod, initial_atoms, nbon_properties_initial)
    index_of_atom = list(original_atom_prop.keys())[0]
    set_pretemplate_properties(nbon_properties_final, nbon_properties_new, original_atom_prop, index_of_atom, steps)

    nbon_section = write_nbon_section(final_template_content, nbon_properties_final)

    bonds = get_bonds(final_template_content)
    initial_bonds = get_bonds(initial_template_content)

    set_bonds(bonds, nbon_properties_new, steps)
    transform_bond_length(original_atom_to_mod[0], heavy_atom, atom_to_transform, initial_atoms, final_atoms,
                          initial_bonds, bonds, step=0, total_steps=steps)  #  step 0 because we are in the i = 0 in the main loop

    bond_section = write_bond_section(final_template_content, bonds)

    write_template(final_template_content, output_template_filename, nbon_section, bond_section, path)


def generate_starting_template(initial_template_file, final_template_file, original_atom_to_mod, heavy_atom, atom_to_transform,
                               output_template_filename="grwz_ref", path="DataLocal/Templates/OPLS2005/HeteroAtoms/",
                               steps=10):
    """
    :param initial_template_file: template file of the initial ligand.
    :param final_template_file: template file of the ligand with the fragment that we want to add.
    :param original_atom_to_mod: PDB atom name of the atom that we want to transform into another (in initial template).
    :param heavy_atom: PDB atom name of the atom that will be transformed (in final template).
    :param output_template_filename: name of the output template.
    :param path: path to templates folder. string.
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
    properties_final = get_atom_properties(atoms_selected_final, final_template)
    new_atoms_properties = get_atom_properties(new_atoms, final_template)
    # Get the bonding data for initial and final dictionaries
    bonds_initial = get_bonds(initial_template)
    bonds_final = get_bonds(final_template)
    # We want to generate a starting template, so we will set the properties correspondent to the first step of growing
    set_properties(properties_final, new_atoms_properties, steps)
    # Now, we will repeat the same process for bonding data
    set_bonds(bonds_final, new_atoms_properties, steps)
    transform_bond_length(original_atom_to_mod[0], heavy_atom, atom_to_transform, atoms_selected_initial, atoms_selected_final,
                          bonds_initial, bonds_final, step=0, total_steps=steps)  #  step 0 because we are in the i = 0 in the main loop

    # Once we have all data in place, we will replace the current content of the final template for the starting
    # values needed to grow
    nbon_section = write_nbon_section(final_template, properties_final)
    bond_section = write_bond_section(final_template, bonds_final)
    # Finally, join everything and write a file with the output template:q!
    write_template(final_template, output_template_filename, nbon_section, bond_section, path)


def grow_parameters_in_template(starting_template_file, initial_template_file, final_template_file,
                                original_atom_to_mod, heavy_atom, atom_to_transform, output_template_filename,
                                path, step, total_steps):
    """
    :param starting_template_file: template file resultant of generate_starting_template() which contain the properties
    of the atoms and bonds modified to start the growing.
    :param initial_template_file: template file of the core ligand.
    :param final_template_file: template file of the ligand with the fragment that we want to add.
    :param original_atom_to_mod: PDB atom name of the atom that we want to transform into another (in initial template).
    :param heavy_atom: PDB atom name of the atom that will be transformed (in final template).
    :param output_template_filename: name of the output template.
    :param path: path to templates folder. string
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
    bonds_final = get_bonds(final_template)
    # We want to add values to a starting template
    modify_properties(properties_starting, new_atoms_properties, step)
    # Now, we will repeat the same process for bonding data
    modify_bonds(bonds_starting, new_atoms_properties, step)
    # Here we will set again the bonding distance of the original bond to the initial template one
    transform_bond_length_grow(original_atom_to_mod[0], heavy_atom, atom_to_transform, atoms_selected_initial,
                               atoms_selected_starting, bonds_initial, bonds_starting, bonds_final, step, total_steps,
                               True)
    # Once we have all data in place, we will replace the current content of the final template for the starting
    # values needed to grow
    nbon_section = write_nbon_section(final_template, properties_starting)
    bond_section = write_bond_section(final_template, bonds_starting)

    # Finally, join everything and write a file with the output template
    write_template(final_template, output_template_filename, nbon_section, bond_section, path)


#grow_parameters_in_template("grwz_ref","0kkz","grwz", ['H15'] ,"C15" ,"G4","grwz_2", "/home/carlespl/project/growing/grow/4DJU_4DJV/DataLocal/Templates/OPLS2005/HeteroAtoms/growing_templates",9,10)
