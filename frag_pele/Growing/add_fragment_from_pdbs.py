import prody
import sys
import logging
import numpy as np
import pickle
from scipy.spatial import distance
import math
import os
import shutil
import re
import Bio.PDB as bio
# Local imports
import frag_pele.constants as c
from frag_pele.Helpers import checker
from frag_pele.Growing.AddingFragHelpers import complex_to_prody, pdb_joiner, atom_constants
try:
    import rdkit
    RDKIT = True
except ImportError:
    RDKIT = False
if RDKIT:
    from lib_prep.FragmentTools import tree_detector

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)
# Lists
LIST_OF_IONS = ["ZN", "MN", "FE", "CO", "NI", "CA", "CD"]


def extract_heteroatoms_pdbs(pdb, create_file=True, ligand_chain="L", get_ligand=False, output_folder="."):
    """
    From a pdb file, it extracts the chain L and checks if the structure has hydrogens. After that, the chain L is
    written in a new PDB file which will have the following format: "{residue name}.pdb".
    :param pdb: pdb file (with a ligand in the chain L).
    :return: Writes a new pdb file "{residue name}.pdb" with the chain L isolated an returns the residue name (string).
    """
    # Parse the complex file and isolate the ligand core and the fragment
    ligand = complex_to_prody.pdb_parser_ligand(pdb, ligand_chain)
    if ligand is None:
        logger.critical("The ligand can not be found. Ensure that the ligand of {} is the chain {}".format(pdb, ligand_chain))
    # Check if the ligand has H
    complex_to_prody.check_protonation(ligand)
    # Save the ligand in a PDB (the name of the file is the name of the residue)
    ligand_name = ligand.getResnames()[0]
    if create_file:
        prody.writePDB(os.path.join(output_folder, ligand_name), ligand)
        print("The ligand of {} has been extracted and saved in '{}.pdb'".format(pdb, os.path.join(output_folder, ligand_name)))
    if get_ligand is True:
        return ligand
    else:
        return ligand_name


def from_pdb_to_bioatomlist(list_of_pdb_names):
    """
    Given a pdb name string (without the extension ".pdb") the function reads it as a Bio.PDB structure and extract the
    atoms found as a list of Bio.PDB.Atom objects.
    :param list_of_pdb_names: list of strings with pdb names.
    :return: list of lists with the Bio.PDB.Atom objects found in each pdb file.
    """
    list_of_lists = []
    for pdb in list_of_pdb_names:
        # Now, read the PDBs as BioPython object
        bio_structure = pdb_joiner.get_ligand_from_PDB("{}.pdb".format(pdb))
        # Get a list of atoms (BioPython Atoms) for each structure
        bioatomlist = pdb_joiner.get_atoms_from_structure(bio_structure)
        list_of_lists.append(bioatomlist)
    return list_of_lists


def extract_heavy_atoms(pdb_atom_names, lists_of_bioatoms):
    """
    Given a heavy atom name (string) and a list of Bio.PDB.Atom objects, it selects this atom of the list and return it
    as a single object.
    :param pdb_atom_names: heavy atom name (string).
    :param lists_of_bioatoms: list of Bio.PDB.Atom objects.
    :return: Bio.PDB.Atom object correspondent to the heavy atom name.
    """
    heavy_atoms = []
    # Select the heavy atoms for each list that we will want to bond together in further steps
    for atom_name, list_of_bioatoms in zip(pdb_atom_names, lists_of_bioatoms):
        atom_heavy = pdb_joiner.select_atoms_from_list(atom_name, list_of_bioatoms)
        heavy_atoms.append(atom_heavy)
    return heavy_atoms


def extract_hydrogens(pdb_atom_names, lists_of_bioatoms, list_of_pdbs, h_core=None, h_frag=None, c_chain="L", f_chain="L"):
    """
    Given a heavy atom name (string), a list of Bio.PDB.Atoms objects and a list of pdb files, it returns the hydrogens
    at bonding distance of the heavy atom. If there is more than one, a checking of contacts with the
    protein will be performed. In case of finding a possible contact between the hydrogen and the protein, we will
    reject this hydrogen and we will repeat the check in another one. If all the hydrogens have contacts with the
    protein, the first of them will be selected and a warning will be printed.
    :param pdb_atom_names: heavy atom name (string).
    :param lists_of_bioatoms: list of Bio.PDB.Atom objects.
    :param list_of_pdbs: list of PDB files.
    :return: Bio.PDB.Atom object correspondent to the hydrogen bonded to the heavy atom
    """

    hydrogens = []
    selected_hydrogens = [h_core, h_frag]
    chains = [c_chain, f_chain]
    for atom_name, pdb, list_of_bioatoms, sel_h, chain in zip(pdb_atom_names, list_of_pdbs, lists_of_bioatoms, selected_hydrogens, chains):
        complex = prody.parsePDB(pdb)
        # Select name of the H atoms bonded to this heavy atom (the place where we will grow)
        atom_name_hydrogens = pdb_joiner.get_H_bonded_to_grow(atom_name, complex, sel_h, chain=chain)
        # Select this hydrogen atoms
        atom_hydrogen = pdb_joiner.select_atoms_from_list(atom_name_hydrogens, list_of_bioatoms)
        hydrogens.append(atom_hydrogen)
    return hydrogens


def transform_coords_from_bio2prody(fragment_prody, bioatom_list):
    """
    Given a fragment (prody molecule object) and a list of Bio.PDB.Atom objects correspondent to this fragment, it
    replace the coordinates of the prody molecule for the ones of the list of Bio.PDB.Atom objects and returns the new
    coordinates.
    :param fragment_prody: prody molecule object.
    :param bioatom_list: list of Bio.PDB.Atom objects.
    :return: array with the new coordinates of the prody molecule object.
    """
    fragment_coords = pdb_joiner.transform_coords(bioatom_list)
    fragment_prody.setCoords(fragment_coords)
    return fragment_prody.getCoords()


def bond(hydrogen_atom_names, molecules):
    """
    Given a list with names of hydrogens (bonded to the heavy atoms that we want to link) and a list of molecules (prody
    molecule object), this function errase this hydrogens and bond the heavy atoms linked to them. In order to create
    this new bond we concatenate the pairs of molecules (prody molecule object).
    Note that the lists need to have pairs of hydrogens and molecules because the bonding is going to happen between
    adjacent pairs of hydrogens atoms and molecules.
    :param hydrogen_atom_names: names of hydrogens (bonded to the heavy atoms that we want to link). Each two elements
    of the list will be considered as hydrogens that will form the bond.
    :param molecules: prody molecule object. Each two elements of the list will be considered as molecules
    that will be bonded.
    :return: list of prody molecule objects as a result of merging the pairs of molecules.
    """
    list_of_pairs = []
    for hydrogen, molecule in zip(hydrogen_atom_names, molecules):
        # Now, we have to remove the hydrogens of the binding, so we will select everything except these H
        mol_no_h = molecule.select("not name {}".format(hydrogen))
        list_of_pairs.append(mol_no_h)
    # Merging both parts into a single one
    i = 0
    bonds = []
    while i < len(list_of_pairs):
        merged = list_of_pairs[i].copy() + list_of_pairs[i+1].copy()
        bonds.append(merged)
        i += 2
    return bonds

def join_structures_to_rotate(core_bond, fragment_bond, list_of_atoms, core_structure, fragment_structure):
    """
    It joins two ProDy structures into a single one, merging both bonds (core bond and fragment bond) creating a unique
    bond between the molecules. In order to do that this function performs a cross superimposition (in BioPython) of
    the whole fragment using as reference (fixed part) the atoms of the bond. Then, it transforms this BioPython object
    into ProDy molecule with the coordinates modified. Once we have all ready, the Hydrogens of the bonds will be
    deleted and both structures will be concatenated into a single ProDy object.
    :param core_bond: Bio.PDB.Atom list with two elements: [heavy atom, hydrogen atom]. These two atoms have to be the
    ones participating in the bond of the core that we would like to use as linking point between core and fragment.
    :param fragment_bond: Bio.PDB.Atom list with two elements: [hydrogen atom, heavy atom]. These two atoms have to be the
    ones participating in the bond of the fragment that we would like to use as linking point between fragment and core.
    :param list_of_atoms: Bio.PDB.Atom list with all the atoms of the fragment. These atoms have to contain
    coordinates.
    :param core_structure: ProDy molecule that contain only the core ligand.
    :param fragment_structure: ProDy molecule that contain only the fragment ligand.
    :return: ProDy molecule with the core_structure and the fragment_structure (with the coordinates modified)
    concatenated.
    """
    # Superimpose atoms of the fragment to the core bond
    pdb_joiner.superimpose(core_bond, fragment_bond, list_of_atoms)
    # Get the new coords and change them in prody
    transform_coords_from_bio2prody(fragment_structure, list_of_atoms)
    # Now, we have to remove the hydrogens of the binding
    h_atom_names = [core_bond[1].name, fragment_bond[0].name]
    merged_structure = bond(h_atom_names, [core_structure, fragment_structure])
    return merged_structure


def join_structures(core_bond, fragment_bond, core_structure, fragment_structure, pdb_complex,
                    pdb_fragment, chain_complex, chain_fragment, output_path, only_grow=False):
    """
    It joins two ProDy structures into a single one, merging both bonds (core bond and fragment bond) creating a unique
    bond between the molecules. In order to do that this function performs a cross superimposition (in BioPython) of
    the whole fragment using as reference (fixed part) the atoms of the bond. Then, it transforms this BioPython object
    into ProDy molecule with the coordinates modified. Once we have all ready, the Hydrogens of the bonds will be
    deleted accordingly with the bond type, and both structures will be concatenated into a single ProDy object.

    :param core_bond: Bio.PDB.Atom list with two elements: [heavy atom, hydrogen atom]. These two atoms have to be the
    ones participating in the bond of the core that we would like to use as linking point between core and fragment.
    :param fragment_bond: Bio.PDB.Atom list with two elements: [hydrogen atom, heavy atom]. These two atoms have to be the
    ones participating in the bond of the fragment that we would like to use as linking point between fragment and core.
    :param list_of_atoms: Bio.PDB.Atom list with all the atoms of the fragment. These atoms have to contain
    coordinates.
    :param core_structure: ProDy molecule that contain only the core ligand.
    :param fragment_structure: ProDy molecule that contain only the fragment ligand.
    :param pdb_complex: path to the file in PDB with the protein-ligand complex.
    :param pdb_fragment: path to the file in PDB with the fragment.
    :param chain_complex: label of the ligand chain in the pdb_complex.
    :param chain_fragment: label of the ligand chain in the pdb_fragment.
    :return: ProDy molecule with the core_structure and the fragment_structure (with the coordinates modified)
    concatenated.
    """
    name_to_replace_core = core_bond[1].name
    name_to_replace_fragment = fragment_bond[0].name
    if only_grow:
        return 0, name_to_replace_core, name_to_replace_fragment
    if RDKIT:
        atoms_to_delete_core = tree_detector.main(pdb_complex, (core_bond[0].name, name_to_replace_core),
                                                  chain_ligand=chain_complex)
        atoms_to_delete_fragment = tree_detector.main(pdb_fragment, (fragment_bond[1].name, name_to_replace_fragment),
                                                      chain_ligand=chain_fragment)
    else:
        print("WARNING: YOU CAN NOT REPLACE HEAVY ATOMS FOR HYDROGENS WITHOUT RDKIT!")
    fragment_rdkit = rdkit.Chem.MolFromPDBFile(pdb_fragment, removeHs=False)
    bond_type = detect_bond_type(fragment_rdkit, fragment_bond[0], fragment_bond[1])
    if RDKIT:
        if core_bond[1].element != "H":
            atom_replaced_idx = replace_heavy_by_hydrogen(core_bond[1], core_structure)
            new_coords, new_dist = correct_hydrogen_position(hydrogen_atom=core_structure[atom_replaced_idx],
                                                             atom_to_bond_with=core_bond[0],
                                                             structure=core_structure[atom_replaced_idx])
            core_structure[atom_replaced_idx].setCoords(new_coords)
            names_to_keep = list(
                set(core_structure.getNames()) ^ set(atoms_to_delete_core))  # Compare two sets and get the common items
            names_to_keep.remove(name_to_replace_core)
            core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))
            prody.writePDB(os.path.join(output_path, "{}.pdb".format(core_structure.getResnames()[0])),
                           core_structure)  # Overwrite the initial structure
            name_to_replace_core = core_structure[atom_replaced_idx].getName()
            core_bond[1].coord = new_coords
        if fragment_bond[0].element != "H":
            atom_replaced_idx = replace_heavy_by_hydrogen(fragment_bond[0], fragment_structure)
            new_coords, new_dist = correct_hydrogen_position(hydrogen_atom=fragment_structure[atom_replaced_idx],
                                                             atom_to_bond_with=fragment_bond[1],
                                                             structure=fragment_structure[atom_replaced_idx])
            fragment_structure[atom_replaced_idx].setCoords(new_coords)
            names_to_keep = list(
                set(fragment_structure.getNames()) ^ set(atoms_to_delete_fragment))  # Compare two sets and get the common items
            names_to_keep.remove(name_to_replace_fragment)
            fragment_structure = fragment_structure.select("name {}".format(" ".join(names_to_keep)))
            prody.writePDB(os.path.join(output_path, "{}.pdb".format(fragment_structure.getResnames()[0])),
                           fragment_structure)  # Overwrite the initial structure
            name_to_replace_fragment = fragment_structure[atom_replaced_idx].getName()
            fragment_bond[0].coord = new_coords
    bio_list = from_pdb_to_bioatomlist([os.path.join(output_path, "{}".format(fragment_structure.getResnames()[0]))])[0] # Its a list, so we keep only the unique element that is inside
    # Superimpose atoms of the fragment to the core bond
    pdb_joiner.superimpose(core_bond, fragment_bond, bio_list)
    # Get the new coords and change them in prody
    transform_coords_from_bio2prody(fragment_structure, bio_list)
    # Now, we have to remove the hydrogens of the binding
    h_atom_names = [name_to_replace_core, name_to_replace_fragment]
    # Correcting linking distance
    new_coords, new_distance = correct_bonding_distance(atom_reference=core_bond[0], atom_to_correct=fragment_bond[1],
                                                        reference_structure=core_structure, movil_structure=fragment_structure,
                                                        bond_type=bond_type)
    if bond_type == "double":
        hydrogen_to_delete = pdb_joiner.get_H_bonded_to_atom(core_bond[0].name, core_structure,
                                                             bond_dist=new_distance,
                                                 	     banned_hydrogen=name_to_replace_core,
                                                             chain=chain_complex)
        names_to_keep = list(core_structure.getNames())
        names_to_keep.remove(hydrogen_to_delete)
        core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))

    if bond_type == "triple":
        for n in range(2):
            hydrogen_to_delete = pdb_joiner.get_H_bonded_to_atom(core_bond[0].name, core_structure,
                                                                 bond_dist=new_distance,
                                                                 banned_hydrogen=name_to_replace_core,
                                                                 chain=chain_complex)
            names_to_keep = list(core_structure.getNames())
            names_to_keep.remove(hydrogen_to_delete)
            core_structure = core_structure.select("name {}".format(" ".join(names_to_keep)))

    fragment_structure.setCoords(new_coords)
    merged_structure = bond(h_atom_names, [core_structure, fragment_structure])
    return merged_structure, name_to_replace_core, name_to_replace_fragment, new_distance # new_distance to modify the clash treshold
    # in check_collisions()


def detect_fragment_bond_type(fragment_structure, atom1, atom2, dictionary_of_distances):
    coords_atom1 = find_coords_of_atom(atom1.name, fragment_structure)
    coords_atom2 = find_coords_of_atom(atom2.name, fragment_structure)
    distance = compute_distance_between_atoms(coords_atom1, coords_atom2)
    for bond, bond_distance in dictionary_of_distances.items():
        if bond[0] == atom1.element and bond[1] == atom2.element and bond_distance+0.06 >= distance > bond_distance-0.06:
            return bond[0], bond[1], bond[2]

def convert_pdbnames_to_rdkit_idx(molecule):
   names_relations = {}
   for atom in molecule.GetAtoms():
       names_relations[atom.GetPDBResidueInfo().GetName().strip()] = atom.GetIdx()
   return names_relations

def get_bond_types_from_mol(molecule):
    bonds = {}
    for bond in molecule.GetBonds():
        bonds[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())] = str(bond.GetBondType()).lower()
        bonds[(bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())] = str(bond.GetBondType()).lower()
    return bonds

def detect_bond_type(molecule, atom1, atom2):
    names = convert_pdbnames_to_rdkit_idx(molecule)
    bonds = get_bond_types_from_mol(molecule)
    bond_names = (names[atom1.name], names[atom2.name])
    bond_type = bonds[bond_names]
    return bond_type

def correct_bonding_distance(atom_reference, atom_to_correct, reference_structure, movil_structure, bond_type="single"):
    hvy_core_coords = find_coords_of_atom(atom_reference.name, reference_structure)
    hvy_frag_coords = find_coords_of_atom(atom_to_correct.name, movil_structure)
    try:
        new_distance = atom_constants.BONDING_DISTANCES[atom_reference.element.upper(), atom_to_correct.element.upper(), bond_type]
    except KeyError:
        raise KeyError(f"Not defined bond {atom_reference.element.upper()}-{ atom_to_correct.element.upper()}. \
please define the maestro's bond distance into the file: frag_pele/frag_pele/Growing/AddingFragHelpers/atom_constants.py")
    new_coords = modify_distance_between_structures(hvy_core_coords, hvy_frag_coords, movil_structure.getCoords(),
                                                    new_distance)
    return new_coords, new_distance


def correct_hydrogen_position(hydrogen_atom, atom_to_bond_with, structure):
    coords_atom_to_bond = atom_to_bond_with.coord
    coords_H_to_move = hydrogen_atom.getCoords()
    new_distance = atom_constants.BONDING_DISTANCES["H", atom_to_bond_with.element, "single"]
    new_coords = modify_distance_between_structures(coords_atom_to_bond, coords_H_to_move, structure.getCoords(),
                                                    new_distance)
    return new_coords, new_distance


def replace_heavy_by_hydrogen(heavy_atom, structure):
    print("HEAVY ATOM TO REPLACE DETECTED!")
    index_of_atom = find_index_by_name(heavy_atom.name, structure)
    structure[index_of_atom].setElement("H")
    hydrogen_name = autoname_atoms(structure.getNames(), "H")
    structure[index_of_atom].setName(hydrogen_name)
    print("{} has been replaced for {}.".format(heavy_atom.name, structure.getNames()[index_of_atom]))
    return index_of_atom


def find_index_by_name(name, structure):
    for n, atom in enumerate(structure.getNames()):
        if atom == name:
            return n
    print("Atom not found!")


def find_coords_of_atom(atom_to_find, structure_prody):
    for atomname, coords in zip(structure_prody.getNames(), structure_prody.getCoords()):
        if atomname == atom_to_find:
            return coords


def autoname_atoms(list_of_atom_names, element):
    dictionary_of_elements = {}
    for name in list_of_atom_names:
        element_of_name = " ".join(re.split("[^a-zA-Z]*", name)).strip()
        if element_of_name in dictionary_of_elements.keys():
            dictionary_of_elements[element_of_name] += 1
        else:
            dictionary_of_elements[element_of_name] = 1
    while True:
        final_name = "{}{}".format(element, dictionary_of_elements[element])
        if final_name not in list_of_atom_names:
            return final_name
        else:
            dictionary_of_elements[element] += 1


def find_coords_of_atom(atom_to_find, structure_prody):
    try:
        for atomname, coords in zip(structure_prody.getNames(), structure_prody.getCoords()):
            if atomname == atom_to_find:
                return coords
    except AttributeError:  # If we get the attribute error probably the structure has a single atom.
        coords = structure_prody.getCoords()
        return coords


def compute_vector_between_atoms(coords_atom_1, coords_atom_2):
    vector = coords_atom_2 - coords_atom_1
    return vector  # Vector from point 1 to point 2


def compute_distance_between_atoms(coords_atom_1, coords_atom_2):
    vector = compute_vector_between_atoms(coords_atom_1, coords_atom_2)
    module = np.linalg.norm(vector)
    return module


def compute_unit_vector_between_atoms(coords_atom_1, coords_atom_2):
    vector = compute_vector_between_atoms(coords_atom_1, coords_atom_2)
    module = np.linalg.norm(vector)
    unit_vector = vector / module
    return unit_vector


def modify_distance_between_structures(coords_core, coords_fragment, coords_to_move, new_distance):
    unit_vector = compute_unit_vector_between_atoms(coords_core, coords_fragment)
    new_point = coords_core + (unit_vector*new_distance)
    vector_to_add = compute_vector_between_atoms(coords_fragment, new_point)
    new_coords = coords_to_move + vector_to_add
    return new_coords


def rotation_thought_axis(bond, theta, core_bond, list_of_atoms, fragment_bond, core_structure, fragment_structure,
                          pdb_complex, pdb_fragment, chain_complex, chain_fragment, output_path, only_grow=False):
    """
    Given a core molecule and a fragment, this function rotates the fragment atoms a certain theta angle around an axis
    (set by the bond).
    :param bond: Bio.PDB.Atom list composed by two elements: [heavy atom of the core, heavy atom of the fragment]
    :param theta: Rotation angle in rads.
    :param core_bond: Bio.PDB.Atom list with two elements: [heavy atom, hydrogen atom]. These two atoms have to be the
    ones participating in the bond of the core that we would like to use as linking point between core and fragment.
    :param list_of_atoms: list_of_atoms: Bio.PDB.Atom list with all the atoms of the fragment. These atoms have to contain
    coordinates.
    :param fragment_bond: Bio.PDB.Atom list with two elements: [hydrogen atom, heavy atom]. These two atoms have to be the
    ones participating in the bond of the fragment that we would like to use as linking point between fragment and core.
    :param core_structure: ProDy molecule that contain only the core ligand.
    :param fragment_structure: ProDy molecule that contain only the fragment ligand.
    :return: ProDy molecule with the core_structure and the fragment_structure rotated around the axis of the bond.
    """
    # Obtain the axis that we want to use as reference for the rotation
    vector = bond[1].get_vector() - bond[0].get_vector()
    # Obtain the rotation matrix for the vector (axis) and the angle (theta)
    rot_mat = bio.rotaxis(theta, vector)
    for atom in list_of_atoms:
        # Multiply the matrix of coordinates for the transpose of the rotation matrix to get the coordinates rotated
        atom.transform(rot_mat, (0, 0, 0))
        transform_coords_from_bio2prody(fragment_structure, list_of_atoms)
    rotated_structure = join_structures_to_rotate(core_bond, fragment_bond, list_of_atoms, core_structure, fragment_structure)
    return rotated_structure


def rotate_throught_bond(bond, angle, rotated_atoms, atoms_fixed):
    # Obtain the axis that we want to use as reference for the rotation
    vector = bond.getCoords()[0] - bond.getCoords()[1]
    vector = bio.Vector(vector)
    # Obtain the rotation matrix for the vector (axis) and the angle (theta)
    rot_mat = bio.rotaxis(angle, vector)
    new_coords = []
    for coords in rotated_atoms.getCoords():
        new_coord = np.dot(coords, rot_mat)
        new_coords.append(new_coord)
    rotated_atoms.setCoords(new_coords)
    structure_result = atoms_fixed + rotated_atoms
    return structure_result


def check_collision(merged_structure, bond, theta, theta_interval, core_bond, list_of_atoms, fragment_bond,
                    core_structure, fragment_structure, pdb_complex, pdb_fragment, chain_complex, chain_fragment,
                    output_path, threshold_clash=1.70, only_grow=False):
    """
    Given a structure composed by a core and a fragment, it checks that there is not collisions between the atoms of
    both. If it finds a collision, the molecule will be rotated "theta_interval" radians and the checking will be
    repeated. If it is not possible to find a conformation without atom collisions, it will print a warning.
    :param merged_structure: ProDy molecule with the core_structure and the fragment_structure concatenated.
    :param bond: Bio.PDB.Atom list composed by two elements: [heavy atom of the core, heavy atom of the fragment]
    :param theta: Initial rotation angle in rads.
    :param theta_interval: Rotation angle that will be added to theta.
    :param core_bond: Bio.PDB.Atom list with two elements: [heavy atom, hydrogen atom]. These two atoms have to be the
    ones participating in the bond of the core that we would like to use as linking point between core and fragment.
    :param list_of_atoms: list_of_atoms: Bio.PDB.Atom list with all the atoms of the fragment. These atoms have to contain
    coordinates.
    :param fragment_bond: Bio.PDB.Atom list with two elements: [hydrogen atom, heavy atom]. These two atoms have to be the
    ones participating in the bond of the fragment that we would like to use as linking point between fragment and core.
    :param core_structure: ProDy molecule that contain only the core ligand.
    :param fragment_structure: ProDy molecule that contain only the fragment ligand.
    :return: ProDy molecule with the core_structure and the fragment_structure (rotated and without intra-molecular
    clashes) around the axis of the bond.
    """
    core_resname = bond[0].get_parent().get_resname()
    frag_resname = bond[1].get_parent().get_resname()
    if core_resname is frag_resname:
        logger.critical("The resname of the core and the fragment is the same. Please, change one of both")
    print("resname {} and within {} of resname {}".format(core_resname,
                                                                                                       threshold_clash,
                                                                                                        frag_resname))
    check_possible_collision = merged_structure.select("resname {} and within {} of resname {}".format(core_resname,
                                                                                                       threshold_clash,
                                                                                                        frag_resname))

    # This list only should have the atom of the fragment that will be bonded to the core, so if not we will have to
    # solve it
    if len(check_possible_collision.getNames()) > 1 or check_possible_collision.getNames()[0] != bond[0].name:
        print("We have a collision between atoms of the fragment {} and the core {} within {} A!" 
              " Rotating the fragment to solve it...".format(frag_resname, core_resname, threshold_clash))
        theta = theta + theta_interval
        if theta >= math.pi*2:
            print("Not possible solution, decreasing the angle of rotation...")
        else:
            rotated_structure = rotation_thought_axis(bond, theta, core_bond, list_of_atoms, fragment_bond, core_structure,
                                                      fragment_structure, pdb_complex, pdb_fragment, chain_complex, chain_fragment,
                                                      output_path=output_path, only_grow=only_grow)
            #prody.writePDB("testing_{}.pdb".format(theta), rotated_structure[0])
            recall = check_collision(merged_structure=rotated_structure[0], bond=bond, theta=theta, theta_interval=theta_interval, 
                                     core_bond=core_bond, list_of_atoms=list_of_atoms, fragment_bond=fragment_bond, core_structure=core_structure,
                                     fragment_structure=fragment_structure, pdb_complex=pdb_complex, pdb_fragment=pdb_fragment, 
                                     chain_complex=chain_complex, chain_fragment=chain_fragment,
                                     output_path=output_path, threshold_clash=threshold_clash, only_grow=only_grow)
            return recall
    else:
        return merged_structure


def get_previous_bond(structure, core_atom, core_resname):
    bond_selection = structure.select("resname {} and within 1.75 of name {}".format(core_resname, core_atom))
    return bond_selection


def finishing_joining(molecule, chain):
    """
    Given a ProDy molecule this function change the Resname of the atoms to "GRW" and the Resnum to "1". Following this
    process it is possible to transform a ProDy object with more than one element with different resnums and resnames
    into a single molecule.
    :param molecule: ProDy molecule.
    :return: ProDy molecule with Resname "GRW" and Resnum "1".
    """
    molecule.setResnames("GRW")
    molecule.setResnums(1)
    molecule.setChids(chain)


def compute_centroid(molecule):
    """
    Given a ProDy molecule, the function extract the coordinates of their atoms and compute the centroid of the
    molecule.
    :param molecule: ProDy molecule object.
    :return: centroid of the molecule, tuple(X,Y,Z).
    """
    coords = molecule.getCoords()
    x = []
    y = []
    z = []
    for coord in coords:
        x.append(float(coord[0]))
        y.append(float(coord[1]))
        z.append(float(coord[2]))
    centroid = (np.mean(x),  np.mean(y), np.mean(z))
    return centroid


def move_atom_along_vector(initial_coord, final_coord, position_proportion):
    """
    Given two points (atom coordinates), this function moves the initial point a distance of "length of the vector
    formed by the two coordinates" * "position_proportion" on the vector's direction.
    :param initial_coord: initial 3D coordinates (X, Y, Z). numpy.ndarray
    :param final_coord: final 3D coordinates (X, Y, Z). numpy.ndarray
    :param position_proportion: proportion of movement that we would like to apply on the initial atom on the vector's
    direction. float(generally between 0 and 1)
    :return:
    """
    vector = final_coord - initial_coord
    new_coords = initial_coord + (position_proportion * vector)
    return new_coords


def reduce_molecule_size(molecule, residue, steps):
    """
    This function performs a reduction of the size of a given residue of a ProDy molecule object.
    :param molecule: ProDy molecule object.
    :param residue: Resname of the residue of the molecule that we want to reduce. string
    :param proportion: proportion of reduction of the size that we want to apply to the selected residue (between 0 and
    1). float
    :return: modify the coordinates of the selected residue for the result of the reduction.
    """
    proportion = 1-(1/(steps+1))
    if proportion >= 0 and proportion <= 1:
        selection = molecule.select("resname {}".format(residue))
        centroid = compute_centroid(selection)
        for atom in selection:
            atom_coords = atom.getCoords()
            new_coords = move_atom_along_vector(atom_coords, centroid, proportion)
            atom.setCoords(new_coords)
    else:
        logger.critical("Sorry, reduce_molecule_size() needs a proportion value between 0 and 1!")


def translate_to_position(initial_pos, final_pos, molecule):
    """
    This function applies a translation of a whole molecule, using the vector from the initial_pos to the final_pos.
    :param initial_pos: initial position in 3D coordinates. Generally we use the coordinates of an atom. numpy. ndarray
    :param final_pos: final position in 3D coordinates. Generally we use the coordinates of an atom. numpy. ndarray
    :param molecule: ProDy molecule object.
    :return: modify the coordinates of the molecule.
    """
    translation = initial_pos - final_pos
    coords_to_move = molecule.getCoords()
    list_of_new_coords = []
    for coords in coords_to_move:
        new_coords = coords + translation
        list_of_new_coords.append(new_coords[0])
    molecule.setCoords(list_of_new_coords)


def extract_protein_from_complex(pdb_file):
    """
    Given a pdb file containing a complex (ligand + protein) it returns only the protein.
    :param pdb_file: pdb file with a complex. string.
    :return: ProDy molecule with only the protein.
    """
    complex = prody.parsePDB(pdb_file)
    protein = complex.select("protein")
    return protein


def get_waters_or_ions_in_pdb(pdb_input):
    """
    Given a pdb file return a string with the waters contained in the file.
    :param pdb_input: pdb input file
    :return: part of the pdb that contains the waters
    """
    water_lines = []
    ion_lines = []
    with open(pdb_input) as pdb:
        for line in pdb:
            if line.startswith("HETATM"):
                if line[21] == "A" and line[17:20].split()[0] == "HOH":
                    water_lines.append(line)
                elif line[21] == "A" and line[17:20].split()[0] in LIST_OF_IONS:
                    ion_lines.append(line)
    if len(water_lines) > 0 and len(ion_lines) == 0:
        water = "".join(water + "TER\n" * (n % 3 == 2) for n, water in enumerate(water_lines))
        return water
    if len(ion_lines) > 0 and len(water_lines) == 0:
        ions = "".join(ion + "TER\n" for ion in ion_lines)
        return ions
    if len(ion_lines) > 0 and len(water_lines) > 0:
        water = "".join(water + "TER\n" * (n % 3 == 2) for n, water in enumerate(water_lines))
        ions = "".join(ion + "TER\n" for ion in ion_lines)
        return "".join(ions+water)


def get_everything_except_ligand(pdb_input, ligand_chain):
    pdb_content = []
    with open(pdb_input) as pdb:
        for line in pdb:
            if (line.startswith("ATOM") or line.startswith("HETATM") and line[21] != ligand_chain) or line.startswith("TER"):
                pdb_content.append(line)
    return "".join(pdb_content)


def check_water(pdb_input):
    """
    Given a pdb file checks if it contains water molecules.
    :param pdb_input: pdb input file
    :return: True or False
    """
    checker = False
    with open(pdb_input) as pdb:
        for line in pdb:
            if "HETATM" in line:
                if line.split()[3] == "HOH":
                    print("Your pdb file contains water molecules")
                    checker = True
                    break
    return checker


def lignames_replacer(pdb_file, original_ligname, new_ligname):
    """
    Given a PDB file it replace the name of a ligand for a new one.
    :param pdb_file: file in PDB format
    :param original_ligname: original name of the ligand
    :param new_ligname: new name of the ligand that will replace the original name
    :return:
    """
    with open(pdb_file) as pdb:
        content = pdb.readlines()
    for index, line in enumerate(content):
        if line.startswith("HETATM"):
            line = line.replace(original_ligname, new_ligname)
            content[index] = line
    pdb_modified = "".join(content)
    with open(pdb_file, "w") as overwrite_pdb:
        overwrite_pdb.write(pdb_modified)


def check_and_fix_repeated_lignames(pdb1, pdb2, ligand_chain_1="L", ligand_chain_2="L"):
    """
    It checks if two pdbs have the same ligand name or if the pdb file 1 has as ligand name "GRW" and it is replaced
    by "LIG".
    :param pdb1: pdb file 1
    :param pdb2: pdb file 2
    :return:
    """
    ligname_1 = extract_heteroatoms_pdbs(pdb1, create_file=False, ligand_chain=ligand_chain_1)
    ligname_2 = extract_heteroatoms_pdbs(pdb2, create_file=False, ligand_chain=ligand_chain_2)
    if ligname_1 == ligname_2 or ligname_1 == "GRW":
        logging.warning("REPEATED NAMES IN LIGANDS FOR THE FILES: '{}' and '{}'. {} replaced by LIG ".format(pdb1, pdb2, ligname_1))
        lignames_replacer(pdb1, ligname_1, "LIG")


def main(pdb_complex_core, pdb_fragment, pdb_atom_core_name, pdb_atom_fragment_name, steps, core_chain="L",
         fragment_chain="L", output_file_to_tmpl="growing_result.pdb", output_file_to_grow="initialization_grow.pdb",
         h_core=None, h_frag=None, rename=False, threshold_clash=1.70, output_path=None, only_grow=False):
    """
    From a core (protein + ligand core = core_chain) and fragment (fragment_chain) pdb files, given the heavy atoms
    names that we want to connect, this function add the fragment to the core structure. We will get three PDB files:
    (1) the ligand core of the complex isolated, that will be used in further steps to generate the template of the
    initial structure; (2) the ligand completed with the core and the fragment added, also prepared to generate the
    template of the final structure; (3) the pdb file that will be used to initialise PELE simulations. Here we have
    the core structure with the fragment added, but this fragment has been size-reduced in order to get small bond
    lengths between its atoms. (During PELE simulations this distance will increase linearly until it reaches the
    bond length given by the template of the final structure)
    :param pdb_complex_core: pdb file with a complex (protein + ligand) that will be used as core to perform the
    addition of the fragment. The chain of the ligand needs to be named as "L". We will also use the information of the
    protein to perform calculations of contacts with the ligand.
    :param pdb_fragment: pdb file, normally with only the ligand (please, put "L" as name of the chain that contain the
    ligand), that will be added to the core.
    :param pdb_atom_core_name: heavy atom name (string) of the ligand core where we want to add the fragment and
    form a new bond.
    :param pdb_atom_fragment_name: heavy atom name (string) of the ligand fragment where we want to perform the
    connection to form a new bond with the core.
    :param core_chain: name of the chain which contains the ligand in the pdb of the core complex. string
    :param fragment_chain: name of the chain which contains the ligand in the pdb of the fragment. string
    :param output_file_to_tmpl: name of the pdb file with the result of the connection between the core and the fragment
    (single ligand). string. The resname of the molecule will be "GRW" and the resnum "1". "growing_result.pdb" by
    default.
    :param output_file_to_grow: name of the pdb file that will be used to initialise PELE simulations. string.
    "initialization_grow.pdb" by default.
    :param h_core: if the user wants to select an specific hydrogen atom of the core to create the new bond, its name
    must be specified here.
    :param h_frag: if the user wants to select an specific hydrogen atom of the fragment to create the new bond, its name
    must be specified here.
    :param rename: if set, the names of the pdb atom names will be replaced with "G+atom_number_fragment".
    :param threshold_clash: distance that will be used to identity which atoms are doing clashes between atoms of the
    fragment and the core.
    :returns: [changing_names_dictionary, hydrogen_atoms, "{}.pdb".format(core_residue_name), output_file_to_tmpl,
    output_file_to_grow, core_original_atom, fragment_original_atom]

    """
    WORK_PATH = os.path.join(output_path, c.PRE_WORKING_DIR)
    if not os.path.exists(WORK_PATH):
        os.mkdir(WORK_PATH)
    # Check that ligand names are not repeated
    check_and_fix_repeated_lignames(pdb_complex_core, pdb_fragment, core_chain, fragment_chain)
    for pdb_file in (pdb_complex_core, pdb_fragment):
        logging.info("Checking {} ...".format(pdb_file))
        checker.check_and_fix_pdbatomnames(pdb_file)
    # Get the selected chain from the core and the fragment and convert them into ProDy molecules.
    ligand_core = complex_to_prody.pdb_parser_ligand(pdb_complex_core, core_chain)
    fragment = complex_to_prody.pdb_parser_ligand(pdb_fragment, fragment_chain)
    # We will check that the structures are protonated. We will also create a new PDB file for each one and we will get
    # the residue name of each ligand.
    core_residue_name = extract_heteroatoms_pdbs(pdb_complex_core, True, core_chain, output_folder=WORK_PATH)
    frag_residue_name = extract_heteroatoms_pdbs(pdb_fragment, True, fragment_chain, output_folder=WORK_PATH)
    # We will use the PDBs previously generated to get a list of Bio.PDB.Atoms for each structure
    bioatoms_core_and_frag = from_pdb_to_bioatomlist([os.path.join(WORK_PATH, core_residue_name),
                                                     os.path.join(WORK_PATH, frag_residue_name)])
    # Then, we will have to transform the atom names of the core and the fragment to a list object
    # (format required by functions)
    pdb_atom_names = [pdb_atom_core_name, pdb_atom_fragment_name]
    # Using the Bio.PDB.Atoms lists and this names we will get the heavy atoms that we will use later to do the bonding
    heavy_atoms = extract_heavy_atoms(pdb_atom_names, bioatoms_core_and_frag)
    # Once we have the heavy atoms, for each structure we will obtain the hydrogens bonded to each heavy atom.
    # We will need pdbs because we will use the information of the protein to select the hydrogens properly.
    hydrogen_atoms = extract_hydrogens(pdb_atom_names, bioatoms_core_and_frag, [pdb_complex_core, pdb_fragment], h_core,
                                       h_frag, core_chain, fragment_chain)
    # Create a list with the atoms that form a bond in core and fragment.
    core_bond = [heavy_atoms[0], hydrogen_atoms[0]]
    fragment_bond = [hydrogen_atoms[1], heavy_atoms[1]]  # This has to be in inverted order to do correctly the superimposition
    logger.info("Performing a superimposition of bond {} of the fragment on bond {} of the core..."
                .format(fragment_bond, core_bond))
    # Using the previous information we will superimpose the whole fragment on the bond of the core in order to place
    # the fragment in the correct position, deleting the H.
    merged_structure, core_original_atom, fragment_original_atom, new_dist = join_structures(core_bond, fragment_bond,
                                                                                             ligand_core, fragment,
                                                                                             pdb_complex_core, pdb_fragment,
                                                                                             core_chain, fragment_chain,
                                                                                             output_path=WORK_PATH, only_grow=only_grow)
    prody.writePDB("merged.pdb", merged_structure[0])
    if not only_grow:
        # It is possible to create intramolecular clashes after placing the fragment on the bond of the core, so we will
        # check if this is happening, and if it is, we will perform rotations of 10º until avoid the clash.
        # check_results = check_collision(merged_structure[0], heavy_atoms, 0, math.pi/18, core_bond,
        #                            bioatoms_core_and_frag[1], fragment_bond, ligand_core, fragment)
        # check_collision(merged_structure, bond, theta, theta_interval, core_bond, list_of_atoms, fragment_bond,
        #            core_structure, fragment_structure)
        check_results = check_collision(merged_structure=merged_structure[0], bond=heavy_atoms, theta=0,
                                        theta_interval=math.pi/18, core_bond=core_bond, list_of_atoms=bioatoms_core_and_frag[1],
                                        fragment_bond=fragment_bond, core_structure=ligand_core, fragment_structure=fragment,
                                        pdb_complex=pdb_complex_core, pdb_fragment=pdb_fragment, chain_complex=core_chain,
                                        chain_fragment=fragment_chain, output_path=WORK_PATH, threshold_clash=new_dist+0.01,
                                        only_grow=only_grow)
        # If we do not find a solution in the previous step, we will repeat the rotations applying only increments of 1º
        if not check_results:
            check_results = check_collision(merged_structure=merged_structure[0], bond=heavy_atoms, theta=0,
                                            theta_interval=math.pi/180, core_bond=core_bond, list_of_atoms=bioatoms_core_and_frag[1],
                                            fragment_bond=fragment_bond, core_structure=ligand_core, fragment_structure=fragment,
                                            pdb_complex=pdb_complex_core, pdb_fragment=pdb_fragment, chain_complex=core_chain,
                                            chain_fragment=fragment_chain, output_path=WORK_PATH,
                                            threshold_clash=new_dist+0.01, only_grow=only_grow)
        # Now, we want to extract this structure in a PDB to create the template file after the growing. We will do a copy
        # of the structure because then we will need to resize the fragment part, so be need to keep it as two different
        # residues.
        try:
            structure_to_template = check_results.copy()
        except AttributeError:
            raise AttributeError("Frag cannot superimpose the fragment onto the core's hydrogen.  \
In order to create space for the fragment \
manually rotate the hydrogen bond of the core where the fragment will be attached to.   \
We are currently working to fix this automatically")

        # Once we have all the atom names unique, we will rename the resname and the resnum of both, core and fragment, to
        # GRW and 1. Doing this, the molecule composed by two parts will be transformed into a single one.
        changing_names = pdb_joiner.extract_and_change_atomnames(structure_to_template, fragment.getResnames()[0],
                                                                 core_residue_name, rename=rename)
        molecule_names_changed, changing_names_dictionary = changing_names

        # Check if there is still overlapping names
        if pdb_joiner.check_overlapping_names(molecule_names_changed):
            logger.critical("{} is repeated in the fragment and the core. Please, change this atom name of the core by"
                            " another one.".format(pdb_joiner.check_overlapping_names(molecule_names_changed)))
        logger.info("The following names of the fragment have been changed:")
        for transformation in changing_names_dictionary:
            logger.info("{} --> {}".format(transformation, changing_names_dictionary[transformation]))
        finishing_joining(molecule_names_changed, core_chain)
        # Extract a PDB file to do the templates
        prody.writePDB(os.path.join(WORK_PATH, output_file_to_tmpl), molecule_names_changed)
        logger.info("The result of core + fragment has been saved in '{}'. This will be used to create the template file."
                    .format(os.path.join(WORK_PATH, output_file_to_tmpl)))
        # Now, we will use the original molecule to do the resizing of the fragment.
        reduce_molecule_size(check_results, frag_residue_name, steps)
        point_reference = check_results.select("name {} and resname {}".format(pdb_atom_fragment_name, frag_residue_name))
        fragment_segment = check_results.select("resname {}".format(frag_residue_name))
        translate_to_position(hydrogen_atoms[0].get_coord(), point_reference.getCoords(), fragment_segment)

        # Repeat all the preparation process to finish the writing of the molecule.
        changing_names = pdb_joiner.extract_and_change_atomnames(check_results, fragment.getResnames()[0], core_residue_name, rename=rename)
        molecule_names_changed, changing_names_dictionary = changing_names
        finishing_joining(molecule_names_changed, core_chain)
        logger.info("The result of core + fragment(small) has been saved in '{}'. This will be used to initialise the growing."
                    .format(os.path.join(WORK_PATH, output_file_to_grow)))
        # Add the protein to the ligand
        output_ligand_grown_path = os.path.join(WORK_PATH, "ligand_grown.pdb")
        prody.writePDB(output_ligand_grown_path, molecule_names_changed)

        with open(output_ligand_grown_path) as lig:
            content_lig = lig.readlines()
            content_lig = content_lig[1:]
            content_lig = "".join(content_lig)

        # Join all parts of the PDB
        output_file = []
        chain_not_lig = get_everything_except_ligand(pdb_complex_core, core_chain)
        output_file.append(chain_not_lig)
        output_file.append("{}TER".format(content_lig))
        out_joined = "".join(output_file)
        with open(os.path.join(WORK_PATH, output_file_to_grow), "w") as output: # Save the file in the pregrow folder
            output.write(out_joined)
        # Make a copy of output files in the main directory
        shutil.copy(os.path.join(WORK_PATH, output_file_to_grow), ".")  # We assume that the user will be running FrAG in PELE's main folder...
        # In further steps we will probably need to recover the names of the atoms for the fragment, so for this reason we
        # are returning this dictionary in the function.
        with open(os.path.join(WORK_PATH, "changingatoms.dict"), "wb") as pkl:
            pickle.dump(changing_names_dictionary, pkl)
    else:
        with open( os.path.join(WORK_PATH, "changingatoms.dict"), "rb") as pkl:
            changing_names_dictionary = pickle.load(pkl)


    return changing_names_dictionary, hydrogen_atoms, "{}.pdb".format(core_residue_name), \
           os.path.join(WORK_PATH, output_file_to_tmpl), \
           os.path.join(WORK_PATH, output_file_to_grow), \
           core_original_atom, fragment_original_atom


