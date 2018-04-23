import prody
import logging
import numpy as np
import math
import os
import shutil
import sys
import Bio.PDB as bio
# Local imports
import Growing.AddingFragHelpers.complex_to_prody as c2p
import Growing.AddingFragHelpers.pdb_joiner as pj



# Getting the name of the module for the log system
logger = logging.getLogger(__name__)
# Working directory variable
PRE_WORKING_DIR = "pregrow"


def extract_heteroatoms_pdbs(pdb, create_file = True, ligand_chain="L"):
    """
    From a pdb file, it extracts the chain L and checks if the structure has hydrogens. After that, the chain L is
    written in a new PDB file which will have the following format: "{residue name}.pdb".
    :param pdb: pdb file (with a ligand in the chain L).
    :return: Writes a new pdb file "{residue name}.pdb" with the chain L isolated an returns the residue name (string).
    """
    # Parse the complex file and isolate the ligand core and the fragment
    ligand = c2p.pdb_parser_ligand(pdb, ligand_chain)
    if ligand is None:
        logger.critical("The ligand can not be found. Ensure that the ligand of {} is the chain L".format(pdb))
    # Check if the ligand has H
    c2p.check_protonation(ligand)
    # Save the ligand in a PDB (the name of the file is the name of the residue)
    ligand_name = ligand.getResnames()[0]
    if create_file:
        prody.writePDB(ligand_name, ligand)
        logger.info("The ligand of {} has been extracted and saved in '{}.pdb'".format(pdb, ligand_name))
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
        bio_structure = pj.get_ligand_from_PDB("{}.pdb".format(pdb))
        # Get a list of atoms (BioPython Atoms) for each structure
        bioatomlist = pj.get_atoms_from_structure(bio_structure)
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
        atom_heavy = pj.select_atoms_from_list(atom_name, list_of_bioatoms)
        heavy_atoms.append(atom_heavy)
    return heavy_atoms


def extract_hydrogens(pdb_atom_names, lists_of_bioatoms, list_of_pdbs):
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
    for atom_name, pdb, list_of_bioatoms in zip(pdb_atom_names, list_of_pdbs, lists_of_bioatoms):
        complex = prody.parsePDB(pdb)
        # Select name of the H atoms bonded to this heavy atom (the place where we will grow)
        atom_name_hydrogens = pj.get_H_bonded_to_grow(atom_name, complex)
        # Select this hydrogen atoms
        atom_hydrogen = pj.select_atoms_from_list(atom_name_hydrogens, list_of_bioatoms)
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
    fragment_coords = pj.transform_coords(bioatom_list)
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


def join_structures(core_bond, fragment_bond, list_of_atoms, core_structure, fragment_structure):
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
    pj.superimpose(core_bond, fragment_bond, list_of_atoms)
    # Get the new coords and change them in prody
    transform_coords_from_bio2prody(fragment_structure, list_of_atoms)
    # Now, we have to remove the hydrogens of the binding
    h_atom_names = [core_bond[1].name, fragment_bond[0].name]
    merged_structure = bond(h_atom_names, [core_structure, fragment_structure])
    return merged_structure


def rotation_thought_axis(bond, theta, core_bond, list_of_atoms, fragment_bond, core_structure, fragment_structure):
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
    rotated_structure = join_structures(core_bond, fragment_bond, list_of_atoms, core_structure, fragment_structure)
    return rotated_structure


def check_collision(merged_structure, bond, theta, theta_interval, core_bond, list_of_atoms, fragment_bond,
                    core_structure, fragment_structure):
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
    check_possible_collision = merged_structure.select("resname {} and within 1.7 of resname {}".format(core_resname,
                                                                                                        frag_resname))
    # This list only should have the atom of the fragment that will be bonded to the core, so if not we will have to
    # solve it
    if len(check_possible_collision.getNames()) > 1 or check_possible_collision.getNames()[0] != bond[0].name:
        logger.info("We have a collision between atoms of the fragment and the core! Rotating the fragment to solve it...")
        theta = theta + theta_interval
        if theta >= math.pi*2:
            logger.warning("Not possible solution, decreasing the angle of rotation...")
        else:
            rotated_structure = rotation_thought_axis(bond, theta, core_bond, list_of_atoms, fragment_bond, core_structure,
                                                      fragment_structure)
            recall = check_collision(rotated_structure[0], bond, theta, theta_interval, core_bond, list_of_atoms,
                                     fragment_bond, core_structure, fragment_structure)
            return recall
    else:
        return merged_structure


def finishing_joining(molecule):
    """
    Given a ProDy molecule this function change the Resname of the atoms to "GRW" and the Resnum to "1". Following this
    process it is possible to transform a ProDy object with more than one element with different resnums and resnames
    into a single molecule.
    :param molecule: ProDy molecule.
    :return: ProDy molecule with Resname "GRW" and Resnum "1".
    """
    molecule.setResnames("GRW")
    molecule.setResnums(1)


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


def get_waters_in_pdb(pdb_input):
    """
    Given a pdb file return a string with the waters contained in the file.
    :param pdb_input: pdb input file
    :return: part of the pdb that contains the waters
    """
    water_lines = []
    with open(pdb_input) as pdb:
        for line in pdb:
            if line.split()[0] == "HETATM" and line.split()[3] == "HOH":
                water_lines.append(line)
    return "".join(water + "TER\n" * (n % 3 == 2) for n, water in enumerate(water_lines))


def main(pdb_complex_core, pdb_fragment, pdb_atom_core_name, pdb_atom_fragment_name, steps, core_chain="L",
         fragment_chain="L", output_file_to_tmpl="growing_result.pdb", output_file_to_grow="initialization_grow.pdb"):
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
    """
    if not PRE_WORKING_DIR:
        os.mkdir(PRE_WORKING_DIR)
        os.chdir(PRE_WORKING_DIR)
    else:
        os.chdir(PRE_WORKING_DIR)
    # Get the selected chain from the core and the fragment and convert them into ProDy molecules.
    ligand_core = c2p.pdb_parser_ligand(pdb_complex_core, core_chain)
    fragment = c2p.pdb_parser_ligand(pdb_fragment, fragment_chain)
    # We will check that the structures are protonated. We will also create a new PDB file for each one and we will get
    # the residue name of each ligand.
    core_residue_name = extract_heteroatoms_pdbs(pdb_complex_core, True, core_chain)
    frag_residue_name = extract_heteroatoms_pdbs(pdb_fragment, True, fragment_chain)
    # We will use the PDBs previously generated to get a list of Bio.PDB.Atoms for each structure
    bioatoms_core_and_frag = from_pdb_to_bioatomlist([core_residue_name, frag_residue_name])
    # Then, we will have to transform the atom names of the core and the fragment to a list object
    # (format required by functions)
    pdb_atom_names = [pdb_atom_core_name, pdb_atom_fragment_name]
    # Using the Bio.PDB.Atoms lists and this names we will get the heavy atoms that we will use later to do the bonding
    heavy_atoms = extract_heavy_atoms(pdb_atom_names, bioatoms_core_and_frag)
    # Once we have the heavy atoms, for each structure we will obtain the hydrogens bonded to each heavy atom.
    # We will need pdbs because we will use the information of the protein to select the hydrogens properly.
    hydrogen_atoms = extract_hydrogens(pdb_atom_names, bioatoms_core_and_frag, [pdb_complex_core, pdb_fragment])
    # Create a list with the atoms that form a bond in core and fragment.
    core_bond = [heavy_atoms[0], hydrogen_atoms[0]]
    fragment_bond = [hydrogen_atoms[1], heavy_atoms[1]]  # This has to be in inverted order to do correctly the superimposition
    logger.info("Performing a superimposition of bond {} of the fragment on bond {} of the core..."
                .format(fragment_bond, core_bond))
    # Using the previous information we will superimpose the whole fragment on the bond of the core in order to place
    # the fragment in the correct position, deleting the H.
    merged_structure = join_structures(core_bond, fragment_bond, bioatoms_core_and_frag[1], ligand_core, fragment)
    # It is possible to create intramolecular clashes after placing the fragment on the bond of the core, so we will
    # check if this is happening, and if it is, we will perform rotations of 10º until avoid the clash.
    check_results = check_collision(merged_structure[0], heavy_atoms, 0, math.pi/18, core_bond,
                                    bioatoms_core_and_frag[1], fragment_bond, ligand_core, fragment)
    # If we do not find a solution in the previous step, we will repeat the rotations applying only increments of 1º
    if not check_results:
        check_results = check_collision(merged_structure[0], heavy_atoms, 0, math.pi/180, core_bond,
                                        bioatoms_core_and_frag[1], fragment_bond, ligand_core, fragment)
    # Now, we want to extract this structure in a PDB to create the template file after the growing. We will do a copy
    # of the structure because then we will need to resize the fragment part, so be need to keep it as two different
    # residues.
    structure_to_template = check_results.copy()
    # Once we have all the atom names unique, we will rename the resname and the resnum of both, core and fragment, to
    # GRW and 1. Doing this, the molecule composed by two parts will be transformed into a single one.
    changing_names = pj.extract_and_change_atomnames(structure_to_template, fragment.getResnames()[0])
    molecule_names_changed, changing_names_dictionary = changing_names

    # Check if there is still overlapping names
    if pj.check_overlapping_names(molecule_names_changed):
        logger.critical("{} is repeated in the fragment and the core. Please, change this atom name of the core by"
                        " another one.".format(pj.check_overlapping_names(molecule_names_changed)))
    logger.info("The following names of the fragment have been changed:")
    for transformation in changing_names_dictionary:
        logger.info("{} --> {}".format(transformation, changing_names_dictionary[transformation]))
    finishing_joining(molecule_names_changed)
    # Extract a PDB file to do the templates
    prody.writePDB(output_file_to_tmpl, molecule_names_changed)
    logger.info("The result of core + fragment has been saved in '{}'. This will be used to create the template file."
                .format(output_file_to_tmpl))
    # Now, we will use the original molecule to do the resizing of the fragment.
    reduce_molecule_size(check_results, frag_residue_name, steps)
    point_reference = check_results.select("name {} and resname {}".format(pdb_atom_fragment_name, frag_residue_name))
    fragment_segment = check_results.select("resname {}".format(frag_residue_name))
    translate_to_position(hydrogen_atoms[0].get_coord(), point_reference.getCoords(), fragment_segment)

    # Repeat all the preparation process to finish the writing of the molecule.
    changing_names = pj.extract_and_change_atomnames(check_results, fragment.getResnames()[0])
    molecule_names_changed, changing_names_dictionary = changing_names
    finishing_joining(molecule_names_changed)
    logger.info("The result of core + fragment(small) has been saved in '{}'. This will be used to initialise the growing."
                .format(output_file_to_grow))
    # Add the protein to the ligand
    prot = extract_protein_from_complex(pdb_complex_core)
    prody.writePDB("protein.pdb", prot)
    prody.writePDB("ligand_grown.pdb", molecule_names_changed)
    with open("protein.pdb") as protein:
        content_prot = protein.readlines()
        content_prot = content_prot[1:]
        content_prot = "".join(content_prot)
    with open("ligand_grown.pdb") as lig:
        content_lig = lig.readlines()
        content_lig = content_lig[1:]
        content_lig = "".join(content_lig)
    output_file = []
    output_file.append("{}TER\n".format(content_prot))
    waters = get_waters_in_pdb(pdb_complex_core)
    output_file.append(waters)
    output_file.append("{}TER".format(content_lig))
    out_joined = "".join(output_file)
    with open(output_file_to_grow, "w") as output :
        output.write(out_joined)
    # Make a copy of output files in the main directory
    shutil.copy(output_file_to_grow, "../")
    os.chdir("../")
    # In further steps we will probably need to recover the names of the atoms for the fragment, so for this reason we
    # are returning this dictionary in the function.
    return changing_names_dictionary, hydrogen_atoms, "{}.pdb".format(core_residue_name), output_file_to_tmpl, output_file_to_grow
