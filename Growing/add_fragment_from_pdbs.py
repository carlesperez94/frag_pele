import prody
import logging
import numpy as np
import sys
# Local imports
import complex_to_smile as c2s
import pdb_joiner as pj
import Bio.PDB as bio


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def extract_heteroatoms_pdbs(pdb):
    """
    From a pdb file, it extracts the chain L and checks if the structure has hydrogens. After that, the chain L is
    written in a new PDB file which will have the following format: "{residue name}.pdb".
    :param pdb: pdb file (with a ligand in the chain L).
    :return: Writes a new pdb file "{residue name}.pdb" with the chain L isolated an returns the residue name (string).
    """
    # Parse the complex file and isolate the ligand core and the fragment
    ligand = c2s.pdb_parser_ligand(pdb, ligand_chain="L")
    # Check if the ligand has H
    c2s.check_protonation(ligand)
    # Save the ligand in a PDB (the name of the file is the name of the residue)
    ligand_name = ligand.getResnames()[0]
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


def main(pdb_complex_core, pdb_fragment, pdb_atom_core_name, pdb_atom_fragment_name,
                          output_file="growing_result.pdb"):
    """
    From a core (protein + ligand core = chain L) and fragment (chain L) pdb files, given the heavy atoms names that we
    want to connect, this function writes a new PDB with the fragment added to the core structure.
    :param pdb_complex_core: pdb file with a complex (protein + ligand) that will be used as core to perform the
    addition of the fragment. The chain of the ligand needs to be named as "L". We will also use the information of the
    protein to perform calculations of contacts with the ligand.
    :param pdb_fragment: pdb file, normally with only the ligand (please, put "L" as name of the chain that contain the
    ligand), that will be added to the core.
    :param pdb_atom_core_name: heavy atom name (string) of the ligand core where we want to add the fragment and
    form a new bond.
    :param pdb_atom_fragment_name: heavy atom name (string) of the ligand fragment where we want to perform the
    connection to form a new bond with the core.
    :param output_file: pdb file with the result of the connection between the core and the fragment (single ligand).
    The resname of the molecule will be "GRW" and the resnum "1". "growing_result.pdb" by default.
    """
    ligand_core = c2s.pdb_parser_ligand(pdb_complex_core, ligand_chain="L")
    fragment = c2s.pdb_parser_ligand(pdb_fragment, ligand_chain="L")
    core_residue_name = extract_heteroatoms_pdbs(pdb_complex_core)
    frag_residue_name = extract_heteroatoms_pdbs(pdb_fragment)
    bioatoms_core_and_frag = from_pdb_to_bioatomlist([core_residue_name, frag_residue_name])
    pdb_atom_names = [pdb_atom_core_name, pdb_atom_fragment_name]
    heavy_atoms = extract_heavy_atoms(pdb_atom_names, bioatoms_core_and_frag)
    hydrogen_atoms = extract_hydrogens(pdb_atom_names, bioatoms_core_and_frag, [pdb_complex_core, pdb_fragment])
    # Create a list with the atoms that form a bond in core and fragment
    core_bond = [heavy_atoms[0], hydrogen_atoms[0]]
    fragment_bond = [hydrogen_atoms[1], heavy_atoms[1]]  # This has to be in inverted order
    logger.info("Performing a superimposition of bond {} of the fragment on bond {} of the core..."
                .format(fragment_bond, core_bond))
    # Superimpose atoms of the fragment to the core bond
    pj.superimpose(core_bond, fragment_bond, bioatoms_core_and_frag[1])
    # Get the new coords and change them in prody
    transform_coords_from_bio2prody(fragment, bioatoms_core_and_frag[1])
    # Now, we have to remove the hydrogens of the binding
    h_atom_names = []
    for h in hydrogen_atoms:
        name = h.name
        h_atom_names.append(name)
    merged_structure = bond(h_atom_names, [ligand_core, fragment])
    # Change the resnames and resnums to set everything as a single molecule
    # Now its only one element of a list, but we keep it as list if we would like to increase the amount of bonds later
    merged_structure[0].setResnames("GRW")
    merged_structure[0].setResnums(1)
    prody.writePDB(output_file, merged_structure[0])
    logger.info("The result of core + fragment has been saved in '{}'".format(output_file))


