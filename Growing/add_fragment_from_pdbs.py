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


def extract_heteroatoms_pdbs(pdb_complex_core, pdb_fragment):
    # Parse the complex file and isolate the ligand core and the fragment
    ligand_core = c2s.pdb_parser_ligand(pdb_complex_core, ligand_chain="L")
    fragment = c2s.pdb_parser_ligand(pdb_fragment, ligand_chain="L")
    # Check if the ligand and the fragment has H
    c2s.check_protonation(ligand_core)
    c2s.check_protonation(fragment)
    # Save the ligand and the fragment in a PDB (the name of the file is the name of the residue)
    ligand_core_name = ligand_core.getResnames()[0]
    fragment_name = fragment.getResnames()[0]
    prody.writePDB(ligand_core_name, ligand_core)
    logger.info("The ligand of {} has been extracted and saved in '{}.pdb'".format(pdb_complex_core, ligand_core_name))
    prody.writePDB(fragment_name, fragment)
    return ligand_core_name, fragment_name


def from_pdb_to_bioatomlist(list_of_pdb_names):
    list_of_lists = []
    for pdb in list_of_pdb_names:
        # Now, read the PDBs as BioPython object
        bio_structure = pj.get_ligand_from_PDB("{}.pdb".format(pdb))
        # Get a list of atoms (BioPython Atoms) for each structure
        bioatomlist = pj.get_atoms_from_structure(bio_structure)
        list_of_lists.append(bioatomlist)
    return list_of_lists


def extract_heavy_atoms(pdb_atom_names, lists_of_bioatoms):
    heavy_atoms = []
    # Select the heavy atoms for each list that we will want to bond together in further steps
    for atom_name, list_of_bioatoms in zip(pdb_atom_names, lists_of_bioatoms):
        atom_heavy = pj.select_atoms_from_list(atom_name, list_of_bioatoms)
        heavy_atoms.append(atom_heavy)
    return heavy_atoms


def extract_hydrogens(pdb_atom_names, lists_of_bioatoms, list_of_pdbs):
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
    fragment_coords = pj.transform_coords(bioatom_list)
    fragment_prody.setCoords(fragment_coords)
    return fragment_prody.getCoords()


def bond(hydrogen_atom_names, molecules):
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
    ligand_core = c2s.pdb_parser_ligand(pdb_complex_core, ligand_chain="L")
    fragment = c2s.pdb_parser_ligand(pdb_fragment, ligand_chain="L")
    residue_names = extract_heteroatoms_pdbs(pdb_complex_core, pdb_fragment)
    bioatoms_core_and_frag = from_pdb_to_bioatomlist(residue_names)
    pdb_atom_names = [pdb_atom_core_name, pdb_atom_fragment_name]
    heavy_atoms = extract_heavy_atoms(pdb_atom_names, bioatoms_core_and_frag)
    hydrogen_atoms = extract_hydrogens(pdb_atom_names, bioatoms_core_and_frag, [pdb_complex_core, pdb_fragment])
    # Create a list with the atoms that form a bond in core and fragment
    core_bond = [heavy_atoms[0], hydrogen_atoms[0]]
    print(core_bond)
    fragment_bond = [hydrogen_atoms[1], heavy_atoms[1]]  # This has to be in inverted order
    print(fragment_bond)
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


