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


def put_fragment_in_place(pdb_complex_core, pdb_fragment, pdb_atom_core_name, pdb_atom_fragment_name,
                          output_file="growing_result.pdb"):
    protein_complex = prody.parsePDB(pdb_complex_core)
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
    # Now, read the PDBs as BioPython object
    core_bio_structure = pj.get_ligand_from_PDB("{}.pdb".format(ligand_core_name))
    frag_bio_structure = pj.get_ligand_from_PDB("{}.pdb".format(fragment_name))
    # Get a list of atoms (BioPython Atoms) for each structure
    atoms_core_bio = pj.get_atoms_from_structure(core_bio_structure)
    atoms_frag_bio = pj.get_atoms_from_structure(frag_bio_structure)
    # Select the heavy atoms for each list that we will want to bond together in further steps
    atom_core_heavy = pj.select_atoms_from_list(pdb_atom_core_name, atoms_core_bio)
    atom_frag_heavy = pj.select_atoms_from_list(pdb_atom_fragment_name, atoms_frag_bio)
    # Select name of the H atoms bonded to this heavy atom (the place where we will grow)
    atom_name_core_hydrogen = pj.get_H_bonded_to_grow(pdb_atom_core_name, protein_complex)
    atom_name_frag_hydrogen = pj.get_H_bonded_to_grow(pdb_atom_fragment_name, fragment)
    # Select this hydrogen atoms
    atom_core_hydrogen = pj.select_atoms_from_list(atom_name_core_hydrogen, atoms_core_bio)
    atom_frag_hydrogen = pj.select_atoms_from_list(atom_name_frag_hydrogen, atoms_frag_bio)
    # Create a list with the atoms that form a bond in core and fragment
    core_bond = [atom_core_heavy, atom_core_hydrogen]
    fragment_bond = [atom_frag_hydrogen, atom_frag_heavy]  # This has to be in inverted order
    logger.info("Performing a superimposition of bond {} of the fragment on bond {} of the core..."
                .format(fragment_bond, core_bond))
    # Do the superimposition with BioPython
    sup = bio.Superimposer()
    # Set the vectors: first element is the fix vector (bond of the core) and second is the moving (bond of the fragment)
    sup.set_atoms(core_bond, fragment_bond)
    # Apply the transformation to the atoms of the fragment (translate and rotate)
    sup.apply(atoms_frag_bio)
    # Get the new coords and change them
    fragment_coords = pj.transform_coords(atoms_frag_bio)
    fragment.setCoords(fragment_coords)
    # Now, we have to remove the hydrogens of the binding, so we will select everything except these H
    core_no_H = ligand_core.select("not name {}".format(atom_name_core_hydrogen))
    fragment_no_H = fragment.select("not name {}".format(atom_name_frag_hydrogen))
    # Merging both parts into a single one
    merged = core_no_H.copy() + fragment_no_H.copy()
    # Change the resnames and resnums to set everything as a single molecule
    merged.setResnames("GRW")
    merged.setResnums(1)
    prody.writePDB(output_file, merged)
    logger.info("The result of core + fragment has been saved in '{}'".format(output_file))


