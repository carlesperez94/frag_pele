import prody
import pybel
import logging
import sys
# Local imports
import complex_to_smile as c2s
import smiles_joiner as sj
import smile_to_pdb as s2p


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def add_fragment_to_pdb(pdb_file, fragment_smi, atom_type_where_grow, atom_index_where_grow, output_filename,
                        ligand_chain="L"):
    """
    This function uses a fragment in smile format to add it to a core ligand of a PDB file.
    :param pdb_file: input pdb with the ligand that will be used as core.
    :param fragment_smi: fragment in smile format that will be added to the core. Take into account that you can not use
    canonical smiles because the first atom will be bonded to the core (probably solved in further updates).
    :param atom_type_where_grow: string with the type of heavy atom (C, O, N, S, etc) that you want to select to do the
    growing.
    :param atom_index_where_grow: number of atom the specific atom where you want to grow once you have selected the
    type (p.ex: I want to grow the Carbon 9, so in atom_type_where_grow I will put a "C" and here I will put a "9".
    :param output_filename: name of the output pdb file.
    :param ligand_chain: chain where the ligand is in the PDB input. "L" by default.
    :return: PDB file of the molecule with the fragment added inplace.
    """
    # Read PDB and select the ligand chain
    core_from_pdb = c2s.pdb_parser_ligand(pdb_file, ligand_chain)
    # Creating a PDB with only the ligand
    core_pdb_file = c2s.selection_to_pdb(core_from_pdb)
    logger.info("{} has been created".format(core_pdb_file))
    # Transform this PDB to a SMILE code
    core_pdb_smile = c2s.pdb_to_smile(core_pdb_file)
    logger.info("CORE SMILE: {}".format(core_pdb_smile))
    # Insert this smile to get a pybel object
    core_pybel_smile = sj.read_smi(core_pdb_smile)

    # Obtain the atomic number of the atom type (C, N, O, S, etc.) where we want to grow
    atomicnum_of_selected_atom = sj.get_atomicnum(atom_type_where_grow)
    # Transform the smile into a list of (index, atomic number) in order to find the correct place to grow
    core_pybel_listed = sj.list_of_atomicnum(core_pybel_smile)
    # Select of the previous list which atom has the type and the index that we want in order to find the place to grow
    selected_atom = sj.find_atomic_index(core_pybel_listed, atomicnum_of_selected_atom, atom_index_where_grow)

    # Transform the smile of the core to a list in order to add the fragment properly
    core_smile_as_list = sj.smile_to_list(core_pdb_smile)
    # Check that the previous list has the same length than the list of atomic numbers (that we know that is right)
    sj.check_smile_len(core_pybel_listed, core_smile_as_list)

    # Add the fragment in smile format to the core smile in the selected place
    new_smile = sj.add_fragment(fragment_smi, core_smile_as_list, selected_atom)
    logger.info("RESULTANT SMILE: {}".format(new_smile))
    # Now, transform it in pybel object
    new_molecule = sj.read_smi(new_smile)
    # Add hydrogens and 3D parameters
    s2p.molecule_to_h_and_3d(new_molecule)
    logger.info("Adding Hydrogens and 3D parameters...")
    # Get the pdb of the resultant structure
    s2p.molecule_to_pdb(new_molecule, output_filename)
    logger.info("The final molecule has been saved in {}".format(output_filename))
