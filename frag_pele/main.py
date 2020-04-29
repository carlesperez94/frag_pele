# Python Imports
import os  # Todo: Try to change this to increment portability
import sys
import logging
import traceback
from logging.config import fileConfig

# Third-Party Imports

# Project Imports
from frag_pele.Helpers import correct_fragment_names
from frag_pele import serie_handler
import frag_pele.constants as c
from frag_pele.Helpers.argument_parser import extract_arguments
from frag_pele.frag.frag_runner import Frag
from frag_pele.frag.Parameters.parameter_instanciator import create_all_parameters_objects

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "AdaptivePELE_repo"))

# Calling configuration file for log system
FilePath = os.path.abspath(__file__)
PackagePath = os.path.dirname(FilePath)
LogPath = os.path.join(PackagePath, c.CONFIG_PATH)
fileConfig(LogPath)

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Get current path
curr_dir = os.path.abspath(os.path.curdir)

if __name__ == '__main__':
    # Create Frag object
    frag_object = Frag()
    # HOT FIX!! Fix it properly
    original_dir = os.path.abspath(os.getcwd())

    # Sacamos los argumentos
    required_arguments_namespace, frag_standard_arguments_namespace, plop_arguments_namespace, \
    pele_arguments_namespace, clustering_arguments_namespace, mae, rename = extract_arguments()

    # Se sacan instrucciones
    list_of_instructions = serie_handler.read_instructions_from_file(required_arguments_namespace.serie_file)
    print("READING INSTRUCTIONS... You will perform the growing of {} fragments. GOOD LUCK and ENJOY the "
          "trip :)".format(len(list_of_instructions)))
    dict_traceback = correct_fragment_names.main(required_arguments_namespace.complex_pdb)

    # We will iterate trough all individual instructions of file.
    for instruction in list_of_instructions:
        # SUCCESSIVE GROWING
        # If in the individual instruction we have more than one command means successive growing.
        if type(instruction) == list:
            # Doing so we will determinate how many successive growings the user wants to do.
            growing_counter = len(instruction)
            atomname_mappig = []
            for i in range(int(growing_counter)):
                core_from_previous_fragment = instruction[i][-1]
                fragment_pdb, core_atom, fragment_atom = instruction[i][0], instruction[i][1], instruction[i][2]
                atoms_if_bond = serie_handler.extract_hydrogens_from_instructions(
                    [fragment_pdb, core_atom, fragment_atom])

                if atoms_if_bond:
                    core_atom = atoms_if_bond[0]
                    if core_from_previous_fragment and i != 0:
                        previous_fragment_atomnames_map = atomname_mappig[int(instruction[i][4]) - 1]
                        core_atom = previous_fragment_atomnames_map[core_atom]
                        h_core = previous_fragment_atomnames_map[atoms_if_bond[1]]
                    else:
                        h_core = atoms_if_bond[1]
                    fragment_atom = atoms_if_bond[2]
                    h_frag = atoms_if_bond[3]
                else:
                    if core_from_previous_fragment and i != 0:
                        previous_fragment_atomnames_map = atomname_mappig[int(instruction[i][4]) - 1]
                        core_atom = previous_fragment_atomnames_map[core_atom]
                    h_core = None
                    h_frag = None

                # In the first iteration we will use the complex_pdb as input.
                if i == 0:
                    ID = instruction[i][3]
                else:  # If is not the first we will use as input the output of the previous iteration
                    pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
                    if "/" in pdb_basename:
                        pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name
                    working_dir = "{}_{}".format(pdb_basename, ID)
                    complex_pdb = os.path.join(working_dir, c.PRE_WORKING_DIR)
                    dict_traceback = correct_fragment_names.main(complex_pdb)
                    ID_completed = []
                    for id in instruction[0:i + 1]:
                        ID_completed.append(id[3])
                    ID = "".join(ID_completed)
                try:
                    ID = ID.split("/")[-1]
                except Exception:
                    traceback.print_exc()
                try:
                    serie_handler.check_instructions(instruction[i], required_arguments_namespace.complex_pdb,
                                                     frag_standard_arguments_namespace.c_chain,
                                                     frag_standard_arguments_namespace.f_chain)
                    # Create parameters objects
                    # todo quit so many arguments
                    pele_parameters, cluster_parameters, plop_parameters, \
                    frag_parameters = create_all_parameters_objects(pele_arguments_namespace,
                                                                    clustering_arguments_namespace,
                                                                    plop_arguments_namespace, complex_pdb, fragment_pdb,
                                                                    core_atom, fragment_atom, h_core, h_frag,
                                                                    frag_standard_arguments_namespace, ID)



                    print("PERFORMING SUCCESSIVE GROWING...")
                    print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                    atomname_map = frag_object.run_frag(pele_parameters, cluster_parameters, plop_parameters,
                                                        frag_parameters, mae, rename)
                    atomname_mappig.append(atomname_map)
                except Exception:
                    os.chdir(original_dir)
                    traceback.print_exc()
        # INDIVIDUAL GROWING
        else:
            # Initialize the growing for each line in the file
            fragment_pdb, core_atom, fragment_atom, ID = instruction[0], instruction[1], instruction[2], instruction[3]
            atoms_if_bond = serie_handler.extract_hydrogens_from_instructions([fragment_pdb, core_atom, fragment_atom])
            try:
                ID = ID.split("/")[-1]
            except Exception:
                traceback.print_exc()

            if atoms_if_bond:
                core_atom = atoms_if_bond[0]
                h_core = atoms_if_bond[1]
                fragment_atom = atoms_if_bond[2]
                h_frag = atoms_if_bond[3]
            else:
                h_core = None
                h_frag = None
            try:
                serie_handler.check_instructions(instruction[i], complex_pdb, frag_standard_arguments_namespace.c_chain,
                                                 frag_standard_arguments_namespace.f_chain)
                # Create parameters objects
                pele_parameters, cluster_parameters, plop_parameters, \
                frag_parameters = create_all_parameters_objects(pele_arguments_namespace,
                                                                clustering_arguments_namespace,
                                                                plop_arguments_namespace, complex_pdb, fragment_pdb,
                                                                core_atom, fragment_atom, h_core, h_frag,
                                                                frag_standard_arguments_namespace, ID)
                print("PERFORMING INDIVIDUAL GROWING...")
                print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                frag_object.run_frag(pele_parameters, cluster_parameters, plop_parameters,
                                     frag_parameters)
            except Exception:
                os.chdir(original_dir)
                traceback.print_exc()
