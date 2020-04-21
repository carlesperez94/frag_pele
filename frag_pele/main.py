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
from frag_pele.Helpers.argument_parser import parse_arguments
from frag_pele.frag.PeleParameters.pele_parameter_archives import PeleParameterArchives
from frag_pele.frag.PeleParameters.pele_parameter_paths import PeleParameterPaths
from frag_pele.frag.PeleParameters.pele_parameter_sim_values import PeleParameterSimulationValues
from frag_pele.frag.PeleParameters.pele_parameters import PeleParameters
from frag_pele.frag.frag_runner import Frag

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
    #HOT FIX!! Fix it properly
    original_dir = os.path.abspath(os.getcwd())

    # Get arguments
    complex_pdb, serie_file, no_check, growing_steps, criteria, restart, c_chain, f_chain, threshold_clash, sampling_control,\
    only_prepare, only_grow, plop_path, sch_python, rotamers, pele_dir, contrl, license, resfold, report, traject, \
    cpus, steps, pele_eq_steps, min_overlap, max_overlap, temperature, seed,steering, translation_high, rotation_high, \
    translation_low, rotation_low, radius_box, data, documents, distcont, threshold, epsilon, condition, metricweights,\
    nclusters, pdbout, banned, limit, explorative, mae, rename = parse_arguments()

    list_of_instructions = serie_handler.read_instructions_from_file(serie_file)
    print("READING INSTRUCTIONS... You will perform the growing of {} fragments. GOOD LUCK and ENJOY the "
          "trip :)".format(len(list_of_instructions)))
    dict_traceback = correct_fragment_names.main(complex_pdb)
    for instruction in list_of_instructions:
        # We will iterate trough all individual instructions of file.
        # SUCCESSIVE GROWING
        if type(instruction) == list:  #  If in the individual instruction we have more than one command means successive growing.
            growing_counter = len(instruction)  #  Doing so we will determinate how many successive growings the user wants to do.
            atomname_mappig = []
            for i in range(int(growing_counter)):
                core_from_previous_fragment = instruction[i][-1]
                fragment_pdb, core_atom, fragment_atom = instruction[i][0], instruction[i][1], instruction[i][2]
                atoms_if_bond = serie_handler.extract_hydrogens_from_instructions([fragment_pdb, core_atom, fragment_atom])
                if atoms_if_bond:
                    core_atom = atoms_if_bond[0]
                    if core_from_previous_fragment and i!=0:
                        previous_fragment_atomnames_map = atomname_mappig[int(instruction[i][4])-1]
                        core_atom = previous_fragment_atomnames_map[core_atom]
                        h_core = previous_fragment_atomnames_map[atoms_if_bond[1]]
                    else:
                        h_core = atoms_if_bond[1]
                    fragment_atom = atoms_if_bond[2]
                    h_frag = atoms_if_bond[3]
                else:
                    if core_from_previous_fragment and i!=0:
                        previous_fragment_atomnames_map = atomname_mappig[int(instruction[i][4])-1]
                        core_atom = previous_fragment_atomnames_map[core_atom]
                    h_core = None
                    h_frag = None
                if i == 0:  # In the first iteration we will use the complex_pdb as input.
                    ID = instruction[i][3]
                else:  # If is not the first we will use as input the output of the previous iteration
                    pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
                    if "/" in pdb_basename:
                        pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name
                    working_dir = "{}_{}".format(pdb_basename, ID)
                    complex_pdb = os.path.join(working_dir, c.PRE_WORKING_DIR)
                    dict_traceback = correct_fragment_names.main(complex_pdb)
                    ID_completed = []
                    for id in instruction[0:i+1]:
                        ID_completed.append(id[3])
                    ID = "".join(ID_completed)
                try:
                    ID = ID.split("/")[-1]
                except Exception:
                    traceback.print_exc()
                try:
                    serie_handler.check_instructions(instruction[i], complex_pdb, c_chain, f_chain)
                    # Create PELE objects
                    pele_params_paths = PeleParameterPaths(pele_dir, license, data, documents)
                    pele_params_archives = PeleParameterArchives(contrl, resfold, report, traject)
                    pele_params_sim_values = PeleParameterSimulationValues(cpus, steps, pele_eq_steps, min_overlap,
                                                                           max_overlap, temperature, seed, steering, translation_high,
                                                                           rotation_high, translation_low, rotation_low, radius_box)
                    pele_parameters = PeleParameters(pele_params_paths, pele_params_archives, pele_params_sim_values)

                    print("PERFORMING SUCCESSIVE GROWING...")
                    print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                    atomname_map = frag_object.main(pele_parameters, complex_pdb, fragment_pdb, core_atom, fragment_atom, growing_steps, criteria,
                                                    plop_path, sch_python, pdbout, distcont, threshold, epsilon, condition, metricweights,
                                                    nclusters, restart, ID, h_core, h_frag, c_chain, f_chain, rotamers,
                                                    banned, limit, mae, rename, threshold_clash, explorative,
                                                    sampling_control, only_prepare, only_grow, no_check)
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
                # Create PELE objects
                pele_params_paths = PeleParameterPaths(pele_dir, license, data, documents)
                pele_params_archives = PeleParameterArchives(contrl, resfold, report, traject)
                pele_params_sim_values = PeleParameterSimulationValues(cpus, steps, pele_eq_steps, min_overlap,
                                                                       max_overlap, temperature, seed, steering,
                                                                       translation_high,
                                                                       rotation_high, translation_low, rotation_low,
                                                                       radius_box)
                pele_parameters = PeleParameters(pele_params_paths, pele_params_archives, pele_params_sim_values)

                print("PERFORMING INDIVIDUAL GROWING...")
                print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                frag_object.main(pele_parameters, complex_pdb, fragment_pdb, core_atom, fragment_atom, growing_steps, criteria,
                                 plop_path, sch_python, pdbout, distcont, threshold, epsilon,
                     condition, metricweights, nclusters, restart, ID, h_core,
                     h_frag, c_chain, f_chain, rotamers, banned, limit, mae, rename,
                     threshold_clash, explorative, sampling_control, only_prepare, only_grow, no_check)
            except Exception:
                os.chdir(original_dir)
                traceback.print_exc()
