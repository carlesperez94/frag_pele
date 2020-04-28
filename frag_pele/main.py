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
from frag_pele.frag.ClusterParameters.cluster_parameters import ClusterParameters
from frag_pele.frag.FragParameters.frag_configuration_parameters import FragConfigurationParameters
from frag_pele.frag.FragParameters.frag_identification_parameters import FragIdentificationParameters
from frag_pele.frag.FragParameters.frag_parameters import FragParameters
from frag_pele.frag.FragParameters.frag_protocol_parameters import FragProtocolParameters
from frag_pele.frag.FragParameters.frag_running_modes_parameters import FragRunningModesParameters
from frag_pele.frag.FragParameters.frag_structural_configuration_parameters import FragStructuralConfigurationParameters
from frag_pele.frag.FragParameters.frag_structural_files_parameters import FragStructuralFilesParameters
from frag_pele.frag.PeleParameters.pele_parameter_archives import PeleParameterArchives
from frag_pele.frag.PeleParameters.pele_parameter_paths import PeleParameterPaths
from frag_pele.frag.PeleParameters.pele_parameter_sim_values import PeleParameterSimulationValues
from frag_pele.frag.PeleParameters.pele_parameters import PeleParameters
from frag_pele.frag.PlopParameters.plop_parameters import PlopParameters
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


def create_pele_parameter_object(pele_arguments_namespace):

    pele_params_paths = PeleParameterPaths(pele_arguments_namespace.pele_dir, pele_arguments_namespace.license,
                                           pele_arguments_namespace.data, pele_arguments_namespace.documents)
    pele_params_archives = PeleParameterArchives(pele_arguments_namespace.contrl, pele_arguments_namespace.resfold,
                                                 pele_arguments_namespace.report, pele_arguments_namespace.traject)
    pele_params_sim_values = PeleParameterSimulationValues(pele_arguments_namespace.cpus, pele_arguments_namespace.steps
                                                           , pele_arguments_namespace.pele_eq_steps,
                                                           pele_arguments_namespace.min_overlap,
                                                           pele_arguments_namespace.max_overlap,
                                                           pele_arguments_namespace.temperature,
                                                           pele_arguments_namespace.seed,
                                                           pele_arguments_namespace.steering,
                                                           pele_arguments_namespace.translation_high,
                                                           pele_arguments_namespace.rotation_high,
                                                           pele_arguments_namespace.translation_low,
                                                           pele_arguments_namespace.rotation_low,
                                                           pele_arguments_namespace.radius_box)
    pele_parameters = PeleParameters(pele_params_paths, pele_params_archives, pele_params_sim_values)

    return pele_parameters


def create_cluster_parameters_object(clustering_arguments_namespace):
    cluster_parameters = ClusterParameters(clustering_arguments_namespace.distcont,
                                           clustering_arguments_namespace.cluster_threshold,
                                           clustering_arguments_namespace.epsilon,
                                           clustering_arguments_namespace.condition,
                                           clustering_arguments_namespace.metricweights,
                                           clustering_arguments_namespace.nclusters,
                                           clustering_arguments_namespace.pdbout,
                                           clustering_arguments_namespace.banned,
                                           clustering_arguments_namespace.limit)

    return cluster_parameters


def create_plop_parameters_object(plop_arguments_namespace):
    plop_parameters = PlopParameters(plop_arguments_namespace.plop_path,
                                     plop_arguments_namespace.sch_python,
                                     plop_arguments_namespace.rotamers)

    return plop_parameters


def create_frag_parameters_object(complex_pdb, fragment_pdb, core_atom, fragment_atom, h_core, h_frag,
                                  frag_standard_arguments_namespace):
    structural_files_parameters = FragStructuralFilesParameters(complex_pdb, fragment_pdb)
    structural_configuration_parameters = FragStructuralConfigurationParameters(core_atom, fragment_atom,
                                                                                h_core, h_frag,
                                                                                frag_standard_arguments_namespace.c_chain,
                                                                                frag_standard_arguments_namespace.f_chain,
                                                                                frag_standard_arguments_namespace.threshold_clash)
    identification_parameters = FragIdentificationParameters(ID)
    configuration_parameters = FragConfigurationParameters(frag_standard_arguments_namespace.growing_steps,
                                                           frag_standard_arguments_namespace.criteria,
                                                           frag_standard_arguments_namespace.sampling_control)
    protocol_parameters = FragProtocolParameters(frag_standard_arguments_namespace.only_prepare,
                                                 frag_standard_arguments_namespace.only_grow,
                                                 frag_standard_arguments_namespace.explorative)
    running_modes_parameters = FragRunningModesParameters(frag_standard_arguments_namespace.restart,
                                                          frag_standard_arguments_namespace.no_check)
    frag_parameters = FragParameters(structural_files_parameters, structural_configuration_parameters,
                                     identification_parameters, configuration_parameters,
                                     protocol_parameters, running_modes_parameters)

    return frag_parameters


if __name__ == '__main__':
    # Create Frag object
    frag_object = Frag()
    # HOT FIX!! Fix it properly
    original_dir = os.path.abspath(os.getcwd())

    # Get arguments
    arguments_dict = parse_arguments()
    required_arguments_namespace = arguments_dict['required named arguments']  # todo change this keys to constants
    frag_standard_arguments_namespace = arguments_dict['Frag Standard Arguments']
    plop_arguments_namespace = arguments_dict['Plop Related Arguments']
    pele_arguments_namespace = arguments_dict['PELE Related Arguments']
    clustering_arguments_namespace = arguments_dict['Clustering Related Arguments']
    mae = arguments_dict['optional arguments'].mae
    rename = arguments_dict['optional arguments'].rename

    list_of_instructions = serie_handler.read_instructions_from_file(required_arguments_namespace.serie_file)
    print("READING INSTRUCTIONS... You will perform the growing of {} fragments. GOOD LUCK and ENJOY the "
          "trip :)".format(len(list_of_instructions)))
    dict_traceback = correct_fragment_names.main(required_arguments_namespace.complex_pdb)
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
                    serie_handler.check_instructions(instruction[i], required_arguments_namespace.complex_pdb,
                                                     frag_standard_arguments_namespace.c_chain,
                                                     frag_standard_arguments_namespace.f_chain)
                    # Create PELE parameters objects
                    pele_parameters = create_pele_parameter_object(pele_arguments_namespace)
                    # Create Cluster parameters object
                    cluster_parameters = create_cluster_parameters_object(clustering_arguments_namespace)
                    # Create Plop parameters object
                    plop_parameters = create_plop_parameters_object(plop_arguments_namespace)
                    # Create Frag parameters object
                    frag_parameters = create_frag_parameters_object(complex_pdb, fragment_pdb, core_atom, fragment_atom,
                                                                    h_core, h_frag, frag_standard_arguments_namespace)

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
                serie_handler.check_instructions(instruction[i], complex_pdb, c_chain, f_chain)
                # Create PELE parameters objects
                pele_parameters = create_pele_parameter_object(pele_arguments_namespace)
                # Create Cluster parameters object
                cluster_parameters = create_cluster_parameters_object(clustering_arguments_namespace)
                # Create Plop parameters object
                plop_parameters = create_plop_parameters_object(plop_arguments_namespace)
                # Create Frag parameters object
                frag_parameters = create_frag_parameters_object(complex_pdb, fragment_pdb, core_atom, fragment_atom,
                                                                h_core, h_frag, frag_standard_arguments_namespace)
                print("PERFORMING INDIVIDUAL GROWING...")
                print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                frag_object.run_frag(pele_parameters, cluster_parameters, plop_parameters,
                                     frag_parameters)
            except Exception:
                os.chdir(original_dir)
                traceback.print_exc()
