# Python Imports

# Third-Party Imports

# Project Imports
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
                                  frag_standard_arguments_namespace, ID):
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


def create_all_parameters_objects(pele_arguments_namespace, clustering_arguments_namespace, plop_arguments_namespace,
                                  complex_pdb, fragment_pdb, core_atom, fragment_atom, h_core, h_frag,
                                  frag_standard_arguments_namespace):  # todo quit so many arguments
    # Create PELE parameters objects
    pele_parameters = create_pele_parameter_object(pele_arguments_namespace)
    # Create Cluster parameters object
    cluster_parameters = create_cluster_parameters_object(clustering_arguments_namespace)
    # Create Plop parameters object
    plop_parameters = create_plop_parameters_object(plop_arguments_namespace)
    # Create Frag parameters object
    frag_parameters = create_frag_parameters_object(complex_pdb, fragment_pdb, core_atom, fragment_atom,
                                                    h_core, h_frag, frag_standard_arguments_namespace)

    return pele_parameters, cluster_parameters, plop_parameters, frag_parameters
