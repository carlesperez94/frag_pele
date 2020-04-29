# Python Imports
import argparse

# Third-Party Imports

# Project Imports
import frag_pele.constants as const


def _create_parser():
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    text_1 = 'Description: FrAG is a Fragment-based ligand growing software which performs automatically '
    text_2 = 'the addition of several fragments to a core structure of the ligand in a protein-ligand complex.'
    parser = argparse.ArgumentParser(description=text_1 + text_2)

    _add_all_arguments(parser)

    return parser


def _add_all_arguments(parser):
    _add_frag_required_named_arguments(parser)
    _add_frag_standard_arguments(parser)
    _add_plop_arguments(parser)
    _add_pele_conf_arguments(parser)
    _add_clustering_arguments(parser)
    _add_protocol_arguments(parser)
    _add_output_format_arguments(parser)
    _add_others_arguments(parser)


def _add_frag_required_named_arguments(parser):
    required_named = parser.add_argument_group(const.KEY_REQUIRED_ARGUMENTS)

    # FrAG related arguments
    required_named.add_argument("-cp", "--complex_pdb", required=True,
                                help="""Path to the PDB file which must contain a protein-ligand complex. Its ligand 
                                will be used as the core structure. Remember to rename the ligand chain with a 
                                different character in order to detect it.""")

    required_named.add_argument("-sef", "--serie_file", required=True,
                                help=""" Name of the tabular file which must contain the instructions required to 
                                perform several successive growings, using different fragments or different growing 
                                positions.

                        To do simple growings:
                        col1   col2    col3
                        To do successive growings:
                        (col1   col2    col3) x n_growings

                        Where col1 is the path to the PDB file of the fragment that will be added to the core structure
                        (name the chain of the fragment "L" by default).
                        Col2 is a string with the PDB atom name of the heavy atom of the core (the ligand contained in
                        the complex) where you would like to start the growing and create a new bond with the fragment.
                        And col3 is a string with the PDB atom name of the heavy atom of the fragment that will be used
                        to perform the bonding with the core.
                        """)


def _add_frag_standard_arguments(parser):
    frag_group = parser.add_argument_group(const.KEY_FRAG_ARGUMENTS)

    frag_group.add_argument("--core", type=str, default=None)
    frag_group.add_argument("-nc", "--no_check", action="store_true",
                            help="Don't perform the environment variables check")
    frag_group.add_argument("-x", "--growing_steps", type=int, default=const.GROWING_STEPS,
                            help="""Number of Growing Steps (GS). By default = {}.""".format(const.GROWING_STEPS))
    frag_group.add_argument("-cr", "--criteria", default=const.SELECTION_CRITERIA,
                            help="""Name of the column of the report file to select the structures that will spawn in the 
                        next GS. Additionally, this parameter will be the selection criteria to extract the best 
                        structure after completing the growing. By default = {}.""".format(const.SELECTION_CRITERIA))
    frag_group.add_argument("-rst", "--restart", action="store_true",
                            help="If set FrAG will continue from the last GS detected. If all GS are finished it will"
                                 "restart in the equilibration phase.")
    frag_group.add_argument("-cc", "--c_chain", default=const.C_CHAIN, help="Chain name of the core. By default = 'L'")
    frag_group.add_argument("-fc", "--f_chain", default=const.F_CHAIN,
                            help="Chain name of the fragment. By default = 'L'")
    frag_group.add_argument("-tc", "--clash_thr", default=const.CLASH_THRESHOLD,
                            help="Threshold distance that would to classify intramolecular clashes.")
    frag_group.add_argument("-sc", "--sampling_control", default=None,
                            help="If set, templatized control file to use in the sampling simulation.")
    frag_group.add_argument("-op", "--only_prepare", action="store_true",
                            help="If set, all files to run growing are prepared, it stops before running PELE.")
    frag_group.add_argument("-og", "--only_grow", action="store_true",
                            help="If set, it runs all growings of folders already prepared.")
    frag_group.add_argument("-EX", "--explorative", action="store_true",
                            help="Run frag pele explorative mode: sampling simulation with high movement of the ligand.")


def _add_plop_arguments(parser):
    plop_group = parser.add_argument_group(const.KEY_PLOP_ARGUMENTS)

    plop_group.add_argument("-pl", "--plop_path", default=const.PLOP_PATH,
                            help="Absolute path to PlopRotTemp.py. By default = {}".format(const.PLOP_PATH))
    plop_group.add_argument("-sp", "--sch_python", default=const.SCHRODINGER_PY_PATH,
                            help="""Absolute path to Schrodinger's python. 
                        By default = {}""".format(const.SCHRODINGER_PY_PATH))
    plop_group.add_argument("-rot", "--rotamers", default=const.ROTRES, type=str,
                            help="""Rotamers threshold used in the rotamers' library. 
                        By default = {}""".format(const.ROTRES))


def _add_pele_conf_arguments(parser):
    pele_group = parser.add_argument_group(const.KEY_PELE_ARGUMENTS)

    pele_group.add_argument("-d", "--pele_dir", default=const.PATH_TO_PELE,
                            help="Complete path to Pele_serial. "
                                 "By default = {}".format(const.PATH_TO_PELE))
    pele_group.add_argument("-c", "--contrl", default=const.CONTROL_TEMPLATE,
                            help="Path to PELE's control file templatized. "
                                 "By default = {}".format(const.CONTROL_TEMPLATE))
    pele_group.add_argument("-l", "--license", default=const.PATH_TO_LICENSE,
                            help="Absolute path to PELE's licenses folder. "
                                 " By default = {}".format(const.PATH_TO_LICENSE))
    pele_group.add_argument("-r", "--resfold", default=const.RESULTS_FOLDER,
                            help="Name for PELE's results folder. "
                                 "By default = {}".format(const.RESULTS_FOLDER))
    pele_group.add_argument("-rp", "--report", default=const.REPORT_NAME,
                            help="Suffix name of the report file from PELE. "
                                 "By default = {}".format(const.REPORT_NAME))
    pele_group.add_argument("-tj", "--traject", default=const.TRAJECTORY_NAME,
                            help="Suffix name of the trajectory file from PELE. "
                                 "By default = {}".format(const.TRAJECTORY_NAME))
    pele_group.add_argument("-cs", "--cpus", type=int, default=const.CPUS,
                            help="Number of cores (computational power) to paralellize PELE's simulations."
                                 "By default = {}".format(const.CPUS))
    pele_group.add_argument("-stp", "--steps", type=int, default=const.STEPS,
                            help="""Number of simulation steps inside each GS. 
                        By default = {}""".format(const.STEPS))
    pele_group.add_argument("-es", "--pele_eq_steps", default=const.PELE_EQ_STEPS,
                            help="Number of PELE steps in equilibration. "
                                 "By default = {}".format(const.PELE_EQ_STEPS))
    pele_group.add_argument("-miov", "--min_overlap", default=const.MIN_OVERLAP,
                            help="Minimum value of overlapping factor used in the control_file of PELE. "
                                 "By default = {}".format(const.MIN_OVERLAP))
    pele_group.add_argument("-maov", "--max_overlap", default=const.MAX_OVERLAP,
                            help="Maximum value of overlapping factor used in the control_file of PELE."
                                 " By default = {}".format(const.MAX_OVERLAP))
    pele_group.add_argument("-tmp", "--temperature", default=const.TEMPERATURE,
                            help="Temperature value to add in the control file. If the temperature is high more steps of "
                                 "PELE will be accepted when applying the Metropolis Criteria. "
                                 "By default = {}".format(const.TEMPERATURE))
    pele_group.add_argument("-sd", "--seed", default=const.SEED,
                            help="Seed to get the random numbers that will be used in the PELE simulation"
                                 "By default = {}".format(const.SEED))
    pele_group.add_argument("-st", "--steering", default=const.STEERING,
                            help="Steering that will be used in the PELE simulation"
                                 "By default = {}".format(const.STEERING))
    pele_group.add_argument("-trh", "--translation_high", default=const.TRANSLATION_HIGH,
                            help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                                 "By default = {}".format(const.TRANSLATION_HIGH))
    pele_group.add_argument("-roth", "--rotation_high", default=const.ROTATION_HIGH,
                            help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                                 "Low value. By default = {}".format(const.ROTATION_HIGH))
    pele_group.add_argument("-trl", "--translation_low", default=const.TRANSLATION_LOW,
                            help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                                 "By default = {}".format(const.TRANSLATION_LOW))
    pele_group.add_argument("-rotl", "--rotation_low", default=const.ROTATION_LOW,
                            help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                                 "Low value. By default = {}".format(const.ROTATION_LOW))
    pele_group.add_argument("-rad", "--radius_box", default=const.RADIUS_BOX,
                            help="Size of the radius to define the box in the PELE simulation where the ligand will be"
                                 "perturbed. By default = {}".format(const.RADIUS_BOX))
    pele_group.add_argument("-dat", "--data", default=const.PATH_TO_PELE_DATA,
                            help="Path to PELE Data folder.")
    pele_group.add_argument("-doc", "--documents", default=const.PATH_TO_PELE_DOCUMENTS,
                            help="Path to PELE Documents folder.")


def _add_clustering_arguments(parser):
    cluster_group = parser.add_argument_group(const.KEY_CLUSTERING_ARGUMENTS)

    cluster_group.add_argument("-dis", "--distcont", default=const.DISTANCE_COUNTER,
                               help="""Distance used to determine which amino acids are in contact with the ligand to generate 
                        different clusters of structures to initialize the next GS. 
                        By default = {}""".format(const.DISTANCE_COUNTER))
    cluster_group.add_argument("-ct", "--threshold", default=const.CONTACT_THRESHOLD,
                               help="Threshold distance used in the clustering. "
                                    "By default = {}".format(const.CONTACT_THRESHOLD))
    cluster_group.add_argument("-e", "--epsilon", default=const.EPSILON,
                               help="An epsilon fraction of processors are distributed proportionally to the value of a metric"
                                    "and the rest are inverselyProportional distributed. A param n can be specified to only "
                                    "consider the n clusters with best metric. "
                                    "By default = {}".format(const.EPSILON))
    cluster_group.add_argument("-cn", "--condition", default=const.CONDITION,
                               help="Selects wether to take into account maximum or minimum values in epsilon related spawning"
                                    "values are min or max. "
                                    "By default = {}".format(const.CONDITION))
    cluster_group.add_argument("-mw", "--metricweights", default=const.METRICS_WEIGHTS,
                               help="Selects how to distribute the weights of the cluster according to its metric, "
                                    "two options: linear (proportional to metric) or Boltzmann weigths (proportional "
                                    "to exp(-metric/T). Needs to define the temperature T. "
                                    "By default = {}".format(const.METRICS_WEIGHTS))
    cluster_group.add_argument("-ncl", "--nclusters", default=const.NUM_CLUSTERS,
                               help="Number of initial structures that we want to use in each new GS. "
                                    "By default = {}".format(const.NUM_CLUSTERS))
    cluster_group.add_argument("-pdbf", "--pdbout", default=const.PDBS_OUTPUT_FOLDER,
                               help="Folder where PDBs selected to spawn in the next GS will be stored."
                                    "By default = {}".format(const.PDBS_OUTPUT_FOLDER))
    cluster_group.add_argument("-ban", "--banned", default=const.BANNED_DIHEDRALS_ATOMS, type=str, nargs='+',
                               action='append',
                               help="List of lists of quartets of PDB atom names that form the banned dihedrals."
                                    "By default = {}".format(const.BANNED_DIHEDRALS_ATOMS))
    cluster_group.add_argument("-lim", "--limit", default=const.BANNED_ANGLE_THRESHOLD, type=int,
                               help="Maximum degrees that can accept the banned dihedral."
                                    "By default = {}".format(const.BANNED_ANGLE_THRESHOLD))


def _add_protocol_arguments(parser):
    parser.add_argument("-HT", "--highthroughput", action="store_true", help="Run frag pele high-throughput mode")
    parser.add_argument("--test", action="store_true", help="run test config")


def _add_output_format_arguments(parser):
    parser.add_argument("--mae", action="store_true", help="Retrieve .mae files intead of pdbs")


def _add_others_arguments(parser):
    parser.add_argument("--rename", action="store_true", help="Avoid core renaming")


def _check_highthroughput_in_args(args):
    if args.highthroughput:
        args.growing_steps = 1
        args.steps = 3
        args.pele_eq_steps = 10


def _check_test_in_args(args):
    if args.test:
        args.growing_steps = 1
        args.steps = 1
        args.pele_eq_steps = 1
        args.temperature = 1000000


def _extract_groups(args, parser):
    arg_groups = {}

    for group in parser._action_groups:
        group_dict = {a.dest: getattr(args, a.dest, None) for a in group._group_actions}
        arg_groups[group.title] = argparse.Namespace(**group_dict)

    return arg_groups


def parse_arguments():
    """
    Parse user arguments
    Output: dict with all the user arguments
    """
    parser = _create_parser()

    args = parser.parse_args()

    _check_highthroughput_in_args(args)
    _check_test_in_args(args)

    arg_groups = _extract_groups(args, parser)

    return arg_groups


def extract_arguments():
    arguments_dict = parse_arguments()
    required_arguments_namespace = arguments_dict[const.KEY_REQUIRED_ARGUMENTS]
    frag_standard_arguments_namespace = arguments_dict[const.KEY_FRAG_ARGUMENTS]
    plop_arguments_namespace = arguments_dict[const.KEY_PLOP_ARGUMENTS]
    pele_arguments_namespace = arguments_dict[const.KEY_PELE_ARGUMENTS]
    clustering_arguments_namespace = arguments_dict[const.KEY_CLUSTERING_ARGUMENTS]
    mae = arguments_dict[const.KEY_OPTIONAL_ARGUMENTS].mae
    rename = arguments_dict[const.KEY_OPTIONAL_ARGUMENTS].rename

    return required_arguments_namespace, frag_standard_arguments_namespace, plop_arguments_namespace, \
           pele_arguments_namespace, clustering_arguments_namespace, mae, rename
