# General imports
import sys
import time
import glob
import argparse
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "AdaptivePELE_repo"))
import logging
from logging.config import fileConfig
import shutil
import subprocess
import traceback
# Local imports
from frag_pele.Helpers import clusterizer, checker, folder_handler, runner, constraints, check_constants
from frag_pele.Helpers import helpers, correct_fragment_names, center_of_mass
from frag_pele.Growing import template_fragmenter, simulations_linker
from frag_pele.Growing import add_fragment_from_pdbs, bestStructs
from frag_pele.Analysis import analyser
from frag_pele.Banner import Detector as dt
from frag_pele import serie_handler
import frag_pele.constants as c

# Calling configuration file for log system
FilePath = os.path.abspath(__file__)
PackagePath = os.path.dirname(FilePath)
LogPath = os.path.join(PackagePath, c.CONFIG_PATH)
fileConfig(LogPath)

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Get current path
curr_dir = os.path.abspath(os.path.curdir)


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """
    # All the docstrings are very provisional and some of them are old, they would be changed in further steps!!
    parser = argparse.ArgumentParser(description="""Description: FrAG is a Fragment-based ligand growing software which
    performs automatically the addition of several fragments to a core structure of the ligand in a protein-ligand
    complex.""")
    required_named = parser.add_argument_group('required named arguments')

    # FrAG related arguments
    required_named.add_argument("-cp", "--complex_pdb", required=True,
                        help="""Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
                        used as the core structure. Remember to rename the ligand chain with a different character in
                        order to detect it.""")
    required_named.add_argument("-sef", "--serie_file", required=True,
                        help=""" Name of the tabular file which must contain the instructions required to perform several
                        successive growings, using different fragments or different growing positions.
                        
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
    parser.add_argument("--debug", action="store_true", help="Run Frag without launching PELE simulation")
    parser.add_argument("-nc", "--no_check", action="store_true", help="Don't perform the environment variables check")
    parser.add_argument("-x", "--growing_steps", type=int, default=c.GROWING_STEPS,
                        help="""Number of Growing Steps (GS). By default = {}.""".format(c.GROWING_STEPS))
    parser.add_argument("-cr", "--criteria", default=c.SELECTION_CRITERIA,
                        help="""Name of the column of the report file to select the structures that will spawn in the 
                        next GS. Additionally, this parameter will be the selection criteria to extract the best 
                        structure after completing the growing. By default = {}.""".format(c.SELECTION_CRITERIA))
    parser.add_argument("-rst", "--restart", action="store_true",
                        help="If set FrAG will continue from the last GS detected. If all GS are finished it will"
                             "restart in the equilibration phase.")
    parser.add_argument("-cc", "--c_chain", default="L", help="Chain name of the core. By default = 'L'")
    parser.add_argument("-fc", "--f_chain", default="L", help="Chain name of the fragment. By default = 'L'")
    parser.add_argument("-tc", "--clash_thr", default=1.7, help="Threshold distance that would to classify intramolecular"
                                                                "clashes.")
    parser.add_argument("-sc",  "--sampling_control", default=None, help="If set, templatized control file to use in the"
                                                                         " sampling simulation.")
    parser.add_argument("-op",  "--only_prepare", action="store_true", help="If set, all files to run growing are"
                                                                            " prepared, it stops before running PELE.")
    parser.add_argument("-og", "--only_grow", action="store_true", help="If set, it runs all growings of folders "
                                                                        "already prepared.")

    # Plop related arguments
    parser.add_argument("-pl", "--plop_path", default=c.PLOP_PATH,
                        help="Absolute path to PlopRotTemp.py. By default = {}".format(c.PLOP_PATH))
    parser.add_argument("-sp", "--sch_python", default=c.SCHRODINGER_PY_PATH,
                        help="""Absolute path to Schrodinger's python. 
                        By default = {}""".format(c.SCHRODINGER_PY_PATH))
    parser.add_argument("-rot", "--rotamers", default=c.ROTRES, type=str,
                        help="""Rotamers threshold used in the rotamers' library. 
                            By default = {}""".format(c.ROTRES))


    # PELE configuration arguments
    parser.add_argument("-d", "--pele_dir", default=c.PATH_TO_PELE,
                        help="Complete path to Pele_serial. "
                             "By default = {}".format(c.PATH_TO_PELE))
    parser.add_argument("-c", "--contrl", default=c.CONTROL_TEMPLATE,
                        help="Path to PELE's control file templatized. By default = {}".format(c.CONTROL_TEMPLATE))
    parser.add_argument("-l", "--license", default=c.PATH_TO_LICENSE,
                        help="Absolute path to PELE's licenses folder. "
                             " By default = {}".format(c.PATH_TO_LICENSE))
    parser.add_argument("-r", "--resfold", default=c.RESULTS_FOLDER,
                        help="Name for PELE's results folder. By default = {}".format(c.RESULTS_FOLDER))
    parser.add_argument("-rp", "--report", default=c.REPORT_NAME,
                        help="Suffix name of the report file from PELE. By default = {}".format(c.REPORT_NAME))
    parser.add_argument("-tj", "--traject", default=c.TRAJECTORY_NAME,
                        help="Suffix name of the trajectory file from PELE. By default = {}".format(c.TRAJECTORY_NAME))
    parser.add_argument("-cs", "--cpus", type=int, default=c.CPUS,
                        help="Number of cores (computational power) to paralellize PELE's simulations."
                             "By default = {}".format(c.CPUS))
    parser.add_argument("-stp", "--steps", type=int, default=c.STEPS,
                        help="""Number of simulation steps inside each GS. By default = {}""".format(c.STEPS))
    parser.add_argument("-es", "--pele_eq_steps", default=c.PELE_EQ_STEPS,
                        help="Number of PELE steps in equilibration. By default = {}".format(c.PELE_EQ_STEPS))
    parser.add_argument("-miov", "--min_overlap", default=c.MIN_OVERLAP,
                        help="Minimum value of overlapping factor used in the control_file of PELE. "
                             "By default = {}".format(c.MIN_OVERLAP))
    parser.add_argument("-maov", "--max_overlap", default=c.MAX_OVERLAP,
                        help="Maximum value of overlapping factor used in the control_file of PELE."
                             " By default = {}".format(c.MAX_OVERLAP))
    parser.add_argument("-tmp", "--temperature", default=c.TEMPERATURE,
                        help="Temperature value to add in the control file. If the temperature is high more steps of "
                             "PELE will be accepted when applying the Metropolis Criteria. "
                             "By default = {}".format(c.TEMPERATURE))
    parser.add_argument("-sd", "--seed", default=c.SEED,
                        help="Seed to get the random numbers that will be used in the PELE simulation"
                             "By default = {}".format(c.SEED))
    parser.add_argument("-st", "--steering", default=c.STEERING,
                        help="Steering that will be used in the PELE simulation"
                             "By default = {}".format(c.STEERING))
    parser.add_argument("-trh", "--translation_high", default=c.TRANSLATION_HIGH,
                        help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                             "By default = {}".format(c.TRANSLATION_HIGH))
    parser.add_argument("-roth", "--rotation_high", default=c.ROTATION_HIGH,
                        help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                             "Low value. By default = {}".format(c.ROTATION_HIGH))
    parser.add_argument("-trl", "--translation_low", default=c.TRANSLATION_LOW,
                        help="Translation that will be used in the PELE simulation to displace the ligand. Low value."
                             "By default = {}".format(c.TRANSLATION_LOW))
    parser.add_argument("-rotl", "--rotation_low", default=c.ROTATION_LOW,
                        help="Rad of rotation that will be used in the PELE simulation to perturb the ligand. "
                             "Low value. By default = {}".format(c.ROTATION_LOW))
    parser.add_argument("-rad", "--radius_box", default=c.RADIUS_BOX,
                        help="Size of the radius to define the box in the PELE simulation where the ligand will be"
                             "perturbed. By default = {}".format(c.RADIUS_BOX))
    parser.add_argument("-dat", "--data", default=c.PATH_TO_PELE_DATA,
                        help="Path to PELE Data folder.")
    parser.add_argument("-doc", "--documents", default=c.PATH_TO_PELE_DOCUMENTS,
                        help="Path to PELE Documents folder.")

    # Clustering related arguments
    parser.add_argument("-dis", "--distcont", default=c.DISTANCE_COUNTER,
                        help="""Distance used to determine which amino acids are in contact with the ligand to generate 
                        different clusters of structures to initialize the next GS.
                        By default = {}""".format(c.DISTANCE_COUNTER))
    parser.add_argument("-ct", "--threshold", default=c.CONTACT_THRESHOLD,
                        help="Threshold distance used in the clustering. By default = {}".format(c.CONTACT_THRESHOLD))
    parser.add_argument("-e", "--epsilon", default=c.EPSILON,
                        help="An epsilon fraction of processors are distributed proportionally to the value of a metric,"
                             "and the rest are inverselyProportional distributed. A param n can be specified to only "
                             "consider the n clusters with best metric. By default = {}".format(c.EPSILON))
    parser.add_argument("-cn", "--condition", default=c.CONDITION,
                        help="Selects wether to take into account maximum or minimum values in epsilon related spawning,"
                             "values are min or max. By default = {}".format(c.CONDITION))
    parser.add_argument("-mw", "--metricweights", default=c.METRICS_WEIGHTS,
                        help="Selects how to distribute the weights of the cluster according to its metric, "
                             "two options: linear (proportional to metric) or Boltzmann weigths (proportional "
                             "to exp(-metric/T). Needs to define the temperature T. "
                             "By default = {}".format(c.METRICS_WEIGHTS))
    parser.add_argument("-ncl", "--nclusters", default=c.NUM_CLUSTERS,
                        help="Number of initial structures that we want to use in each new GS. "
                             "By default = {}".format(c.NUM_CLUSTERS))

    parser.add_argument("-pdbf", "--pdbout", default=c.PDBS_OUTPUT_FOLDER,
                        help="Folder where PDBs selected to spawn in the next GS will be stored."
                             "By default = {}".format(c.PDBS_OUTPUT_FOLDER))
    parser.add_argument("-ban", "--banned", default=c.BANNED_DIHEDRALS_ATOMS, type=str, nargs='+', action='append',
                        help="List of lists of quartets of PDB atom names that form the banned dihedrals."
                             "By default = {}".format(c.BANNED_DIHEDRALS_ATOMS))
    parser.add_argument("-lim", "--limit", default=c.BANNED_ANGLE_THRESHOLD, type=int,
                        help="Maximum degrees that can accept the banned dihedral."
                             "By default = {}".format(c.BANNED_ANGLE_THRESHOLD))

    #Protocol argument
    parser.add_argument("-HT", "--highthroughput", action="store_true",
                        help="Run frag pele high-throughput mode")

    parser.add_argument("-EX", "--explorative", action="store_true",
                        help="Run frag pele explorative mode: sampling simulation with high movement of the ligand.")

    parser.add_argument("--test", action="store_true", help="run test config")

    #Output format option
    parser.add_argument("--mae", action="store_true",
                        help="Retrieve .mae files intead of pdbs")

    #Others
    parser.add_argument("--rename", action="store_true",
                        help="Avoid core renaming")

    args = parser.parse_args()


    return args.complex_pdb, args.growing_steps, \
           args.criteria, args.plop_path, args.sch_python, args.pele_dir, args.contrl, args.license, \
           args.resfold, args.report, args.traject, args.pdbout, args.cpus, \
           args.distcont, args.threshold, args.epsilon, args.condition, args.metricweights, args.nclusters, \
           args.pele_eq_steps, args.restart, args.min_overlap, args.max_overlap, args.serie_file, \
           args.c_chain, args.f_chain, args.steps, args.temperature, args.seed, args.rotamers, \
           args.banned, args.limit, args.mae, args.rename, args.clash_thr, args.steering, \
           args.translation_high, args.rotation_high, args.translation_low, args.rotation_low, args.explorative, \
           args.radius_box, args.sampling_control, args.data, args.documents, args.only_prepare, args.only_grow, \
           args.no_check, args.debug, args.highthroughput, args.test


def grow_fragment(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
                  pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distance_contact, clusterThreshold,
                  epsilon, condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, ID,
                  h_core=None, h_frag=None, c_chain="L", f_chain="L", steps=6, temperature=1000, seed=1279183, rotamers="30.0",
                  banned=None, limit=None, mae=False, rename=False, threshold_clash=1.7, steering=0,
                  translation_high=0.05, rotation_high=0.10, translation_low=0.02, rotation_low=0.05, explorative=False,
                  radius_box=4, sampling_control=None, data=None, documents=None, only_prepare=False, only_grow=False, 
                  no_check=False, debug=False):
    """
    Description: FrAG is a Fragment-based ligand growing software which performs automatically the addition of several
    fragments to a core structure of the ligand in a protein-ligand complex.
    :param complex_pdb: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
    used as the core structure. Remember to rename the ligand chain with a different character in
    order to detect it.
    :type complex_pdb: str
    :param fragment_pdb: Path to the PDB file containing the fragment.
    :type fragment_pdb: str
    :param core_atom: PDB-atom-name of the atom of the core that will be used as starting point to grow the fragment.
    :type core_atom: str (max length: 4)
    :param fragment_atom: PDB-atom-name of the atom of the fragment that will bond the core.
    :type fragment_atom: str (max length: 4)
    :param iterations: Number of Growing Steps (GS).
    :type iterations: int
    :param criteria: Name of the column of the report file to select the structures that will spawn in the
    next GS. Additionally, this parameter will be the selection criteria to extract the best
    structure after completing the growing.
    :type criteria: str
    :param plop_path: Absolute path to PlopRotTemp.py.
    :type plop_path: str
    :param sch_python: Absolute path to Schrodinger's python.
    :type sch_python: str
    :param pele_dir: Absolute path to PELE executables.
    :type pele_dir: str
    :param contrl: Name of the PELE's control file templatized.
    :type contrl: str
    :param license: Absolute path to PELE's license file.
    :type license: str
    :param resfold: Name of the results folder.
    :type resfold: str
    :param report: Prefix name of PELE's output file.
    :type report: str
    :param traject: Prefix name of the trajectory file generated by PELE.
    :type traject: str
    :param pdbout: Prefix name of the output PDB extracted from FrAG in each iteration.
    :type pdbout: str
    :param cpus: Number of CPUs that will be used to perform the simulation.
    :type cpus: int
    :param distance_contact: This distance will be used to determine which amino acids are in contact with the ligand.
    Then, this information will be useful to generate different clusters of structures to initialize the next iteration.
    :type distance_contact: float
    :param clusterThreshold: Threshold that will be used in the clustering step.
    :type clusterThreshold: float (from 0 to 1)
    :param epsilon: Epsilon parameter used in the clustering.
    :type epsilon: float (from 0 to 1)
    :param condition: Condition to select best values to cluster: minimum or maximum.
    :type condition: str
    :param metricweights: Selects how to distribute the weights of the cluster according to its metric,
    two options: linear (proportional to metric) or Boltzmann weigths (proportional to exp(-metric/T).
    Needs to define the temperature T.
    :type metricweights: str
    :param nclusters: Number of cluster that we would like to generate.
    :type nclusters: int
    :param pele_eq_steps: Number of PELE steps that we want to perform in the equilibration.
    :type pele_eq_steps: int
    :param restart: If set FrAG will find your last iteration completed and will restart from this point.
    :type restart: bool
    :param min_overlap: Minimun value of overlapping factor that will be replaced in the control file.
    :type min_overlap: float (from 0 to 1)
    :param max_overlap: Maximum value of overlapping factor that will be replaced in the control file.
    :type max_overlap: float (from 0 to 1)
    :param ID: Name to identify the new ligand in certain folder.
    :type ID: str
    :param h_core: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the core that will be replaced for
     the fragment, being the initial point of the growing.
    :type h_core: str (max length: 4)
    :param h_frag: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the fragment that will removed to
    bond the fragment with the core.
    :type h_frag: str (max length: 4)
    :param c_chain: Chain name for the ligand in the complex_pdb.
    :type c_chain: str (max length: 1)
    :param f_chain: Chain name for the ligand in the fragment_pdb.
    :type f_chain: str (max length: 1)
    :param steps: PELE steps to do in each GS.
    :type steps: int
    :param temperature: Temperature to add in the control file.
    :type temperature: int
    :param seed: Seed in the PELE's control file.
    :type seed: int
    :param rotamers: Degrees of rotation used to select the rotamer's library. Higher value, faster but less accuracy.
    :type rotamers: str
    :param banned: If set, list of tuples with the atom names of the angles that will not be modified.
    :type banned: list
    :param limit: Limit angle in degrees of the banned angles. Any structure with a higher value will be discarted.
    :type limit: int
    :param mae: If set, the output structures will be saved in "mae." format.
    :type mae: bool
    :param rename: If set, the pdb atom names of the fragment will be renamed with "G+atom_number".
    :type rename: bool
    :param threshold_clash: Distance that will be used to detect contacts between atoms of the fragment and atoms
    of the core.
    :type threshold_clash: float
    :param steering: Steering parameter of the PELE's control file.
    :type steering: int
    :param translation_high: Translation (high value) of the PELE's control file. In anstrongs.
    :type translation_high: float
    :param translation_low: Translation (low value) of the PELE's control file. In anstrongs.
    :type translation_low: float
    :param rotation_high: Rotation (high value) of PELE's control file. In rad.
    :type rotation_high: float
    :param rotation_low: Rotation (low value) of PELE's control file. In rad.
    :type rotation_low: float
    :param explorative: Set all parameters to perform a more explorative simulation after the growing, in the sampling
    simulation.
    :type explorative: bool
    :param radius_box: Radius box size of PELE's control file. In anstrongs.
    :type radius_box: float
    :param sampling_control: templatized control file to be used in the sampling simulation.
    :type sampling_control: str
    :return:
    """
    #Check harcoded path in constants.py
    if not no_check:
        check_constants.check()
    # Time computations
    start_time = time.time()
    # Global variable to keep info
    simulation_info = []
    # Path definition
    plop_relative_path = os.path.join(PackagePath, plop_path)
    current_path = os.path.abspath(".")

    pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
    if "/" in pdb_basename:
        pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name

    working_dir = os.path.join(current_path, "{}_{}".format(pdb_basename, ID))
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)  # Creating a working directory for each PDB-fragment combination
    pdbout_folder = os.path.join(working_dir, pdbout)
    path_to_templates_generated = os.path.join(working_dir,
                                               "DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated")
    path_to_templates = os.path.join(working_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")
    path_to_lib = os.path.join(working_dir, "DataLocal/LigandRotamerLibs")
    # Creation of output folder
    folder_handler.check_and_create_DataLocal(working_dir=working_dir)
    # Creating constraints
    const = "\n".join(constraints.retrieve_constraints(complex_pdb, {}, {}, 5, 5, 10))
    # Creating symbolic links
    helpers.create_symlinks(data, os.path.join(working_dir, 'Data'))
    helpers.create_symlinks(documents, os.path.join(working_dir, 'Documents'))
    #  ---------------------------------------Pre-growing part - PREPARATION -------------------------------------------
    fragment_names_dict, hydrogen_atoms, pdb_to_initial_template, pdb_to_final_template, pdb_initialize, \
    core_original_atom, fragment_original_atom = add_fragment_from_pdbs.main(complex_pdb, fragment_pdb, core_atom,
                                                                             fragment_atom, iterations, h_core=h_core,
                                                                             h_frag=h_frag, core_chain=c_chain,
                                                                             fragment_chain=f_chain, rename=rename,
                                                                             threshold_clash=threshold_clash,
                                                                             output_path=working_dir,
                                                                             only_grow=only_grow)

    # Create the templates for the initial and final structures
    template_resnames = []
    for pdb_to_template in [pdb_to_initial_template, pdb_to_final_template]:
        if not only_grow and not restart:
            cmd = "{} {} {} {} {} {}".format(sch_python, plop_relative_path, os.path.join(working_dir,
                                             add_fragment_from_pdbs.c.PRE_WORKING_DIR, pdb_to_template), rotamers,
                                             path_to_templates_generated, path_to_lib)

            try:
                subprocess.call(cmd.split())
            except OSError:
                raise OSError("Path {} not foud. Change schrodinger path under frag_pele/constants.py".format(sch_python))
        template_resname = add_fragment_from_pdbs.extract_heteroatoms_pdbs(os.path.join(working_dir, add_fragment_from_pdbs.
                                                                           c.PRE_WORKING_DIR, pdb_to_template),
                                                                           False, c_chain, f_chain)
        template_resnames.append(template_resname)

    # Set box center from ligand COM
    resname_core = template_resnames[0]
    center = center_of_mass.center_of_mass(os.path.join(working_dir, c.PRE_WORKING_DIR, "{}.pdb".format(resname_core)))

    # Get template filenames
    template_initial, template_final = ["{}z".format(resname.lower()) for resname in template_resnames]

    if only_prepare:
        print("Files of {} prepared".format(ID))
        return

    # --------------------------------------------GROWING SECTION-------------------------------------------------------
    # Lists definitions

    templates = ["{}_{}".format(os.path.join(path_to_templates_generated, template_final), n) for n in range(0, iterations+1)]

    results = [os.path.join(working_dir, c.OUTPUT_FOLDER, str(n)) for n in range(0, iterations+1)]

    pdbs = [pdb_initialize if n == 0 else "{}_{}".format(n, pdb_initialize) for n in range(0, iterations+1)]

    pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, cpus-1)]

    # Generate starting templates
    template_fragmenter.main(template_initial_path=os.path.join(path_to_templates_generated, template_initial),
                             template_grown_path=os.path.join(path_to_templates_generated, template_final),
                             step=1, total_steps=iterations, hydrogen_to_replace=core_original_atom,
                             core_atom_linker=core_atom,
                             tmpl_out_path=os.path.join(path_to_templates_generated, "{}_0".format(template_final)))

    # Make a copy in the main folder of Templates in order to use it as template for the simulation
    shutil.copy(os.path.join(path_to_templates_generated, "{}_0".format(template_final)),
                os.path.join(path_to_templates, template_final))  # Replace the original template in the folder

    # Clear PDBs folder
    if not restart:
        list_of_subfolders = glob.glob("{}*".format(pdbout_folder))
        for subfolder in list_of_subfolders:
            shutil.rmtree(subfolder)

    # Simulation loop - LOOP CORE
    for i, (template, pdb_file, result) in enumerate(zip(templates, pdbs, results)):

        # Only if reset
        if restart:
            if os.path.exists(os.path.join(pdbout_folder, "{}".format(i))) and os.path.exists(os.path.join(pdbout_folder,
                                                                                                    "{}".format(i),
                                                                                                    "initial_0_0.pdb")):
                print("STEP {} ALREADY DONE, JUMPING TO THE NEXT STEP...".format(i))
                continue
        # Otherwise start from the beggining
        pdb_input_paths = ["{}".format(os.path.join(pdbout_folder, str(i-1), pdb_file)) for pdb_file in pdb_selected_names]
        # Banned dihedrals will be checked here
        if banned:
            pdbs_with_banned_dihedrals = Detector.check_folder(folder=os.path.join(pdbout_folder, str(i-1)),
                                                               threshold=limit,
                                                               dihedrals=banned,
                                                               lig_chain=c_chain,
                                                               processors=cpus)
            pdb_input_paths = [pdb_file for pdb_file, flag in pdbs_with_banned_dihedrals.items() if flag]

        # Control file modification
        overlapping_factor = float(min_overlap) + (((float(max_overlap) - float(min_overlap))*i) / iterations)
        overlapping_factor = "{0:.2f}".format(overlapping_factor)

        if i != 0:
            # Check atom overlapping
            pdbs_with_overlapping = clusterizer.check_atom_overlapping(pdb_input_paths)
            pdb_input_paths_checked = []
            for pdb in pdb_input_paths:
                if pdb not in pdbs_with_overlapping:
                    pdb_input_paths_checked.append(pdb)
            simulation_file = simulations_linker.control_file_modifier(contrl, pdb=pdb_input_paths_checked, step=i,
                                                                       license=license,
                                                                       working_dir=working_dir,
                                                                       overlap=overlapping_factor, results_path=result,
                                                                       steps=steps,
                                                                       chain=c_chain, constraints=const, center=center,
                                                                       temperature=temperature, seed=seed,
                                                                       steering=steering,
                                                                       translation_high=translation_high,
                                                                       translation_low=translation_low,
                                                                       rotation_high=rotation_high,
                                                                       rotation_low=rotation_low,
                                                                       radius=radius_box)
        else:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_initialize, result, i))
            simulation_file = simulations_linker.control_file_modifier(contrl, pdb=[pdb_initialize], step=i,
                                                                       license=license,
                                                                       working_dir=working_dir,
                                                                       overlap=overlapping_factor, results_path=result,
                                                                       steps=steps,
                                                                       chain=c_chain, constraints=const, center=center,
                                                                       temperature=temperature, seed=seed,
                                                                       steering=steering,
                                                                       translation_high=translation_high,
                                                                       translation_low=translation_low,
                                                                       rotation_high=rotation_high,
                                                                       rotation_low=rotation_low,
                                                                       radius=radius_box)

        logger.info(c.LINES_MESSAGE)
        if i != 0:
            template_fragmenter.main(template_initial_path=os.path.join(path_to_templates_generated, template_initial),
                                     template_grown_path=os.path.join(path_to_templates_generated, template_final),
                                     step=i+1, total_steps=iterations, hydrogen_to_replace=core_original_atom,
                                     core_atom_linker=core_atom,
                                     tmpl_out_path=os.path.join(path_to_templates, template_final))

        # Make a copy of the template file in growing_templates folder
        shutil.copy(os.path.join(path_to_templates, template_final), template)

        # Creating results folder
        folder_handler.check_and_create_results_folder(result, working_dir)
        # ------SIMULATION PART------
        # Change directory to the working one
        os.chdir(working_dir)
        if debug:
            return 
        else:
            simulations_linker.simulation_runner(pele_dir, simulation_file, cpus)
        logger.info(c.LINES_MESSAGE)
        logger.info(c.FINISH_SIM_MESSAGE.format(result))
        # Before selecting a step from a trajectory we will save the input PDB file in a folder
        folder_handler.check_and_create_pdb_clusters_folder(pdbout_folder, i)

        # ---------------------------------------------------CLUSTERING-------------------------------------------------
        # Transform column name of the criteria to column number
        result_abs = os.path.abspath(result)
        logger.info("Looking structures to cluster in '{}'".format(result_abs))
        column_number = clusterizer.get_column_num(result_abs, criteria, report)
        # Selection of the trajectory used as new input
        clusterizer.cluster_traject(str(template_resnames[1]), cpus-1, column_number, distance_contact,
                                        clusterThreshold, "{}*".format(os.path.join(result_abs, traject)),
                                        os.path.join(pdbout_folder, str(i)), os.path.join(result_abs),
                                        epsilon, report, condition, metricweights, nclusters)
    # ----------------------------------------------------EQUILIBRATION-------------------------------------------------
    # Set input PDBs
    pdb_inputs = ["{}".format(os.path.join(pdbout_folder, str(iterations), pdb_file)) for pdb_file in pdb_selected_names]
    if banned:
        pdbs_with_banned_dihedrals = Detector.check_folder(folder=os.path.join(pdbout_folder, str(iterations)),
                                                           threshold=limit,
                                                           dihedrals=banned,
                                                           lig_chain=c_chain,
                                                           processors=cpus)
        pdb_inputs = [pdb_file for pdb_file, flag in pdbs_with_banned_dihedrals.items() if flag]
    if not os.path.exists(os.path.join(working_dir, "sampling_result")):  # Create the folder if it does not exist
        os.mkdir(os.path.join(working_dir, "sampling_result"))
    # Modify the control file to increase the steps TO THE SAMPLING SIMULATION
    if sampling_control:
        simulation_file = simulations_linker.control_file_modifier(sampling_control, pdb=pdb_inputs, step=iterations,
                                                                   license=license, working_dir=working_dir,
                                                                   overlap=max_overlap,
                                                                   results_path=os.path.join(working_dir, "sampling_result"),
                                                                   steps=pele_eq_steps, chain=c_chain,
                                                                   constraints=const, center=center,
                                                                   temperature=temperature, seed=seed, steering=steering,
                                                                   translation_high=translation_high,
                                                                   translation_low=translation_low,
                                                                   rotation_high=rotation_high, rotation_low=rotation_low,
                                                                   radius=radius_box)
    elif explorative and not sampling_control:
        simulation_file = simulations_linker.control_file_modifier(contrl, pdb=pdb_inputs, license=license,
                                                                   working_dir=working_dir, step=iterations,
                                                                   overlap=max_overlap,
                                                                   results_path=os.path.join(working_dir, "sampling_result"),
                                                                   steps=pele_eq_steps,
                                                                   chain=c_chain, constraints=const, center=center,
                                                                   temperature=temperature, seed=seed,
                                                                   steering=2,
                                                                   translation_high=0.5,
                                                                   translation_low=0.3,
                                                                   rotation_high=0.4,
                                                                   rotation_low=0.15,
                                                                   radius=25)
    else:
        simulation_file = simulations_linker.control_file_modifier(contrl, pdb=pdb_inputs, step=iterations,
                                                                   license=license, overlap=max_overlap,
                                                                   working_dir=working_dir,
                                                                   results_path=os.path.join(working_dir, "sampling_result"),
                                                                   steps=pele_eq_steps, chain=c_chain,
                                                                   constraints=const, center=center,
                                                                   temperature=temperature, seed=seed, steering=steering,
                                                                   translation_high=translation_high,
                                                                   translation_low=translation_low,
                                                                   rotation_high=rotation_high, rotation_low=rotation_low,
                                                                   radius=radius_box)

    # EQUILIBRATION SIMULATION
    # Change directory to the working one
    os.chdir(working_dir)
    shutil.copy(os.path.join(path_to_templates_generated, template_final), path_to_templates)
    if not (restart and os.path.exists("top_result")):
        logger.info(".....STARTING EQUILIBRATION.....")
        simulations_linker.simulation_runner(pele_dir, simulation_file, cpus)
    os.chdir(curr_dir)
    equilibration_path = os.path.join(working_dir, "sampling_result")
    # SELECTION OF BEST STRUCTURES
    selected_results_path = os.path.join(working_dir, "top_result")
    if not os.path.exists(selected_results_path):  # Create the folder if it does not exist
        os.mkdir(selected_results_path)
    best_structure_file, all_output_files = bestStructs.main(criteria, selected_results_path, path=equilibration_path,
                                                             n_structs=50)

    shutil.copy(os.path.join(selected_results_path, best_structure_file), os.path.join(working_dir, c.PRE_WORKING_DIR,
                                                                                       selected_results_path + ".pdb"))
    # COMPUTE AND SAVE THE SCORE
    analyser.analyse_at_epoch(report_prefix=report, path_to_equilibration=equilibration_path, execution_dir=curr_dir,
                              column=criteria, quantile_value=0.25)

    
    #MOVE FROM PDB TO MAE
    if mae:
        if sch_python.endswith("python"):
            schrodinger_path = os.path.dirname(os.path.dirname(sch_python))
        elif sch_python.endswith("run"):
            schrodinger_path = os.path.dirname(sch_python)
        python_file = os.path.join(os.path.dirname(FilePath), "Analysis/output_files.py")
        for outputfile in all_output_files:
            filename = os.path.join(selected_results_path, outputfile)
            command = "{} {} {} --schr {} {}".format(sch_python, python_file, filename, schrodinger_path, "--remove") 
            subprocess.call(command.split())

    # COMPUTE TIME
    end_time = time.time()
    total_time = (end_time - start_time) / 60
    logging.info("Growing of {} in {} min".format(fragment_pdb, total_time))

    return fragment_names_dict
    

def main(complex_pdb, serie_file, iterations=c.GROWING_STEPS, criteria=c.SELECTION_CRITERIA, plop_path=c.PLOP_PATH, sch_python=c.SCHRODINGER_PY_PATH, pele_dir=c.PATH_TO_PELE, contrl=c.CONTROL_TEMPLATE, license=c.PATH_TO_LICENSE, resfold=c.RESULTS_FOLDER, 
    report=c.REPORT_NAME, traject=c.TRAJECTORY_NAME, pdbout=c.PDBS_OUTPUT_FOLDER, cpus=c.CPUS, distcont=c.DISTANCE_COUNTER, threshold=c.CONTACT_THRESHOLD, epsilon=c.EPSILON, condition=c.CONDITION, metricweights=c.METRICS_WEIGHTS, 
    nclusters=c.NUM_CLUSTERS, pele_eq_steps=c.PELE_EQ_STEPS, restart=False, min_overlap=c.MIN_OVERLAP, max_overlap=c.MAX_OVERLAP,
    c_chain="L", f_chain="L", steps=c.STEPS, temperature=c.TEMPERATURE, seed=c.SEED, rotamers=c.ROTRES, banned=c.BANNED_DIHEDRALS_ATOMS, limit=c.BANNED_ANGLE_THRESHOLD, mae=False,
    rename=None, threshold_clash=1.7, steering=c.STEERING, translation_high=c.TRANSLATION_HIGH, rotation_high=c.ROTATION_HIGH, 
    translation_low=c.TRANSLATION_LOW, rotation_low=c.ROTATION_LOW, explorative=False, radius_box=c.RADIUS_BOX, sampling_control=None, data=c.PATH_TO_PELE_DATA, documents=c.PATH_TO_PELE_DOCUMENTS, 
    only_prepare=False, only_grow=False, no_check=False, debug=False, protocol=False, test=False):

    if protocol == "HT":
        iteration = 1
        steps = 3
        pele_eq_steps = 10

    if test:
        iteration = 1
        steps = 1
        steps = 1
        pele_eq_steps = 1
        temp = 1000000

    #HOT FIX!! Fix it properly
    original_dir = os.path.abspath(os.getcwd())
    list_of_instructions = serie_handler.read_instructions_from_file(serie_file)
    print("READING INSTRUCTIONS... You will perform the growing of {} fragments. GOOD LUCK and ENJOY the "
          "trip :)".format(len(list_of_instructions)))
    dict_traceback = correct_fragment_names.main(complex_pdb)
    for instruction in list_of_instructions:
        # We will iterate trough all individual instructions of file.
        os.chdir(original_dir)
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
                    complex_sequential_pdb = complex_pdb
                    ID = instruction[i][3]
                else:  # If is not the first we will use as input the output of the previous iteration
                    pdb_basename = complex_pdb.split(".pdb")[0]  # Get the name of the pdb without extension
                    if "/" in pdb_basename:
                        pdb_basename = pdb_basename.split("/")[-1]  # And if it is a path, get only the name
                    working_dir = "{}_{}".format(pdb_basename, ID)
                    complex_sequential_pdb = os.path.join(working_dir, "top_result.pdb")
                    dict_traceback = correct_fragment_names.main(complex_sequential_pdb)
                    ID_completed = []
                    for id in instruction[0:i+1]:
                        ID_completed.append(id[3])
                    ID = "".join(ID_completed)
                try:
                    ID = ID.split("/")[-1]
                except Exception:
                    os.chdir(original_dir)
                    traceback.print_exc()
                    if debug: raise Exception()
                try:
                    ID = ID.split("/")[-1]
                except Exception:
                    os.chdir(original_dir)
                    traceback.print_exc()
                    if debug: raise Exception()
                try:
                    serie_handler.check_instructions(instruction[i], complex_sequential_pdb, c_chain, f_chain)
                    print("PERFORMING SUCCESSIVE GROWING...")
                    print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                    atomname_map = grow_fragment(complex_sequential_pdb, fragment_pdb, core_atom, fragment_atom, 
                                   iterations, criteria, plop_path, sch_python,pele_dir, contrl, license, resfold, 
                                   report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, 
                                   metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, 
                                   ID, h_core, h_frag, c_chain, f_chain, steps, temperature, seed, rotamers, banned,
                                   limit, mae, rename, threshold_clash, steering, translation_high, rotation_high,
                                   translation_low, rotation_low, explorative, radius_box, sampling_control, data, documents,
                                   only_prepare, only_grow, no_check, debug)
                    atomname_mappig.append(atomname_map)
 
                except Exception:
                    os.chdir(original_dir)
                    traceback.print_exc()
                    if debug: raise Exception()
        # INDIVIDUAL GROWING
        else:
            # Initialize the growing for each line in the file
            fragment_pdb, core_atom, fragment_atom, ID = instruction[0], instruction[1], instruction[2], instruction[3]
            atoms_if_bond = serie_handler.extract_hydrogens_from_instructions([fragment_pdb, core_atom, fragment_atom])
            try:
                complex_pdb = complex_pdb
                ID = ID.split("/")[-1]
            except Exception:
                os.chdir(original_dir)
                traceback.print_exc()
                if debug: raise Exception()
 
            if atoms_if_bond:
                core_atom = atoms_if_bond[0]
                h_core = atoms_if_bond[1]
                fragment_atom = atoms_if_bond[2]
                h_frag = atoms_if_bond[3]
            else:
                h_core = None
                h_frag = None
            try:
                print("PERFORMING INDIVIDUAL GROWING...")
                print("HYDROGEN ATOMS IN INSTRUCTIONS:  {}    {}".format(h_core, h_frag))
                grow_fragment(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
                     pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon,
                     condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, ID, h_core,
                     h_frag, c_chain, f_chain, steps, temperature, seed, rotamers, banned, limit, mae, rename,
                     threshold_clash, steering, translation_high, rotation_high,
                     translation_low, rotation_low, explorative, radius_box, sampling_control, data, documents,
                     only_prepare, only_grow, no_check, debug)
            except Exception:
                os.chdir(original_dir)
                traceback.print_exc()
                if debug: raise Exception()
        os.chdir(original_dir) 

if __name__ == '__main__':
    complex_pdb, iterations, criteria, plop_path, sch_python, pele_dir, \
    contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, metricweights, \
    nclusters, pele_eq_steps, restart, min_overlap, max_overlap, serie_file, \
    c_chain, f_chain, steps, temperature, seed, rotamers, banned, limit, mae, \
    rename, threshold_clash, steering, translation_high, rotation_high, \
    translation_low, rotation_low, explorative, radius_box, sampling_control, data, documents, \
    only_prepare, only_grow, no_check, debug, protocol, test = parse_arguments()
    
    main(complex_pdb, serie_file, iterations, criteria, plop_path, sch_python, pele_dir, contrl, license, resfold,
             report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, metricweights,
             nclusters, pele_eq_steps, restart, min_overlap, max_overlap,
             c_chain, f_chain, steps, temperature, seed, rotamers, banned, limit, mae,
             rename, threshold_clash, steering, translation_high, rotation_high,
             translation_low, rotation_low, explorative, radius_box, sampling_control, data, documents,
             only_prepare, only_grow, no_check, debug, protocol, test)

