# General imports
import sys
import time
import glob
import argparse
import os
import logging
from logging.config import fileConfig
import shutil
import subprocess
import traceback
# Local imports
import Helpers.clusterizer
import Helpers.checker
import Helpers.modify_rotamers
import Helpers.folder_handler
import Helpers.runner
import Helpers.constraints as cs
import Helpers.helpers as hp
import Helpers.center_of_mass as cm
import Growing.template_fragmenter
import Growing.simulations_linker
import Growing.add_fragment_from_pdbs
import Growing.AddingFragHelpers
import Growing.bestStructs
import constants as c
import serie_handler as sh

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
                        successive growings, using different fragments or different growing positions.\n

                        To do simple growings:\n
                        col1   col2    col3\n
                        To do successive growings:\n
                        (col1   col2    col3) x n_growings\n

                        Where col1 is the path to the PDB file of the fragment that will be added to the core structure
                        (name the chain of the fragment "L" by default).
                        Col2 is a string with the PDB atom name of the heavy atom of the core (the ligand contained in
                        the complex) where you would like to start the growing and create a new bond with the fragment.
                        And col3 is a string with the PDB atom name of the heavy atom of the fragment that will be used
                        to perform the bonding with the core.
                        """)
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
    parser.add_argument("-docon", "--docontrolsim", action="store_true",
                        help="""When it is set, FrAG runs a negative-control simulation (without growing).""")


    # Plop related arguments
    parser.add_argument("-pl", "--plop_path", default=c.PLOP_PATH,
                        help="Absolute path to PlopRotTemp.py. By default = {}".format(c.PLOP_PATH))
    parser.add_argument("-sp", "--sch_python", default=c.SCHRODINGER_PY_PATH,
                        help="""Absolute path to Schrödinger's python. 
                        By default = {}""".format(c.SCHRODINGER_PY_PATH))

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

    args = parser.parse_args()

    return args.complex_pdb, args.growing_steps, \
           args.criteria, args.plop_path, args.sch_python, args.pele_dir, args.contrl, args.license, \
           args.resfold, args.report, args.traject, args.pdbout, args.cpus, \
           args.distcont, args.threshold, args.epsilon, args.condition, args.metricweights, args.nclusters, \
           args.pele_eq_steps, args.restart, args.min_overlap, args.max_overlap, args.serie_file, \
           args.c_chain, args.f_chain, args.docontrolsim, args.steps, args.temperature


def main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
         pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distance_contact, clusterThreshold,
         epsilon, condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, ID,
         h_core=None, h_frag=None, c_chain="L", f_chain="L", steps=6, temperature=1000):
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
    :return:
    """
    # Time computations
    start_time = time.time()
    # Path definition
    plop_relative_path = os.path.join(PackagePath, plop_path)
    pdbout_folder = "{}_{}".format(pdbout, ID)
    path_to_templates_generated = "DataLocal/Templates/OPLS2005/HeteroAtoms/templates_generated"
    path_to_templates = "DataLocal/Templates/OPLS2005/HeteroAtoms"
    # Creation of output folder
    Helpers.folder_handler.check_and_create_DataLocal()
    # Creating constraints
    const = "\n".join(cs.retrieve_constraints(complex_pdb, {}, {}, 5, 5, 10))
    # Creating symbolic links
    hp.create_symlinks(c.PATH_TO_PELE_DATA, 'Data')
    hp.create_symlinks(c.PATH_TO_PELE_DOCUMENTS, 'Documents')

    #  ---------------------------------------Pre-growing part - PREPARATION -------------------------------------------
    fragment_names_dict, hydrogen_atoms, pdb_to_initial_template, pdb_to_final_template, pdb_initialize = Growing.\
        add_fragment_from_pdbs.main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, h_core=h_core,
                                    h_frag=h_frag, core_chain=c_chain, fragment_chain=f_chain)

    # Create the templates for the initial and final structures
    template_resnames = []
    for pdb_to_template in [pdb_to_initial_template, pdb_to_final_template]:
        cmd = "{} {} {}".format(sch_python, plop_relative_path, os.path.join(curr_dir,
                                  Growing.add_fragment_from_pdbs.c.PRE_WORKING_DIR, pdb_to_template))
        new_env = os.environ.copy()
        new_env["PYTHONPATH"] = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/internal/lib/python2.7/site-packages/"
        subprocess.call(cmd.split(), env=new_env)
        template_resname = Growing.add_fragment_from_pdbs.extract_heteroatoms_pdbs(os.path.join(
                                                                                   Growing.add_fragment_from_pdbs.
                                                                                   c.PRE_WORKING_DIR, pdb_to_template),
                                                                                   False)
        template_resnames.append(template_resname)

    # Set box center from ligand COM
    resname_core = template_resnames[0]
    center = cm.center_of_mass(os.path.join("pregrow", "{}.pdb".format(resname_core)))

    # Now, move the templates to their respective folders
    template_initial, template_final = ["{}z".format(resname.lower()) for resname in template_resnames]

    # --------------------------------------------GROWING SECTION-------------------------------------------------------
    # Lists definitions

    templates = ["{}_{}".format(os.path.join(path_to_templates_generated, template_final), n) for n in range(0, iterations+1)]

    results = ["{}{}_{}{}".format(c.OUTPUT_FOLDER, ID, resfold, n) for n in range(0, iterations+1)]

    pdbs = [pdb_initialize if n == 0 else "{}_{}".format(n, pdb_initialize) for n in range(0, iterations+1)]

    pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, cpus-1)]

    original_atom = hydrogen_atoms[0].get_name()  # Hydrogen of the core that we will use as growing point
    # Generate starting templates
    replacement_atom = fragment_names_dict[fragment_atom]
    Growing.template_fragmenter.create_initial_template(initial_template=os.path.join(path_to_templates_generated,
                                                                                      template_initial),
                                                        final_template=os.path.join(path_to_templates_generated,
                                                                                    template_final),
                                                        original_atom_to_mod=[original_atom],
                                                        heavy_atom=core_atom,
                                                        atom_to_transform=replacement_atom,
                                                        output_template_filename="{}_0".format(os.path.join(
                                                                                        path_to_templates_generated,
                                                                                        template_final)),
                                                        path="", steps=iterations)
    Growing.template_fragmenter.generate_starting_template(initial_template_file=os.path.join(path_to_templates_generated,
                                                                                          template_initial),
                                                           final_template_file=os.path.join(path_to_templates_generated,
                                                                                          template_final),
                                                           original_atom_to_mod=[original_atom],
                                                           heavy_atom=core_atom,
                                                           atom_to_transform=replacement_atom,
                                                           output_template_filename="{}_ref".format(os.path.join(
                                                                                        path_to_templates_generated,
                                                                                        template_final)),
                                                           path="", steps=iterations)

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

        # Control file modification
        overlapping_factor = float(min_overlap) + (((float(max_overlap) - float(min_overlap))*i) / iterations)
        overlapping_factor = "{0:.2f}".format(overlapping_factor)

        if i != 0:
            simulation_file = Growing.simulations_linker.control_file_modifier(contrl, pdb_input_paths, i, license, overlapping_factor,
                                                             result, steps=steps, chain=c_chain, constraints=const, center=center,
                                                                               temperature=temperature)
        else:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_initialize, result, i))
            simulation_file = Growing.simulations_linker.control_file_modifier(contrl, [pdb_initialize], i, license, overlapping_factor,
                                                             result, steps=steps, chain=c_chain, constraints=const, center=center,
                                                                               temperature=temperature)
                                                                     #  We have put [] in pdb_initialize
                                                                     #  because by default we have to use
                                                                     #  a list as input
        logger.info(c.LINES_MESSAGE)
        if i != 0 and i != iterations:
            Growing.template_fragmenter.grow_parameters_in_template(starting_template_file="{}_ref".format(
                                                                    os.path.join(path_to_templates_generated,
                                                                                 template_final)),
                                                                    initial_template_file=os.path.join(
                                                                                           path_to_templates_generated,
                                                                                           template_initial),
                                                                    final_template_file=os.path.join(
                                                                                           path_to_templates_generated,
                                                                                           template_final),
                                                                    original_atom_to_mod=[original_atom],
                                                                    heavy_atom=core_atom,
                                                                    atom_to_transform=replacement_atom,
                                                                    output_template_filename=os.path.join(
                                                                        path_to_templates, template_final),
                                                                    path="", step=i, total_steps=iterations)
        elif i == iterations:
            shutil.copy(os.path.join(path_to_templates_generated, template_final), os.path.join(path_to_templates,
            template_final))

        # Make a copy of the template file in growing_templates folder
        shutil.copy(os.path.join(path_to_templates, template_final), template)

        # Creating results folder
        Helpers.folder_handler.check_and_create_results_folder(result)
        Growing.simulations_linker.simulation_runner(pele_dir, simulation_file, cpus)
        logger.info(c.LINES_MESSAGE)
        logger.info(c.FINISH_SIM_MESSAGE.format(result))
        # Before selecting a step from a trajectory we will save the input PDB file in a folder
        Helpers.folder_handler.check_and_create_pdb_clusters_folder(pdbout_folder, i)

        # ---------------------------------------------------CLUSTERING-------------------------------------------------
        # Transform column name of the criteria to column number
        result_abs = os.path.abspath(result)
        logger.info("Looking structures to cluster in '{}'".format(result_abs))
        column_number = Helpers.clusterizer.get_column_num(result_abs, criteria, report)
        # Selection of the trajectory used as new input
        Helpers.clusterizer.cluster_traject(str(template_resnames[1]), cpus-1, column_number, distance_contact,
                                            clusterThreshold, "{}*".format(os.path.join(result_abs, traject)),
                                            os.path.join(pdbout_folder, str(i)), epsilon, report, condition, metricweights,
                                            nclusters)
    # ----------------------------------------------------EQUILIBRATION-------------------------------------------------
    # Set input PDBs
    pdb_inputs = ["{}".format(os.path.join(pdbout_folder, str(iterations), pdb_file)) for pdb_file in pdb_selected_names]
    if not os.path.exists("equilibration_result_{}".format(ID)):  # Create the folder if it does not exist
        os.mkdir("equilibration_result_{}".format(ID))
    # Modify the control file to increase the steps to 20 and change the output path
    simulation_file = Growing.simulations_linker.control_file_modifier(contrl, pdb_inputs, iterations, license, max_overlap,
                                                     "equilibration_result_{}".format(ID), steps=pele_eq_steps, chain=c_chain,
                                                     constraints=const, center=center, temperature=temperature)
    # Call PELE to run the simulation
    if not (restart and os.path.exists("selected_result_{}".format(ID))):
        logger.info(".....STARTING EQUILIBRATION.....")
        Growing.simulations_linker.simulation_runner(pele_dir, simulation_file, cpus)
    equilibration_path = os.path.join(os.path.abspath(os.path.curdir), "equilibration_result_{}".format(ID))
    selected_results_path = "selected_result_{}".format(ID)
    if not os.path.exists(selected_results_path):  # Create the folder if it does not exist
        os.mkdir(selected_results_path)
    best_structure_file = Growing.bestStructs.main(criteria, selected_results_path, path=equilibration_path,
                             n_structs=10)
    shutil.copy(os.path.join(selected_results_path, best_structure_file), os.path.join(c.PRE_WORKING_DIR, selected_results_path + ".pdb"))
    end_time = time.time()
    total_time = (end_time - start_time) / 60
    logging.info("Growing of {} in {} min".format(fragment_pdb, total_time))


if __name__ == '__main__':
    complex_pdb, iterations, criteria, plop_path, sch_python, pele_dir, \
    contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, metricweights, \
    nclusters, pele_eq_steps, restart, min_overlap, max_overlap, serie_file, \
    c_chain, f_chain, docontrolsim, steps, temperature = parse_arguments()

    list_of_instructions = sh.read_instructions_from_file(serie_file)
    print("READING INSTRUCTIONS... You will perform the growing of {} fragments. GOOD LUCK and ENJOY the trip :)".format(len(list_of_instructions)))
    sh.check_instructions(list_of_instructions, complex_pdb, c_chain, f_chain)
    for instruction in list_of_instructions:
        # We will iterate trough all individual instructions of file.
        if type(instruction) == list:  #  If in the individual instruction we have more than one command means successive growing.
            growing_counter = len(instruction)  #  Doing so we will determinate how many successive growings the user wants to do.
            for i in range(int(growing_counter)):
                fragment_pdb, core_atom, fragment_atom = instruction[i][0], instruction[i][1], instruction[i][2]
                atoms_if_bond = sh.extract_hydrogens_from_instructions([fragment_pdb, core_atom, fragment_atom])
                if atoms_if_bond:
                    core_atom = atoms_if_bond[0]
                    h_core = atoms_if_bond[1]
                    fragment_atom = atoms_if_bond[2]
                    h_frag = atoms_if_bond[3]
                else:
                    pass
                if i == 0:  # In the first iteration we will use the complex_pdb as input.
                    ID = instruction[i][3]
                else:  # If is not the first we will use as input the output of the previous iteration
                    complex_pdb = "pregrow/selected_result_{}.pdb".format(ID)
                    ID_completed = []
                    for id in instruction:
                        ID_completed.append(id[3])
                    ID = "".join(ID_completed)

                try:
                    if "/" in ID:
                        ID = ID.split("/")[-1]
                    main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path,
                     sch_python,pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distcont,
                     threshold, epsilon, condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap,
                     max_overlap, ID, h_core, h_frag, c_chain, f_chain, steps, temperature)

                except Exception:
                    traceback.print_exc()

        else:
            # Initialize the growing for each line in the file
            fragment_pdb, core_atom, fragment_atom, ID = instruction[0], instruction[1], instruction[2], instruction[3]
            if "/" in ID:
                ID = ID.split("/")[-1]
            atoms_if_bond = sh.extract_hydrogens_from_instructions([fragment_pdb, core_atom, fragment_atom])
            if atoms_if_bond:
                core_atom = atoms_if_bond[0]
                h_core = atoms_if_bond[1]
                fragment_atom = atoms_if_bond[2]
                h_frag = atoms_if_bond[3]

            try:
                main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
                 pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon,
                 condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, ID, h_core,
                 h_frag, c_chain, f_chain, steps, temperature)
            except Exception:
                traceback.print_exc()
    if docontrolsim:
        logging.info("INITIALIZING NON-GROWING SIMULATION :)")
        ID = "{}_cntrl".format(complex_pdb)
        Helpers.runner.run_no_grow(complex_pdb, sch_python, resfold, iterations, cpus, ID, pele_dir,
                                   pdbout, restart, min_overlap, max_overlap, contrl, criteria, report,
                                   epsilon, threshold, distcont, condition, traject,
                                   metricweights, nclusters, pele_eq_steps, license)
