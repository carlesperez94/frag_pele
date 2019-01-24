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
    # Growing related arguments
    required_named.add_argument("-cp", "--complex_pdb", required=True,
                        help="""PDB file which contains a protein-ligand complex that we will use as 
                        core for growing (name the chain of the ligand "L")""")
    required_named.add_argument("-sef", "--serie_file", required=True,
                        help="""
                        Name of the tabular file which contains the information required to perform several 
                        successive growings using different fragments or different growing positions.
                        
                        To do simple growings:
                        col1   col2    col3
                        To do successive growings:
                        (col1   col2    col3) x n_growings 
                        
                        Where col1 is PDB file which contains the fragment that will be added to the core structure 
                        (name the chain of the fragment "L" by default).
                        Col2 is a string with the PDB atom name of the heavy atom of the core (the ligand contained in 
                        the complex) where you would like to start the growing and create a new bond with the fragment.
                        And col3 is a string with the PDB atom name of the heavy atom of the fragment that will be used 
                        to perform the bonding with the core.
                        
                        """)
    parser.add_argument("-x", "--iterations", type=int, default=c.ITERATIONS,
                        help="""Number of intermediate templates that you want to generate""")
    parser.add_argument("-cr", "--criteria", default=c.SELECTION_CRITERIA,
                        help="""Name of the column used as criteria in order to select the template used as input for 
                                 successive simulations.""")
    parser.add_argument("-rst", "--restart", action="store_true",
                        help="If restart is true FrAG will continue from the last epoch detected.")

    parser.add_argument("-hc", "--h_core", default=None,
                        help="""String with the PDB atom name of the hydrogen atom of the core (the ligand contained in 
                        the complex) where you would like to start the growing and create a new bond with the fragment.""")
    parser.add_argument("-hf", "--h_frag", default=None,
                        help="""String with the PDB atom name of the hydrogen atom of the fragment
                         where you would like create a new bond with the core.""")

    parser.add_argument("-cc", "--c_chain", default="L", help="Chain name of the core")

    parser.add_argument("-fc", "--f_chain", default="L", help="Chain name of the fragment")

    parser.add_argument("-docon", "--docontrolsim", default=False,
                        help="""When it is true FrAG runs a control simulation (without growing).""")

    # Plop related arguments
    parser.add_argument("-pl", "--plop_path", default=c.PLOP_PATH,
                        help="Absolute path to PlopRotTemp.py.")
    parser.add_argument("-sp", "--sch_python", default=c.SCHRODINGER_PY_PATH,
                        help="Absolute path to Schrödinger's python.")
    # PELE configuration arguments
    parser.add_argument("-d", "--pele_dir", default=c.PATH_TO_PELE,
                        help="Complete path to Pele_serial")
    parser.add_argument("-c", "--contrl", default=c.CONTROL_TEMPLATE,
                        help="Control file templatized.")
    parser.add_argument("-l", "--license", default=c.PATH_TO_LICENSE,
                        help="Absolute path to PELE's licenses folder.")
    parser.add_argument("-r", "--resfold", default=c.RESULTS_FOLDER,
                        help="Name for results folder")
    parser.add_argument("-rp", "--report", default=c.REPORT_NAME,
                        help="Suffix name of the report file from PELE.")
    parser.add_argument("-tj", "--traject", default=c.TRAJECTORY_NAME,
                        help="Suffix name of the trajectory file from PELE.")
    parser.add_argument("-pdbf", "--pdbout", default=c.PDBS_OUTPUT_FOLDER,
                        help="PDBs output folder")
    parser.add_argument("-cs", "--cpus", default=c.CPUS,
                        help="Amount of CPU's that you want to use in mpirun of PELE")
    parser.add_argument("-es", "--pele_eq_steps", default=c.PELE_EQ_STEPS,
                        help="Number of PELE steps in equilibration")
    parser.add_argument("-miov", "--min_overlap", default=c.MIN_OVERLAP,
                        help="Minimum value of overlapping factor used in the control_file of PELE.")
    parser.add_argument("-maov", "--max_overlap", default=c.MAX_OVERLAP,
                        help="Maximum value of overlapping factor used in the control_file of PELE.")

    # Clustering related arguments
    parser.add_argument("-dis", "--distcont", default=c.DISTANCE_COUNTER,
                        help="""This distance will be used to determine which amino acids are in contact with the ligand.
                        Then, this information will be useful to generate different clusters of structures to initialize
                         the next iteration.""")
    parser.add_argument("-ct", "--threshold", default=c.CONTACT_THRESHOLD,
                        help="Threshold that will be used in the clustering step.")
    parser.add_argument("-e", "--epsilon", default=c.EPSILON,
                        help="Epsilon parameter used in the clustering.")
    parser.add_argument("-cn", "--condition", default=c.CONDITION,
                        help="Condition to select best values to cluster: minimum or maximum.")
    parser.add_argument("-mw", "--metricweights", default=c.METRICS_WEIGHTS,
                        help="")
    parser.add_argument("-ncl", "--nclusters", default=c.NUM_CLUSTERS,
                        help="Number of initial structures that we want to use in each simulation.")

    args = parser.parse_args()

    return args.complex_pdb, args.iterations, \
           args.criteria, args.plop_path, args.sch_python, args.pele_dir, args.contrl, args.license, \
           args.resfold, args.report, args.traject, args.pdbout, args.cpus, \
           args.distcont, args.threshold, args.epsilon, args.condition, args.metricweights, args.nclusters, \
           args.pele_eq_steps, args.restart, args.min_overlap, args.max_overlap, args.serie_file, \
           args.h_core, args.h_frag, args.c_chain, args.f_chain, args.docontrolsim


def main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
         pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distance_contact, clusterThreshold,
         epsilon, condition, metricweights, nclusters, pele_eq_steps, restart, min_overlap, max_overlap, ID, h_core,
         h_frag, c_chain, f_chain):
    """
    Description: FrAG is a Fragment-based ligand growing software which performs automatically the addition of several
    fragments to a core structure of the ligand in a protein-ligand complex.
    :param complex_pdb: name of the PDB file containing the protein and the core structure. Str
    :param fragment_pdb: name of the PDB file containing the fragment. Str
    :param core_atom: PDB-atom-name of the atom of the core that will be used as starting point to grow the fragment.
    Str (len 4)
    :param fragment_atom: PDB-atom-name of the atom of the fragment that will bond the core. Str (len 4)
    :param iterations: number of growing steps. Int
    :param criteria: Name of the report's column used as selection criteria to spawn the structures that will be used in
    the next growing step. Str
    :param plop_path: Absolute path to PlopRotTemp.py. Str
    :param sch_python: Absolute path to Schrodinger's python. Str
    :param pele_dir: Absolute path to PELE executables. Str
    :param contrl: Name of the PELE's control file templatized. Str
    :param license: Absolute path to PELE's license file. Str
    :param resfold: Name of the results folder. Str
    :param report: Prefix name of PELE's output file. Str
    :param traject: Prefix name of the trajectory file generated by PELE. Str
    :param pdbout: Prefix name of the output PDB extracted from FrAG in each iteration. Str
    :param cpus: Number of CPUs that will be used to perform the simulation. Int
    :param distance_contact: This distance will be used to determine which amino acids are in contact with the ligand.
    Then, this information will be useful to generate different clusters of structures to initialize the next iteration.
    :param clusterThreshold: Threshold that will be used in the clustering step. Float (0.0-1.0)
    :param epsilon: Epsilon parameter used in the clustering. Float (0.0-1.0)
    :param condition: Condition to select best values to cluster: minimum or maximum. Str "min" or "max".
    :param metricweights: Str "linear"
    :param nclusters: Number of cluster that we would like to generate. Int
    :param pele_eq_steps: Number of PELE steps that we want to perform in the equilibration. Int
    :param restart: If set FrAG will find your last iteration completed and will restart from this point. Bolean
    :param min_overlap: Minimun value of overlapping factor that will be replaced in the control file. Float (0.0-1.0)
    :param max_overlap: Maximum value of overlapping factor that will be replaced in the control file. Float (0.0-1.0)
    :param ID: Name that will be generated to identify the new ligand in certain folder.
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
    template_names = []
    for resname in template_resnames:
        template_name = "{}z".format(resname.lower())
        template_names.append(template_name)
        Helpers.modify_rotamers.change_angle(os.path.join(curr_dir, c.ROTAMERS_PATH, "{}.rot.assign".format(resname)), c.ROTRES)
    template_initial, template_final = template_names

    # --------------------------------------------GROWING SECTION-------------------------------------------------------
    # Lists definitions

    templates = ["{}_{}".format(os.path.join(path_to_templates_generated, template_final), n) for n in range(0, iterations+1)]

    results = ["{}{}_{}{}".format(c.OUTPUT_FOLDER, ID, resfold, n) for n in range(0, iterations+1)]

    pdbs = [pdb_initialize if n == 0 else "{}_{}".format(n, pdb_initialize) for n in range(0, iterations+1)]

    pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, int(cpus)-1)]

    original_atom = hydrogen_atoms[0].get_name()  # Hydrogen of the core that we will use as growing point
    # Generate starting templates
    replacement_atom = fragment_names_dict[fragment_atom]
    Growing.template_fragmenter.create_initial_template(template_initial, template_final, [original_atom], core_atom,
                                                        replacement_atom, "{}_0".format(template_final),
                                                        path_to_templates_generated,
                                                        iterations)
    Growing.template_fragmenter.generate_starting_template(template_initial, template_final, [original_atom],core_atom,
                                                           replacement_atom, "{}_ref".format(template_final),
                                                           path_to_templates_generated,
                                                           iterations)

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
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_selected_names, result, i))
            Growing.simulations_linker.control_file_modifier(contrl, pdb_input_paths, i, license, overlapping_factor,
                                                             result, chain=c_chain, constraints=const, center=center)
        else:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_initialize, result, i))
            Growing.simulations_linker.control_file_modifier(contrl, [pdb_initialize], i, license, overlapping_factor,
                                                             result, chain=c_chain, constraints=const, center=center) 
                                                                     #  We have put [] in pdb_initialize
                                                                     #  because by default we have to use
                                                                     #  a list as input
        logger.info(c.LINES_MESSAGE)
        if i != 0 and i != iterations:
            Growing.template_fragmenter.grow_parameters_in_template("{}_ref".format(template_final),
                                                                    os.path.join(path_to_templates_generated, template_initial),
                                                                    os.path.join(path_to_templates_generated, template_final),
                                                                    [original_atom], core_atom, replacement_atom, template_final,
                                                                    os.path.join(path_to_templates, i, iterations))
        elif i == iterations:
            shutil.copy(os.path.join(path_to_templates, template_final), path_to_templates_generated, template_final)

        # Make a copy of the template file in growing_templates folder
        shutil.copy(os.path.join(path_to_templates, template_final), template)

        # Creating results folder
        Helpers.folder_handler.check_and_create_results_folder(result)
        # Create symbolic links
        if not os.path.exists("Data"):
            os.symlink(c.PATH_TO_DATA, "Data")
        if not os.path.exists("Documents"):
            os.symlink(c.PATH_TO_DATA, "Documents")
        # Start simulation
        Growing.simulations_linker.simulation_runner(pele_dir, contrl, int(cpus))
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

        Helpers.clusterizer.cluster_traject(str(template_resnames[1]), int(cpus-1), column_number, distance_contact,
                                            clusterThreshold, "{}*".format(os.path.join(result_abs, traject)),
                                            os.path.join(pdbout_folder, str(i)), epsilon, report, condition, metricweights,
                                            nclusters)
    # ----------------------------------------------------EQUILIBRATION-------------------------------------------------
    # Set input PDBs
    pdb_inputs = ["{}".format(os.path.join(pdbout_folder, str(iterations), pdb_file)) for pdb_file in pdb_selected_names]
    logger.info(".....STARTING EQUILIBRATION.....")
    if not os.path.exists("equilibration_result_{}".format(ID)):  # Create the folder if it does not exist
        os.mkdir("equilibration_result_{}".format(ID))
    # Modify the control file to increase the steps to 20 and change the output path
    Growing.simulations_linker.control_file_modifier(contrl, pdb_inputs, iterations, license, max_overlap,
                                                     "equilibration_result_{}".format(ID), pele_eq_steps)
    # Call PELE to run the simulation
    Growing.simulations_linker.simulation_runner(pele_dir, contrl, int(cpus))
    equilibration_path = os.path.join(os.path.abspath(os.path.curdir), "equilibration_result_{}".format(ID))
    if not os.path.exists("selected_result_{}".format(ID)):  # Create the folder if it does not exist
        os.mkdir("selected_result_{}".format(ID))
    os.chdir("selected_result_{}".format(ID))
    Growing.bestStructs.main(criteria, "best_structure.pdb", path=equilibration_path,
                             n_structs=10)
    shutil.copy("sel_0_best_structure.pdb", "../pregrow/selection_{}.pdb".format(ID))
    os.chdir("../")
    end_time = time.time()
    total_time = (end_time - start_time) / 60
    logging.info("Growing of {} in {} min".format(fragment_pdb, total_time))


if __name__ == '__main__':
    complex_pdb, iterations, criteria, plop_path, sch_python, pele_dir, \
    contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, metricweights, \
    nclusters, pele_eq_steps, restart, min_overlap, max_overlap, serie_file, h_core, h_frag, \
    c_chain, f_chain, docontrolsim = parse_arguments()

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
                    complex_pdb = "selection_{}.pdb".format(ID)
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
                     max_overlap, ID, h_core, h_frag, c_chain, f_chain)

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
                 h_frag, c_chain, f_chain)
            except Exception:
                traceback.print_exc()
    if docontrolsim:
        logging.info("INITIALIZING NON-GROWING SIMULATION :)")
        ID = "{}_cntrl".format(complex_pdb)
        Helpers.runner.run_no_grow(complex_pdb, sch_python, resfold, iterations, cpus, ID, pele_dir,
                                   pdbout, restart, min_overlap, max_overlap, contrl, criteria, report,
                                   epsilon, threshold, distcont, condition, traject,
                                   metricweights, nclusters, pele_eq_steps, license)
