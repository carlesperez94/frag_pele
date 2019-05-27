# General imports
import sys
import os
import time
import subprocess
import glob
import shutil
import logging
# Local imports
import frag_pele.constants as c
import frag_pele.Helpers.folder_handler
import frag_pele.Helpers.clusterizer
import frag_pele.Growing.add_fragment_from_pdbs
import frag_pele.Growing.simulations_linker
import frag_pele.Growing.bestStructs


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)
curr_dir = os.path.abspath(os.path.curdir)
FilePath = os.path.relpath(__file__)


def run_no_grow(complex_pdb, sch_python, resfold, iterations, cpus, ID, pele_dir,
                pdbout, restart, min_overlap, max_overlap, contrl, criteria, report,
                epsilon, clusterThreshold, distance_contact, condition, traject,
                metricweights, nclusters, pele_eq_steps, license):
    # Time computations
    start_time = time.time()
    pdbout_folder = "{}_{}".format(pdbout, ID)

    # Path definition
    plop_relative_path = os.path.join("/".join(FilePath.split("/")[:-2]), c.PLOP_PATH)
    Helpers.folder_handler.check_and_create_DataLocal()

    # Get ligand and extract templates
    pdb_input = os.path.join(c.PRE_WORKING_DIR, complex_pdb)
    ligand_name = Growing.add_fragment_from_pdbs.extract_heteroatoms_pdbs(pdb_input)
    pdb_to_template = "{}.pdb".format(ligand_name)
    shutil.copy(pdb_to_template, c.TEMPLATES_PATH)
    cmd = "{} {} {}".format(sch_python, plop_relative_path, os.path.join(curr_dir, c.TEMPLATES_PATH, pdb_to_template))
    subprocess.call(cmd.split())
    template_name = "{}z".format(ligand_name.lower())
    shutil.copy(template_name, os.path.join(curr_dir, c.TEMPLATES_PATH))
    shutil.copy("{}.rot.assign".format(ligand_name), os.path.join(curr_dir, c.ROTAMERS_PATH))

    # Preparation of results and output folders
    results = ["{}{}_{}{}".format(c.OUTPUT_FOLDER, ID, resfold, n) for n in range(0, iterations + 1)]
    pdbs = [complex_pdb if n == 0 else "{}_{}".format(n, complex_pdb) for n in range(0, iterations + 1)]
    pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, int(cpus) - 1)]

    if not restart:
        list_of_subfolders = glob.glob("{}*".format(pdbout_folder))
        for subfolder in list_of_subfolders:
            shutil.rmtree(subfolder)
    # Simulation loop - LOOP CORE
    for i, (pdb_file, result) in enumerate(zip(pdbs, results)):
        # Only if reset
        if restart:
            if os.path.exists(os.path.join(pdbout_folder, "{}".format(i))) and os.path.exists(os.path.join(pdbout_folder,
                                                                                                    "{}".format(i),
                                                                                                    "initial_0_0.pdb")):
                print("STEP {} ALREADY DONE, JUMPING TO THE NEXT STEP...".format(i))
                continue
        # Otherwise start from the beggining
        pdb_input_paths = ["{}".format(os.path.join(pdbout_folder, str(i - 1), pdb_file)) for pdb_file in
                           pdb_selected_names]

        # Control file modification
        overlapping_factor = float(min_overlap) + (((float(max_overlap) - float(min_overlap))*i) / iterations)
        overlapping_factor = "{0:.2f}".format(overlapping_factor)

        if i != 0:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_selected_names, result, i))
            Growing.simulations_linker.control_file_modifier(contrl, pdb_input_paths, i, license, overlapping_factor,
                                                             result)
        else:
            logger.info(c.SELECTED_MESSAGE.format(contrl, os.path.join(c.PRE_WORKING_DIR, complex_pdb), result, i))
            Growing.simulations_linker.control_file_modifier(contrl, [os.path.join(c.PRE_WORKING_DIR, complex_pdb)],
                                                             i, license, overlapping_factor,
                                                             result) #  We have put [] in pdb_initialize
                                                                     #  because by default we have to use
                                                                     #  a list as input
        logger.info(c.LINES_MESSAGE)

        # Running PELE simulation
        Helpers.folder_handler.check_and_create_results_folder(result)
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

        Helpers.clusterizer.cluster_traject(str(ligand_name), int(cpus - 1), column_number, distance_contact,
                                            clusterThreshold, "{}*".format(os.path.join(result_abs, traject)),
                                            os.path.join(pdbout_folder, str(i)), epsilon, report, condition,
                                            metricweights, nclusters)

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
    Growing.bestStructs.main(criteria, "best_structure.pdb", path=equilibration_path, n_structs=10)
    shutil.copy("sel_0_best_structure.pdb", "../pregrow/selection_{}.pdb".format(ID))
    os.chdir("../")
    end_time = time.time()
    total_time = (end_time - start_time) / 60
    logging.info("Simulation of {} performed in {} min".format(complex_pdb, total_time))
