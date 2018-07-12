# General imports
import sys
import argparse
import os
import logging
from logging.config import fileConfig
import shutil
import subprocess
# Local imports
import Helpers.clusterizer
import Growing.template_fragmenter
import Growing.simulations_linker
import Growing.add_fragment_from_pdbs
import Growing.AddingFragHelpers
import Growing.bestStructs
import constants as c

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
    parser = argparse.ArgumentParser(description="""From an input file, correspondent to the template of the initial 
    structure of the ligand,this function will generate "x_fragments" intermediate templates until reach the final 
    template,modifying Vann Der Waals, bond lengths and deleting charges.""")
    required_named = parser.add_argument_group('required named arguments')
    # Growing related arguments
    required_named.add_argument("-cp", "--complex_pdb", required=True,
                        help="""PDB file which contains a protein-ligand complex that we will use as 
                        core for growing (name the chain of the ligand "L")""")
    required_named.add_argument("-fp", "--fragment_pdb", required=True,
                        help="""PDB file which contains the fragment that will be added to the core structure (name the 
                        chain of the fragment "L")""")
    required_named.add_argument("-ca", "--core_atom", required=True,
                        help="""String with the PDB atom name of the heavy atom of the core (the ligand contained in 
                        the complex) where you would like to start the growing and create a new bond with the fragment.
                        """)
    required_named.add_argument("-fa", "--fragment_atom", required=True,
                        help="""String with the PDB atom name of the heavy atom of the fragment that will be used to 
                        perform the bonding with the core.""")

    parser.add_argument("-x", "--iterations", type=int, default=c.ITERATIONS,
                        help="""Number of intermediate templates that you want to generate""")
    parser.add_argument("-cr", "--criteria", default=c.SELECTION_CRITERIA,
                        help="""Name of the column used as criteria in order to select the template used as input for 
                                 successive simulations.""")
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
    # Clustering related arguments
    parser.add_argument("-dis", "--distcont", default=c.DISTANCE_COUNTER,
                        help="Name for results folder")
    parser.add_argument("-ct", "--threshold", default=c.CONTACT_THRESHOLD,
                        help="Suffix name of the report file from PELE.")
    parser.add_argument("-e", "--epsilon", default=c.EPSILON,
                        help="Suffix name of the trajectory file from PELE.")
    parser.add_argument("-cn", "--condition", default=c.CONDITION,
                        help="PDBs output folder")
    parser.add_argument("-mw", "--metricweights", default=c.METRICS_WEIGHTS,
                        help="Amount of CPU's that you want to use in mpirun of PELE")
    parser.add_argument("-ncl", "--nclusters", default=c.NUM_CLUSTERS,
                        help="Number of initial structures that we want to use in each simulation.")
    parser.add_argument("-rr", "--rotres", default=c.ROTRES,
                        help="")
    args = parser.parse_args()

    return args.complex_pdb, args.fragment_pdb, args.core_atom, args.fragment_atom, args.iterations, \
           args.criteria, args.plop_path, args.sch_python, args.pele_dir, args.contrl, args.license, \
           args.resfold, args.report, args.traject, args.pdbout, args.cpus, \
           args.distcont, args.threshold, args.epsilon, args.condition, args.metricweights, args.nclusters, args.rotres


def main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python,
         pele_dir, contrl, license, resfold, report, traject, pdbout, cpus, distance_contact, clusterThreshold,
         epsilon, condition, metricweights, nclusters, rotres):
    """
        Description: This function is the main core of the program. It creates N intermediate templates
        and control files for PELE. Then, it perform N successive simulations automatically.

        Input:

        :param: template_initial: Name of the input file correspondent to the initial template for the ligand that you
        want to grow.
        :param template_final: Name of the input file correspondent to the final template for the ligand that you want
        to get.
        :param n_files: Number of intermediate templates that you want to generate
        :param original_atom: When an atom is transformed into another one we want to conserve the properties.
        The PDB atom name of the initial template that we want to transform into another of the final template has to be
        specified here.
       :param final_atom: When an atom is transformed into another one we want to conserve the properties.
        The PDB atom name of the final template that will be used as target to be transformed into the original atom.
        :param control_template: Template control file used to generate intermediates control files.
        :param pdb: Initial pdb file which already contain the ligand with the fragment that we want to grow
        but with bond lengths correspondent to the initial ligand (dummy-like).
        :param results_f_name: Name for results folder.
        :param criteria: Name of the column of the report file used as criteria in order to select the template
        used as input for successive simulations.
        :param path_pele: Complete path to Pele_serial.
        :param report: Name of the report file of PELE simulation's results.
        :param traject: Name of the trajectory file of PELE simulation's results.
        :return The algorithm is formed by two main parts:
        First, it prepare the files needed in a PELE simulation: control file and ligand templates.
        Second, after analyzing the report file it selects the best structure in the trajectory file in order to be used
        as input for the next growing step.
        This process will be repeated until get the final grown ligand (protein + ligand).
        """
    # Path definition
    plop_relative_path = os.path.join(PackagePath, plop_path)

    #  ---------------------------------------Pre-growing part - PREPARATION -------------------------------------------

    fragment_names_dict, hydrogen_atoms, pdb_to_initial_template, pdb_to_final_template, pdb_initialize = Growing.\
        add_fragment_from_pdbs.main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations)
    print(hydrogen_atoms)
    # Create the templates for the initial and final structures
    template_resnames = []
    for pdb_to_template in [pdb_to_initial_template, pdb_to_final_template]:
        cmd = "{} {} {}".format(sch_python, plop_relative_path, os.path.join(curr_dir,
                                  Growing.add_fragment_from_pdbs.PRE_WORKING_DIR, pdb_to_template))
		                  #,"--gridres {}".format(rotres))
        subprocess.call(cmd.split())
        template_resname = Growing.add_fragment_from_pdbs.extract_heteroatoms_pdbs(os.path.join(
                                                                                   Growing.add_fragment_from_pdbs.
                                                                                   PRE_WORKING_DIR, pdb_to_template),
                                                                                   False)
        template_resnames.append(template_resname)
    # Now, move the templates to their respective folders
    template_names = []
    for resname in template_resnames:
        template_name = "{}z".format(resname.lower())
        template_names.append(template_name)
        shutil.copy(template_name, os.path.join(curr_dir, c.TEMPLATES_PATH))
        #shutil.copy("{}.rot.assign".format(resname), os.path.join(curr_dir, c.ROTAMERS_PATH))
    template_initial, template_final = template_names

    # --------------------------------------------GROWING SECTION-------------------------------------------------------
    # Lists definitions

    templates = ["{}_{}".format(template_final, n) for n in range(0, iterations+1)]

    results = ["{}{}_{}".format(c.OUTPUT_FOLDER, resfold, n) for n in range(0, iterations+1)]

    pdbs = [pdb_initialize if n == 0 else "{}_{}".format(n, pdb_initialize) for n in range(0, iterations+1)]

    pdb_selected_names = ["initial_0_{}.pdb".format(n) for n in range(0, cpus-1)]


    # Create a copy of the original templates in growing_templates folder
    shutil.copy(os.path.join(curr_dir, c.TEMPLATES_PATH, template_initial),
                os.path.join(os.path.join(curr_dir, c.TEMPLATES_PATH, c.TEMPLATES_FOLDER), template_initial))

    shutil.copy(os.path.join(curr_dir, c.TEMPLATES_PATH, template_final),
                os.path.join(os.path.join(curr_dir, c.TEMPLATES_PATH, c.TEMPLATES_FOLDER), template_final))

    original_atom = hydrogen_atoms[0].get_name()  # Hydrogen of the core that we will use as growing point
    print(original_atom)
    # Generate starting templates
    replacement_atom = fragment_names_dict[fragment_atom]
    Growing.template_fragmenter.create_initial_template(template_initial, template_final, [original_atom], core_atom,
                                                        replacement_atom, "{}_0".format(template_final),
                                                        os.path.join(curr_dir, c.TEMPLATES_PATH),
                                                        iterations)
    Growing.template_fragmenter.generate_starting_template(template_initial, template_final, [original_atom],core_atom,
                                                           replacement_atom, "{}_ref".format(template_final),
                                                           os.path.join(curr_dir, c.TEMPLATES_PATH),
                                                           iterations)

    # Make a copy in the main folder of Templates in order to use it as template for the simulation
    shutil.copy(os.path.join(curr_dir, c.TEMPLATES_PATH, "{}_0".format(template_final)),
                os.path.join(curr_dir, c.TEMPLATES_PATH, template_final))  # Replace the original template in the folder

    # Simulation loop - LOOP CORE
    for i, (template, pdb_file, result) in enumerate(zip(templates, pdbs, results)):

        pdb_input_paths = ["{}".format(os.path.join(pdbout, str(i-1), pdb_file)) for pdb_file in pdb_selected_names]

        # Control file modification

        if i != 0:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_selected_names, result, i))
            Growing.simulations_linker.control_file_modifier(contrl, pdb_input_paths, i, license, result)
        else:
            logger.info(c.SELECTED_MESSAGE.format(contrl, pdb_initialize, result, i))
            Growing.simulations_linker.control_file_modifier(contrl, [pdb_initialize], i, license, result) #  We have put [] in pdb_initialize
                                                                                                           #  because by default we have to use
                                                                                                           #  a list as input
        logger.info(c.LINES_MESSAGE)
        if i != 0 and i != iterations:
            Growing.template_fragmenter.grow_parameters_in_template("{}_ref".format(template_final),
                                                                    os.path.join(curr_dir, c.TEMPLATES_PATH
                                                                    , c.TEMPLATES_FOLDER, template_initial),
                                                                    os.path.join(curr_dir, c.TEMPLATES_PATH
                                                                    , c.TEMPLATES_FOLDER, template_final),
                                                                    [original_atom], core_atom, replacement_atom,
                                                                    template_final,
                                                                    os.path.join(curr_dir, c.TEMPLATES_PATH),
                                                                    i, iterations)
        elif i == iterations:
            shutil.copy(os.path.join(os.path.join(curr_dir, c.TEMPLATES_PATH, c.TEMPLATES_FOLDER), template_final),
                        os.path.join(os.path.join(curr_dir, c.TEMPLATES_PATH, template_final)))

        # Make a copy of the template file in growing_templates folder
        shutil.copy(os.path.join(curr_dir, c.TEMPLATES_PATH, template_final),
                    os.path.join(os.path.join(curr_dir, c.TEMPLATES_PATH, c.TEMPLATES_FOLDER), template))

        # Running PELE simulation
        if not os.path.exists(c.OUTPUT_FOLDER):
            os.mkdir(c.OUTPUT_FOLDER)
        if not os.path.exists(result):
            os.chdir(c.OUTPUT_FOLDER)
            os.mkdir("{}_{}".format(resfold, (int(i))))
            os.chdir("../")

        Growing.simulations_linker.simulation_runner(pele_dir, contrl, cpus)
        logger.info(c.LINES_MESSAGE)
        logger.info(c.FINISH_SIM_MESSAGE.format(result))
        # Before selecting a step from a trajectory we will save the input PDB file in a folder
        if i == 0:
            if not os.path.exists(pdbout):  # Create the folder if it does not exist
                os.mkdir(pdbout)
            if not os.path.exists(os.path.join(pdbout, "{}".format(i))):  # The same if the subfolder does not exist
                os.chdir(pdbout)
                os.mkdir("{}".format(i))  # Put the name of the iteration
                os.chdir("../../")
        else:
            if not os.path.exists(os.path.join(pdbout, "{}".format(i))):
                os.chdir(pdbout)
                os.mkdir("{}".format(i))
                os.chdir("../")

        # ---------------------------------------------------CLUSTERING-------------------------------------------------
        # Transform column name of the criteria to column number
        result = os.path.abspath(result)
        logger.info("Looking structures to cluster in '{}'".format(result))
        column_number = Helpers.clusterizer.get_column_num(result, criteria, report)
        # Selection of the trajectory used as new input

        Helpers.clusterizer.cluster_traject(str(template_resnames[1]), cpus, column_number, distance_contact,
                                            clusterThreshold, "{}*".format(os.path.join(result, traject)),
                                            os.path.join(pdbout, str(i)), epsilon, report, condition, metricweights,
                                            nclusters)
    # ----------------------------------------------------EQUILIBRATION-------------------------------------------------
    # Set input PDBs
    pdb_inputs = ["{}".format(os.path.join(pdbout, str(iterations), pdb_file)) for pdb_file in pdb_selected_names]
    logger.info(".....STARTING EQUILIBRATION.....")
    if not os.path.exists("equilibration_result"):  # Create the folder if it does not exist
        os.mkdir("equilibration_result")
    # Modify the control file to increase the steps to 20 and change the output path
    Growing.simulations_linker.control_file_modifier(contrl, pdb_inputs, iterations, license,
                                                     "equilibration_result", steps=20)
    # Call PELE to run the simulation
    Growing.simulations_linker.simulation_runner(pele_dir, contrl, cpus)
    equilibration_path = os.path.join(os.path.curdir, "equilibration_result")
    Growing.bestStructs.main(criteria, "best_structure.pdb", path=equilibration_path,
                             n_structs=10)


if __name__ == '__main__':
    complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python, pele_dir, \
    contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon, condition, metricweights, \
    nclusters, rotres = parse_arguments()

    main(complex_pdb, fragment_pdb, core_atom, fragment_atom, iterations, criteria, plop_path, sch_python, pele_dir,
         contrl, license, resfold, report, traject, pdbout, cpus, distcont, threshold, epsilon, condition,
         metricweights, nclusters, rotres)
