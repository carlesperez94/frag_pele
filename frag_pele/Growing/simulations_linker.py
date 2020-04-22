# General imports
import os
import shutil
import subprocess
import logging

# Local imports
import frag_pele.Helpers.templatize as tp

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def control_file_modifier(control_template, pdb, license, working_dir, overlap=0.7, step=0,
                          results_path="/growing_output", steps=6, chain="L", constraints=" ", center="",
                          temperature=1000, seed=1279183, steering=0, translation_high=0.05, translation_low=0.02,
                          rotation_high=0.10, rotation_low=0.05, radius=4):
    """
    This function creates n control files for each intermediate template created in order to change
    the logPath, reportPath and trajectoryPath to have all control files prepared for PELE simulations.
    """

    ctrl_fold_name = os.path.join(working_dir, "control_folder")
    controlfile_path = control_template
    control_template = os.path.basename(control_template)

    # Then, in the main loop we will do a copy of control files, so we will print this in the logger
    logger.info("Intermediate control files created will be stored in '{}'".format(ctrl_fold_name))

    # As pdb is a list of complexes, we have to create one line per complex in pdb and then put them all in the control
    list_of_lines_complex = []
    for complex in pdb:
        control_file_complex = '{"files" : [{"path": "%s" }] }' % (complex)
        list_of_lines_complex.append(control_file_complex)
    lines_complex = ",\n".join(list_of_lines_complex)

    # Definition of the keywords that we are going to substitute from the template
    keywords = {"LICENSE": license,
                "RESULTS_PATH": results_path,
                "CHAIN": chain,
                "CONSTRAINTS": constraints,
                "CENTER": center, 
                "PDB": lines_complex,
                "STEPS": steps,
                "OVERLAP": overlap,
                "TEMPERATURE": temperature,
                "SEED": seed,
                "STEERING": steering,
                "TRANSLATION_HIGH": translation_high,
                "TRANSLATION_LOW": translation_low,
                "ROTATION_HIGH": rotation_high,
                "ROTATION_LOW": rotation_low,
                "RADIUS": radius
                }

    # Creation of a folder where we are going to contain our control files, just if needed
    if not os.path.exists(ctrl_fold_name):
        os.mkdir(ctrl_fold_name)

    # Create a copy of the control template in the control folder, because templatize.py replace the original template
    if not os.path.exists(os.path.join(ctrl_fold_name, control_template)):
        shutil.copyfile(controlfile_path, os.path.join(ctrl_fold_name, control_template))
        shutil.copyfile(os.path.join(ctrl_fold_name, control_template), os.path.join(working_dir, control_template))

    # Else, if has been created this means that we already have a template in this folder, so we will need a copy of the
    # file in the main folder to then replace the template for a real control file
    else:
        shutil.copyfile(os.path.join(ctrl_fold_name, control_template), os.path.join(working_dir, control_template))

    # Modifying the control file template
    tp.TemplateBuilder(os.path.join(working_dir, control_template), keywords)
    # Make a copy in the control files folder
    simulation_file = os.path.join(ctrl_fold_name, "{}_{}".format(step, control_template))
    shutil.copyfile(os.path.join(working_dir, control_template), simulation_file)
    logger.info("{}_{} has been created successfully!".format(step, control_template))

    return simulation_file


def simulation_runner(path_to_pele, control_in, cpus=4):
    """
    Runs a PELE simulation with the parameters described in the input control file.

    Input:

    path_to_pele --> Complete path to PELE folder

    control_in --> Name of the control file with the parameters to run PELE
    """
    if cpus:
        cpus = int(cpus)
        if cpus < 2:
            logger.critical("Sorry, to run mpi PELE you need at least 2 CPUs!")
        else:
            logger.info("Starting PELE simulation. You will run mpi PELE with {} cores.".format(cpus))
            cmd = "mpirun -np {} {} {}".format(cpus, path_to_pele, control_in)
            print(cmd)
            logger.info("Running {}".format(cmd))
            subprocess.call(cmd.split())
    else:
        logger.info("Starting PELE simulation. You will run serial PELE.")
        cmd = "{} {}".format(path_to_pele, control_in)
        logger.info("Running {}".format(cmd))
        subprocess.call(cmd.split())
