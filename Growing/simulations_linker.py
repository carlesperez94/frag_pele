# General imports
import os
import shutil
import subprocess
import logging

# Local imports
import Helpers.templatize

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def control_file_modifier(control_template, pdb, step, results_path="/growing_output",
                          license="/opt/PELErev12492/licenses"):
    """
    This function creates n control files for each intermediate template created in order to change
    the logPath, reportPath and trajectoryPath to have all control files prepared for PELE simulations.
    """
    # Obtain the path of the directory
    PATH = os.getcwd()
    print(step)
    # Then, in the main loop we will do a copy of control files, so we will print this in the logger
    control_path = os.path.join(PATH, "control_folder")
    logger.info("Intermediate control files created will be stored in '{}'".format(control_path))

    # Definition of the keywords that we are going to substitute from the template
    keywords = {"LICENSE": license,
                "RESULTS_PATH": results_path,
                "PDB": pdb
                }
    # Creation of a folder where we are going to contain our control files, just if needed
    if not os.path.exists(os.path.join(PATH, "control_folder")):
        os.mkdir(os.path.join(PATH, "control_folder"))
    else:
        pass
    # Changing the name of the results folder in the dictionary
    keywords.update({"RESULTS_PATH": "{}_{}".format(os.path.join(PATH, results_path), step)})

    # Create a copy of the control template in the control folder, because templatize.py replace the original template
    if not os.path.exists(os.path.join(control_path, control_template)):
        shutil.copyfile(control_template, os.path.join(control_path, control_template))
    # Else, if has been created this means that we already have a template in this folder, so we will need a copy of the
    # file in the main folder to then replace the template for a real control file
    else:
        shutil.copyfile(os.path.join(control_path, control_template), control_template)
    # Modifying the control file template
    Helpers.templatize.TemplateBuilder(control_template, keywords)
    # Make a copy in the control files folder
    shutil.copyfile(control_template, os.path.join(control_path, "{}_{}".format(step, control_template)))
    logger.info("{}_{} has been created successfully!".format(step, control_template))

def simulation_runner(path_to_pele, control_in):
    """
    Runs a PELE simulation with the parameters described in the input control file.

    Input:

    path_to_pele --> Complete path to Pele_serial

    control_in --> Name of the control file with the parameters to run PELE
    """

    subprocess.call(["{}".format(path_to_pele), "{}".format(control_in)])




    




