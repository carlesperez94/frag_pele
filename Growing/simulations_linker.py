# General imports
import os
import string
import logging
import shutil
import subprocess
# Local imports
import Helpers.templatize

logging.basicConfig(filename="output.log",format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)


def control_file_modifier(control_file, pdb, results_path="/growing_output", n_files = 10, license="/opt/PELErev12492/licenses"):
    """
    This function creates n control files for each intermediate template created in order to change
    the logPath, reportPath and trajectoryPath to have all control files prepared for PELE simulations.
    """
    # Obtain the path of the directory

    PATH = os.getcwd()

    # Definition of the keywords that we are going to substitute from the template
    keywords = {"LICENSE": license,
                "RESULTS_PATH": results_path,
                "PDB": pdb
                }
    # Creation of a folder where we are going to contain our control files
    if not os.path.exists(os.path.join(PATH, "control_folder")):
        os.mkdir(os.path.join(PATH, "control_folder"))
    else:
       pass
    # Then, we are going to put all files inside this folder, so we will create another path variable
    control_path = os.path.join(PATH, "control_folder")

    # As we have several templates that uses as input different PDBs we have to do a loop to create this
    # different templates

    for n in range(0, n_files):

        # First, we are going to create n_files copies of the template in order to keep the original unmodified.
        control_renamed = os.path.join(control_path, "control_growing_{}.conf".format(string.ascii_lowercase[n]))
        shutil.copyfile(control_file, control_renamed)

        # Changing the name of the pdb for each control file
        keywords.update({"RESULTS_PATH": "{}_{}".format(os.path.join(PATH, results_path), string.ascii_lowercase[n]),
                         "PDB": "{}_{}.pdb".format(pdb, string.ascii_lowercase[n])})

        # Modifying the control file template for each step
        Helpers.templatize.TemplateBuilder(control_renamed, keywords)


def simulation_runner(path_to_pele, control_in):
    """
    Runs a PELE simulation with the parameters described in the input control file.

    Input:

    path_to_pele --> Complete path to Pele_serial

    control_in --> Name of the control file with the parameters to run PELE
    """

    subprocess.call(["{}".format(path_to_pele), "{}".format(control_in)])




    




