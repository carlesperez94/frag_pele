# General imports
import os
import string
import logging
import shutil
# Local imports
import Ligand_growing.Helpers.templatize

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
    #if not os.path.exists(os.path.join(PATH, "control_folder")):
    #    os.mkdir(os.path.join(PATH, "control_folder"))
    #else:
    #    pass
    # Then, we are going to put all files inside this folder, so we will create another path variable
    #working_path = os.path.join(PATH, "control_folder")

    # As we have several templates that uses as input different PDBs we have to do a loop to create this
    # different templates

    for n in range(0, n_files):

        # First, we are going to create n_files copies of the template in order to keep the original unmodified.
        control_renamed = os.path.join(PATH, "control_growing_{}.conf".format(string.ascii_lowercase[n]))
        shutil.copyfile(control_file, control_renamed)

        # Changing the name of the pdb for each control file
        keywords.update({"RESULTS_PATH": "{}_{}".format(results_path, string.ascii_lowercase[n]),
                         "PDB": "{}_{}.pdb".format(pdb, string.ascii_lowercase[n])})

        # Modifying the control file template for each step
        Helpers.templatize.TemplateBuilder(control_renamed, keywords)


def simulation_runner(control_in):
    #counter=1
    #quantity will depend on the amount of simulations that we will want to do
    #while counter<quantity:
    print(control_in)
    os.system("/opt/PELErev12492/bin/Pele_serial {}".format(control_in))

    #print("\n\n----------____SIMULATION {} FINISHED____-------------\n\n".format(counter))
    #counter=1+counter



    




