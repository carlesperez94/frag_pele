import os
import logging
import constants as c


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_and_create_DataLocal():
    # Get current path
    curr_dir = os.path.abspath(os.path.curdir)
    if not os.path.exists(os.path.join(curr_dir, "DataLocal")):
        os.mkdir(os.path.join(curr_dir, "DataLocal"))
    if not os.path.exists(os.path.join(curr_dir, "DataLocal/Templates")):
        os.mkdir(os.path.join(curr_dir, "DataLocal/Templates"))
    if not os.path.exists(os.path.join(curr_dir, "DataLocal/LigandRotamerLibs")):
        os.mkdir(os.path.join(curr_dir, "DataLocal/LigandRotamerLibs"))
    if not os.path.exists(os.path.join(curr_dir, "DataLocal/LigandRotamerLibs")):
        os.mkdir(os.path.join(curr_dir, "DataLocal/LigandRotamerLibs"))
    if not os.path.exists(os.path.join(curr_dir, "DataLocal/Templates/OPLS2005")):
        os.mkdir(os.path.join(curr_dir, "DataLocal/Templates/OPLS2005"))
    if not os.path.exists(os.path.join(curr_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")):
        os.mkdir(os.path.join(curr_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms"))


def check_and_create_results_folder(results_folder):
    if not os.path.exists(c.OUTPUT_FOLDER):
        os.mkdir(c.OUTPUT_FOLDER)
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)


def check_and_create_pdb_clusters_folder(pdbout_folder, iteration):
    if iteration == 0:
        if not os.path.exists(pdbout_folder):  # Create the folder if it does not exist
            os.mkdir(pdbout_folder)
        if not os.path.exists(os.path.join(pdbout_folder, "{}".format(iteration))):  # Do the same if the subfolder does not exist
            os.mkdir((os.path.join(pdbout_folder, "{}".format(iteration))))  # Put the name of the iteration
    else:
        if not os.path.exists(os.path.join(pdbout_folder, "{}".format(iteration))):
            os.mkdir(os.path.join(pdbout_folder, "{}".format(iteration)))
