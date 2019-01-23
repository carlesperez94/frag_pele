import os
import logging
import constants as c


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_and_create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)


def check_and_create_DataLocal():
    # Get current path
    curr_dir = os.path.abspath(os.path.curdir)
    check_and_create_folder((os.path.join(curr_dir, "DataLocal")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/LigandRotamerLibs")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates/OPLS2005")))
    check_and_create_folder((os.path.join(curr_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")))


def check_and_create_results_folder(results_folder):
    check_and_create_folder(c.OUTPUT_FOLDER)
    check_and_create_folder(results_folder)


def check_and_create_pdb_clusters_folder(pdbout_folder, iteration):
    if iteration == 0:
        check_and_create_folder(pdbout_folder)
        check_and_create_folder(os.path.join(pdbout_folder, "{}".format(iteration)))
    else:
        check_and_create_folder(os.path.join(pdbout_folder, "{}".format(iteration)))
