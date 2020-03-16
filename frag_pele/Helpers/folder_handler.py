import os
import logging
# Local import
import frag_pele.constants as c


# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def check_and_create_folder(path):
    if not os.path.exists(path):
        os.mkdir(path)


def check_and_create_DataLocal(working_dir):
    # Get current path
    check_and_create_folder((os.path.join(working_dir, "DataLocal")))
    check_and_create_folder((os.path.join(working_dir, "DataLocal/LigandRotamerLibs")))
    check_and_create_folder((os.path.join(working_dir, "DataLocal/Templates")))
    check_and_create_folder((os.path.join(working_dir, "DataLocal/Templates/OPLS2005")))
    check_and_create_folder((os.path.join(working_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms")))
    check_and_create_folder((os.path.join(working_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms", "templates_generated")))


def check_and_create_results_folder(results_folder, workdir):
    check_and_create_folder(os.path.join(workdir, c.OUTPUT_FOLDER))
    check_and_create_folder(results_folder)


def check_and_create_pdb_clusters_folder(pdbout_folder, iteration):
    if iteration == 0:
        check_and_create_folder(pdbout_folder)
        check_and_create_folder(os.path.join(pdbout_folder, "{}".format(iteration)))
    else:
        check_and_create_folder(os.path.join(pdbout_folder, "{}".format(iteration)))
