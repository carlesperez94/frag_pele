import sys
import os
import socket

DIR = os.path.dirname(__file__)

##############################################

# PUBLIC CONSTANTS (to change by the user)
# Preparation inputs to grow

machine = socket.getfqdn()
if "bsc.mn" in machine:
    # PELE parameters
    PATH_TO_PELE = "/gpfs/projects/bsc72/PELE++/mniv/rev12536/bin/Pele_mpi"
    PATH_TO_PELE_DATA = "/gpfs/projects/bsc72/PELE++/data/rev12360/Data"
    PATH_TO_PELE_DOCUMENTS = "/gpfs/projects/bsc72/PELE++/Documents/rev12360"
    PATH_TO_LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    # PlopRotTemp parameters
    SCHRODINGER_PY_PATH = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/utilities/python"
else:
    # PELE parameters
    PATH_TO_PELE = "$PELE_BIN"
    PATH_TO_PELE_DATA = os.path.join("$PELE", "Data")
    PATH_TO_PELE_DOCUMENTS = os.path.join("$PELE", "Documents")
    PATH_TO_LICENSE = "$LICENSE"
    SCHRODINGER_PY_PATH = os.path.join("$SCHRODINGER", "utilities/python")

CONTROL_TEMPLATE = os.path.join(DIR, "Templates/control_template.conf")
RESULTS_FOLDER = "growing_output"
GROWING_STEPS = 10
SELECTION_CRITERIA = "Binding Energy"
SERIE_FILE = False
REPORT_NAME = "report"
TRAJECTORY_NAME = "trajectory"
CPUS = 48
PELE_EQ_STEPS = 20
RESTART = False
STEPS = 6
TEMPERATURE = 1000
MAX_OVERLAP = 0.70
MIN_OVERLAP = 0.50
TEMPERATURE = 1000
# Clustering parameters
DISTANCE_COUNTER = 4
CONTACT_THRESHOLD = 0.3
EPSILON = 0.5


##############################################

# PRIVATE CONSTANTS (not to change)
PRE_WORKING_DIR = "pregrow"
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
ROTAMERS_PATH = "DataLocal/LigandRotamerLibs/"
PDBS_OUTPUT_FOLDER = "PDBs_growing"
OUTPUT_FOLDER = "growing_results/"
TEMPLATES_FOLDER = "growing_templates"
CONFIG_PATH = "log_configure.ini"
PLOP_PATH = "PlopRotTemp_S_2017/ligand_prep.py"
ROTRES = 30
# Clustering constants
CONDITION = "min"   #   min or max
METRICS_WEIGHTS = "linear"
NUM_CLUSTERS = 5
# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
LINES_MESSAGE = "\n•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.••*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•\n"
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\nStep: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "

##############################################
