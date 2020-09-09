import sys
import os
import socket

DIR = os.path.dirname(__file__)

##############################################

# PUBLIC CONSTANTS (to change by the user)
# Preparation inputs to grow

machine = socket.getfqdn()
print(machine)

# Paths definitions (IMPORTANT!)
if "bsc.mn" in machine:
    # PELE parameters
    PATH_TO_PELE = "/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-SideChainPert/bin/PELE-1.6.2_mpi"
    PATH_TO_PELE_DATA = "/gpfs/projects/bsc72/PELE++/mniv/V1.6.2-SideChainPert/Data"
    PATH_TO_PELE_DOCUMENTS = None # Data and documents is added automatically from PELE >= 1.6
    PATH_TO_LICENSE = "/gpfs/projects/bsc72/PELE++/license"
    # PlopRotTemp parameters
    SCHRODINGER_PY_PATH = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/utilities/python"
    ENV_PYTHON = "/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC/internal/lib/python2.7/site-packages/"
elif "bsccv" in machine:
    # PELE parameters
    PATH_TO_PELE = "/data/EAPM/PELE/PELE++/bin/rev12360/Pele_rev12360_mpi"
    PATH_TO_PELE_DATA = "/data/EAPM/PELE/PELE++/data/rev12360/Data"
    PATH_TO_PELE_DOCUMENTS = "/data/EAPM/PELE/PELE++/Documents/rev12360"
    PATH_TO_LICENSE = "/data/EAPM/PELE/PELE++/license"
    SCHRODINGER_PY_PATH = "/data2/bsc72/SCHRODINGER_ACADEMIC/utilities/python"
    ENV_PYTHON = "/data2/bsc72/SCHRODINGER_ACADEMIC/internal/lib/python2.7/site-packages/"
else:
    # PELE parameters
    PATH_TO_PELE = "/home/user/pelepath/PELE_mpi"
    PATH_TO_PELE_DATA = None
    PATH_TO_PELE_DOCUMENTS = None
    PATH_TO_LICENSE = "/home/user/pelepath/licenses"
    SCHRODINGER_PY_PATH = "/home/user/schrodingerVVVV-VV/run"
    ENV_PYTHON = ""

# FragPELE configuration
CONTROL_TEMPLATE = os.path.join(DIR, "Templates/control_template.conf")
GROWING_STEPS = 10
SELECTION_CRITERIA = "Binding Energy"
CPUS = 48
PELE_EQ_STEPS = 20  # PELE steps in the equilibration
RESTART = False
STEPS = 6  # PELE steps for growing step
BANNED_DIHEDRALS_ATOMS = None
BANNED_ANGLE_THRESHOLD = None

# PELE control file configuration
REPORT_NAME = "report"
TRAJECTORY_NAME = "trajectory"
RESULTS_FOLDER = "growing_output"
TEMPERATURE = 1000
MAX_OVERLAP = 0.70
MIN_OVERLAP = 0.50
SEED = 1279183
STEERING = 0
TRANSLATION_HIGH = 0.05 #Amstrongs
ROTATION_HIGH = 0.10 # Radiants
TRANSLATION_LOW = 0.02
ROTATION_LOW = 0.05 
RADIUS_BOX = 4  # Amstrongs
# PlopRotTemp configuration
ROTRES = "10.0"

# Clustering parameters
DISTANCE_COUNTER = 4
CONTACT_THRESHOLD = 0.3
EPSILON = 0.5
CONDITION = "min"   #   min or max
METRICS_WEIGHTS = "linear"
NUM_CLUSTERS = 5
##############################################

# PRIVATE CONSTANTS (not to change)
PRE_WORKING_DIR = "pregrow"
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
ROTAMERS_PATH = "DataLocal/LigandRotamerLibs/"
PDBS_OUTPUT_FOLDER = "clustering_PDBs"
OUTPUT_FOLDER = "growing_steps/"
TEMPLATES_FOLDER = "growing_templates"
CONFIG_PATH = "log_configure.ini"
PLOP_PATH = "PlopRotTemp_S_2017/ligand_prep.py"

# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
LINES_MESSAGE = "\n•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.••*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•\n"
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\nStep: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "

# LISTS
AA_LIST = ["ALA", "ARG", "ASP", "CYS", "GLN", "GLY", "GLU", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYS", 
           "MET", "PHE", "PRO", "TRP", "TYR", "THR", "VAL"]
##############################################
