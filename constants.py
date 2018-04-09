import sys

##############################################

# PUBLIC CONSTANTS (to change by the user)
# Preparation inputs to grow
ITERATIONS = 10
SELECTION_CRITERIA = "Binding Energy"
# PELE parameters
CONTROL_TEMPLATE = "control_template.conf"
RESULTS_FOLDER = "growing_output"
PATH_TO_PELE = "/opt/PELErev12492/bin/Pele_mpi"
PATH_TO_LICENSE = "/opt/PELErev12492/licenses"
REPORT_NAME = "report"
TRAJECTORY_NAME = "trajectory.pdb"
CPUS = 4
N_INI_STRUCTURES = 2
# PlopRotTemp parameters
SCHRODINGER_PY_PATH = "/opt/schrodinger2016-4/utilities/python"

##############################################

# PRIVATE CONSTANTS (not to change)
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
ROTAMERS_PATH = "DataLocal/LigandRotamerLibs/"
PDBS_OUTPUT_FOLDER = "PDBs_growing"
OUTPUT_FOLDER = "growing_results/"
TEMPLATES_FOLDER = "growing_templates"
CONFIG_PATH = "log_configure.ini"
PLOP_PATH = "PlopRotTemp_S_2017/ligand_prep.py"
# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
LINES_MESSAGE = "•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.••*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•*´¨`*•.¸¸.•"
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\nStep: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "

##############################################
