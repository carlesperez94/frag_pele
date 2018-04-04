import sys

##############################################

# PUBLIC CONSTANTS (to change by the user)
# Preparation inputs to grow
COMPLEX_PDB = "4e20_prep.pdb"
FRAGMENT_PDB = "frag.pdb"
CORE_ATOM = "N3"
FRAGMENT_ATOM = "C7"
ITERATIONS = 10
SELECTION_CRITERIA = "Binding Energy"
CONFIG_PATH = "/home/carlespl/project/growing/Ligand_growing/log_configure.ini"
# PELE parameters
CONTROL_TEMPLATE = "control_template.conf"
RESULTS_FOLDER = "growing_output"
PATH_TO_PELE = "/opt/PELErev12492"
REPORT_NAME = "report"
TRAJECTORY_NAME = "trajectory.pdb"
CPUS = 4
# PlopRotTemp parameters
PLOP_PATH = "/home/carlespl/project/growing/Ligand_growing/PlopRotTemp_S_2017/ligand_prep.py"
SCHRODINGER_PY_PATH = "/opt/schrodinger2016-4/utilities/python"

##############################################

# PRIVATE CONSTANTS (not to change)
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
ROTAMERS_PATH = "DataLocal/LigandRotamerLibs/"
PDBS_OUTPUT_FOLDER = "PDBs_growing"
OUTPUT_FOLDER = "growing_results/"
TEMPLATES_FOLDER = "growing_templates"
# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\nStep: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "

##############################################
