import sys

# PUBLIC CONSTANTS (to change by the user)
IN_TEMPLATE = "mbez"
FIN_TEMPLATE = "pyjz"
CONTROL_TEMPLATE = "control_template.conf"
PDB = "181l_prepared_met.pdb"
ITERATIONS = 10
RESULTS_FOLDER = "growing_output"
SELECTION_CRITERIA = "Binding Energy"
PATH_TO_PELE = "/opt/PELErev12492"
REPORT_NAME = "report"
TRAJECTORY_NAME = "trajectory.pdb"
CONFIG_PATH = "/home/carlespl/project/growing/Ligand_growing/log_configure.ini"

# PRIVATE CONSTANTS (not to change)
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
ORIGINAL_ATOM = "_H8_"
FINAL_ATOM = "_C8_"
PDBS_OUTPUT_FOLDER = "PDBs_growing"
# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "