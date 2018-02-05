# General imports
import sys
import argparse
import os
import logging
from logging.config import fileConfig
import shutil
# Local imports
import Growing.template_selector
import Growing.template_fragmenter
import Growing.simulations_linker

# Calling configuration file for log system
fileConfig("/home/carlespl/project/growing/Ligand_growing/log_configure.ini")

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Path variables
TEMPLATES_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"

# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION {} COMPLETED!!! "


def parse_arguments():
    """
        Parse user arguments

        Output: list with all the user arguments
    """
    parser = argparse.ArgumentParser(description="""From an input file, correspondent to the template of the initial 
    structure of the ligand,this function will generate "x_fragments" intermediate templates until reach the final 
    template,modifying Vann Der Waals, bond lengths and deleting charges.""")

    required_named = parser.add_argument_group('required named arguments')  # We have to do this in order
    # to print properly the required arguments when using the defined arguments method

    # Required arguments
    required_named.add_argument("-i", "--initial", required=True,
                                help="""Input file correspondent to the 
                                initial template for the ligand that you 
                                want to grow.""")
    required_named.add_argument("-oa", "--original_atom", required=True,
                                help="""When an atom is transformed into another
                                one we want to conserve the properties
                                of the first one until being changed in the 
                                last template. The PDB atom name of the first atom
                                that we want to transform """)
    required_named.add_argument("-fa", "--final_atom", required=True,
                                help="""When an atom is transformed into another
                                one we want to conserve the properties
                                of the first one until being changed in the 
                                last template. The PDB atom name of the final atom
                                that will have their properties replaced.
                                """)
    required_named.add_argument("-c", "--contrl", required=True,
                                help='Initial control file.')
    required_named.add_argument("-p", "--pdb", required=True,
                                help="""Initial pdb file which already contain the ligand with 
                                the fragment that we want to grow but with bond lengths correspondent 
                                to the initial ligand (dummy-like).""")

    # In the future we will remove this argument
    required_named.add_argument("-f", "--final", required=True,
                                help="""Input file correspondent to the
                                final template for the ligand that you 
                                want to get.""")

    # Optional arguments
    parser.add_argument("-x", "--frag", type=int, default=10,
                        help="""Number of intermediate templates that you want 
                             to generate""")
    parser.add_argument("-r", "--resfold", default="growing_output",
                        help="Name for results folder")
    parser.add_argument("-cr", "--criteria", default="Binding Energy",
                        help="""Name of the column used as criteria in order
                             to select the template used as input for 
                             successive simulations.""")
    parser.add_argument("-d", "--dir", default="/opt/PELErev12492",
                        help="Complete path to Pele_serial")
    parser.add_argument("-rp", "--report", default="report",
                        help="Name of the report file from PELE.")
    parser.add_argument("-tj", "--traject", default="trajectory.pdb",
                        help="Name of the trajectory file from PELE.")
    args = parser.parse_args()

    return args.initial, args.final, args.frag, args.contrl, args.original_atom, args.final_atom, args.pdb, args.resfold, args.criteria, args.dir, args.report, args.traject


def main(template_initial, template_final, n_files, control_template, original_atom, final_atom, pdb, results_f_name,
         criteria, path_pele, report, traject):
    """
        Description: This function is the main core of the program. It creates N intermediate templates
        and control files for PELE. Then, it perform N successive simulations automatically.

        Input:

        :param: template_initial: Name of the input file correspondent to the initial template for the ligand that you
        want to grow.
        :param template_final: Name of the input file correspondent to the final template for the ligand that you want
        to get.
        :param n_files: Number of intermediate templates that you want to generate
        :param original_atom: When an atom is transformed into another one we want to conserve the properties.
        The PDB atom name of the initial template that we want to transform into another of the final template has to be
        specified here.
       :param final_atom: When an atom is transformed into another one we want to conserve the properties.
        The PDB atom name of the final template that will be used as target to be transformed into the original atom.
        :param control_template: Template control file used to generate intermediates control files.
        :param pdb: Initial pdb file which already contain the ligand with the fragment that we want to grow
        but with bond lengths correspondent to the initial ligand (dummy-like).
        :param results_f_name: Name for results folder.
        :param criteria: Name of the column of the report file used as criteria in order to select the template
        used as input for successive simulations.
        :param path_pele: Complete path to Pele_serial.
        :param report: Name of the report file of PELE simulation's results.
        :param traject: Name of the trajectory file of PELE simulation's results.
        :return The algorithm is formed by two main parts:
        First, it prepare the files needed in a PELE simulation: control file and ligand templates.
        Second, after analyzing the report file it selects the best structure in the trajectory file in order to be used
        as input for the next growing step.
        This process will be repeated until get the final grown ligand (protein + ligand).
        """
    # Main loop

    templates = [template_final if n == (n_files-1) else "{}_{}".format(template_final, n) for n in range(0, n_files)]

    results = ["{}_{}".format(results_f_name, n) for n in range(0, n_files) ]

    pdbs = [pdb if n == 0 else "{}_{}".format(n, pdb) for n in range(0, n_files)]

    if not os.path.exists("PDBs_growing"):
        os.mkdir("PDBs_growing")

    # Create a copy of the original templates in growing_templates folder
    shutil.copy(os.path.join(TEMPLATES_PATH, template_initial),
                os.path.join(os.path.join(TEMPLATES_PATH, "growing_templates"), template_initial))
    shutil.copy(os.path.join(TEMPLATES_PATH, template_final),
                os.path.join(os.path.join(TEMPLATES_PATH, "growing_templates"), template_final))

    Growing.template_fragmenter.generate_starting_template(template_initial, template_final, original_atom,
                                                           final_atom, "{}_0".format(template_final),
                                                           n_files)

    for i, (template, pdb_file, result) in enumerate(zip(templates, pdbs, results)):

        # Control file modification
        logger.info(SELECTED_MESSAGE.format(control_template, pdb, result))
        Growing.simulations_linker.control_file_modifier(control_template, pdb, result, path_pele)

        if i != 0:
            Growing.template_fragmenter.grow_parameters_in_template("{}_0".format(template_final),
                                                                    os.path.join("growing_templates", template_initial),
                                                                    os.path.join("growing_templates", template_final),
                                                                    original_atom, final_atom, template_final, i)

            # Make a copy of the template file in growing_templates folder
        shutil.copy(os.path.join(TEMPLATES_PATH, template_final),
                    os.path.join(os.path.join(TEMPLATES_PATH, "growing_templates"), template))


        # Running PELE simulation
        if not os.path.exists(result):
            os.mkdir(result)

        logger.info(FINISH_SIM_MESSAGE.format(result))
        Growing.simulations_linker.simulation_runner(path_pele, control_template)

        # Before selecting a step from a trajectory we will save the input PDB file in a folder
        shutil.copy(pdb, os.path.join("PDBs_growing", pdb_file))

        # Selection of the trajectory used as new input
        Growing.template_selector.trajectory_selector(pdb, result, report, traject, criteria)


if __name__ == '__main__':
    init, final, frag, control, original_atom, final_atom, pdb, res_fold, criteria, path_pele, report, traject = parse_arguments()
    main(init, final, frag, control, original_atom, final_atom, pdb, res_fold, criteria, path_pele, report, traject)
