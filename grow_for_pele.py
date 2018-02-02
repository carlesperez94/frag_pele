# General imports
import sys
import argparse
import string
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
CONTROL_PATH = "control_folder/control_growing_{}.conf"

# Messages constants
TEMPLATE_MESSAGE = "We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting..."
SELECTED_MESSAGE = "\n============ Files selected ============\nControl file: {}\nPDB file: {}\nResults folder name: {}\n"
FINISH_SIM_MESSAGE = "SIMULATION FOR control_growing_{} COMPLETED!!! "


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
    parser.add_argument("-d", "--dir", default="/opt/PELErev12492/bin/Pele_serial",
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

        "template_initial" --> Name of the input file correspondent to the initial template for the ligand that you
                               want to grow.

        "template_final" --> Name of the input file correspondent to the final template for the ligand that you want
                             to get.

        "n_files" --> Number of intermediate templates that you want to generate

        "transformation" --> When an atom is transformed into another one we want to conserve the properties.
                             The atom name of the initial template that we want to transform into another of
                             the final template has to be specified in a text file separated by whitespaces.

        "control_file" --> Template control file used to generate intermediates control files.

        "pdb" --> Initial pdb file which already contain the ligand with the fragment that we want to grow
                  but with bond lengths correspondent to the initial ligand (dummy-like).

        "results_f_name" --> Name for results folder

        "criteria" --> Name of the column of the report file used as criteria in order to select the template
                       used as input for successive simulations.

        "path_pele" --> Complete path to Pele_serial

        Output:

        First, intermediate control files and templates. Then, the results for each simulation and a pdb file for
        each selected trajectory (the last selected trajectory is the final structure with the ligand grown).

        """
    # Main loop

    for n in range(0, n_files):
        print(n)

        # Template creation part
        if n == 0:
            # If we are starting the main loop we have to generate a starting template
            Growing.template_fragmenter.generate_starting_template(template_initial, template_final, original_atom,
                                                                   final_atom, template_final, n_files)
            # Save a copy of the file in growing_templates
            shutil.copy(os.path.join(TEMPLATES_PATH, template_final),
                                     os.path.join(os.path.join(TEMPLATES_PATH, "growing_templates"),
                                     "{}_{}".format(template_final,n)))
        if n > 0 and n < n_files:
            # Now, we have already create this starting template so we just need to increase their values
            Growing.template_fragmenter.grow_parameters_in_template(template_initial, template_final,
                                                                    original_atom, final_atom,
                                                                    template_final, n)
        # We have to end up with the final template, so in the last step we will copy it into the templates folder
        if n == n_files:
            shutil.copy(os.path.join(os.path.join(TEMPLATES_PATH, "growing_templates"),
                        template_final), os.path.join(TEMPLATES_PATH, template_final))
        # Control file modification
        Growing.simulations_linker.control_file_modifier(control_template, pdb, n, results_f_name)
        logger.info(SELECTED_MESSAGE.format(control_template, pdb, results_f_name))

        # Running PELE simulation
        if not os.path.exists("{}_{}".format(results_f_name, n)):
            os.mkdir("{}_{}".format(results_f_name, n))
        else:
            pass
        Growing.simulations_linker.simulation_runner(path_pele, CONTROL_PATH)
        logger.info(FINISH_SIM_MESSAGE.format(n))

        # Selection of the trajectory used as new input
        Growing.template_selector.trajectory_selector(pdb, "{}".format(results_f_name), report, traject, criteria)

        logger.info("Step of the Trajectory selected in {}_{}.pdb".format(pdb, string.ascii_lowercase[n + 1]))


if __name__ == '__main__':
    init, final, frag, control, original_atom, final_atom, pdb, res_fold, criteria, path_pele, report, traject = parse_arguments()
    main(init, final, frag, control, original_atom, final_atom, pdb, res_fold, criteria, path_pele, report, traject)
