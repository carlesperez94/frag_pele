# General imports
import sys
import argparse
import string
import os
import logging
from logging.config import fileConfig
# Local imports
import Growing.template_selector
import Growing.template_fragmenter
import Growing.simulations_linker

# Calling configuration file for log system
fileConfig("/home/carlespl/project/growing/Ligand_growing/log_configure.ini")

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Path variables
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
                                help="""input file correspondent to the 
                                initial template for the ligand that you 
                                want to grow.""")
    required_named.add_argument("-t", "--trans", required=True,
                                help="""When an atom is transformed into another
                                one we want to conserve the properties
                                of the first one until being changed in the 
                                last template. The atom name of the initial template
                                that we want to transform into another of the final
                                template has to be specified in 
                                a text file separated by whitespaces.""")
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

    return args.initial, args.final, args.frag, args.trans, args.contrl, args.pdb, args.resfold, args.criteria, args.dir, args.report, args.traject


def main(template_initial, template_final, n_files, transformation, control_file, pdb, results_f_name, criteria, path_pele, report, traject):
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
    # Creating template files
    logger.info((TEMPLATE_MESSAGE.format(template_initial, template_final, n_files)))
    templates = Growing.template_fragmenter.fragmenter(template_initial, template_final, transformation, n_files)
    # Creating control files
    logger.info(SELECTED_MESSAGE.format(control_file, pdb, results_f_name, n_files))
    control_files = Growing.simulations_linker.control_file_modifier(control_file, pdb, results_f_name, n_files)

    # Run Pele for each control file

    # for template, control_file in zip(templates, control_files):
    for n in range(0, n_files):
        # Run Pele
        if not os.path.exists("{}_{}".format(results_f_name, string.ascii_lowercase[n])):
            os.mkdir("{}_{}".format(results_f_name, string.ascii_lowercase[n]))
        else:
            pass
        Growing.simulations_linker.simulation_runner(path_pele, CONTROL_PATH.format(string.ascii_lowercase[n]))
        logger.info(FINISH_SIM_MESSAGE.format(string.ascii_lowercase[n]))

        # Choose the best trajectory
        Growing.template_selector.trajectory_selector("{}_{}_tmp.pdb".format(pdb, string.ascii_lowercase[n+1]),
                                                      "{}_{}".format(results_f_name, string.ascii_lowercase[n]),
                                                      report, traject, criteria)
        Growing.template_selector.change_ligandname("{}_{}_tmp.pdb".format(pdb, string.ascii_lowercase[n + 1]),
                                                    "{}_{}.pdb".format(pdb, string.ascii_lowercase[n + 1]))

        if not os.path.isfile("{}_{}.pdb".format(pdb, string.ascii_lowercase[n + 1])):
            logger.critical("We could not create {}_{}.pdb".format(pdb, string.ascii_lowercase[n + 1]))
            exit()
        else:
            logger.info("Step of the Trajectory selected in {}_{}.pdb".format(pdb, string.ascii_lowercase[n + 1]))


if __name__ == '__main__':
    init, final, frag, trans, control, pdb, res_fold, criteria, path_pele, report, traject = parse_arguments()
    main(init, final, frag, trans, control, pdb, res_fold, criteria, path_pele, report, traject)
