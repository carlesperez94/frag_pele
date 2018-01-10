# General imports
import argparse
import string
import os
import logging
# Local imports
import template_selector
import template_fragmenter_2
import simulations_linker

# Logging constants
LOG_FILENAME = "output.out"
LOG_FORMAT = "%(asctime)s:%(levelname)s:%(message)s"

# Logging definition block
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter(LOG_FORMAT)

file_handler = logging.FileHandler(LOG_FILENAME)
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)


def parse_arguments():
    """
        Parse user arguments

        Output: list with all the user arguments
    """
    parser = argparse.ArgumentParser(description="""From an input file, correspondent to the template of the initial structure of the ligand,
    throught this function we will generate "x_fragments" intermediate templates until reach the final template,
    modifing Vann Der Waals, bond lengths and deleting charges""")

    required_named = parser.add_argument_group('required named arguments')  # We have to do this in order
    # to print properly the required arguments when using the defined arguments method

    # Required arguments
    required_named.add_argument("-i", "--initial",
                                required=True,
                                help="""input file correspondent to the 
                                initial template for the ligand that you 
                                want to grow""")
    required_named.add_argument("-t", "--trans", required=True,
                                help="""When an atom is transformed into another
                                one we want to conserve the properties
                                of the first one until being changed in the 
                                last template. This has to be specified in 
                                the program.""")
    required_named.add_argument("-c", "--contrl", required=True,
                                help='Initial control file')
    required_named.add_argument("-p", "--pdb", required=True,
                                help='Initial pdb file')
    required_named.add_argument("-cr", "--criteria",
                                required=True,
                                help="""Name of the column used as criteria in order
                                to select the template used as input for 
                                succesive simulations.""")

    # In the future we will remove this argument
    required_named.add_argument("-f", "--final", required=True,
                                help="""Input file correspondent to the
                                final template for the ligand that you 
                                want to reach""")

    # Optional arguments
    parser.add_argument("-x", "--frag", type=int,
                        default=10,
                        help="""Number of intermediate templates that you want 
                             to generate""")

    parser.add_argument("-r", "--resfold",
                        default="growing_output",
                        help="Name for results folder")
    args = parser.parse_args()

    return args.initial, args.final, args.frag, args.trans, args.contrl, args.pdb, args.resfold, args.criteria


def main(template_initial, template_final,
         n_files, transformation,
         control_file, pdb,
         results_f_name, criteria):
    """
        Description:

        Input:

        Output:

    """

    # Creating template files
    logger.info(("We are going to transform the template _{}_ \
                 into _{}_ in _{}_ steps! Starting...".format(
        template_initial, template_final,
        n_files)))

    templates = template_fragmenter_2.fragmenter(template_initial,
                                                 template_final,
                                                 transformation,
                                                 n_files)
    # Creating control files
    logger.info(
        """=================== Files selected ================\n
        Control file: {} \n
        PDB file: {} \n
        Results folder name: {}\n""".format(
            control_file,
            pdb,
            results_f_name)
                )

    control_files = simulations_linker.control_file_modifier(
        control_file,
        pdb,
        results_f_name,
        n_files
    )

    # Run Pele for each control file

    # for template, control_file in zip(templates, control_files):
    for n in range(0, n_files):
        # Run Pele
        if not os.path.exists("{}_{}".format(results_f_name,
                                             string.ascii_lowercase[n]
                                             )):
            os.mkdir("{}_{}".format(results_f_name,
                                    string.ascii_lowercase[n])
                     )
            simulations_linker.simulation_runner(
                "control_file_grw_{}".format(string.ascii_lowercase[n])
            )
        logger.info("SIMULATION FOR control_file_grw_{} COMPLETED!!!!!!!".format(
            string.ascii_lowercase[n])
                    )
        # Choose the best trajectory
        template_selector.trajectory_selector("{}_{}".format(results_f_name,
                                                             string.ascii_lowercase[n]
                                                             ),
                                              "{}_{}_tmp.pdb".format(pdb,
                                                                     string.ascii_lowercase[n + 1]),
                                              "{}".format(criteria)
                                              )
        if os.path.exists("{}_{}".format(pdb,
                                         string.ascii_lowercase[n + 1])
                          ):
            logger.info("Step of the Trajectory selected in {}_{}.pdb".format(
                pdb,
                string.ascii_lowercase[n + 1])
                        )
        else:
            logger.critical("We could not create {}_{}.pdb".format(pdb,
                                                                   string.ascii_lowercase[n + 1])
                            )
            exit()


if __name__ == '__main__':
    init, final, frag, trans, control, pdb, res_fold, criteria = parse_arguments()
    main(init, final, frag, trans, control, pdb, res_fold, criteria)
