#import template_selector
#import template_fragmenter_2
#import simulations_linker
import argparse
import string
import os
import logging

# if os.path.exists("output.log"):
#    os.remove("output.log")
logging.basicConfig(filename="output.log", format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)


# The main block of code
def main():
    templates_path = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
    # Definition of arguments
    parser = argparse.ArgumentParser(description='''From an input file, correspondent to the template of the initial structure of the ligand,
    throught this function we will generate "x_fragments" intermediate templates until reach the final template,
    modifing Vann Der Waals, bond lengths and deleting charges''')

    parser.add_argument("-i", "--initial", dest="template_initial", required=True,
                        help='input file correspondent to the initial template for the ligand that you want to grow')
    parser.add_argument("-f", "--final", dest="template_final", required=True,
                        help='input file correspondent to the final template for the ligand that you want to reach')
    parser.add_argument("-x", "--frag", type=int, dest="n_files", required=True,
                        help='number of intermediate templates that you want to generate')
    parser.add_argument("-t", "--trans", dest="transformation", required=True,
                        help='When an atom is transformed in another one we want to conserve the properties of the '
                             'first one until being changed in the last template. This has to be specified in the program.')
    parser.add_argument("-c", "--contrl", dest="control_file", required=True,
                        help='Initial control file')
    parser.add_argument("-p", "--pdb", dest="pdb", required=True,
                        help='Initial pdb file')
    parser.add_argument("-r", "--resfold", dest="results_f_name", required=True,
                        help='Name for results folder')
    parser.add_argument("-cr", "--criteria", dest="criteria", required=True,
                        help='Name of the column used as criteria in order to select the template used as input for succesive simulations')
    args = parser.parse_args()
    # Main algorithm
    print("We are going to transform the template _{}_ into _{}_ in _{}_ steps!\nStarting...\n".format(args.template_initial, args.template_final, args.n_files))
    logging.info(("We are going to transform the template _{}_ into _{}_ in _{}_ steps! Starting...".format(args.template_initial, args.template_final, args.n_files)))
    template_fragmenter_2.fragmenter("{}".format(args.template_initial), "{}".format(args.template_final), "{}".format(args.transformation), args.n_files)
    print("========== Files selected ========\nControl file: {}\nPDB file: {}\nResults folder name: {}\n=====================".format(args.control_file, args.pdb, args.results_f_name))
    logging.info("===== Files selected ====== Control file: {} PDB file: {} Results folder name: {}".format(args.control_file, args.pdb, args.results_f_name))
    simulations_linker.control_file_modifier(args.control_file, args.pdb, args.results_f_name, args.n_files)
    for n in range(0, args.n_files):
        if os.path.exists("{}_{}".format(args.results_f_name, string.ascii_lowercase[n])) == False:
            os.mkdir("{}_{}".format(args.results_f_name, string.ascii_lowercase[n]))
        simulations_linker.simulation_runner("control_file_grw_{}".format(string.ascii_lowercase[n]))
        print("SIMULATION FOR control_file_grw_{} COMPLETED!!!!!!!".format(string.ascii_lowercase[n]))
        logging.info("SIMULATION FOR control_file_grw_{} COMPLETED!!!!!!!".format(string.ascii_lowercase[n]))
        template_selector.trajectory_selector("{}_{}".format(args.results_f_name, string.ascii_lowercase[n]), "{}_{}_tmp.pdb".format(args.pdb, string.ascii_lowercase[n + 1]), "{}".format(args.criteria))
        template_selector.change_ligandname("{}_{}_tmp.pdb".format(args.pdb, string.ascii_lowercase[n + 1]), "{}_{}.pdb".format(args.pdb, string.ascii_lowercase[n + 1]))
        if os.path.exists("{}_{}".format(args.pdb, string.ascii_lowercase[n + 1])):
            print("Step of the Trajectory selected in {}_{}.pdb".format(args.pdb, string.ascii_lowercase[n + 1]))
            logging.info("Step of the Trajectory selected in {}_{}.pdb".format(args.pdb, string.ascii_lowercase[n + 1]))
        else:
            logging.critical("We could not create {}_{}.pdb".format(args.pdb, string.ascii_lowercase[n + 1]))
            exit("CRITICAL ERROR!!! Check the log file for more information.")


if __name__ == '__main__':
    main()
