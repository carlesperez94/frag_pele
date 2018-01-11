# General imports
import os
import string
import logging
# Local imports
import Helpers.templatize

logging.basicConfig(filename="output.log",format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)


def control_file_modifier(control_file, pdb, n_files, license="/opt/PELErev12492/licenses", results_path="growing_output"):
    """
    This function creates n control files for each intermediate template created in order to change
    the logPath, reportPath and trajectoryPath to have all control files prepared for PELE simulations.
    """
    # Definition of the keywords that we are going to substitute from the template

    keywords = {"LICENSE": license,
                "RESULTS_PATH": results_path,
                "PDB": pdb
                }

    # As we have several templates that uses as input different PDBs we have to do a loop to create this
    # different templates
    for n in range(0, n_files):

        # Changing the name of the pdb for each control file
        keywords.update({"RESULTS_PATH": "{}_{}".format(results_path, string.ascii_lowercase[n]),
                         "PDB": "{}_{}".format(pdb, string.ascii_lowercase[n])})
        Helpers.templatize.TemplateBuilder(control_file, keywords)

#control_file_modifier("control_template.conf", "181l_prepared_met" , 10, license="/opt/PELErev12492/licenses",
#                          report_path="growing_output")
#    # Variables
#    PATH = os.getcwd()
#
#    # Defining regex pattern's
#    LOG_PATH = re.compile(r'"simulationLogPath"\s+:\s+"(.*)/logFile.txt"')
#    REPORT_PATH = re.compile(r'"trajectoryPath"\s+:\s+"(.*)/trajectory.pdb"')
#    TRAJECTORY_PATH = re.compile(r'"reportPath"\s+:\s+"(.*)/report"')
#    PDB_NAME = re.compile(r'\"Complex\"\s*:\s*{\s*\"files\"\s*:\s*\[\s*{\s*\"path\":\s*\"(.*)\"\s*}\s*\]\s*},')
#
#    if not os.path.exists("control_files_folder"):
#        os.mkdir("control_files_folder")
#
#    for n in range (x_fragments):
#        #we are going to call the control files as control_file_"n", where n is an integer.
#        with open("control_files_folder/control_file_grw_%s" %string.ascii_lowercase[n],'w') as ouput_f:
#            with open(control_file) as control_i:
#                content=control_i.read()
#            #searching the patterns on the file
#            lo= LOG_PATH.search(content)
#            if not lo:
#                logging.critical("We could not find the LogPath in {}!!!".format(control_file))
#                exit("CRITICAL ERROR!!! Check the log file for more information.")
#            report = REPORT_PATH.search(content)
#            if not report:
#                logging.critical("We could not find the reportPath in {}!!!".format(control_file))
#                exit("CRITICAL ERROR!!! Check the log file for more information.")
#            trajectory = TRAJECTORY_PATH.search(content)
#            if not trajectory:
#                logging.critical("We could not find the trajectoryPath in {}!!!".format(control_file))
#                exit("CRITICAL ERROR!!! Check the log file for more information.")
#            pdb_n= PDB_NAME.search(content)
#            if not pdb_n:
#                logging.critical("We could not find the Path of the Complex in {}!!!".format(control_file))
#                exit("CRITICAL ERROR!!! Check the log file for more information.")
#            #modifying the content
#            content = content.replace(lo.group(1),'{}/{}_{}'.format(PATH,results_f_name,string.ascii_lowercase[n]))
#            content = content.replace(trajectory.group(1),'{}/{}_{}'.format(PATH,results_f_name,string.ascii_lowercase[n]))
#            content = content.replace(report.group(1),'{}/{}_{}'.format(PATH,results_f_name, string.ascii_lowercase[n]))
#            #here we are using the pdb sufix for the pdb that we want. For example: if we use tyrosine.pdb to generate n intermediates pdb, we will have as name tyrosine_n.pdb
#            content = content.replace(pdb_n.group(1),'{}/{}_{}.pdb'.format(PATH,pdb, string.ascii_lowercase[n]))
#            #rewriting
#            ouput_f.write(content)

def simulation_runner(control_in):
    #counter=1
    #quantity will depend on the amount of simulations that we will want to do
    #while counter<quantity:
    os.system("/opt/PELErev12492/bin/Pele_serial control_files_folder/{}".format (control_in))
    #print("\n\n----------____SIMULATION {} FINISHED____-------------\n\n".format(counter))
    #counter=1+counter



    




