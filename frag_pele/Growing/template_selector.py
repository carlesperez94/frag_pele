import sys
import re
import os
import string
import pandas as pd
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


def trajectory_selector(output, path_to_file="/growing_output", report="report", trajectory="trajectory.pdb",
                        criteria="Binding Energy"):
    """
    This function select the step of a trajectory of PELE
    with the minimum value of the criteria selected
    and extract it as a single pdb file.
    """
    # Storing the trajectory as a string
    with open(os.path.join(path_to_file, trajectory), 'r') as input_file:
        file_content = input_file.read()

    # Storing the report file as pandas data-frame
    data = pd.read_csv(os.path.join(path_to_file, report), sep='    ', engine='python')

    # Now, select only the columns correspondent to the numberOfAcceptedPeleSteps
    # (if you sum 1 to this number you can obtain the step in the trajectory) and the criteria
    selected_data = data.loc[1:, ['numberOfAcceptedPeleSteps', criteria]]

    # Find the minimum value of the criteria
    min_energy = selected_data.min(axis=0)[1]

    # Once we have this value, find the step that correspond to this minimum value
    min_trajectory = selected_data.loc[selected_data[criteria] == min_energy, 'numberOfAcceptedPeleSteps'].iloc[0]

    # Using this step number we will select the MODEL (step of the trajectory) that correspond to this value.
    # We do min_trajectory+1 because the difference between accepted steps and models numbers
    trajectory_selected = re.search('MODEL\s+%d(.*?)ENDMDL' %int(min_trajectory+1), file_content, re.DOTALL)

    # Writing this step into a PDB with a single structure. We will do it with a list in order to reduce
    # computing time
    output_file_content = []
    output_file_content.append("MODEL     {}".format(int(min_trajectory+1)))
    output_file_content.append("{}".format(trajectory_selected.group(1)))
    output_file_content.append("ENDMDL")
    output_file_content = "\n".join(output_file_content)
    with open(output, 'w') as output_file:
        output_file.write(output_file_content)
        logger.info(":::.MODEL     {} has been selected.:::".format(int(min_trajectory + 1)))


# Old function. We are not using it now.
def change_ligandname(input_file, output):
    """
    From an input pdb file this function replace the first character of the ligand name
    string to the next one in alphabetic order
    """
    # Creating a list of capital letters
    letters = list(string.ascii_uppercase)
    with open(output, 'w') as output_f:
        with open(input_file) as input_f:
            for line in input_f:
                # We only want to read lines that contain information about the ligand
                if line.startswith("HETATM"):
                    # This is the ligandname of the original file, column 4.
                    ligandname_old = line.split()[3]

                    # We will transform this name into a list of letters (strings in python are unmutable)
                    ligandname_new = list(ligandname_old)

                    # Now, the index that we will select is the first letter of the ligandname
                    index = letters.index(ligandname_old[0])

                    # Doing this we are going to transform the first letter of the ligand name into
                    # the next in the alphabet. F.ex: CMB -> DMB
                    new_letter = letters[index+1]
                    ligandname_new[0] = new_letter

                    # Join again all the letters into a string
                    ligandname_new = "".join(ligandname_new)

                    # Write it in the file
                    output_f.write(line.replace(ligandname_old, ligandname_new))
                else:
                    # The rest of the file that don't have information about the ligand will stay unchanged
                    output_f.write(line)
    # We would like to delete the temporary input files
    os.remove(input_file)
