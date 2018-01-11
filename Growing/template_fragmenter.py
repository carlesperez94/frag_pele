import re
import os
import string
from shutil import copyfile
import logging

logging.basicConfig(filename="output.log",format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG)

def section_selector(template_input,beg_pattern,end_pattern):
    '''Using this function you can select a section of the template file which is between beg_patter and end_pattern'''
    templates_path = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
    #Open the template file and save it as a string
    if not os.path.exists("{}{}".format(templates_path,template_input)):
        logging.critical("Template {} does not exist!!".format(template_input))
    with open("{}{}".format(templates_path,template_input), 'r') as input_file:
        file_content = input_file.read()
        if file_content == "":
            logging.critical("Template file {} is empty!".format(template_input))
    #Select and return everything between two patterns
    section_selected = re.search('{}\n(.*?){}'.format(beg_pattern,end_pattern), file_content, re.DOTALL)
    if section_selected == None:
        logging.critical("Template file contain ERRORS!!!".format(beg_pattern,end_pattern))
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    return section_selected.group(1)

def fragmenter(template_initial,template_final,transformation,n_files):
    templates_path = "DataLocal/Templates/OPLS2005/HeteroAtoms/"
    #Select the section of the templates where is placed the information related with which atoms are there
    atoms_i_section = section_selector(template_initial,"\*", "NBON")
    atoms_f_section = section_selector(template_final, "\*", "NBON")
    #Definition of a patter with regex in order to detect the important information
    atoms_initial = re.findall('\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)', atoms_i_section)
    print("\n\nThe following atoms have been found in the initial template:\n")
    logging.info("The following atoms have been found in the initial template:")
    for atom in atoms_initial:
        print("{}".format(atom[4]))
        logging.info("{}".format(atom[4]))
    atoms_final = re.findall('\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)', atoms_f_section)
    print("\nThe following atoms have been found in the final template:\n")
    logging.info("The following atoms have been found in the final template:")
    for atom in atoms_final:
        print("{}".format(atom[4]))
        logging.info("{}".format(atom[4]))
    #Searching for new atoms in final template
    #First make two lists with all the symbols for each atom
    atoms_initial_symbol=[]
    atoms_final_symbol=[]
    for atom in atoms_initial:
        atoms_initial_symbol.append(atom[4])
    for atom in atoms_final:
        atoms_final_symbol.append(atom[4])
    #Now, find with symbols are different between initial and final template
    new_atoms_symbol=[]
    for atom in atoms_final_symbol:
        if atom not in atoms_initial_symbol:
            new_atoms_symbol.append(atom)
    print("The following new atoms have been detected:")
    logging.info("The following new atoms have been detected:")
    for atom in new_atoms_symbol:
        print(atom)
        logging.info(atom)
    if new_atoms_symbol == []:
        logging.critical("Something went wrong... We could not find new atoms!!! Check templates and try again!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    #Now we have the atoms that have been added in the list new_atoms_symbols. We want to obtain their index:
    print("\n\nChecking indexing for new atoms...\n")
    logging.info("Checking indexing for new atoms...")
    new_atoms_indexes=[]
    for atom_new in new_atoms_symbol:
        for atom in atoms_final:
            if atom[4] == atom_new:
                new_atoms_indexes.append(atom[0])
    #Creation of the intermediate templates and starting to fill them
    for n in range(n_files):
        f=open("{}{}{}z".format(templates_path,string.ascii_lowercase[n], template_final[0:2]), 'w')
        find_ligname = re.search('^({})\s+'.format(template_final[0:3].upper()),atoms_f_section,re.DOTALL)
        #Replacing the name of the ligand and writing the files:
        f.write("* LIGAND DATABASE FILE (OPLS2005)\n*\n{}".format(atoms_f_section.replace("{}".format(find_ligname.group(1)),"{}{}".format(string.ascii_uppercase[n],template_final[0:2].upper()))))
    #Reading the atoms that we are going to transform. F.ex: an H6 into C7.
    with open(transformation, 'r') as transf:
        transformations = []
        index_transformation=[]
        for line in transf:
            transformations.append((line.split()[0], line.split()[1]))
        for transformation in transformations:
            print("\n\nTransforming {} to {}...\n".format(transformation[0], transformation[1]))
            logging.info("Transforming {} to {}...".format(transformation[0], transformation[1]))
            if transformation == []:
                logging.warning("There are not atoms to be transformed!")
            if transformation[1] in new_atoms_symbol and transformation[0] in atoms_initial_symbol:
                index_transformation.append((int(atoms_initial_symbol.index(transformation[0]))+1,new_atoms_indexes[new_atoms_symbol.index(transformation[1])]))
    #Reading NBON parameters
    print("\nReading NBON parameters...\n")
    logging.info("Reading NBON parameters...")
    nbon_section_i=section_selector(template_initial,"NBON","BOND")
    nbon_section_f=section_selector(template_final,"NBON","BOND")
    #Definition of a patter with regex in order to detect the important information
    nbon_initial = re.findall(
        '\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',
        nbon_section_i)
    if nbon_initial == []:
        logging.critical("Something went wrong... We cannot read NBON section of the initial template!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    nbon_final = re.findall(
        '\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)',
        nbon_section_f)
    if nbon_final == []:
        logging.critical("Something went wrong... We cannot read NBON section of the final template!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    print("\nReading VDW...\n")
    logging.info("Reading VDW...")
    #We are going to modify nbon_final, so we want to transform it into a list of lists.
    nbon_final = [list(elem) for elem in nbon_final]
    #Reading Vann Der Waals and Charges
    nbond_to_modify=[]
    for element in nbon_final:
        if element[0] in new_atoms_indexes:
            nbond_to_modify.append((element[0],element[1], element[3]))
    print("\nThe following NBON parameters are going to be modified:\n")
    logging.info("The following NBON parameters are going to be modified:")
    print("Index      VDW   Charge")
    logging.info("Index      VDW   Charge")
    if nbond_to_modify == []:
        logging.critical("Something went wrong... We cannot modify NBON parameters!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    for element in nbond_to_modify:
        print("{:>2}      {}     {: 3.9f}".format(element[0],element[1],float(element[2])))
        logging.info("{:>2}      {}     {: 3.9f}".format(element[0],element[1],float(element[2])))
    #Tranforming atoms
    banned_index=[]
    for something in nbond_to_modify:
        for index in index_transformation:
            if int(something[0]) == int(index[1]):
                banned_index.append(index[1])
                #We are doing substitutions of the lines that contain atoms that need to be transformed
                nbon_final[int(index[1])-1][1] = nbon_initial[int(index[0])-1][1]
                nbon_final[int(index[1]) - 1][3] = nbon_initial[int(index[0]) - 1][3]
    #Writting NBON section:
    os.chdir("{}".format(templates_path, template_final))
    for n in range(1,n_files):
        f = open("{}{}z".format(string.ascii_lowercase[n-1], template_final[0:2]), 'a')
        f.write("NBON\n")
        for element in nbon_final:
            #We are selecting the data of the atoms that are new but they do not belong to transformed group
            if element[0] in new_atoms_indexes and element[0] not in banned_index:
                f.write(" {:5d}   {:3.4f}   {:3.4f}  {: 3.6f}   {:3.4f}   {:3.4f}   {:3.9f}  {: 3.9f}\n".format(
                int(element[0]), float((n*float(element[1]))/int(n_files)), float(element[2]), float(0.0), float(element[4]),
                float(element[5]), float(element[6]), float(element[7])))
            #Now, the transformed ones and the not new atoms.
            else:
                f.write(" {:5d}   {:3.4f}   {:3.4f}  {: 3.6f}   {:3.4f}   {:3.4f}   {:3.9f}  {: 3.9f}\n".format(
                    int(element[0]), float(element[1]), float(element[2]), float(element[3]), float(element[4]),
                    float(element[5]), float(element[6]), float(element[7])))
    os.chdir("../../../../")
    #Reading BOND parameters
    print("\nReading BOND parameters...\n")
    logging.info("Reading BOND parameters...")
    bond_section_i = section_selector(template_initial, "BOND", "THET")
    bond_section_f = section_selector(template_final, "BOND", "THET")
    #Definition of a patter with regex in order to detect the important information
    bond_initial = re.findall('\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)',bond_section_i)
    if bond_initial == []:
        logging.critical("Something went wrong... We cannot read BOND section of the initial template!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    bond_final = re.findall('\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)',bond_section_f)
    if bond_final == []:
        logging.critical("Something went wrong... We cannot read BOND section of the final template!")
        exit("CRITICAL ERROR!!! Check the log file for more information.")
    print("\nModifying BOND parameters...\n")
    logging.info("Modifying BOND parameters...")
    #Writting BOND section:
    os.chdir("{}".format(templates_path, template_final))
    for n in range(1, n_files):
        f = open("{}{}z".format(string.ascii_lowercase[n - 1], template_final[0:2]), 'a')
        f.write("BOND\n")
        for bond in bond_final:
            if bond[1] in new_atoms_indexes:
                for index in index_transformation:
                    if bond[1] != index[1]:
                        f.write(" {:5d} {:5d}   {:5.3f} {: 2.3f}\n".format(
                            int(bond[0]), int(bond[1]), float(bond[2]), float(float(bond[3])*n/n_files)))
                    else:
                        for element in bond_initial:
                            if bond[1] == element[1]:
                                f.write(" {:5d} {:5d}   {:5.3f} {: 2.3f}\n".format(
                                    int(bond[0]), int(bond[1]), float(bond[2]), float(element[3])))
            else:
                f.write(" {:5d} {:5d}   {:5.3f} {: 2.3f}\n".format(
                    int(bond[0]), int(bond[1]), float(bond[2]), float(bond[3])))
    os.chdir("../../../../")
    #Write the final part of the file:
    print("\nFilling the rest of the file...\n")
    logging.info("Filling the rest of the file...")
    end_section_f = section_selector(template_final, "THET", "END")
    os.chdir("{}".format(templates_path, template_final))
    for n in range(1, n_files+1):
        f = open("{}{}z".format(string.ascii_lowercase[n - 1], template_final[0:2]), 'a')
        f.write("THET\n")
        f.write("{}".format(end_section_f))
        f.write("END\n")
        f.close()
    os.chdir("../../../../")
    #Add the final template:
    copyfile("{}{}".format(templates_path,template_final),"{}{}{}z".format(templates_path,string.ascii_lowercase[n_files-1], template_final[0:2]))
    fr = open("{}{}{}z".format(templates_path,string.ascii_lowercase[n_files-1], template_final[0:2]),'r')
    content = fr.read()
    find_ligname = re.search('({})\s+'.format(template_final[0:3].upper()), content)
    fw = open(
        "{}/{}{}z".format(templates_path, string.ascii_lowercase[n_files - 1], template_final[0:2]),
        'w')
    fw.write("{}".format(
            content.replace("{}".format(find_ligname.group(1)),
                                    "{}{}".format(string.ascii_uppercase[n_files-1], template_final[0:2].upper()))))
    fr.close()
    fw.close()
    print("=========== Process finished ===========\n")
    logging.info("=========== Process finished ===========")