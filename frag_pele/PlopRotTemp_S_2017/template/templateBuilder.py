import subprocess
import __future__
import os
import shutil
import re
import warnings
from schrodinger import structure
from frag_pele.PlopRotTemp_S_2017.template.chargeHandler import ChargeHandler
from itertools import groupby

try:
    from PlopRotTemp.PlopRotTemp import preproces_file_lines, find_bonds_in_mae, find_resnames_in_mae, find_names_in_mae, xyz2int
except ImportError:
    from PlopRotTemp import preproces_file_lines, find_resnames_in_mae, find_names_in_mae, xyz2int, find_bonds_in_mae

CONVERSION = {'CT1': 'CT',
              'O2Z':'O2',
              'CT1G':'CT',
              'NS': 'N',
              'NG': 'NY',
              'CG1': 'CA8',
              'NE':'NT2',
              'CB':'CA2',
              'NNC': 'NY2',
              'NND': 'NY',
              'NNB': 'NY',
              'NI': 'NY',
              'NIP': 'NP'
              }

DEFAULT_INT = 6
DEFAULT_SPRING_K = '268.0'
DEFAULT_EQ_DIST = '1.529'
DEFAULT_ATOMTYPE = [0.000, 0.000, 1.500, 1.250, 0.005000000, 0.000000000] 

FILE_DIR_PATH = os.path.dirname(os.path.dirname(__file__))
SIMILARITY_PATH = os.path.join(FILE_DIR_PATH, 'param/similarity.param') ##impact??
PARAM_PATH = os.path.join(FILE_DIR_PATH, 'param/f14_sgbnp.param') ##impact??
OPLS_CONVERSION_FILE = 'param.dat'
OPLS_VERSION = '14'


class TemplateBuilder:

    """
    :Description: base builder class for ligand template creation

    :Attributes:
        - input_file: input .mae file
        - output_file : output template file

    :Author: Daniel Soler
    """

    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        #self.build_template()
        

    def build_template(self, charges_from_file=None):
      """
        :Description: Build ligand template from mae file

        :Input:
          - charges_from_file: charge to import the charges from

        :Output:
          - output: ligand template
          - resname: residue name

        :Explanation:

        Preparation of param:
          1- Use the flld_server schrodinger command line server
             to create the 2005OPLS paramaters for the ligand
          2- Parse the previous parameters extracting 
             atomtypes, bond info, non bonded info, torsions,
             dihedrals and imporpers rotations.
          3- From this parameters create the internal coordinates
             and get some extra parameters as the nonbonded
             sgb_radius, gammas, alphas from our proper database
             which is in PARAM_PATH in case it needs to be extended
             by the user.
          4- With the previous parsed parameters create the several
             sections of the template
          5- Write everything to a RESIDUENAME.hetgrp file

        :Author: Daniel Soler
      """

      #Create ligand params with ffld_sever from schrodinger
      ffld_server_command = os.path.join(os.environ['SCHRODINGER'], 'utilities/ffld_server')
      subprocess.call([ffld_server_command, "-imae", self.input_file, "-version",
        OPLS_VERSION, "-print_parameters", "-out_file", OPLS_CONVERSION_FILE])

      #Retrieve all the useful information from that params
      atom_names_param = self.retrieve_atom_names(OPLS_CONVERSION_FILE)
      #shutil.copyfile(OPLS_CONVERSION_FILE, "a.txt")
      self.search_and_replace(OPLS_CONVERSION_FILE, atom_names_param)
      atom_types, parents, charges, sigmas, epsilons, stretchings, tors, phis, impropers = self.parse_param(OPLS_CONVERSION_FILE, atom_names_param)


      #Connectivity information from Mae
      res_name = find_resnames_in_mae(self.input_file)[0] #Ligand must be defined as a whole residue
      atom_names = find_names_in_mae(self.input_file, undersc=True)
      bonds = find_bonds_in_mae(self.input_file)

      #Fix parents
      parents = self.fix_parents_rings(parents, atom_names)

      
      zmat = self.create_zmatrix(parents)
      number_bonds = len(stretchings)
      number_torsions = len(tors)
      number_atoms = len(atom_names)
      number_phis = len(phis)
      number_improper = len(impropers)

      #Charges from file
      if charges_from_file: charges = ChargeHandler(self.input_file).get_charges()
      #fix C=O amides
      new_atom_types = self.fix_atomtype('O', 'N', 'OCN1', 2, atom_types)
      #fix H
      new_atom_types = self.fix_atomtype('H', 'N', 'HN', 1, new_atom_types)
      #fix aromatics
      new_atom_types = self.fix_aromatics(new_atom_types)

      #Missing NBND information in param
      sgb_radius, vdw_radius_param, gammas, alphas = self.SGB_paramaters(new_atom_types)

      #Choosign param files for hidrogens vdwr or ffldserver params for the rest.
      vdw_radius = []
      for vdw_params, sigma, atom_type in zip(vdw_radius_param, sigmas, atom_types):
        vdw_radius.append(vdw_params)


      ################################Template Creation########################
      
      header = ["* LIGAND DATABASE FILE (OPLS2005)",
                "*",
                '{0:>0}  {1:>6} {2:>5} {3:>5} {4:>7} {5:>7}'.format(res_name, number_atoms, number_bonds,
                                                          number_torsions, (number_phis + number_improper), 0)
                ]

      connectivity_section = []
      for i, (atom_name, atom_type, parent, zmat_row) in enumerate(
        zip(atom_names, atom_types, parents, zmat), 1):
        connectivity_line = '{0:>5} {1:>5} {2:>0} {3:>5} {4:>5} {5:>5} {6:>11.5f} {7:>11.5f} {8:>11.5f} '.format(
          i, int(parent)+1, 'S', atom_type, atom_name,
          DEFAULT_INT, zmat_row[0], zmat_row[1], zmat_row[2])
        connectivity_section.append(connectivity_line)

      NBOND_section = []
      for i, (atom_type, charge, sigma, epsilon, radius, vdw, alpha, gamma) in enumerate(
        zip(atom_types, charges, sigmas, epsilons, sgb_radius, vdw_radius, alphas, gammas), 1):
            NBOND_section.append('{0:>5} {1:>8.4f} {2:>8.4f} {3:>10.6f} {4:>8.4f} {5:>8.4f} {6:>13.9f} {7:>13.9f}'.format(
              i, float(sigma), float(epsilon), float(charge), float(radius), float(vdw), float(gamma), float(alpha)))
     
      strech_section = []
      for stretching in stretchings:
            strech_section.append('{0:>5} {1:>5} {2:>9.3f}  {3:>5.3f}'.format(
              stretching[0], stretching[1], float(stretching[2]), float(stretching[3])))
       
      tors_section = []
      for tor in tors:
        tors_section.append('{0:>5} {1:>5} {2:>5} {3:>11.5f} {4:>10.5f}'.format(
          tor[0], tor[1], tor[2], float(tor[3]), float(tor[4])))

      phi_section = []
      phis = self.descompose_dihedrals(phis)
      for phi in phis:
        phi_section.append('{0:>5} {1:>5} {2:>5} {3:>5} {4:>9.5f} {5:>4.1f} {6:>3.1f}'.format(
          phi[0], phi[1], phi[2], phi[3], (float(phi[4])/2.0), phi[5], abs(float(phi[6]))))

      iphi_section = []
      for improper in impropers:
        iphi_section.append('{0:>5} {1:>5} {2:>5} {3:>5} {4:>9.5f} {5:>3.1f} {6:>3.1f}'.format(
          improper[0], improper[1], improper[2], improper[3], float(improper[4])/2.0, -1, 2))

      file_content = []
      file_content.extend( header +
                           connectivity_section +
                           ['NBON'] +
                           NBOND_section +
                           ['BOND'] +
                           strech_section +
                           ['THET'] +
                           tors_section +
                           ['PHI'] +
                           phi_section +
                           ['IPHI'] +
                           iphi_section +
                           ['END'])
      #Write to file
      with open(res_name, 'w') as f:
        f.write('\n'.join(file_content))
      
      # #Remove param.dat file
      try:
        os.remove(OPLS_CONVERSION_FILE)
      except OSError:
        print("Error, param.dat not created. Make sure $SCHRODINGER/utilities/ffld_server is up and running in your computer.")

      #stdout
      print("Template {} generated successfully".format(self.output_file))

      return res_name, self.output_file, self.input_file, self.output_file, res_name


    def build_triangular_matrix(self, stretchings, tors, phis, atom_names):
        """
            Build triangular interaction matrix.
            You can find more information in PELE's docs
        """
        bonds = stretchings

        counts = []
        connections = []
        for i, atom_name in enumerate(atom_names):
            count = 0
            connected = []
            for stretching in bonds[:]:
                if i in stretching:
                    count += 1
                    if i == stretching[0]:
                        connected.append(str(stretching[1]+1))
                    else:
                        connected.append(str(stretching[0]+1))
                    bonds.remove(stretching)
            for tor in tors[:]:
                if str(i+1) in tor[0:3]:
                    count += 1
                    tors.remove(tor)
            for phi in phis[:]:
                if str(i+1) in phi[0:4]:
                    count += 1
                    phis.remove(phi)
            counts.append(count)
            connections.append(connected)
        final_counts = " ".join([str(count) for count in counts if count != 0 ])
        connections_repited = [ " ".join(connection) if connection else "0" for connection in connections ]
        final_connections = [x[0] for x in groupby(connections_repited)]
        triangular_matrix = final_counts+"\n"+"\n".join(final_connections)







    def fix_parents_rings(self, parents, atom_names):
        """
            For every ring in the structure assign as parent 
            of each atom the previous. To close the ring assign
            as parent of the initial atom the last.
        """
        str1 = next(structure.StructureReader(self.input_file))
        rings = str1.ring
        for ring in rings:
            ring_atoms = ring.getAtomList()
            initial_atom = ring_atoms[0]-1
            last_atom= ring_atoms[-1]-1

            start = True
            current_atom = initial_atom
            next_atom = ring_atoms[1] -1
            counter=1
            while (current_atom != initial_atom  or start is True):
                try:
                    parents[next_atom] = current_atom
                    counter+=1
                    current_atom=next_atom
                    next_atom=ring_atoms[counter]-1
                    start=False
                except IndexError:
                    parents[initial_atom] = last_atom
                    break
        return parents


    def create_zmatrix(self, parents):
      """
        :Description: Retrieve the internal coordinates (zmatrix)
        from the complex cooridnates & topology:

        :Input:
          - str1: structure.structure class object
          - parents: atom's parent list

        :Ouput:
          - zmatrix: structure internal coord
      """
      str1 = next(structure.StructureReader(self.input_file))
      order = [i for i in range(len(parents))]
      coordinates = [atom.xyz for atom in str1.atom]
      zmat = xyz2int(coordinates, order, parents)
      return zmat


    def fix_atomtype(self, target_atom, neighbour, new_atom_type, bonds_dist, atom_types):
        """
            :Description: For each target_atom check whether or not
            there is an atom bonded to the atom neighbour
            bonded at a certain bonds_dist, and then,
            if there is, set the target_atom to new_atom_type

            :Input:
                - target_atom :Atom targeted to reassigned its atom type
                - neighbour: Atom bonded to the target_atom at a certain bonds_dist
                - new_atom_type: Atom type to assign to target_atom if the neighbour atom
                                 is connected to him at a certain bonds_dist
                -bonds_dist: number of bonds from the target_atom to look for the nighbour_atom
                - atom_types: Atom types of the system

            :Output: 
                - new_atom_types: New set of atom type where we may have changed the 
                target_atom atom type if the condition we looked for where accomplished.

            :Example: 

                - new_atom_types = fix_atomtype('O', 'N', 'OCN1', 2, atom_types)
                   
               Look for an carbonyl oxigen atom (O) bonded to a nitrogen (N) at 
               a bond distance of (2). If its found the oxigen atom will be reassign
               with the new atom type (OCN1). And finally the new atom types will be
               returned.

            :Author: Daniel Soler

        """
        new_atom_types = atom_types[:]
        st1 = next(structure.StructureReader(self.input_file))

        atoms_to_study = []
        indexes=[]
        for i, atom in enumerate(st1.atom):
            if atom_types[i] == target_atom:
                atoms_to_study.append(atom)
                indexes.append(i)

        for i, atom in enumerate(atoms_to_study):
            current_atom = atom
            for j in range(bonds_dist):
                atoms_connected = current_atom.bonded_atoms
                for atom_connected in atoms_connected:
                    if(atom_connected == atom):
                        continue
                    else:
                        if(atom_connected.element==neighbour):
                                atom_type_index = indexes[i]
                                new_atom_types[atom_type_index] = new_atom_type
                        current_atom = atom_connected
        return new_atom_types

    def fix_aromatics(self, atom_types):
        """
          :Description: Look for the atom types ['CA', 'C5A', 'C5B', 'CN', 'CB']
          and convert them to 'CA2' if they have 3 bonds, in other words,
          if they are sp2 carbons.

          :Author: Daniel Soler
        """
        new_atom_types = atom_types[:]
        st1 = next(structure.StructureReader(self.input_file))
        target_atoms = ['CA', 'C5A','CA5', 'C5B', 'C5BC', 'CN', 'CB', 'C56A', 'C56B', 'CT4', 'CRA', 'CN56', 'C5X', 'C5BB', 'CR3']
        new_atom_type = 'CA2'

        atoms_to_study = []
        indexes=[]
        for i, atom in enumerate(st1.atom):
            if atom_types[i] in target_atoms:
                atoms_to_study.append(atom)
                indexes.append(i)

        for i, atom in enumerate(atoms_to_study):
            atoms_connected = list(atom.bonded_atoms)
            if(len(atoms_connected)==3):
              atom_type_index = indexes[i]
              new_atom_types[atom_type_index] = new_atom_type
                      
        return new_atom_types


    @classmethod
    def SGB_paramaters(cls, atom_types, tried = []):
      """
        :Description: Parse the param/f14_sgbnp.param file
        to obtain all the SGB Nonbonded parameters
        to build the template file.

        If some atom_type is not found the param/similarity.param
        is parsed to look for the next most similar atom type.
        If that is not found as well, we append the tried
        atom_types to tried [] for efficiency and keep looking for
        other similar atom types.

        :Arguments: 
            - atom_types: atom_types of the system
            - tried: atom_types that we try to obtain its parameters
                     but they are not contained in the f14_sgbnp.param file.

        :Returns:
            - radius: Atoms SGB radius
            - vdw_r: Atoms Van der Waals
            - gammas: Atoms gamma solvatation factor
            - alphas: Atoms alpha solvatation factor

            For more information on this parameters refer
            to the PELE documentation.

        Example with one atom_type (ONN):

            1- We look whether or not its in the conversion atom dictionary
                1.1- If it is we change the atom type to the one in the dict
            2- We search the SGB parameters in the PARAM_PATH file
                2.1- If found return them separately
                2.2 -If not found call find_similar_atomtype_params to look
                     for similar atom types parameters
            3- We return the similar SGB parameters or a default
               if they were not found

        :Author: Daniel Soler
      """

      radius = []
      vdw_r = []
      gammas = []
      alphas = [] 
      lines = preproces_file_lines(PARAM_PATH)
      for atom_type in atom_types:
        if atom_type in CONVERSION:
            atom_type = CONVERSION[atom_type]

        found = False
        for line in lines:
          if not line.startswith('#'):
              line = line.split()
              atom_type_file = line[1]
              if(atom_type == atom_type_file):
                radius.append(line[4])
                vdw_r.append(line[5])
                gammas.append(line[6])
                alphas.append(line[7])
                found = True
                break

        if not found:
          new_params = cls.find_similar_atomtype_params(atom_type, tried=[])
          if(new_params):
            radius.append(new_params[4])
            vdw_r.append(new_params[5])
            gammas.append(new_params[6])
            alphas.append(new_params[7])
            # warnings.warn("Paramaters of {} used for {}.".format(
            # new_params[1], atom_type))
          else:
            # warnings.warn("Defaults Paramaters used for{}.".format(atom_type))
            radius.append(DEFAULT_ATOMTYPE[2])
            vdw_r.append(DEFAULT_ATOMTYPE[3])
            gammas.append(DEFAULT_ATOMTYPE[4])
            alphas.append(DEFAULT_ATOMTYPE[5])

      return(radius, vdw_r, gammas, alphas)


    @classmethod
    def find_similar_atomtype_params(cls, atom_type, tried):
      """
        :Description: If some atom_type is not found in the SGB_paramaters() func
        the param/similarity.param is parsed to look for the next most
        similar atom type.
        If that is not found as well, we append the tried
        atom_types to tried [] for efficiency and keep looking for
        other similar atom types.

        :Input: 
            -atom_type: atom_type of the atom not found.
            - tried: similar atom_types that we try to obtain its parameters
                     but they are not contained in the f14_sgbnp.param file.
        :Output:
            - new_params: [
                          similar_atom_type,
                          SGB_radius,
                          Van der Waals radius,
                          Solvatation gamma factor,
                          Solvatation alpha factor
                          ]

      """
      new_atom_type=False
      similarity = 0

      #Search for the most similar atom
      lines = preproces_file_lines(SIMILARITY_PATH)
      for line in lines:
        line = line.split()
        if(line[0]==atom_type and float(line[2])>float(similarity) and line[1] not in tried):
          similarity = line[2]
          new_atom_type = line[1]
        elif(line[1]==atom_type and float(line[2])>float(similarity)  and line[0] not in tried):
          similarity = line[2]
          new_atom_type = line[0]
      if new_atom_type:
        tried.append(new_atom_type)

        lines = preproces_file_lines(PARAM_PATH)
        for line in lines:
          if not line.startswith('#'):
              line = line.split()
              atom_type_file = line[1]
              if(atom_type == atom_type_file):
                return(line)
        new_params = cls.find_similar_atomtype_params(new_atom_type, tried)
        return new_params

      else:
        return []

    @staticmethod  
    def retrieve_atom_names(OPLS_CONVERSION_FILE):
      """
        :Description: parse the param.data and retrieve
        the atom names of the molecule

        :Input: OPLS_CONVERSION_FILE: ffldserver output file
                   cointaining the system parameters

        :Output: atom_names: Atom_names of the systems

        :Author: Daniel Soler
      """
      atom_names = []
      lines = preproces_file_lines(OPLS_CONVERSION_FILE)
      start_found_NBND = False
      for line in lines:
          if not line.startswith("----------"):
            if(line.startswith("atom type vdw")):
              start_found_NBND = True
            elif(start_found_NBND):
              try:
                line = line.split()
                atom_names.append(line[0])
              except IndexError:
                return atom_names

    @staticmethod
    def search_and_replace(file, to_search):
      """
        Search and replace atom_names for numbers
      """
      to_replace = range(1, len(to_search)+1)

      with open(file, "r+") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
          lines[i] = ' ' + line.strip('\n')


      with open(file, "w") as f:
        f.write('\n'.join(lines))
      
      with open(file, "r") as f:
        filedata = f.read()
     

      for (item_to_search, item_to_replace) in zip(to_search, to_replace):
          filedata = filedata.replace(' ' + item_to_search.strip('_')+' ', ' '+str(item_to_replace)+' ') #atom types are like _O1_ (strip)



      f = open(file,'w')
      f.write(filedata)
      f.close()

    def parse_param(self, param_file, atom_names):
      """
        :Description: Parse the OPLS conversion param file
        and get the atomtypes.

        Go through the param.dat through several stages
        while jumping al the ----- lines:


        1. Look for the BSCI section for each atom
           and its parent appending both like
           bond [A,B] and parent [B].
        2. Set end_connectivity_found=True to start
           parsing the NBND section
        
        In the NBND section:
        2. Look for the atom type vd2w line
        3. Parse each line of the non bonding section
        4. set NBND_finished= True to start parsing the bonded section
        

        In the BOND section
        5. For each subsection ([stretchings, bendings, proper_tors, improper_tors]):
           1. Search for the respective title indicating the beggining of the section
           (["Stretch", "Bending", "proper Torsion", "improper Torsion"])
           2. Then parse the section obtaining the values from the columns indicated in columns to take
           3.when all sections are parsed return all values

        :Author: Daniel Soler

      """

      atom_types = []
      parents = [-1] * len(atom_names)
      charges = []
      sigmas = []
      epsilons = []
      vdws = []
      stretchings = []
      bendings = []
      proper_tors = []
      improper_tors = []
      lists = [stretchings, bendings, proper_tors, improper_tors]
      keywords = ["Stretch", "Bending", "proper Torsion", "improper Torsion"]
      columns_to_take = [[0,1, 2, 3], [0, 1, 2, 3,4], [0, 1 ,2 , 3, 4, 5, 6, 7], [0, 1, 2, 3, 4]]
      start_connectivity_found = False
      end_connectivity_found = False
      start_found_NBND = False
      NBND_finished = False
      start_found_BND = False


      with open(param_file) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
          #prepared lines
          line = re.sub(' +',' ',line)
          line = line.strip('\n').strip()

          if not end_connectivity_found:
            if(line.startswith("BCI's")):
              start_connectivity_found = True
            elif(start_connectivity_found):
              try:
                line = line.split()
                #parents[int(line[1])-1]-->Atom names start at 1 but list at 0
                if(parents[int(line[1])-1]==-1):
                    parents[int(line[1])-1] = int(line[0])-1
                elif(parents[int(line[0])-1]==-1):
                    parents[int(line[0])-1] = int(line[1])-1
              except IndexError:
                end_connectivity_found = True


          #NBND section
          elif (line.startswith("----------") is False) and (NBND_finished is False) and (end_connectivity_found is True):
            if(line.startswith("atom type vdw")):
              start_found_NBND = True
            elif(start_found_NBND):
              try:
                line = line.split()
                if not line[3].isdigit():
                  atom_types.append(line[3])
                else:
                  atom_types.append(atom_names[int(line[3])-1])
                charges.append(line[4])
                sigmas.append(line[5])
                epsilons.append(line[6])
              except IndexError:
                NBND_finished=True

          elif(NBND_finished):
            for (keyword, indexes, List) in zip(keywords, columns_to_take, lists):
              while(line.startswith(keyword) is False):
                try:
                  line, i = move_line_forward(lines, i)
                except IndexError:
                  #There are no improper torsions!!
                  return atom_types, parents, charges, sigmas, epsilons, stretchings, bendings, proper_tors, improper_tors
              line, i = move_line_forward(lines, i)
              while(line):
                    line = re.sub(' +',' ',line)
                    line = line.strip('\n')
                    line = line.split()
                    values = [line[index] for index in indexes]
                    values = self.amide_trans_cis_hotfix(line, values)
                    List.append(values)
                    line, i = move_line_forward(lines, i)
            return atom_types, parents, charges, sigmas, epsilons, stretchings, bendings, proper_tors, improper_tors

    @staticmethod
    def amide_trans_cis_hotfix(line, values):
        """
            OPLS 2005 have the same probability for cis and trans amide
            but as you can see on http://www.cryst.bbk.ac.uk/PPS2/course/section7/os_torg.html
            there is a dirst term that favours the trans conformation. This patch
            has that into account by adding the first v1 dihedral term
        """
        amide = False
        try:
            if line[12] == "H-N-C-O":
                amide = True
        except IndexError:
            pass
        finally:
            if amide:
                # first dihedral component
                values[4] = "1.300"
        return values

    @staticmethod
    def descompose_dihedrals(phis):
        """
        :Description: For each dihedral line as:
        atom1 atom2 atom3 atom4   V1    V2   V3    V4
          1     2     3      4   0.00 1.000 5.000 3.000

        Separate all components as next:
        atom1 atom2 atom3 atom4   value  which VX???
         1     2     3      4     1.000      2
         1     2     3      4     5.000      3
        """
        new_phis = []
        new_phi = []
        for phi in phis:
            atoms = phi[0:4]
            phis_components = ["{0:.3f}".format(abs(float(component))) for component in phi[4:8]]
            if(phis_components == ['0.000','0.000','0.000','0.000']):
              new_phis.append([atoms[0], atoms[1], atoms[2], atoms[3], phi[4], 1, 1])
              continue
            else:
                for index, component in enumerate(phi[4:8], 1):
                    new_phis.append([atoms[0], atoms[1], atoms[2], atoms[3], component, 1, index])
        return new_phis



def move_line_forward(lines, i):
  i+=1
  line = lines[i].strip()
  return line, i
