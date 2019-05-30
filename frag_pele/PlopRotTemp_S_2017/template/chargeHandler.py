# -*- coding: utf-8 -*-

import sys
import re
from frag_pele.PlopRotTemp_S_2017.template.tmp_helpers import Helper



class ChargeHandler(Helper):

    """
        Handler class in charge of retrieving
        atom charges from file to TemplateBuilder
    """

    def __init__(self, file):
        self.file = file
        

    def get_charges(self):
        """
            Parse self.file for charges and retrieve.
        """
        charges = []

        with open(self.file, "r") as f:
            while f:  # Find Bond Section
                line = f.readline()
                if line.startswith(" m_atom"):
                    break
            keywords = []
            while f:  # Read keywords                
                line = f.readline()
                if (re.search(':::', line)):
                    break
                else:
                    keyword = "index" if "index" in line else line.strip()
                    keywords.append(keyword)
            while f:  # get charges
                line = f.readline()
                if not re.search(':::', line):
                    line = re.sub(' +',' ',line).strip('\n').strip().split()
                    for i, keyword in enumerate(keywords):
                        if keyword == "r_m_charge1":          
                            charge = line[i]
                            charges.append(charge)
                else:
                    if charges: return charges
                    else: raise ValueError("NO CHARGES IN MAE")
            

                    




    
