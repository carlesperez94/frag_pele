class CovalentPDB:

    content = None
    ligand = None
    conects = []

    def __init__(self, pdb_file, lig_chain="L", lig_name="LIG", res_num=None, 
                 res_chain=None, lig_link=None, res_link=None):
        self.pdb_file = pdb_file
        self.lig_chain = lig_chain
        self.lig_name = lig_name
        self.res_num = res_num
        self.res_chain = res_chain
        self.lig_link = lig_link
        self.res_link = res_link
        self._read_pdb_content()
        if not self._check_if_exists_conects() and (res_num==None or res_chain==None):
            raise ValueError("PDB file does not contain CONECTs. \nPlease, fill the number and chain of the residue covalently linked to the ligand")
        self._read_ligand()
        self._read_conects()
        self._get_residue_info_linked()
    
    def _read_pdb_content(self):
        with open(self.pdb_file, "r") as infile:
            self.content = infile.read()    
 
    def _check_if_exists_conects(self):
        pdb_lines = self.content.split("\n")
        for l in pdb_lines:
            if l.startswith("CONECT"):
                return True
        return False

    def _read_ligand(self):
        lig_lines = []
        for l in self.content.split("\n"):
            if self.lig_name == l[17:20] and l.startswith("HETATM"):
                lig_lines.append(l)
        self.ligand = "\n".join(lig_lines)

    def _read_conects(self):
        for l in self.content.split("\n"):
            if l.startswith("CONECT"):
                group = l.split("CONECT")[1].split()
                group = list(map(int, group))
                self.conects.append(group)
    
    def _get_ligand_residue_connection(self):
        conections = []
        all_indexes = get_indexes(self.content)
        ligand_indexes = get_indexes(self.ligand)
        residue_indexes = list(set(all_indexes) - set(ligand_indexes))
        if self.conects:
            for c_group in self.conects:
                # Check if any conection group contains atoms of the ligand and residues at the same time
                if set(c_group).intersection(set(ligand_indexes)) and set(c_group).intersection(set(residue_indexes)):
                    conections.append(c_group)
        return conections
   
    def _get_residue_info_linked(self):
        conections = self._get_ligand_residue_connection()[0]
        for line in self.content.split("\n"):
            if line.startswith("HETATM") and line[6:11].strip() == str(conections[0]):
                self.lig_link = line[12:15].strip()
            if line.startswith("ATOM") and line[6:11].strip() == str(conections[-1]):
                self.res_num = line[23:26].strip()
                self.res_chain = line[21]
                self.res_link = line[12:15].strip()

def get_indexes(pdb_fragment):
    indexes = []
    for l in pdb_fragment.split("\n"):
        if l.startswith("HETATM") or l.startswith("ATOM"):
            indexes.append(int(l[7:11].strip()))
    return indexes

#def main(pdb_file, lig_chain="L", lig_name="LIG", res_num=None, res_chain=None):
pdb = CovalentPDB("/home/bsc72/bsc72292/projects/covid/COVALENT_Frag/covdock_sar26lu7_15.pdb")
print(pdb.res_num, pdb.res_chain, pdb.res_link, pdb.lig_link)
