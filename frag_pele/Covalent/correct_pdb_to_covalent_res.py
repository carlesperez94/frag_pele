import lib_prep.pdb_modifier as pdm

LIST_OF_IONS = ["ZN", "MN", "FE", "CO", "NI", "CA", "CD"]

class Corrector(pdm.PDB):
    def __init__(self, pdb_file, original_reschain, original_resnum, new_resname, ligand_chain=None):
        pdm.PDB.__init__(self, in_pdb=pdb_file)
        self.original_reschain = original_reschain
        self.original_resnum = original_resnum
        self.new_resname = new_resname
        self.ligand_chain = ligand_chain
        self.old_res_lines = self.get_residue(self.original_reschain,
                                              self.original_resnum)
        self.new_res_lines = self.get_atoms_of_resname(self.new_resname)
        self.remove_lines(self.new_res_lines)
        self.correct_atom_section_new_res_lines()
        self.correct_chain_and_resnum_new_res_lines()
        self.add_lines_at_position(self.new_res_lines, 
                                   self.find_index_of_line(self.old_res_lines[0]))
        self.remove_lines(self.old_res_lines)
        self.content = "".join(self.lines)


    def correct_atom_section_new_res_lines(self):
        correction = []
        for line in self.new_res_lines:
            if line.startswith("HETATM"):
                line = list(line)
                line[0:6] = "ATOM  "
                line = "".join(line)
            correction.append(line)
        self.new_res_lines = correction
    
    def check_and_rename_old_resname(self):
        correction = []
        for line in self.lines:
            resname = pdm.get_resname_from_line(line)
            reschain = pdm.get_chain_from_line(line)
            resnum = pdm.get_resnum_from_line(line)
            if str(reschain).strip() == str(self.original_reschain).strip() and \
               str(resnum).strip() == str(self.original_resnum).strip() and \
               str(resname).strip() == str(self.new_resname).strip():
                line = pdm.set_resname_to_line(line, "RES")
            correction.append(line)
        self.lines = correction
 
    def clear_heteroatoms(self):
        to_keep = []
        for line in self.lines:
            if line.startswith("HETATM"):
                if pdm.get_chain_from_line(line) != self.ligand_chain:
                    continue
                if pdm.get_resname_from_line(line) != "HOH":
                    continue
                if pdb.get_resname_from_line not in LIST_OF_IONS:
                    continue
                to_keep.append(line)
        self.lines = to_keep

    def correct_chain_and_resnum_new_res_lines(self):
        correction = []
        for line in self.new_res_lines:
            line = pdm.set_chain_to_line(line, self.original_reschain)
            line = pdm.set_resnum_to_line(line, self.original_resnum)
            correction.append(line)
        self.new_res_lines = correction

    def get_indexes_from_list_of_lines(self, list_of_lines):
        indexes = []
        for line in list_of_lines:
            idx = pdm.get_atom_index_from_line(line)
            indexes.append(idx)
        return indexes

    def add_lines_at_position(self, lines, position):
        lines.reverse()
        for line in lines:
            self.lines.insert(position, line)
   
    def find_index_of_line(self, line):
        index = self.lines.index(line)
        return index

    def remove_lines(self, lines):
        for line in lines:
            self.lines.remove(line)
        

def correct_pdb(pdb_file, reschain, resnum, new_resname):
    correct = Corrector(pdb_file, reschain, resnum, new_resname)
    correct.overwrite_pdb()

def delete_atom_ligand(pdb_file, ligchain, atomname):
    if len(atomname) != 4:
        raise ValueError("The atom name of the linking atom must be a string of length 4! Current length: {}".format(len(atomname)))
    pdb = pdm.PDB(pdb_file)
    atomlines = pdb.atom_section.split("\n\n")
    for l in atomlines:
        at_name = l[12:16]
        at_chain = l[21:22]
        if at_name == atomname and at_chain == ligchain:
            print("Deleting {} from {}".format(atomname, ligchain))
            for n, line in enumerate(pdb.lines):
                if l in line:
                    del pdb.lines[n]
            pdb.content = "".join(pdb.lines)
            pdb.overwrite_pdb()
