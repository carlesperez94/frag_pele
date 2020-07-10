import lib_prep.pdb_modifier as pdm


class Corrector(pdm.PDB):
    def __init__(self, pdb_file, original_reschain, original_resnum, new_resname):
        pdm.PDB.__init__(self, in_pdb=pdb_file)
        self.original_reschain = original_reschain
        self.original_resnum = original_resnum
        self.new_resname = new_resname
        self.old_res_lines = self.get_residue(self.original_reschain,
                                              self.original_resnum)
        self.new_res_lines = self.get_atoms_of_resname(self.new_resname)
        self.remove_lines(self.new_res_lines)
        self.correct_atom_section_new_res_lines()
        self.correct_chain_and_resnum_new_res_lines()
        self.add_lines_at_position(self.new_res_lines, 
                                   self.find_index_of_line(self.old_res_lines[0]))
        self.remove_lines(self.old_res_lines)


    def correct_atom_section_new_res_lines(self):
        correction = []
        for line in self.new_res_lines:
            if line.startswith("HETATM"):
                line = list(line)
                line[0:6] = "ATOM  "
                line = "".join(line)
            correction.append(line)
        self.new_res_lines = correction

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
    with open("test.pdb", "w") as out:
        out.write("".join(correct.lines))

correct_pdb("/gpfs/projects/bsc72/FragPELE/FragPELE2.3.0_testing/frag_pele/tests/receptor_145_cys_frag_resSGC4/pregrow/initialization_grow.pdb",
            "A", "145", "GRW")


