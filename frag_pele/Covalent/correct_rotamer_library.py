import numpy as np
class RotamerModifier():
    def __init__(self, input_file):
        self.input_file = input_file
        self.content = None
        self.lines = []
        self.lines_splitted = []
        self.read_file()
        self.read_rotamers()

    def read_file(self):
        with open(self.input_file) as inp:
            self.lines = inp.readlines()
   
    def read_rotamers(self):
        for line in self.lines:
            split = line.split(" ")
            self.lines_splitted.append(split)

    def find_rotamer(self, atom1, atom2):
        indexes = []
        for n, line in enumerate(self.lines_splitted):
            if line[3] == "sidelib":
                if (line[5] == atom1 and line[6] == atom2) or \
                   (line[5] == atom2 and line[6] == atom1):
                    indexes.append(n)
        return indexes # The index of the list

    def delete_rotamer(self, atom1, atom2):
        indexes_to_del = self.find_rotamer(atom1, atom2)
        for index_to_del in indexes_to_del:
            print("Rotamer '{}' will be deleted".format(" ".join(self.lines_splitted[index_to_del]).strip("\n")))
        np.delete(self.lines_splitted, indexes_to_del)
        self.update_content()

    def update_content(self):
        self.lines = []
        for line in self.lines_splitted:
            self.lines.append(" ".join(line))
        self.content = "".join(self.lines)

    def overwrite_file(self):
        with open(self.input_file, "w") as out:
            out.write(self.content)

def delete_atoms_from_rot_lib(input_file, list_of_atom_pairs):
    rot = RotamerModifier(input_file)
    for pair in list_of_atom_pairs:
        rot.delete_rotamer(pair[0], pair[1])
    rot.overwrite_file()

def delete_atoms_from_rot_lib_by_list(input_file, list_of_atoms):
    rot = RotamerModifier(input_file)
    rot.read_rotamers()
    rotamers = rot.lines_splitted
    for r in rotamers:
        if r[3] == "sidelib":
            if r[5] in list_of_atoms and r[6] in list_of_atoms:
                rot.delete_rotamer(r[5], r[6])
    rot.overwrite_file()
