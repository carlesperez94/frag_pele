import os
import time
import subprocess
from frag_pele.Helpers import correct_fragment_names
import frag_pele.constants as c


class CovalentPDB:

    content = None
    ligand = None
    residue_type_linked = None
    conect_str = []
    conects = []

    def __init__(self, pdb_file, lig_chain="L", res_num=None, 
                 res_chain=None, lig_link=None, res_link=None):
        self.pdb_file = pdb_file
        self.lig_chain = lig_chain
        self.res_num = res_num
        self.res_chain = res_chain
        self.lig_link = lig_link
        self.res_link = res_link
        self._read_pdb_content()
        if not self._check_if_exists_conects() and (res_num==None or res_chain==None):
            raise ValueError("PDB file does not contain CONECTs." 
                             "Please, fill the number and chain of the residue covalently linked to the ligand")
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
            if l.startswith("HETATM"):
                if self.lig_chain == l[21]:
                    lig_lines.append(l)
        self.ligand = "\n".join(lig_lines)

    def _read_conects(self):
        for l in self.content.split("\n"):
            if l.startswith("CONECT"):
                group = l.split("CONECT")[1].split()
                group = list(map(int, group))
                self.conects.append(group)
    
    def _read_conects_str(self):
        lines = []
        for l in self.content.split("\n"):
            if l.startswith("CONECT"):
                lines.append(l)
        return "\n".join(lines)

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
                self.residue_type_linked = line[17:20].strip()
                print(line[17:20].strip())

    def get_residue_attached(self):
        lines = []
        for l in self.content.split("\n"):
            if l.startswith("ATOM"):
                if l[23:26].strip() == self.res_num and l[21] == self.res_chain:
                    lines.append(l)
        return "\n".join(lines)
  
    def replace_residue_name(self, new_name):
        lines = []
        for l in self.content.split("\n"):
            if l.startswith("ATOM"):
                if l[23:26].strip() == self.res_num and l[21] == self.res_chain:
                    l = list(l)
                    l[17:20] = new_name
                    l = "".join(l)
            lines.append(l)
        self.content = "\n".join(lines)

    def get_ligand_and_res(self):
        lines = []
        mix = self.get_residue_attached() + "\n" + self.ligand
        for l in mix.split("\n"):
            l = list(l)
            l[0:6] = "HETATM"
            l[17:20] = "MIX"
            l[21:22] = "L"
            l = "".join(l)
            lines.append(l)
        mix = "\n".join(lines)
        con = self._read_conects_str()
        return mix + "\n" + con

    def get_pdb_prepared(self):
        restype = self.residue_type_linked
        new_resname = restype[0:2] + "Y"
        self.replace_residue_name(new_resname) 
        

def get_indexes(pdb_fragment):
    indexes = []
    for l in pdb_fragment.split("\n"):
        if l.startswith("HETATM") or l.startswith("ATOM"):
            indexes.append(int(l[7:11].strip()))
    return indexes

def write_string_to_file(string, filename):
    with open(filename, "w") as outf:
        outf.write(string)

def read_file(filename):
    with open(filename) as inf:
        cont = inf.read()
    return cont

def prepare_pdb(pdb_in, pdb_out, sch_path):
    command = [os.path.join(sch_path, "utilities/prepwizard"), pdb_in, pdb_out, "-noepik", "-noprotassign",
              "-noccd", "-noimpref"]
    subprocess.call(command)

def check_and_solve_duplicate_atomnames(pdb_file, chain="L"):
    content = read_file(pdb_file)
    elements = []
    pdb_out = []
    names = []
    new_names = []
    for line in content.split("\n"):
        if line[0:6] == "HETATM" and line[21:22] == chain:
            element = line[76:78]
            elements.append(element)
            atom_names = line[12:16]
            names.append(atom_names)
            list_of_unique_elements = set(elements)
            counters = {i: elements.count(i) for i in list_of_unique_elements}
            counter = counters[line[76:78]]
            set_to_check = set(names)
            list_to_check = sorted(list(set_to_check))
            sorted_list_names = sorted(names)
            if list_to_check != sorted_list_names:
                new_name = "{}{}".format(line[76:78].strip().upper(), counter)
                while True:
                    new_name = new_name + " "
                    if len(new_name) == 4:
                        break
                line = list(line)
                line[12:16] = new_name
                line = "".join(line)
                if len(new_name) > 4:
                    raise ValueError("Length of the string {} is too long. Only 4 characters accepted.".format(new_name))
                new_names.append(new_name)
        pdb_out.append(line)
    write_string_to_file("\n".join(pdb_out), pdb_file)
    return new_names

def create_template(pdb_file, sch_python=c.SCHRODINGER_PY_PATH, 
                    plop_script_path="../PlopRotTemp_S_2017/ligand_prep.py", rotamers="10.0", 
                    out_templates_path=".", path_to_lib="."):
    cmd = "{} {} {} {} {} {}".format(sch_python, plop_script_path, pdb_file, rotamers,
                                     out_templates_path, path_to_lib)
    print(cmd)
    try:
        subprocess.call(cmd.split())
    except OSError:
        raise OSError("Path {} not foud. Change schrodinger path under frag_pele/constants.py".format(sch_python))

def prepare_complex_and_extract_files_for_covalent(pdb_file, sch_path="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC",
                                                  lig_chain="L", res_num=None, res_chain=None, lig_link=None, 
                                                  res_link=None, rotamers="10.0", out_templates_path=".", 
                                                  path_to_lib="."):
    pdb = CovalentPDB(pdb_file, lig_chain, res_num, res_chain, lig_link, res_link)
    # Modifying PDB content
    pdb.get_pdb_prepared()
    mix = pdb.get_ligand_and_res()
    res = pdb.get_residue_attached()
    # Creating LIGAND, RESIDUE, MIX(COMBINATION) and COMPLEX
    write_string_to_file(pdb.ligand, "LIG.pdb")
    write_string_to_file(res, "{}.pdb".format(pdb.residue_type_linked))
    write_string_to_file(mix, "MIX.pdb")
    write_string_to_file(pdb.content, pdb_file)
    # Preparing them with Schrodinger (Correcting free atoms)
    prepare_pdb("LIG.pdb", "LIG_h.pdb", sch_path=sch_path)
    prepare_pdb("MIX.pdb", "MIX_h.pdb", sch_path=sch_path)
    prepare_pdb("{}.pdb".format(pdb.residue_type_linked), 
                "{}_h.pdb".format(pdb.residue_type_linked), 
                sch_path=sch_path)
    time.sleep(3) # Waiting queue system
    # Solving duplicate names and getting the new H atoms
    h_to_del = check_and_solve_duplicate_atomnames("LIG_h.pdb", chain=lig_chain)
    check_and_solve_duplicate_atomnames("MIX_h.pdb", chain=lig_chain)
    # Creating templates of LIGAND and COMBINATION
    create_template("LIG_h.pdb", rotamers=rotamers, out_templates_path=out_templates_path, 
                    path_to_lib=path_to_lib)
    create_template("MIX_h.pdb", rotamers=rotamers, out_templates_path=out_templates_path, 
                    path_to_lib=path_to_lib)
    print("Files prepared successfully!")
    return pdb.residue_type_linked, pdb.res_link, pdb.lig_link, h_to_del


#def prepare_templates_for_covalent():


def prepare_covalent_growing(pdb_file, sch_path="/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC",
                             lig_chain="L", res_num=None, res_chain=None, lig_link=None,
                             res_link=None, rotamers="10.0", out_templates_path=".",
                             path_to_lib="."):
    prepare_complex_and_extract_files_for_covalent("covdock_sar26lu7_14.pdb")

