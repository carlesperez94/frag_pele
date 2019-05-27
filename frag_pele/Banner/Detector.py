import prody
import logging
import glob
import multiprocessing as mp

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)


class Detector:
    def __init__(self, pdb, threshold, dihedrals, lig_chain="L"):
        self.pdb = pdb
        self.threshold = threshold
        self.lig_chain = lig_chain
        self.dihedrals = dihedrals

    def read_pdb(self):
        pdb = prody.parsePDB(self.pdb)
        return pdb

    def get_ligand(self):
        pdb = self.read_pdb()
        ligand = pdb.select("chain {}".format(self.lig_chain))
        if ligand is None:
            logger.critical("Wrong chain selected!")
        elif ligand.ishetero:
            return ligand
        else:
            logger.critical("The selected chain does not contain heteroatoms!")

    def select_atoms(self):
        ligand = self.get_ligand()
        dihedral_list = []
        for dihedral in self.dihedrals:
            atoms = ligand.select("name {} {} {} {}".format(dihedral[0], dihedral[1], dihedral[2], dihedral[3]))
            new_list = []
            for atom in dihedral:
                for atom2 in atoms:
                    if atom == atom2.getName():
                        new_list.append(atom2)
                        break
            dihedral_list.append(new_list)
        return dihedral_list

    def read_dihedral(self):

        dihedral_list = self.select_atoms()
        dihedral_results = []
        for dihedral in dihedral_list:
            atoms = []
            for atom in dihedral:
                atoms.append(atom)
            dihedral_angle = prody.calcDihedral(atoms[0], atoms[1], atoms[2], atoms[3])
            dihedral_results.append(dihedral_angle)
        return dihedral_results

    def check_threshold_dihedral(self):
        value_list = self.read_dihedral()
        check = True
        for value in value_list:
            if abs(self.threshold) > value:
                check = False
        return check


def checker(pdb, threshold, dihedrals, lig_chain):
    new_detector = Detector(pdb, threshold, dihedrals, lig_chain)
    return new_detector.check_threshold_dihedral(), pdb


def check_folder(folder, threshold, dihedrals, lig_chain="L", processors=4):
    dict_of_checks = {}
    list_of_pdbs = glob.glob("{}/*.pdb".format(folder))
    multi = []
    pool = mp.Pool(processors)
    for pdb in list_of_pdbs:
        multi.append(pool.apply_async(checker, [pdb, threshold, dihedrals, lig_chain]))
    for process in multi:
        flag, pdb = process.get()
        dict_of_checks[pdb] = flag
    return dict_of_checks






