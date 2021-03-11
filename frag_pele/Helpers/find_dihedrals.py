import os
import glob
import numpy as np
import networkx as nx
import peleffy
from frag_pele.constants import SCHRODINGER
from peleffy.topology import Topology
from peleffy.forcefield import OpenForceField, OPLS2005ForceField
from peleffy.topology import Molecule
from peleffy.utils import Logger
from peleffy.topology import molecule
from peleffy.utils.toolkits import RDKitToolkitWrapper

class ComputeDihedrals(object):
    """
    A class to produce a library of dihedral angles.
    """
    def __init__(self, pdb_file, forcefield='OPLS2005'):
        """
        Initializes a ComputeDihedrals object
        Parameters
        ----------
        topology_obj : An peleffy.topology.Topology
            A Topology object that contains the ligand's information
        mode: str
            Whether to extract all dihedrals or only those marked as flexible
        """
        self._pdb_file = pdb_file
        self._forcefield = forcefield
        self._topology, self._molecule = self.load_molecule()
        self.dihedral_library = {}

    def load_molecule(self):
        os.environ['SCHRODINGER'] = SCHRODINGER 
        m = Molecule(self._pdb_file)
        if self._forcefield == 'OPLS2005':
            ff = OPLS2005ForceField()      
        elif self._forcefield == 'OpenForceField':
            ff = OpenForceField('openff_unconstrained-1.2.0.offxml')
        else:
            raise ValueError("Not valid ForceField. Pick between 'OPLS2005' and 'OpenForceField'")
        parameters = ff.parameterize(m)                                        
        top = peleffy.topology.Topology(m, parameters) 
        return top, m

    def list_all_dihedrals(self):
        dihedrals = []
        seen_dihedrals = set()
        for proper in self._topology.propers:
            dihedral_indexes = (proper.atom1_idx, proper.atom2_idx, proper.atom3_idx, proper.atom4_idx)
            if dihedral_indexes in seen_dihedrals:
                continue
            seen_dihedrals.add(dihedral_indexes)
            dihedrals.append(list(dihedral_indexes))
        return dihedrals

    def calculate_cluster_angles(self, dihedral_list):
        """
        Calculate dihedral angles from pdb
        Parameters
        ----------
        pdb_file: str
            Path to the cluster representative conformation
        dihedral_list: list
            List of the tuples containing the atoms that form the dihedrals
        match_indexes: bool
            Whether to use the atom indices from the dihedral list or match to
            the cluster structure before
        """
        rdkit_wrapper = RDKitToolkitWrapper()
        pdb_dihedrals = []
        # use the input molecule as template since the cluster structures
        # probably will not have proper stereochemistry
        mol = molecule.Molecule(self._pdb_file, connectivity_template=self._molecule.rdkit_molecule)
        # we use substructure matching to ensure that the indices in the
        # clusters pdb and the input ligand are the same
        for dihedral in dihedral_list:
            names = [self._topology.atoms[atom].PDB_name for atom in dihedral]
            angle = get_dihedral(mol, *dihedral, units="degrees")
        
            pdb_dihedrals.append(names+[angle])
        self.dihedral_library[self._pdb_file] = pdb_dihedrals
    
    def _calculate_all_dihedrals(self):
        dihedrals = self.list_all_dihedrals()
        self.calculate_cluster_angles(dihedrals)


    def calculate(self):
        """
        Calculate dihedrals library from the bce output
        """
        logger = Logger()
        logger.info(' - Calculating dihedral library')
        self._calculate_all_dihedrals()

def get_dihedral(mol, atom1, atom2, atom3, atom4, units="radians"):
    """
    It calculates the value of the dihedral angle in the specified units
        (default radians)
    Parameters
    ----------
    molecule : an offpele.topology.Molecule
        The offpele's Molecule object
    atom1 : int
        Index of the first atom in the dihedral
    atom2 : int
        Index of the second atom in the dihedral
    atom3 : int
        Index of the third atom in the dihedral
    atom4 : int
        Index of the fourth atom in the dihedral
    units : str
        The units in which to calculate the angle (default is radians, can
        be radians or degrees)
    """
    from rdkit.Chem import rdMolTransforms
    if units == "degrees":
        angle = rdMolTransforms.GetDihedralDeg(mol.rdkit_molecule.GetConformer(), atom1, atom2, atom3, atom4)
    else:
        angle = rdMolTransforms.GetDihedralRad(mol.rdkit_molecule.GetConformer(), atom1, atom2, atom3, atom4)
    return angle


def select_dihedrals(input_dihedrals_list, selected_dihedrals_atoms):
    selected_dihedrals = []
    for dihedrals in input_dihedrals_list:
        atom1, atom2, atom3, atom4, angle = dihedrals
        atoms = [atom1, atom2, atom3, atom4]
        for atoms_of_dih in selected_dihedrals_atoms:
            if atoms_of_dih == atoms or atoms_of_dih[::-1] == atoms: # Try in both reading senses
                print(f"Constraining dihedral: {atoms}")
                selected_dihedrals.append([*atoms, angle])
    if len(selected_dihedrals) == 0:
        raise ValueError(f"The dihedral formed by {selected_dihedrals_atoms} were not found in ligand dihedrals.")
    return selected_dihedrals

