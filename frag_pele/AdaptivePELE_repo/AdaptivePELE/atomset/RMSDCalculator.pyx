from __future__ import unicode_literals
import numpy as np
cimport cython
cimport numpy as np
cimport frag_pele.AdaptivePELE_repo.AdaptivePELE.atomset.atomset as atomset


cdef class RMSDCalculator:
    def __init__(self, symmetries=[]):
        """
            :param symmetries: List of dictionaries with gropus of symmetric atoms atomId:symmetricalAtomId corresponding with the symmetrical atoms
            :type symmetries: list
        """
        self.nonSymmetricalAtomsSet = None
        self.symmetries = symmetries

    def __getstate__(self):
        state = {'nonSymmetricalAtomsSet': self.nonSymmetricalAtomsSet,
                 'symmetries': self.symmetries}
        return state

    def __setstate__(self, state):
        self.nonSymmetricalAtomsSet = state['nonSymmetricalAtomsSet']
        self.symmetries = state['symmetries']

    def computeNonSymmAtoms(self, atomset.PDB PDB):
        cdef set allAtomsSet
        allAtomsSet = set(PDB.atoms.keys())
        for group in self.symmetries:
            symmetriesSet = set(group.keys()).union(set(group.values()))
            allAtomsSet -= symmetriesSet
        self.nonSymmetricalAtomsSet = allAtomsSet

    def computeRMSD2(self, atomset.PDB PDB1, atomset.PDB PDB2):
        """
            Compute the squared RMSD between two PDB

            :param PDB1: First PDB with which the RMSD will be calculated
            :type PDB1: PDB
            :param PDB2: First PDB with which the RMSD will be calculated
            :type PDB2: PDB
            :returns: float -- The squared RMSD between two PDB
        """
        cdef double rmsd, d2, d2sm
        cdef basestring atom1Id, atom2Id, atomId
        cdef atomset.Atom atom11, atom12, atom21, atom22, atom1, atom2
        cdef dict group
        cdef int n
        rmsd = 0.0
        if self.nonSymmetricalAtomsSet is None:
            self.computeNonSymmAtoms(PDB1)

        for group in self.symmetries:
            d2 = 0.0
            d2sm = 0.0
            for atom1Id, atom2Id in group.iteritems():
                try:
                    atom11 = PDB1.getAtom(atom1Id)
                    atom12 = PDB1.getAtom(atom2Id)
                    atom21 = PDB2.getAtom(atom1Id)
                    atom22 = PDB2.getAtom(atom2Id)
                except KeyError as err:
                    raise KeyError("Atom %s not found in PDB" % err.message)
                d2 += atom11.squaredDistance(atom21) + atom12.squaredDistance(atom22)
                d2sm += atom12.squaredDistance(atom21) + atom11.squaredDistance(atom22)
            rmsd += min(d2, d2sm)
        for atomId in self.nonSymmetricalAtomsSet:
            try:
                atom1 = PDB1.getAtom(atomId)
                atom2 = PDB2.getAtom(atomId)
            except KeyError as err:
                raise KeyError("Atom %s not found in PDB" % err.message)
            rmsd += atom1.squaredDistance(atom2)
        n = len(PDB1.atoms.items())
        return rmsd/n

    def computeRMSD(self, atomset.PDB PDB1, atomset.PDB PDB2):
        """
            Compute the RMSD between two PDB

            :param PDB1: First PDB with which the RMSD will be calculated
            :type PDB1: PDB
            :param PDB2: First PDB with which the RMSD will be calculated
            :type PDB2: PDB
            :returns: float -- The squared RMSD between two PDB
        """
        return np.sqrt(self.computeRMSD2(PDB1, PDB2))
