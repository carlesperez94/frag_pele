from __future__ import unicode_literals
import atomset
import numpy as np
cimport cython
cimport numpy as np
cimport AdaptivePELE.atomset.atomset as atomset


cdef class SymmetryContactMapEvaluator:
    def __init__(self, symmetries=[]):
        """
            :param symmetries: List of dictionaries with gropus of symmetric atoms atomId:symmetricalAtomId corresponding with the symmetrical atoms
            :type symmetries: list
        """
        self.symmetricAtoms = set()
        for group in symmetries:
            self.symmetricAtoms.update(set(group.keys()).union(set(group.values())))
        self.symmetries = symmetries
        self.symToRowMap = {}
        self.ligandList = []
        self.proteinList = []


    def __getstate__(self):
        state = {'symmetries': self.symmetries, 'symToRowMap': self.symToRowMap,
                 'symmetricAtoms': self.symmetricAtoms,
                 'ligandList': self.ligandList, 'proteinList': self.proteinList}
        return state

    def __setstate__(self, state):
        self.symmetries = state['symmetries']
        self.symToRowMap = state['symToRowMap']
        self.symmetricAtoms = state['symmetricAtoms']
        self.ligandList = state.get('ligandList', [])
        self.proteinList = state.get('proteinList', [])

    def createContactMap(self, atomset.PDB PDB, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandResChain=""):
        """
            Create the contact map of the protein and ligand. The contact map is
            a boolean matrix that has as many rows as the number of ligand heavy
            atoms and as many columns as the number of alpha carbons. The
            value is one if the ligand atom and the alpha carbons are less than
            contactThresholdDistance Amstrongs away

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param ligandResnum: Integer containing the residue number of the ligand in the pdb
            :type ligandResnum: int
            :param ligandResChain: String containing the chain name of the ligand in the pdb
            :type ligandResChain: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        if self.symToRowMap:
            return self.buildContactMap(PDB, ligandResname, contactThresholdDistance, ligandResnum, ligandResChain)
        else:
            return self.buildContactMapWithRowMap(PDB, ligandResname, contactThresholdDistance, ligandResnum, ligandResChain)

    def buildContactMapWithRowMap(self, atomset.PDB PDBobj, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandResChain=""):
        """
            Create a map that relates symmetric atoms in the ligand to row
            indices in the contactMap matrix and the contact map of the protein
            and ligand. The contact map is a boolean matrix that has as many
            rows as the number of ligand heavy atoms and as many columns as the
            number of alpha carbons. The value is one if the ligand atom and
            the alpha carbons are less than contactThresholdDistance Amstrongs
            away

            :param PDBobj: PDB object that contains the conformation's info
            :type PDBobj: PDB
            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param ligandResnum: Integer containing the residue number of the ligand in the pdb
            :type ligandResnum: int
            :param ligandResChain: String containing the chain name of the ligand in the pdb
            :type ligandResChain: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        cdef int contactThresholdDistance2, rowind, colind
        contactThresholdDistance2 = contactThresholdDistance**2
        cdef atomset.PDB ligandPDB, alphaCarbonsPDB
        cdef np.ndarray[char, ndim=2] contactMap
        cdef basestring ligandAtomID, proteinAtomID
        cdef atomset.Atom proteinAtom, ligandAtom
        cdef set contacts
        cdef double dist2
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(PDBobj.pdb, resname=ligandResname, resnum=ligandResnum,
                             chain=ligandResChain, heavyAtoms=True)
        alphaCarbonsPDB = atomset.PDB()
        alphaCarbonsPDB.initialise(PDBobj.pdb, type="CM")
        if not self.ligandList:
            self.ligandList = ligandPDB.atomList
        if not self.proteinList:
            self.proteinList = alphaCarbonsPDB.atomList
        # empty contact map, rows are atoms of the ligand, columns are protein
        # alpha carbons
        contactMap = np.zeros((len(ligandPDB.atomList),
                               len(alphaCarbonsPDB.atomList)), dtype=np.uint8)
        contacts = set([])
        for rowind in range(len(ligandPDB.atomList)):
            ligandAtomID = ligandPDB.atomList[rowind]
            ligandAtom = ligandPDB.atoms[ligandAtomID]
            if ligandAtomID in self.symmetricAtoms:
                self.symToRowMap[ligandAtomID] = rowind
            for colind in range(len(alphaCarbonsPDB.atomList)):
                proteinAtomID = alphaCarbonsPDB.atomList[colind]
                proteinAtom = alphaCarbonsPDB.atoms[proteinAtomID]
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contactMap[rowind, colind] = True
                if proteinAtom.name == "CA" and dist2 < 64.0:
                    # Contact ratio will be always calculated using a contact
                    # threshold of 8, so that tresholds and densities are
                    # independent of the contact threshold of the contactMap
                    contacts.update([proteinAtomID])
        return contactMap.view(np.bool), len(contacts)

    def buildContactMap(self, atomset.PDB PDBobj, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandResChain=""):
        """
            Create the contact map of the protein and ligand. The contact map is
            a boolean matrix that has as many rows as the number of ligand heavy
            atoms and as many columns as the number of alpha carbons. The
            value is one if the ligand atom and the alpha carbons are less than
            contactThresholdDistance Amstrongs away

            :param PDBobj: PDB object that contains the conformation's info
            :type PDBobj: PDB
            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param ligandResnum: Integer containing the residue number of the ligand in the pdb
            :type ligandResnum: int
            :param ligandResChain: String containing the chain name of the ligand in the pdb
            :type ligandResChain: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: numpy.Array -- The contact map of the ligand and the protein
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        cdef int contactThresholdDistance2, rowind, colind
        contactThresholdDistance2 = contactThresholdDistance**2
        cdef atomset.PDB ligandPDB, alphaCarbonsPDB
        cdef np.ndarray[char, ndim=2] contactMap
        cdef basestring ligandAtomID, proteinAtomID
        cdef atomset.Atom proteinAtom, ligandAtom
        cdef set contacts
        cdef double dist2
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(PDBobj.pdb, resname=ligandResname, resnum=ligandResnum,
                             chain=ligandResChain, heavyAtoms=True)
        alphaCarbonsPDB = atomset.PDB()
        alphaCarbonsPDB.initialise(PDBobj.pdb, type="CM")
        if not self.ligandList:
            self.ligandList = ligandPDB.atomList
        if not self.proteinList:
            self.proteinList = alphaCarbonsPDB.atomList
        # empty contact map, rows are atoms of the ligand, columns are protein
        # alpha carbons
        contactMap = np.zeros((len(ligandPDB.atomList),
                               len(alphaCarbonsPDB.atomList)), dtype=np.uint8)
        contacts = set([])
        for rowind in range(len(ligandPDB.atomList)):
            ligandAtomID = ligandPDB.atomList[rowind]
            ligandAtom = ligandPDB.atoms[ligandAtomID]
            if ligandAtomID in self.symmetricAtoms:
                self.symToRowMap[ligandAtomID] = rowind
            for colind in range(len(alphaCarbonsPDB.atomList)):
                proteinAtomID = alphaCarbonsPDB.atomList[colind]
                proteinAtom = alphaCarbonsPDB.atoms[proteinAtomID]
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contactMap[rowind, colind] = True
                if proteinAtom.name == "CA" and dist2 < 64.0:
                    # Contact ratio will be always calculated using a contact
                    # threshold of 8, so that tresholds and denisities are
                    # independent of the contact threshold of the contactMap
                    contacts.update([proteinAtomID])
        return contactMap.view(np.bool), len(contacts)

    def evaluateJaccard(self, contactMap, clusterContactMap):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                clusterContactMap)
        intersectContactMaps = (permContactMap & clusterContactMap).sum()
        unionContactMaps = permContactMap.sum() + clusterContactMap.sum() - intersectContactMaps
        if unionContactMaps < 1e-7:
            # both contactMaps have zero contacts
            return 0.0
        similarity = float(intersectContactMaps)/unionContactMaps
        distance = 1-similarity
        return distance

    def evaluateCorrelation(self, contactMap, clusterContactMap):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                clusterContactMap)
        similarity = calculateCorrelationContactMaps(permContactMap, clusterContactMap)
        similarity += 1  # Necessary to omit negative correlations
        similarity /= 2.0  # Correlation values need to be higher now
        distance = 1-similarity
        return distance

    def evaluateDifferenceDistance(self, contactMap, clusterContactMap):
        permContactMap = self.buildOptimalPermutationContactMap(contactMap,
                                                                clusterContactMap)
        differenceContactMaps = np.abs(permContactMap^clusterContactMap).sum()
        averageContacts = (0.5*(permContactMap.sum()+clusterContactMap.sum()))
        if not averageContacts:
            # The only way the denominator can be 0 is if both contactMaps are
            # all zeros, thus being equal and belonging to the same cluster
            return 0.0
        else:
            distance = differenceContactMaps/averageContacts
            return distance

    def buildOptimalPermutationContactMap(self, contactMap, clusterContactMap):
        """
            Build a permutated version of the contactMap which maximizes the
            similarity between the two contactMaps
            :param contactMap: contactMap of the conformation to compare
            :type contactMap: numpy.Array
            :param clusterContactMap: cluster contactMap to which we are comparing
            :type cluster: numpy.Array
            :returns: numpy.Array -- The permutated contact map
        """
        if not self.symmetries and not self.symToRowMap:
            return contactMap
        permContactMap = np.array(contactMap)
        for group in self.symmetries:
            for atom1Id, atom2Id in group.iteritems():
                try:
                    atom1Row = self.symToRowMap[atom1Id]
                    atom2Row = self.symToRowMap[atom2Id]
                    line11 = clusterContactMap[atom1Row]
                    line12 = clusterContactMap[atom2Row]
                    line21 = contactMap[atom1Row]
                    line22 = contactMap[atom2Row]
                except KeyError as err:
                    raise KeyError("Atom %s not found in symmetries" % err.message)
                d2 = (line11 == line21).sum() + (line12 == line22).sum()
                d2sm = (line12 == line21).sum() + (line11 == line22).sum()
                if d2sm > d2:
                    # swap rows using advanced slicing from numpy
                    permContactMap[[atom1Row, atom2Row], :] = contactMap[[atom2Row, atom1Row], :]
        return permContactMap


def calculateCorrelationContactMaps(contactMap, clusterContactMap):
    """
        Calculate the correlation of two contactMaps
    """
    contactMap1 = contactMap.reshape((1, -1))
    # Reshape the matrix into a 1D array to use the numpy corrcoef function
    contactMap2 = clusterContactMap.reshape((1, -1))
    total1 = contactMap1.sum()
    total2 = contactMap2.sum()
    if not total1 or not total2:
        # if any array is all zeros the correlation will be NaN
        # if both are zero, correlation is perfect (i.e 1) else it is different
        # and will be in different clusters
        return total1 == total2
    return np.corrcoef(contactMap1, contactMap2)[0, 1]
