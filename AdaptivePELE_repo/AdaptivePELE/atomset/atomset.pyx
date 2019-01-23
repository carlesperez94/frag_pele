from __future__ import unicode_literals
import numpy as np
import re
from io import StringIO, open
cimport cython
cimport numpy as np


cdef class Atom:
    _chargePattern = re.compile("[0-9]|\+|\-")
    _ATOM_WEIGHTS = {"H": 1.00794,
                    "D": 2.01410178,  # deuterium
                    "HE": 4.00,
                    "LI": 6.941,
                    "BE": 9.01,
                    "B": 10.811,
                    "C": 12.0107,
                    "N": 14.0067,
                    "O": 15.9994,
                    "F": 18.998403,
                    "NE": 20.18,
                    "NA": 22.989769,
                    "MG": 24.305,
                    "AL": 26.98,
                    "SI": 28.09,
                    "P": 30.973762,
                    "S": 32.065,
                    "CL": 35.453,
                    "AR": 39.95,
                    "K": 39.0983,
                    "CA": 40.078,
                    "SC": 44.96,
                    "TI": 47.87,
                    "V": 50.94,
                    "CR": 51.9961,
                    "MN": 54.938045,
                    "FE": 55.845,
                    "CO": 58.93,
                    "NI": 58.6934,
                    "CU": 63.546,
                    "ZN": 65.409,
                    "GA": 69.72,
                    "GE": 72.64,
                    "AS": 74.9216,
                    "SE": 78.96,
                    "BR": 79.90,
                    "KR": 83.80,
                    "RB": 85.47,
                    "SR": 87.62,
                    "Y": 88.91,
                    "ZR": 91.22,
                    "NB": 92.91,
                    "W": 95.94,  # Molybdenum?  Not sure why it's not always MO
                    "MO": 95.94,
                    "TC": 98.0,
                    "RU": 101.07,
                    "RH": 102.91,
                    "PD": 106.42,
                    "AG": 107.8682,
                    "CD": 112.411,
                    "IN": 114.82,
                    "SN": 118.71,
                    "SB": 121.76,
                    "TE": 127.60,
                    "I": 126.90447,
                    "XE": 131.29,
                    "CS": 132.91,
                    "BA": 137.33,
                    "PR": 140.91,
                    "EU": 151.96,
                    "GD": 157.25,
                    "TB": 158.93,
                    "IR": 192.22,
                    "PT": 195.084,
                    "AU": 196.96657,
                    "HG": 200.59,
                    "PB": 207.2,
                    "U": 238.03}
    def __init__(self, basestring atomContent=""):
        """ Create an atom from a pdb line

            :param atomContent: Line of the pdb from which the atom will be created
            :type atomContent: basestring
        """
        # Force string attributes to be unicode strings
        self.atomSerial = u""
        self.name = u""
        self.resname = u""
        self.resChain = u""
        self.resnum = u""
        self.type = u""
        # atomContent = atomContent.split()
        if len(atomContent) > 6 and (atomContent[:4] == 'ATOM' or atomContent[:6] == 'HETATM'):
            self.atomSerial = atomContent[6:11].strip()
            self.name = atomContent[12:16].strip()
            self.resname = atomContent[17:20].strip()
            self.resChain = atomContent[21]
            self.resnum = atomContent[22:26].strip()
            self.x = float(atomContent[30:38])
            self.y = float(atomContent[38:46])
            self.z = float(atomContent[46:54])

            self.type = re.sub(self._chargePattern, "", atomContent[76:80]).strip().upper()
            self.mass = self._ATOM_WEIGHTS[self.type]

            if atomContent.startswith('ATOM'):
                self.protein = True
            else:
                self.protein = False

            self.id = self.atomSerial + ":" + self.name + ":" + self.resname
            # self.id = self.atomSerial

    def __getstate__(self):
        # Copy the object's state from
        state = {"atomSerial": self.atomSerial, "name": self.name, "x": self.x,
                 "y": self.y, "z": self.z, "mass": self.mass, "type": self.type,
                 "resname": self.resname, "resChain": self.resChain,
                 "resnum": self.resnum, "protein": self.protein, "id": self.id}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.atomSerial = state['atomSerial']
        self.name = state['name']
        self.resname = state['resname']
        self.resnum = state['resnum']
        self.resChain = state['resChain']
        self.type = state['type']
        self.id = state['id']
        self.mass = state['mass']
        self.x = state['x']
        self.y = state['y']
        self.z = state['z']
        self.protein = state['protein']

    def isHeavyAtom(self):
       """
            Check if Atom is a heavy atom

            :returns: bool -- True if Atom is heavy atom, false otherwise
       """
       return self.type != 'H'

    def isProtein(self):
        """
            Check if Atom is a protein atom

            :returns: bool -- True if Atom is a protein atom, false otherwise
        """
        return self.protein

    def isHeteroAtom(self):
        """
            Check if Atom is an hetero atom

            :returns: bool -- True if Atom is an hetero atom, false otherwise
        """
        return not self.protein

    def printAtom(self):
        """
            Print Atom information
        """
        print self.atomSerial, self.name, self.resname, self.resChain, self.resnum, self.x, self.y, self.z, self.type, self.mass

    def __richcmp__(self, Atom atom2, int op):
        if op == 2:
            #equality
            return self.id == atom2.id
        elif op == 3:
            return self.id != atom2.id
        elif op == 1:
            if self.id == atom2.id:
                return True
            else:
                return self.serial < atom2.serial
        elif op == 5:
            if self.id == atom2.id:
                return True
            else:
                return self.serial > atom2.serial
        elif op == 0:
            return self.serial < atom2.serial
        elif op == 4:
            return self.serial > atom2.serial

    def __str__(self):
        return "%s: %s %s %s [%f, %f, %f] %s %f" % (self.id, self.atomSerial,
                                                    self.resChain, self.resnum,
                                                    self.x, self.y, self.z,
                                                    self.type, self.mass)

    def getAtomCoords(self):
        """
            Get the coordinates of the atom

            :returns: numpy.Array -- Array with the coordinate of the atom
        """
        return np.array([self.x, self.y, self.z])

    def squaredDistance(self, Atom atom2):
        """
            Calculate the squared distance between two atoms

            :param atom2: Second Atom to whom the distance will be calculated
            :type atom2: Atom
            :returns: float -- The distance between the atoms
        """
        return (self.x - atom2.x)**2 + (self.y - atom2.y)**2 + (self.z - atom2.z)**2


cdef class PDB:
    _typeProtein = "PROTEIN"
    _typeHetero = "HETERO"
    _typeAll = "ALL"
    _typeCM = "CM"

    #Atoms to be used in the contact map
    CMAtoms = {"ALA": "empty", "VAL": "empty", "LEU": "empty", "ILE": "empty",
               "MET": "empty", "PRO": "empty", "PHE": "CZ", "TYR": "OH",
               "TRP": "CH2", "SER": "empty", "THR": "empty", "CYS": "empty",
               "ASN": "empty", "GLN": "empty", "LYS": "NZ", "HIS": "CE1",
               "HIE": "CE1", "HID": "CE1", "HIP": "CE1", "ARG": "NE",
               "ASP": "OD1", "GLU": "OE1", "GLY": "empty"}

    def __init__(self):
        """
            Object that will contain the information of a PDB file. Has to call
            the initialise method to load the file
        """
        self.atoms = {}
        # {atomId: atom, ...}
        # Where atomId := serial:atomName:resName
        self.totalMass = 0
        # ensure every string is unicode
        self.pdb = u""
        self.com = None
        self.centroid = None

        # Necessary for contactMaps
        self.atomList = []

    def __richcmp__(self, object other, int op):
        """
            Compare two pdb strings, remark lines should be ignored and only the
            atoms and its information should be compared
        """
        cdef list pdb1, pdb2
        if op == 2:
            pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
            pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
            return pdb1 == pdb2
        elif op == 3:
            pdb1 = [element.strip() for element in self.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
            pdb2 = [element.strip() for element in other.pdb.split('\n') if element.startswith("ATOM") or element.startswith("HETATM")]
            return pdb1 != pdb2
        else:
            print "No boolean operator available for PDB apart from equality"

    def __getstate__(self):
        # Copy the object's state from
        state = {"atoms": self.atoms, "atomList": self.atomList,
                 "com": self.com, "centroid": self.centroid,
                 "totalMass": self.totalMass, "pdb": self.pdb}
        return state


    def __setstate__(self, state):
        # Restore instance attributes
        self.atoms = state['atoms']
        self.atomList = state['atomList']
        self.com = state.get('com')
        self.centroid = state.get('centroid')
        self.totalMass = state['totalMass']
        self.pdb = state['pdb']

    def initialise(self, basestring PDBstr, bint heavyAtoms=True, basestring resname="", basestring atomname="", basestring type="ALL", basestring chain="", int resnum = 0):
        """
            Load the information from a PDB file or a string with the PDB
            contents

            :param PDBstr: may be a path to the PDB file or a string with the contents of the PDB
            :type PDBstr: basestring
            :param heavyAtoms: wether to consider only heavy atoms (True if onl y heavy atoms have to be considered)
            :type heavyAtoms: bool
            :param resname: Residue name to select from the pdb (will only select the residues with that name)
            :type resname: basestring
            :param atomname: Residue name to select from the pdb (will only select the atoms with that name)
            :type atomname: basestring
            :param type: type of atoms to select: may be ALL, PROTEIN or HETERO
            :type type: basestring
            :param chain: Chain name to select from the pdb (will only select the atoms with that name)
            :type chain: basestring
            :param resnum: Residue number to select from the pdb (will only select the atoms with that name)
            :type atomname: int
            :raises: ValueError if the pdb contained no atoms
        """
        cdef object PDBContent
        cdef list stringWithPDBContent
        cdef int atomLineNum
        cdef basestring atomName, resName, atomLine, resnumStr
        cdef Atom atom
        if resnum == 0:
            resnumStr = ""
        else:
            resnumStr = str(resnum)
        PDBContent = StringIO(readPDB(PDBstr))  # Using StringIO
        # creates a buffer that can handle a pdb file or a string containing
        # the PDB
        self.pdb = PDBContent.read()  # in case one wants to write it

        stringWithPDBContent = self.pdb.split('\n')
        for atomLine in stringWithPDBContent:
            if not atomLine.startswith("ATOM") and not atomLine.startswith("HETATM"):
                continue
            if type == self._typeCM:
                atomName = atomLine[12:16].strip()
                resName = atomLine[17:20].strip()
                if resName not in self.CMAtoms:
                    continue
                if atomName != "CA" and atomName != self.CMAtoms[resName]:
                    continue
            else:
                # HUGE optimisation (not to create an atom each time ~ 2e-5 s/atom)
                if resname != "" and not atomLine[17:20].strip() == resname:
                    continue
                if atomname != "" and not atomLine[12:16].strip() == atomname:
                    continue
                if chain != "" and not atomLine[21:22].strip() == chain:
                    continue
                if resnumStr != "" and not atomLine[22:26].strip() == resnumStr:
                    continue

            atom = Atom(atomLine)
            # Here atom will be not null, empty or not.
            # With "try", we prune empty atoms
            try:
                if (not heavyAtoms or atom.isHeavyAtom()) and\
                   (type == self._typeAll or type == self._typeCM or (type == self._typeProtein and atom.isProtein()) or (type == self._typeHetero and atom.isHeteroAtom())):
                        self.atoms.update({atom.id: atom})
                        self.atomList.append(atom.id)
            except:
                pass
        if self.atoms == {}:
            raise ValueError('The input pdb file/string was empty, no atoms loaded!')

    def computeTotalMass(self):
        """
            Calculate the total mass of the PDB
        """
        cdef int atomNum
        self.totalMass = 0.0
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            self.totalMass += atom.mass

    def printAtoms(self):
        """
            Print Atom information for all the atoms in the PDB
        """
        cdef Atom atom
        for atom in self.atoms.values():
            print atom  # atom.printAtom()

    def getNumberOfAtoms(self):
        """
            Get the number of Atoms in the PDB

            :returns: int -- Number of atoms in the PDB
        """
        return len(self.atoms)

    def getAtom(self, atomId):
        """
            Get an Atom in the PDB by its id

            :param atomId: Id of the Atom (in the format "atomserial:atomName:resname")
            :type atomId: basestring
            :returns: int -- Number of atoms in the PDB
            :raises: KeyError if the id is not in the PDB
        """
        return self.atoms[atomId]

    def __len__(self):
        return len(self.atomList)


    def __getitem__(self, atomId):
        return self.atoms[atomId]

    def __setitem__(self, atomId, atom):
        self.atoms[atomId] = atom

    def __delitem__(self, atomId):
        self.atoms.pop(atomId)
        self.atomList.remove(atomId)

    def __iter__(self):
        for atomId in self.atomList:
            yield self.atoms[atomId]

    def extractCOM(self):
        """
            Calculate the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if not self.totalMass:
            self.computeTotalMass()
        cdef list COM
        cdef int atomNum
        COM = [0., 0., 0.]
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            COM[0] += atom.mass * atom.x
            COM[1] += atom.mass * atom.y
            COM[2] += atom.mass * atom.z

        COM[0] /= self.totalMass
        COM[1] /= self.totalMass
        COM[2] /= self.totalMass
        self.com = COM
        return COM

    def getCOM(self):
        """
            Get the PDB's center of mass

            :returns: list -- List with the center of mass coordinates
        """
        if self.com is None:
            return self.extractCOM()
        else:
            return self.com

    def extractCentroid(self):
        """
            Calculate the PDB centroid

            :returns: List -- List with the centroid coordinates
        """
        cdef list centroid
        cdef double n
        cdef int atomNum
        centroid = [0., 0., 0.]
        for atomNum in range(len(self.atomList)):
            atom = self.atoms[self.atomList[atomNum]]
            centroid[0] += atom.x
            centroid[1] += atom.y
            centroid[2] += atom.z

        n = float(len(self.atoms))
        centroid[0] /= n
        centroid[1] /= n
        centroid[2] /= n
        self.centroid = centroid
        return centroid

    def getCentroid(self):
        """
            Get the PDB's centroid

            :returns: list -- List with the centroid coordinates
        """
        if self.centroid is None:
            return self.extractCentroid()
        else:
            return self.centroid

    def writePDB(self, basestring path):
        """
            Write the pdb contents of the file from wich the PDB object was
            created

            :param path: Path of the file where to write the pdb
            :type path: basestring
        """
        cdef object fileHandle
        with open(path, 'w', encoding="utf-8") as fileHandle:
            fileHandle.write(self.pdb)

    def countContacts(self, basestring ligandResname, int contactThresholdDistance, int ligandResnum=0, basestring ligandChain=""):
        """
            Count the number of alpha carbons that are in contact with the
            protein (i.e. less than contactThresholdDistance Amstrogms away)

            :param ligandResname: Residue name of the ligand in the PDB
            :type ligandResname: basestring
            :param contactThresholdDistance: Maximum distance at which two atoms are considered in contanct (in Angstroms)
            :type contactThresholdDistance: int
            :returns: int -- The number of alpha carbons in contact with the ligand
        """
        cdef double contactThresholdDistance2,dist2
        contactThresholdDistance2= contactThresholdDistance**2

        cdef PDB ligandPDB, alphaCarbonsPDB

        ligandPDB = PDB()
        ligandPDB.initialise(self.pdb, resname=ligandResname, resnum=ligandResnum, chain=ligandChain, heavyAtoms=True)

        alphaCarbonsPDB = PDB()
        alphaCarbonsPDB.initialise(self.pdb, type=self._typeProtein,
                                   atomname="CA")
        # count contacts
        cdef set contacts = set([])
        cdef int rowind, colind
        cdef basestring proteinAtomId
        cdef Atom ligandAtom, proteinAtom
        for rowind in range(len(ligandPDB.atomList)):
        # can be optimised with cell list
            ligandAtom = ligandPDB.atoms[ligandPDB.atomList[rowind]]
            for colind in range(len(alphaCarbonsPDB.atomList)):
                proteinAtomId = alphaCarbonsPDB.atomList[colind]
                proteinAtom = alphaCarbonsPDB.atoms[proteinAtomId]
                dist2 = ligandAtom.squaredDistance(proteinAtom)
                if dist2 < contactThresholdDistance2:
                    contacts.update([proteinAtomId])

        return len(contacts)


def computeCOMDifference(PDB1, PDB2):
    """
        Compute the difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The distance between the centers of mass between two PDB
    """
    return np.sqrt(computeCOMSquaredDifference(PDB1, PDB2))


def computeCOMSquaredDifference(PDB PDB1, PDB PDB2):
    """
        Compute the squared difference between the center of mass of two PDB

        :param PDB1: First PDB with which the RMSD will be calculated
        :type PDB1: PDB
        :param PDB2: First PDB with which the RMSD will be calculated
        :type PDB2: PDB
        :returns: float -- The squared distance between the centers of mass between two PDB
    """
    cdef list COM1, COM2
    COM1 = PDB1.getCOM()
    COM2 = PDB2.getCOM()

    dx = COM1[0] - COM2[0]
    dy = COM1[1] - COM2[1]
    dz = COM1[2] - COM2[2]

    return dx*dx + dy*dy + dz*dz


def computeSquaredCentroidDifference(PDB PDB1, PDB PDB2):
    """
        Compute the centroid squared difference between two PDBs

        :param PDB1: First PDB
        :type PDB1: PDB
        :param PDB2: Second PDB
        :type PDB2: PDB
        :returns: float -- The squared centroid distance between two PDB
    """
    cdef list centroid1, centroid2
    centroid1 = PDB1.getCentroid()
    centroid2 = PDB2.getCentroid()

    dx = centroid1[0] - centroid2[0]
    dy = centroid1[1] - centroid2[1]
    dz = centroid1[2] - centroid2[2]

    return dx*dx + dy*dy + dz*dz

def readPDB(pdbfile):
    """
        Helper function, parses a string with PDB content or the path of a pdb file into a string

        :param pdbfile: A string with PDB content or the path of a pdb file
        :type pdbfile: basestring
        :returns: basestring -- A string with PDB content
    """
    try:
        return open(pdbfile, "rt").read()
    except IOError:
        return pdbfile
