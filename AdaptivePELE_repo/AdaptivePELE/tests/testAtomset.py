from __future__ import absolute_import, division, print_function, unicode_literals
from io import open
import os
import unittest
import mdtraj
import numpy as np
import AdaptivePELE.atomset.atomset as atomset
from AdaptivePELE.atomset import RMSDCalculator
from AdaptivePELE.atomset import SymmetryContactMapEvaluator as sym
from AdaptivePELE.clustering import clustering
from AdaptivePELE.utilities import utilities


class atomsetTest(unittest.TestCase):
    """ For the moment the tests include loading from file and string, resname
    and atomname selection. It uses a toy pdb of only 5 lines located in
    tests/data
    """
    def testPDB_from_file(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom1 = atomset.Atom("ATOM      1  N   ASN A   1       7.920  22.268   9.257 1.00 15.18           N1+")
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        atom4 = atomset.Atom("ATOM      4  O   ASN A   1       7.030  20.308   7.700 1.00 16.13           O  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_from_str(self):
        # preparation
        with open("tests/data/pdb_test.pdb", "rt") as pdbfile:
            pdbstring = pdbfile.read()
        pdb = atomset.PDB()

        # function to test
        pdb.initialise(pdbstring)

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom1 = atomset.Atom("ATOM      1  N   ASN A   1       7.920  22.268   9.257 1.00 15.18           N1+")
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        atom4 = atomset.Atom("ATOM      4  O   ASN A   1       7.030  20.308   7.700 1.00 16.13           O  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_resname(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb", resname="CYS")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom5 = atomset.Atom("ATOM      5  CB  CYS A   1       8.108  20.445  11.030 1.00 16.53           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom5.id: atom5}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_atomname(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb", atomname="CA")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS A   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN A   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   1       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   1       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom2 = atomset.Atom("ATOM      2  CA  CYS A   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN A   1       7.870  19.937   8.524 1.00 16.63           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom3.id: atom3}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_protein(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb",
                       type="PROTEIN")

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  ")
        goldenAtomsDict = {atom1.id: atom1}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_hetero(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb", type="HETERO")

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  ")
        atom2 = atomset.Atom("HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_type_heavyAtoms(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test_ligand.pdb", heavyAtoms=False)

        # assertion
        pdbContent = "MODEL        1\n\
ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  \n\
ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  \n\
TER\n\
HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  \n\
TER\n\
HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+\n\
TER\n\
ENDMDL\n\
END   \n"

        atom1 = atomset.Atom("ATOM   1717  H   ASN A 119      25.915   9.925  -7.236  1.00 31.61           H  ")
        atom2 = atomset.Atom("ATOM   1718  CA  ASN A 119      27.159  10.509  -6.736  1.00 33.83           C  ")
        atom3 = atomset.Atom("HETATM 1733  O1  AIN L   1      13.907  16.130   0.624  0.50 28.52           O  ")
        atom4 = atomset.Atom("HETATM 1753 CA    CA B   1      16.636  15.477   0.293  1.00  3.39          Ca2+")
        goldenAtomsDict = {atom1.id: atom1, atom2.id: atom2, atom3.id: atom3, atom4.id: atom4}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_resnum(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test2.pdb", resnum=2)

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS B   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN B   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   2       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   2       8.108  20.445  11.030  1.00 16.53           C  \n"
        atom2 = atomset.Atom("ATOM      4  O   ASN A   2       7.030  20.308   7.700  1.00 16.13           O  ")
        atom3 = atomset.Atom("ATOM      5  CB  CYS A   2       8.108  20.445  11.030  1.00 16.53           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom3.id: atom3}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_sel_chain(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test2.pdb", chain="B")

        # assertion
        pdbContent = "ATOM      1  N   ASN A   1       7.920  22.268   9.257  1.00 15.18           N1+\n\
ATOM      2  CA  CYS B   1       8.394  20.916   9.575  1.00 16.24           C  \n\
ATOM      3  CA  ASN B   1       7.870  19.937   8.524  1.00 16.63           C  \n\
ATOM      4  O   ASN A   2       7.030  20.308   7.700  1.00 16.13           O  \n\
ATOM      5  CB  CYS A   2       8.108  20.445  11.030  1.00 16.53           C  \n"

        atom2 = atomset.Atom("ATOM      2  CA  CYS B   1       8.394  20.916   9.575 1.00 16.24           C  ")
        atom3 = atomset.Atom("ATOM      3  CA  ASN B   1       7.870  19.937   8.524 1.00 16.63           C  ")
        goldenAtomsDict = {atom2.id: atom2, atom3.id: atom3}

        self.assertEqual(pdb.pdb, pdbContent)
        self.assertEqual(pdb.atoms, goldenAtomsDict)

    def testPDB_COM(self):
        # preparation
        pdb = atomset.PDB()

        # function to test
        pdb.initialise("tests/data/pdb_test.pdb")

        # assertion
        total_mass = 66.0382
        COM_array = np.array([516.1336264, 1373.048894, 602.7150822])
        pdb.extractCOM()
        self.assertNotAlmostEqual(pdb.totalMass, 0)
        COM = COM_array/pdb.totalMass

        np.testing.assert_array_almost_equal(COM, pdb.getCOM())
        self.assertAlmostEqual(total_mass, pdb.totalMass)

    def testPDB_write(self):
        # preparation
        pdb = atomset.PDB()
        pdb.initialise("tests/data/pdb_test.pdb")

        # function to test
        pdb.writePDB("tests/data/pdb_test_write.pdb")

        # assertion
        pdbtestfile = open("tests/data/pdb_test_write.pdb", "r")
        pdbtestsstr = pdbtestfile.read()
        pdbtestfile.close()
        self.assertEqual(pdb.pdb, pdbtestsstr)

    def testPDB_RMSD(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/ain_native_fixed.pdb", resname='AIN')
        pdb_traj = atomset.PDB()
        pdb_traj.initialise("tests/data/ain_trajectory.pdb", resname='AIN')
        RMSDCalc = RMSDCalculator.RMSDCalculator()

        # function to test
        RMSD = RMSDCalc.computeRMSD(pdb_native, pdb_traj)
        golden_RMSD = 3.928617
        self.assertAlmostEqual(RMSD, golden_RMSD, 5)

    def testPDB_RMSD_symmetries(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/ain_native_fixed.pdb", resname='AIN')
        pdb_traj = atomset.PDB()
        pdb_traj.initialise("tests/data/ain_trajectory.pdb", resname='AIN')
        symDict = [{"1733:O1:AIN": "1735:O2:AIN"}]
        RMSDCalc = RMSDCalculator.RMSDCalculator(symDict)
        # function to test
        RMSD = RMSDCalc.computeRMSD(pdb_native, pdb_traj)
        reverseRMSD = RMSDCalc.computeRMSD(pdb_traj, pdb_native)
        golden_RMSD = 3.860743
        self.assertAlmostEqual(RMSD, reverseRMSD, 5)
        self.assertAlmostEqual(RMSD, golden_RMSD, 5)

    def test_combination_symmetries(self):
        # preparation
        pdb_0 = atomset.PDB()
        pdb_0.initialise("tests/data/symmetries/cluster_0.pdb", resname='AEN')
        pdb_1 = atomset.PDB()
        pdb_1.initialise("tests/data/symmetries/cluster_1.pdb", resname='AEN')
        pdb_2 = atomset.PDB()
        pdb_2.initialise("tests/data/symmetries/cluster_2.pdb", resname='AEN')
        symmetries3PTB = [{"3225:C3:AEN": "3227:C5:AEN", "3224:C2:AEN": "3228:C6:AEN"},
                          {"3230:N1:AEN": "3231:N2:AEN"}]
        RMSDCalc = RMSDCalculator.RMSDCalculator(symmetries3PTB)
        # funtion to test
        RMSD02 = RMSDCalc.computeRMSD(pdb_0, pdb_2)
        RMSD20 = RMSDCalc.computeRMSD(pdb_2, pdb_0)
        RMSD01 = RMSDCalc.computeRMSD(pdb_0, pdb_1)
        RMSD10 = RMSDCalc.computeRMSD(pdb_1, pdb_0)
        RMSD21 = RMSDCalc.computeRMSD(pdb_2, pdb_1)
        RMSD12 = RMSDCalc.computeRMSD(pdb_1, pdb_2)
        self.assertEqual(RMSD01, RMSD10)
        self.assertEqual(RMSD02, RMSD20)
        self.assertEqual(RMSD21, RMSD12)

    def testPDB_contacts(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/native_ain.pdb")

        # function to test
        contacts = pdb_native.countContacts("AIN", 8)
        golden_contacts = 19
        self.assertEqual(contacts, golden_contacts)

    def testPDB_contactmap(self):
        # preparation
        pdb_native = atomset.PDB()
        pdb_native.initialise("tests/data/pdb_test_contact.pdb")
        symmetryEvaluator = sym.SymmetryContactMapEvaluator([])

        # function to test
        contact_map, contacts = symmetryEvaluator.createContactMap(pdb_native,
                                                                   "AIN", 8)
        golden_contact_map = np.array([[1, 0, 0, 0], [0, 1, 1, 1]])
        golden_contacts = pdb_native.countContacts("AIN", 8)
        np.testing.assert_array_equal(contact_map, golden_contact_map)
        self.assertEqual(golden_contacts, contacts)

    def test_contactMapContacts(self):
        # preparation
        pdb_1 = atomset.PDB()
        pdb_1.initialise("tests/data/pdb_test_contact.pdb", resname='AIN')
        symmetryEvaluator = sym.SymmetryContactMapEvaluator([])

        # function to test
        _, contacts = symmetryEvaluator.createContactMap(pdb_1, "AIN", 16)
        golden_contacts = pdb_1.countContacts("AIN", 8)
        self.assertEqual(golden_contacts, contacts)

    def test_symmetryContactMapJaccard(self):
        pdb_1 = atomset.PDB()
        pdb_1.initialise("tests/data/symmetries/cluster_1.pdb", resname='AEN')
        pdb_1_sym = atomset.PDB()
        pdb_1_sym.initialise("tests/data/symmetries/cluster_1_sym.pdb",
                             resname='AEN')
        symmetries3PTB = [{"3230:N1:AEN": "3231:N2:AEN"}]
        symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries3PTB)
        symmetryEvaluatorEmpty = sym.SymmetryContactMapEvaluator()

        contactMap1, contacts1 = symmetryEvaluator.buildContactMap(pdb_1, 'AEN', 16)
        cluster = clustering.Cluster(pdb_1, contactMap=contactMap1)
        contactMap1Sym, contactsSym = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)
        contactMapNoSym, _ = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)

        goldenJaccard = 0.0
        Jaccard = symmetryEvaluator.evaluateJaccard(contactMap1Sym, cluster.contactMap)
        JaccardNosym = symmetryEvaluatorEmpty.evaluateJaccard(contactMapNoSym, cluster.contactMap)

        self.assertEqual(contacts1, contactsSym)
        self.assertAlmostEqual(goldenJaccard, Jaccard)
        self.assertNotAlmostEqual(Jaccard, JaccardNosym)

    def test_symmetryContactMapCorrelation(self):
        pdb_1 = atomset.PDB()
        pdb_1.initialise("tests/data/symmetries/cluster_1.pdb", resname='AEN')
        pdb_1_sym = atomset.PDB()
        pdb_1_sym.initialise("tests/data/symmetries/cluster_1_sym.pdb",
                             resname='AEN')
        symmetries3PTB = [{"3230:N1:AEN": "3231:N2:AEN"}]
        symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries3PTB)
        symmetryEvaluatorEmpty = sym.SymmetryContactMapEvaluator()

        contactMap1, contacts1 = symmetryEvaluator.buildContactMap(pdb_1, 'AEN', 16)
        cluster = clustering.Cluster(pdb_1, contactMap=contactMap1)
        contactMap1Sym, contactsSym = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)
        contactMapNoSym, _ = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)

        goldenCorrelation = 0.0
        correlationSym = symmetryEvaluator.evaluateCorrelation(contactMap1Sym, cluster.contactMap)
        correlationNosym = symmetryEvaluatorEmpty.evaluateCorrelation(contactMapNoSym, cluster.contactMap)

        self.assertEqual(contacts1, contactsSym)
        self.assertAlmostEqual(goldenCorrelation, correlationSym)
        self.assertNotAlmostEqual(correlationSym, correlationNosym)

    def test_symmetryContactMapDifference(self):
        pdb_1 = atomset.PDB()
        pdb_1.initialise("tests/data/symmetries/cluster_1.pdb", resname='AEN')
        pdb_1_sym = atomset.PDB()
        pdb_1_sym.initialise("tests/data/symmetries/cluster_1_sym.pdb",
                             resname='AEN')
        symmetries3PTB = [{"3230:N1:AEN": "3231:N2:AEN"}]
        symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries3PTB)
        symmetryEvaluatorEmpty = sym.SymmetryContactMapEvaluator()

        contactMap1, contacts1 = symmetryEvaluator.buildContactMap(pdb_1, 'AEN', 16)
        cluster = clustering.Cluster(pdb_1, contactMap=contactMap1)
        contactMap1Sym, contactsSym = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)
        contactMapNoSym, _ = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)

        goldenDifference = 0.0
        DifferenceSym = symmetryEvaluator.evaluateDifferenceDistance(contactMap1Sym, cluster.contactMap)
        DifferenceNosym = symmetryEvaluatorEmpty.evaluateDifferenceDistance(contactMapNoSym, cluster.contactMap)

        self.assertEqual(contacts1, contactsSym)
        self.assertAlmostEqual(goldenDifference, DifferenceSym)
        self.assertNotAlmostEqual(DifferenceSym, DifferenceNosym)

    def test_PDB_interface(self):
        pdb = atomset.PDB()
        pdb.initialise("tests/data/symmetries/cluster_1.pdb", resname='AEN')
        self.assertEqual(len(pdb), 9)
        atomList = [atom for atom in pdb]
        atoms = [pdb.getAtom(a) for a in pdb.atomList]
        self.assertEqual(atomList, atoms)
        atomId = pdb.atomList[0]
        atom = pdb[atomId]
        self.assertEqual(atom, pdb.getAtom(atomId))
        pdb[atomId] = None
        self.assertEqual(None, pdb.getAtom(atomId))

    def test_write_XTC_to_pdb(self):
        golden = "tests/data/ain_native_fixed.pdb"
        output = "xtc_to_pdb.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(xtc_obj.xyz[0], resname="AIN", topology=topology)
        top = utilities.getTopologyFile(golden)
        xtc.writePDB(output)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")
        output_pdb = atomset.PDB()
        output_pdb.initialise(output, resname="AIN")
        os.remove(output)
        self.assertEqual(golden_pdb.atoms, output_pdb.atoms)

    def testPDB_sel_resname_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")
        self.assertEqual(xtc, golden_pdb)

    def testPDB_sel_atomname_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], atomname="CA", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, atomname="CA")
        self.assertEqual(xtc, golden_pdb)

    def testPDB_sel_type_protein_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], type="PROTEIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, type="PROTEIN")
        self.assertEqual(xtc, golden_pdb)

    def testPDB_sel_type_hetero_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], type="HETERO", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, type="HETERO")
        self.assertEqual(xtc, golden_pdb)

    def testPDB_sel_type_heavyAtoms_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], heavyAtoms=False, topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, heavyAtoms=False)
        self.assertEqual(xtc, golden_pdb)

    def testPDB_sel_resnum_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resnum=2, topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resnum=2)
        self.assertEqual(xtc, golden_pdb)

    def testPDB_COM_XTC(self):
        # preparation
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")

        # assertion
        self.assertAlmostEqual(xtc.totalMass, golden_pdb.totalMass, 3)
        np.testing.assert_array_almost_equal(xtc.getCOM(), golden_pdb.getCOM(), decimal=3)

    def testPDB_RMSD_XTC(self):
        # preparation
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")

        # assertion
        RMSDCalc = RMSDCalculator.RMSDCalculator()

        # function to test
        RMSD = RMSDCalc.computeRMSD(golden_pdb, xtc)
        golden_RMSD = 0.0000
        self.assertAlmostEqual(RMSD, golden_RMSD, 2)

    def testPDB_RMSD_symmetries_XTC(self):
        # preparation
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")
        symDict = [{"1733:O1:AIN": "1735:O2:AIN"}]
        RMSDCalc = RMSDCalculator.RMSDCalculator(symDict)
        # function to test
        RMSD = RMSDCalc.computeRMSD(xtc, golden_pdb)
        reverseRMSD = RMSDCalc.computeRMSD(golden_pdb, xtc)
        golden_RMSD = 0.00000
        self.assertAlmostEqual(RMSD, reverseRMSD, 2)
        self.assertAlmostEqual(RMSD, golden_RMSD, 2)

    def testPDB_contacts_XTC(self):
        # preparation
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")

        # function to test
        contacts = golden_pdb.countContacts("AIN", 8)
        contacts_xtc = xtc.countContacts("AIN", 8)
        self.assertEqual(contacts, contacts_xtc)

    def testPDB_contactmap_XTC(self):
        # preparation
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")
        symmetryEvaluator = sym.SymmetryContactMapEvaluator([])

        # function to test
        contact_map, contacts = symmetryEvaluator.createContactMap(golden_pdb,
                                                                   "AIN", 8)
        symmetryEvaluator_xtc = sym.SymmetryContactMapEvaluator([])
        contact_map_xtc, contacts_xtc = symmetryEvaluator_xtc.createContactMap(xtc, "AIN", 8)
        np.testing.assert_array_equal(contact_map, contact_map_xtc)
        self.assertEqual(contacts_xtc, contacts)

    def test_symmetryContactMapJaccard_XTC(self):
        xtc_obj = mdtraj.load("tests/data/symmetries/cluster_1.xtc", top="tests/data/symmetries/cluster_1.pdb")
        topology = utilities.getTopologyFile("tests/data/symmetries/cluster_1.pdb")
        pdb_1 = atomset.PDB()
        pdb_1.initialise(10*xtc_obj.xyz[0], resname='AEN', topology=topology)
        topology = utilities.getTopologyFile("tests/data/symmetries/cluster_1_sym.pdb")
        xtc_obj = mdtraj.load("tests/data/symmetries/cluster_1_sym.xtc", top="tests/data/symmetries/cluster_1_sym.pdb")
        pdb_1_sym = atomset.PDB()
        pdb_1_sym.initialise(10*xtc_obj.xyz[0], resname='AEN', topology=topology)
        symmetries3PTB = [{"3230:N1:AEN": "3231:N2:AEN"}]
        symmetryEvaluator = sym.SymmetryContactMapEvaluator(symmetries3PTB)
        symmetryEvaluatorEmpty = sym.SymmetryContactMapEvaluator()

        contactMap1, contacts1 = symmetryEvaluator.buildContactMap(pdb_1, 'AEN', 16)
        cluster = clustering.Cluster(pdb_1, contactMap=contactMap1)
        contactMap1Sym, contactsSym = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)
        contactMapNoSym, _ = symmetryEvaluator.createContactMap(pdb_1_sym, 'AEN', 16)

        goldenJaccard = 0.0
        Jaccard = symmetryEvaluator.evaluateJaccard(contactMap1Sym, cluster.contactMap)
        JaccardNosym = symmetryEvaluatorEmpty.evaluateJaccard(contactMapNoSym, cluster.contactMap)
        self.assertEqual(contacts1, contactsSym)
        self.assertAlmostEqual(goldenJaccard, Jaccard)
        self.assertNotAlmostEqual(Jaccard, JaccardNosym)

    def test_PDB_interface_XTC(self):
        golden = "tests/data/ain_native_fixed.pdb"
        topology = utilities.getTopologyFile(golden)
        xtc_obj = mdtraj.load("tests/data/ain_native_fixed.xtc", top=golden)
        xtc = atomset.PDB()
        xtc.initialise(10*xtc_obj.xyz[0], resname="AIN", topology=topology)
        golden_pdb = atomset.PDB()
        golden_pdb.initialise(golden, resname="AIN")
        self.assertEqual(len(golden_pdb), len(xtc))
        atomList = [atom for atom in golden_pdb]
        atomList_xtc = [atom for atom in xtc]
        self.assertEqual(atomList, atomList_xtc)
        atomId = xtc.atomList[0]
        atom = xtc[atomId]
        self.assertEqual(atom, xtc.getAtom(atomId))
        xtc[atomId] = None
        self.assertEqual(None, xtc.getAtom(atomId))
