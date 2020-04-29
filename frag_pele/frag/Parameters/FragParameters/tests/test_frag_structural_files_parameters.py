# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.FragParameters import FragStructuralFilesParameters


class TestFragStructuralFilesParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.structural_files = FragStructuralFilesParameters("test_complex_pdb", "test_fragment_pdb")

    def test_structural_files(self):
        self.assertEqual(self.structural_files.fragment_pdb, "test_fragment_pdb")
        self.assertEqual(self.structural_files.complex_pdb, "test_complex_pdb")
