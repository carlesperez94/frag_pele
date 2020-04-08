# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.path_handler.path_handler import PathHandler


class TestPathHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.path_handler = PathHandler("TestComplex", "TestFragment", "TestPlop", "TestSch", "TestPele", "TestLicense")

    def test_path_handler(self):
        self.assertEqual(self.path_handler.complex_pdb, "TestComplex")
        self.assertEqual(self.path_handler.fragment_pdb, "TestFragment")
        self.assertEqual(self.path_handler.plop_path, "TestPlop")
        self.assertEqual(self.path_handler.sch_python, "TestSch")
        self.assertEqual(self.path_handler.pele_dir, "TestPele")
        self.assertEqual(self.path_handler.pele_license, "TestLicense")
