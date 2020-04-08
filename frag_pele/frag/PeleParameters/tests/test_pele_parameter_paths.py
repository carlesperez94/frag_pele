# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.PeleParameters.pele_parameter_paths import PeleParameterPaths


class TestPathHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pele_parameters_paths = PeleParameterPaths("test_pele_dir", "test_pele_license", "test_data",
                                                       "test_documents")

    def test_path_handler(self):
        self.assertEqual(self.pele_parameters_paths.pele_dir, "test_pele_dir")
        self.assertEqual(self.pele_parameters_paths.pele_license, "test_pele_license")
        self.assertEqual(self.pele_parameters_paths.data, "test_data")
        self.assertEqual(self.pele_parameters_paths.documents, "test_documents")
