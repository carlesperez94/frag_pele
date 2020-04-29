# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.FragParameters.frag_structural_configuration_parameters import FragStructuralConfigurationParameters


class TestFragStructuralConfigurationParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.structural_information = FragStructuralConfigurationParameters("test_core_atom", "test_fragment_atom",
                                                                           "test_h_core", "test_h_frag", "test_c_chain",
                                                                           "test_f_chain", 1.0)

    def test_structural_configuration(self):
        self.assertEqual(self.structural_information.core_atom, "test_core_atom")
        self.assertEqual(self.structural_information.fragment_atom, "test_fragment_atom")
        self.assertEqual(self.structural_information.h_core, "test_h_core")
        self.assertEqual(self.structural_information.h_frag, "test_h_frag")
        self.assertEqual(self.structural_information.c_chain, "test_c_chain")
        self.assertEqual(self.structural_information.f_chain, "test_f_chain")
        self.assertEqual(self.structural_information.threshold_clash, 1.0)
