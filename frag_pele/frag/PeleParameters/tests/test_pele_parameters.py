# Python Imports
import unittest
from unittest.mock import Mock

# Third-Party Imports

# Project Imports
from frag_pele.frag.PeleParameters.pele_parameters import PeleParameters


class TestPeleParameters(unittest.TestCase):

    def test_pele_parameters(self):
        mock = Mock()
        pele_parameters = PeleParameters(mock, mock, mock)

        self.assertEqual(len(pele_parameters.__dict__), 3)
        self.assertTrue(hasattr(pele_parameters, "pele_params_path"))
        self.assertTrue(hasattr(pele_parameters, "pele_params_archives"))
        self.assertTrue(hasattr(pele_parameters, "pele_params_sim_values"))
