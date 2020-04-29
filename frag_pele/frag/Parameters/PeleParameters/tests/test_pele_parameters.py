# Python Imports
import unittest
from unittest.mock import Mock, patch

# Third-Party Imports

# Project Imports
from frag_pele.frag.PeleParameters.pele_parameter_archives import PeleParameterArchives
from frag_pele.frag.PeleParameters.pele_parameter_paths import PeleParameterPaths
from frag_pele.frag.PeleParameters.pele_parameter_sim_values import PeleParameterSimulationValues
from frag_pele.frag.PeleParameters.pele_parameters import PeleParameters


class TestPeleParameters(unittest.TestCase):

    def test_pele_parameters(self):
        mock = Mock()
        pele_parameters = PeleParameters(mock, mock, mock)

        self.assertEqual(len(pele_parameters.__dict__), 3)
        self.assertTrue(hasattr(pele_parameters, "pele_params_path"))
        self.assertTrue(hasattr(pele_parameters, "pele_params_archives"))
        self.assertTrue(hasattr(pele_parameters, "pele_params_sim_values"))

    def test_extract_parameters(self):
        mock_pele_params_path = Mock(spec=PeleParameterPaths)
        mock_pele_params_archives = Mock(spec=PeleParameterArchives)
        mock_pele_param_sim = Mock(spec=PeleParameterSimulationValues)

        pele_parameters = PeleParameters(mock_pele_params_path, mock_pele_params_archives, mock_pele_param_sim)

        return_pele_path, return_pele_archives, return_pele_sim = pele_parameters.extract_parameters()

        self.assertIsInstance(return_pele_path, PeleParameterPaths)
        self.assertIsInstance(return_pele_archives, PeleParameterArchives)
        self.assertIsInstance(return_pele_sim, PeleParameterSimulationValues)
