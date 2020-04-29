# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.FragParameters import FragConfigurationParameters


class TestFragConfigurationParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.configuration_parameters = FragConfigurationParameters(1, "test_criteria", "test_sampling_control")

    def test_configuration_parameters(self):
        self.assertEqual(self.configuration_parameters.growing_steps, 1)
        self.assertEqual(self.configuration_parameters.criteria, "test_criteria")
        self.assertEqual(self.configuration_parameters.sampling_control, "test_sampling_control")
