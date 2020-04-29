# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.PeleParameters import PeleParameterSimulationValues


class TestPeleParameterSimValues(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pele_parameter_sim_values = PeleParameterSimulationValues(1, 2, 3, 4.1, 5.2, 6, 7, 8, 9.3, 10.4, 11.5, 12.6,
                                                                      13.7)

    def test_parameter_archives(self):
        self.assertEqual(self.pele_parameter_sim_values.cpus, 1)
        self.assertEqual(self.pele_parameter_sim_values.steps, 2)
        self.assertEqual(self.pele_parameter_sim_values.pele_eq_steps, 3)
        self.assertEqual(self.pele_parameter_sim_values.min_overlap, 4.1)
        self.assertEqual(self.pele_parameter_sim_values.max_overlap, 5.2)
        self.assertEqual(self.pele_parameter_sim_values.temperature, 6)
        self.assertEqual(self.pele_parameter_sim_values.seed, 7)
        self.assertEqual(self.pele_parameter_sim_values.steering, 8)
        self.assertEqual(self.pele_parameter_sim_values.translation_high, 9.3)
        self.assertEqual(self.pele_parameter_sim_values.rotation_high, 10.4)
        self.assertEqual(self.pele_parameter_sim_values.translation_low, 11.5)
        self.assertEqual(self.pele_parameter_sim_values.rotation_low, 12.6)
        self.assertEqual(self.pele_parameter_sim_values.radius_box, 13.7)
