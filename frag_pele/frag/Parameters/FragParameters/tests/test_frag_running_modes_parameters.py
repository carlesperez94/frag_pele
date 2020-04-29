# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.FragParameters.frag_running_modes_parameters import FragRunningModesParameters


class TestFragRunningModesParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.running_modes = FragRunningModesParameters(True, True)

    def test_running_modes(self):
        self.assertEqual(self.running_modes.no_check, True)
        self.assertEqual(self.running_modes.restart, True)
