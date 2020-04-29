# Python Imports

# Third-Party Imports

# Project Imports
import unittest

from frag_pele.frag.PlopParameters.plop_parameters import PlopParameters


class TestPlopParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.plop_parameters = PlopParameters("test_plop_path", "test_sch_python", "test_rotamers")

    def test_plop_parameters(self):
        self.assertEqual(self.plop_parameters.plop_path, "test_plop_path")
        self.assertEqual(self.plop_parameters.sch_python, "test_sch_python")
        self.assertEqual(self.plop_parameters.rotamers, "test_rotamers")