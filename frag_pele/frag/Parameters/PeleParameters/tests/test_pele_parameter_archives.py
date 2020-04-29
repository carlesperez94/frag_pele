# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.PeleParameters.pele_parameter_archives import PeleParameterArchives


class TestPeleParameterArchives(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pele_parameter_archives = PeleParameterArchives("test_control_file", "test_resfold", "test_report",
                                                            "test_traject")

    def test_parameter_archives(self):
        self.assertEqual(self.pele_parameter_archives.control_file, "test_control_file")
        self.assertEqual(self.pele_parameter_archives.resfold, "test_resfold")
        self.assertEqual(self.pele_parameter_archives.report, "test_report")
        self.assertEqual(self.pele_parameter_archives.traject, "test_traject")
