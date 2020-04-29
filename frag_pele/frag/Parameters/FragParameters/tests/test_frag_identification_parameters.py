# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.FragParameters import FragIdentificationParameters


class TestFragIdentificationParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.identification_parameters = FragIdentificationParameters("test_ligand_id")

    def test_identification_parameters(self):
        self.assertEqual(self.identification_parameters.ligand_id, "test_ligand_id")
