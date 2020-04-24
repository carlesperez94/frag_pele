# Python Imports
import unittest

# Third-Party Imports

# Project Imports
from frag_pele.frag.FragParameters.frag_protocol_parameters import FragProtocolParameters


class TestFragProtocolParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.protocol_parameters = FragProtocolParameters(True, True)

    def test_protocol_parameters(self):
        self.assertEqual(self.protocol_parameters.only_grow, True)
        self.assertEqual(self.protocol_parameters.only_prepare, True)
