# Python Imports
import unittest
from unittest.mock import Mock

# Third-Party Imports

# Project Imports

from frag_pele.frag.Parameters.FragParameters import FragParameters


class TestFragParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        mock = Mock()

        cls.frag_parameters = FragParameters(mock, mock, mock, mock, mock, mock)

    def test_frag_parameters(self):
        self.assertTrue(hasattr(self.frag_parameters, "structural_files"), True)
        self.assertTrue(hasattr(self.frag_parameters, "structural_configuration"), True)
        self.assertTrue(hasattr(self.frag_parameters, "identification_parameters"), True)
        self.assertTrue(hasattr(self.frag_parameters, "configuration_parameters"), True)
        self.assertTrue(hasattr(self.frag_parameters, "protocol_parameters"), True)
        self.assertTrue(hasattr(self.frag_parameters, "running_modes"), True)
