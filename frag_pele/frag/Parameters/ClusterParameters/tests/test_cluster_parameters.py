# Python Imports

# Third-Party Imports

# Project Imports
import unittest

from frag_pele.frag.Parameters.ClusterParameters.cluster_parameters import ClusterParameters


class TestClusterParameters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cluster_parameters = ClusterParameters(0.1, 0.2, 0.3, "test_condition", "test_metricweight", 4,
                                                   "test_pdbout", ['test_banned_list'], 5)

    def test_cluster_parameters(self):
        self.assertEqual(self.cluster_parameters.distance_contact, 0.1)
        self.assertEqual(self.cluster_parameters.cluster_threshold, 0.2)
        self.assertEqual(self.cluster_parameters.epsilon, 0.3)
        self.assertEqual(self.cluster_parameters.condition, "test_condition")
        self.assertEqual(self.cluster_parameters.metric_weights, "test_metricweight")
        self.assertEqual(self.cluster_parameters.number_clusters, 4)
        self.assertEqual(self.cluster_parameters.pdnout, "test_pdbout")
        self.assertEqual(self.cluster_parameters.banned_list[0], 'test_banned_list')
        self.assertEqual(self.cluster_parameters.limit, 5)