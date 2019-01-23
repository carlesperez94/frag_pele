from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.clustering import clustering
import unittest


class clusteringTest(unittest.TestCase):
    def testCluster(self):
        # preparation
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringParams = {"type": "rmsd",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8}}
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", 3)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)

        # assertion
        allClusters = clusteringInstance.clusters.clusters
        goldenNumberOfClusters = 2
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].getMetric(), goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].getMetric(), goldenEnergyCluster2, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)

    def test_cluster_accumulative(self):
        # preparation
        clusteringParams = {"type": "contactMap",
                            "params": {"ligandResname": "AIN",
                                       "contactThresholdDistance": 8,
                                       "similarityEvaluator": "correlation"},
                            "thresholdCalculator": {
                                "type": "constant",
                                "params": {
                                    "value": 0.15
                                }
                            }}
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringInstance = clusteringBuilder.buildClustering(clusteringParams,
                                                               "ain_report", 3)

        trajNames = ["tests/data/aspirin_data/traj*"]

        # function to test
        clusteringInstance.cluster(trajNames)
        # assertion
        allClusters = clusteringInstance.clusters.clusters
        # goldenNumberOfClusters = 3
        # goldenEnergyCluster1 = -8421.5
        # goldenEnergyCluster2 = -8424.8
        # goldenEnergyCluster3 = -8453.29
        # goldenElementsCluster1 = 1
        # goldenElementsCluster2 = 1
        # goldenElementsCluster3 = 1

        goldenNumberOfClusters = 2
        goldenEnergyCluster1 = -8424.8
        goldenEnergyCluster2 = -8453.29
        goldenElementsCluster1 = 2
        goldenElementsCluster2 = 1

        self.assertEqual(len(allClusters), goldenNumberOfClusters)
        self.assertAlmostEqual(allClusters[0].getMetric(), goldenEnergyCluster1, 2)
        self.assertAlmostEqual(allClusters[1].getMetric(), goldenEnergyCluster2, 2)
        # self.assertAlmostEqual(allClusters[2].getMetric(), goldenEnergyCluster3, 2)
        self.assertEqual(allClusters[0].elements, goldenElementsCluster1)
        self.assertEqual(allClusters[1].elements, goldenElementsCluster2)
        # self.assertEqual(allClusters[2].elements, goldenElementsCluster3)

    def test_CMInnerLimit(self):
        CMEvaluator = clustering.CMClusteringEvaluator(None, None)
        cluster1_8 = clustering.Cluster(None, contacts=1, contactThreshold=8)
        cluster2_8 = clustering.Cluster(None, contacts=2.2, contactThreshold=8)
        cluster2_6 = clustering.Cluster(None, contacts=2.2, contactThreshold=6)
        cluster1_6 = clustering.Cluster(None, contacts=1.0, contactThreshold=6)
        cluster1_4 = clustering.Cluster(None, contacts=1.0, contactThreshold=4)
        cluster2_4 = clustering.Cluster(None, contacts=2.2, contactThreshold=4)
        # Reduced slope to 8
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_8), 12)
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_8), 4)
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_6), 12)
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_6), 4)
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_4), 12)
        self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_4), 4)
        # Slope 14
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_8), 18)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_8), 4)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_6), 18)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_6), 4)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_4), 18)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_4), 4)
        # Now the limits are always calculated with threshold 8
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_6), 10)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_6), 4)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster1_4), 10)
        # self.assertAlmostEqual(CMEvaluator.getInnerLimit(cluster2_4), 4)
