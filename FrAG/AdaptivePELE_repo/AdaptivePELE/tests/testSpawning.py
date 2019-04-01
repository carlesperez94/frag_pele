from __future__ import absolute_import, division, print_function, unicode_literals
import unittest
import numpy as np
import AdaptivePELE.spawning.spawning as spawning
from AdaptivePELE.clustering import clustering
from AdaptivePELE.spawning import densitycalculator


class TestSpawningCalculator(unittest.TestCase):
    def testDivideTrajAccordingToWeights(self):
        spawningCalculator = spawning.SpawningCalculator()
        weights = [0.5, 0.2, 0.2, 0.1]
        trajToDistribute1 = 12
        degeneracy = spawningCalculator.divideTrajAccordingToWeights(weights, trajToDistribute1)

        golden = [6, 2, 3, 1]
        self.assertEqual(degeneracy, golden)

    def testDivideInverselyProportionalToArray(self):
        spawningCalculator = spawning.SpawningCalculator()
        weights = np.array([0.5, 0.2, 0.2, 0.1])
        trajToDistribute = 12
        degeneracy = spawningCalculator.divideInverselyProportionalToArray(weights, trajToDistribute)

        golden = [1, 3, 3, 5]

        self.assertEqual(degeneracy, golden)

    def testInverselyProportionalToPopulationCalculator(self):
        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        inverselyProp = spawning.InverselyProportionalToPopulationCalculator()
        clusteringParams = None
        trajs = 10
        degeneracy = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
        golden = [1, 2, 2, 5]

        self.assertEqual(degeneracy, golden)

    def testInverselyProportionalToPopulationCalculatorWithNullDensity(self):
        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        densityCalculator = densitycalculator.NullDensityCalculator()

        inverselyProp = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
        clusteringParams = None
        trajs = 10
        degeneracy = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
        golden = [1, 2, 2, 5]

        self.assertEqual(degeneracy, golden)

    def testInverselyProportionalToPopulationCalculatorWithDensity(self):
        # test 3
        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        contacts = [4, 3, 2, 1]
        for size, ncontacts in zip(sizes, contacts):
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            cluster.contacts = ncontacts
            clusters.addCluster(cluster)

        clusteringBlock = {
            "density": {
                "type": "heaviside",
                "params": {
                    "values": [6, 2, 3, 1],
                    "conditions": [3, 2, 1]
                }
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(clusteringBlock)

        inverselyProp = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
        clusteringParams = None
        trajs = 8
        degeneracy = inverselyProp.calculate(clusters.clusters, trajs, clusteringParams)
        golden = [2, 2, 2, 2]

        self.assertEqual(degeneracy, golden)

    def testEpsilonCalculator(self):
        epsilon = spawning.EpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.metricWeights = "linear"
        params.nclusters = 100

        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, -1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            clusters.addCluster(cluster)

        trajs = 20
        degeneracy = epsilon.calculate(clusters.clusters, trajs, params)
        golden = np.array([7, 4, 4, 5])
        np.testing.assert_array_equal(degeneracy, golden)

    def testEpsilonCalculatorWithDensity(self):
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.metricWeights = "linear"
        params.nclusters = 100

        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, 0]
        contacts = [4, 3, 2, 1]
        for size, energy, ncontacts in zip(sizes, energies, contacts):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            cluster.contacts = ncontacts
            clusters.addCluster(cluster)

        clusteringBlock = {
            "density": {
                "type": "heaviside",
                "params": {
                    "values": [6, 2, 3, 1],
                    "conditions": [3, 2, 1]
                }
            }
        }
        densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityCalculatorBuilder.build(clusteringBlock)

        epsilon = spawning.EpsilonDegeneracyCalculator(densityCalculator)
        trajs = 16
        degeneracy = epsilon.calculate(clusters.clusters, trajs, params)
        golden = np.array([6, 4, 4, 2])
        np.testing.assert_array_equal(degeneracy, golden)

    def testEpsilonCalculatorwithMaxValues(self):
        epsilon = spawning.EpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.metricWeights = "linear"
        params.nclusters = 100
        params.condition = "max"

        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, -1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            clusters.addCluster(cluster)

        trajs = 20
        degeneracy = epsilon.calculate(clusters.clusters, trajs, params)
        golden = np.array([1, 5, 5, 9])
        np.testing.assert_array_equal(degeneracy, golden)

    def testSameWeightDegeneracyCalculator(self):
        sameWeightDegCalculator = spawning.SameWeightDegeneracyCalculator()
        params = None

        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        for size in sizes:
            cluster = clustering.Cluster(None, None, None, None)
            cluster.elements = size
            clusters.addCluster(cluster)

        trajs = 10
        degeneracy = sameWeightDegCalculator.calculate(clusters.clusters, trajs, params)
        golden = [1, 1, 1, 1]

        self.assertEqual(degeneracy, golden)

    def testVariableEpsilonCalculator(self):
        variable_epsilon = spawning.VariableEpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.varEpsilonType = "linearVariation"
        params.maxEpsilon = 0.75
        params.minEpsilon = 0.5
        params.variationWindow = 8
        params.maxEpsilonWindow = 2
        params.metricWeights = "linear"
        params_test = {}
        params.period = params_test.get("period", params.variationWindow)
        params.nclusters = 100
        params.period += np.sign(np.abs(params.variationWindow-params.period))
        rateVariation = 0.25/3
        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, -1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            clusters.addCluster(cluster)

        trajs = 20

        degeneracy6 = variable_epsilon.calculate(clusters.clusters, trajs, params, 0)
        golden6 = np.array([7, 4, 4, 5])
        np.testing.assert_array_equal(degeneracy6, golden6)
        self.assertAlmostEqual(params.epsilon, params.minEpsilon)
        degeneracy7 = variable_epsilon.calculate(clusters.clusters, trajs, params, 1)
        # TODO: check degeneracy after next steps
        self.assertAlmostEqual(params.epsilon, params.minEpsilon+rateVariation)
        degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 2)
        self.assertAlmostEqual(params.epsilon, params.minEpsilon+2*rateVariation)
        degeneracy9 = variable_epsilon.calculate(clusters.clusters, trajs, params, 9)
        self.assertAlmostEqual(params.epsilon, params.minEpsilon)

    def testPeriodicVariableEpsilonCalculator(self):
        variable_epsilon = spawning.VariableEpsilonDegeneracyCalculator()
        params = spawning.SpawningParams()
        params.epsilon = 0.5
        params.varEpsilonType = "linearVariation"
        params.maxEpsilon = 0.75
        params.minEpsilon = 0.5
        params.variationWindow = 20
        params.maxEpsilonWindow = 1
        params.metricWeights = "linear"
        params.nclusters = 100
        params_test = {"period": 8}
        params.period = params_test.get("period", params.variationWindow)
        params.period += np.sign(np.abs(params.variationWindow-params.period))
        rateVariation = (params.maxEpsilon-params.minEpsilon)/3
        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, -1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            clusters.addCluster(cluster)

        trajs = 20

        degeneracy6 = variable_epsilon.calculate(clusters.clusters, trajs, params, 0)
        golden6 = np.array([7, 4, 4, 5])
        np.testing.assert_array_equal(degeneracy6, golden6)
        self.assertAlmostEqual(params.epsilon, params.minEpsilon)
        for i in range(1, params.variationWindow):
            degeneracy7 = variable_epsilon.calculate(clusters.clusters, trajs, params, i)
        # TODO: check degeneracy after next steps
        # self.assertAlmostEqual(params.epsilon, params.minEpsilon+rateVariation)
        # degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 2)
        # self.assertAlmostEqual(params.epsilon, params.maxEpsilon)
        # degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 3)
        # self.assertAlmostEqual(params.epsilon,params.maxEpsilon-rateVariation/2)
        # degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 4)
        # self.assertAlmostEqual(params.epsilon, params.minEpsilon)
        # degeneracy8 = variable_epsilon.calculate(clusters.clusters, trajs, params, 5)
        # self.assertAlmostEqual(params.epsilon, params.minEpsilon)
        # degeneracy9 = variable_epsilon.calculate(clusters.clusters, trajs, params, 9)
        # self.assertAlmostEqual(params.epsilon, params.minEpsilon)

    def testUCBCalculator(self):
        UCB = spawning.UCBCalculator()
        params = spawning.SpawningParams()
        params.alpha = 8.0

        clusters = clustering.Clusters()
        sizes = [6, 2, 3, 1]
        energies = [-4, -2, -2, -1]
        for size, energy in zip(sizes, energies):
            cluster = clustering.Cluster(None, None, None, None, metricCol=0)
            cluster.elements = size
            cluster.metrics = [energy]
            clusters.addCluster(cluster)

        trajs = 20
        degeneracy = UCB.calculate(clusters.clusters, trajs, params)
        golden = np.array([3, 5, 3, 9])
        np.testing.assert_array_equal(degeneracy, golden)


def main():
    return unittest.main(exit=False)

if __name__ == '__main__':
    main()
