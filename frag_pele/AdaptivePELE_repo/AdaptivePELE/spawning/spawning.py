from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import math
import sys
import numpy as np
import random
import scipy.optimize as optim
import os
import glob
from AdaptivePELE.constants import blockNames
from AdaptivePELE.constants import constants
from AdaptivePELE.utilities import utilities
from AdaptivePELE.spawning import spawningTypes
from AdaptivePELE.spawning import densitycalculator
from abc import abstractmethod


def reward(x, rews):
    return -(x[:, np.newaxis]*rews).sum()


def return_sign(i, m, n, r):
    """
        Helper function, creates a three-piece step function

        :param i: Value to compare
        :type i: int
        :param m: Middle value
        :type m: int
        :param n: Left extreme value
        :type n: int
        :param r: Right extreme value
        :type r: int

        :returns: int -- Three-piece sign
    """

    if i <= n-m:
        return 1
    elif i <= r:
        return 0
    else:
        return -1


def getSizes(clusters):
    """
        Get the size of the clusters

        :param clusters: Existing clusters
        :type clusters: :py:class:`.Clusters`

        :returns: np.Array -- Array containing the size of the clusters
    """
    sizes = np.zeros(len(clusters))
    for i, cluster in enumerate(clusters):
        sizes[i] = cluster.elements
    return sizes


def calculateContactsVar(deltaR, epsMax):
    """
        Calculate the variation of epsilon according to the contact ratio

        :param deltaR: Change in contact ratio
        :type deltaR: float
        :param epsMax: Maximum value of epsilon
        :type epsMax: float

        :returns: float -- Epsilon variation
    """
    if deltaR < 0.1:
        return 0
    elif deltaR > 1.0:
        return epsMax * 0.09
    else:
        return epsMax * 0.09 * deltaR


class SpawningAlgorithmBuilder:

    def build(self, spawningBlock):
        """
            Build the selected spawning calculator and spawning params objects

            :param spawningBlock: Block of the control file with the spawning
                parameters
            :type spawningBlock: dict

            :returns: :py:class:`.SpawningCalculator`,
                :py:class:`.SpawningParams` -- SpawningCalculator and
                SpawningParams objects
        """
        spawningCalculatorBuilder = SpawningBuilder()
        spawningCalculator = spawningCalculatorBuilder.buildSpawningCalculator(spawningBlock)

        spawningParams = SpawningParams()
        spawningParams.buildSpawningParameters(spawningBlock)

        return spawningCalculator, spawningParams


class SpawningBuilder:

    def buildSpawningCalculator(self, spawningBlock):
        """
            Build the selected spawning calculator object

            :param spawningBlock: Block of the control file with the spawning
                parameters
            :type spawningBlock: dict

            :returns: :py:class:`.SpawningCalculator` -- SpawningCalculator
                object
        """
        densityBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityBuilder.build(spawningBlock)

        spawningTypeString = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningTypeString == blockNames.StringSpawningTypes.sameWeight:
            spawningCalculator = SameWeightDegeneracyCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.independent:
            spawningCalculator = IndependentRunsCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.inverselyProportional:
            spawningCalculator = InverselyProportionalToPopulationCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.epsilon:
            spawningCalculator = EpsilonDegeneracyCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.fast:
            spawningCalculator = FASTDegeneracyCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.variableEpsilon:
            spawningCalculator = VariableEpsilonDegeneracyCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.UCB:
            spawningCalculator = UCBCalculator(densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.REAP:
            spawningCalculator = REAPCalculator()
        elif spawningTypeString == blockNames.StringSpawningTypes.null:
            spawningCalculator = NullSpawningCalculator()
        else:
            sys.exit("Unknown spawning type! Choices are: " + str(spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY.values()))
        return spawningCalculator


class SpawningParams:

    def __init__(self):
        self.epsilon = None
        self.temperature = None
        self.threshold = None
        self.reportFilename = None
        self.reportCol = None
        self.decrement = None
        self.varEpsilonType = None
        self.maxEpsilon = None
        self.minEpsilon = None
        self.variationWindow = None
        self.maxEpsilonWindow = None
        self.metricWeights = None
        self.alpha = None
        self.nclusters = None  # number of clusters to consider in epsilon
        self.period = None
        self.metricInd = None
        self.condition = blockNames.SpawningParams.minValue  # wether to consider min or max values in epsilon

    def buildSpawningParameters(self, spawningBlock):
        """
            Build the selected spawning params objects

            :param spawningBlock: Block of the control file with the spawning
                parameters
            :type spawningBlock: dict

            :returns: :py:class:`.SpawningParams` -- SpawningParams object
        """
        spawningParamsBlock = spawningBlock[blockNames.SpawningParams.params]
        spawningType = spawningBlock[blockNames.StringSpawningTypes.type]

        # Params general to all
        # reportFilename is now mandatory for all spawning
        self.reportFilename = spawningParamsBlock[blockNames.SpawningParams.report_filename]
        # Params specific to epsilon related spawning
        if spawningType == blockNames.StringSpawningTypes.epsilon or \
                spawningType == blockNames.StringSpawningTypes.variableEpsilon:
            self.epsilon = spawningParamsBlock[blockNames.SpawningParams.epsilon]
            self.nclusters = spawningParamsBlock.get(blockNames.SpawningParams.nclusters, 5)
            self.metricWeights = spawningParamsBlock.get(blockNames.SpawningParams.metricWeights,
                                                         blockNames.SpawningParams.linear)
            self.condition = spawningParamsBlock.get(blockNames.SpawningParams.condition,
                                                     blockNames.SpawningParams.minValue)

        if spawningType == blockNames.StringSpawningTypes.epsilon or \
                spawningType == blockNames.StringSpawningTypes.variableEpsilon or\
                spawningType == blockNames.StringSpawningTypes.fast or \
                spawningType == blockNames.StringSpawningTypes.simulatedAnnealing or \
                spawningType == blockNames.StringSpawningTypes.UCB or \
                spawningType == blockNames.StringSpawningTypes.REAP:
            self.temperature = spawningParamsBlock.get(blockNames.SpawningParams.temperature, 1000)
            # Start counting the columns by 1
            self.reportCol = spawningParamsBlock[blockNames.SpawningParams.report_col]-1

        if spawningType == blockNames.StringSpawningTypes.variableEpsilon:
            self.varEpsilonType = spawningParamsBlock[blockNames.SpawningParams.varEpsilonType]
            self.maxEpsilon = spawningParamsBlock[blockNames.SpawningParams.maxEpsilon]
            if self.varEpsilonType == blockNames.VariableEpsilonTypes.linearVariation:
                self.minEpsilon = spawningParamsBlock.get(blockNames.SpawningParams.minEpsilon, self.epsilon)
                self.variationWindow = spawningParamsBlock[blockNames.SpawningParams.variationWindow]
                self.maxEpsilonWindow = spawningParamsBlock[blockNames.SpawningParams.maxEpsilonWindow]
                self.period = spawningParamsBlock.get(blockNames.SpawningParams.period, self.variationWindow)
                self.period += np.sign(np.abs(self.variationWindow-self.period))
                # Add one epoch to the total lenght of the variation in the case of periodic variation to leave a step between variation periods

        if spawningType == blockNames.StringSpawningTypes.UCB:
            self.alpha = spawningParamsBlock.get(blockNames.SpawningParams.alpha, 8.0)

        if spawningType == blockNames.StringSpawningTypes.REAP:
            self.metricInd = spawningParamsBlock.get(blockNames.SpawningParams.metricsInd, -1)


class SpawningCalculator:
    """
        The purpose of this abstract class is to contain the behaviour of the different strategies for the spawning.
        Spawning is the way in which we split the different explorers at the begining of each epoch.
    """

    def __init__(self):
        self.type = "BaseClass"  # change for abstract attribute

    @abstractmethod
    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        pass

    @abstractmethod
    def log(self):
        """
            Log spawning information
        """
        pass

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration):
        """
            Write initial structures for the next iteration

            :param outputPathConstants: Output constants that depend on the path
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param degeneracyOfRepresentatives: List with the degeneracy of
                each cluster (number of processors that will start from that state)
            :type degeneracyOfRepresentatives: list
            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`
            :param iteration: Number of epoch
            :type iteration: int

            :returns: int, list -- number of processors, list with the
                snapshot from which the trajectories will start in the next iteration
        """
        tmpInitialStructuresTemplate = outputPathConstants.tmpInitialStructuresTemplate
        counts = 0
        procMapping = []
        for i, cluster in enumerate(clustering.clusters.clusters):
            for _ in range(int(degeneracyOfRepresentatives[i])):
                outputFilename = tmpInitialStructuresTemplate % (iteration, counts)
                print('Writing to ', outputFilename, 'cluster', i)
                procMapping.append(cluster.writeSpawningStructure(outputFilename))

                counts += 1

        print("counts & cluster centers", counts, np.where(np.array(degeneracyOfRepresentatives) > 0)[0].size)
        return counts, procMapping

    def divideTrajAccordingToWeights(self, weights, trajToDistribute):
        """
            Distribute the trajectories among the clusters according to their
            weight. Weights must be normalized (i.e. sum(weights) = 1)

            :param weights: Weight of each cluster
            :type weights: np.Array
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int

            :returns: list -- List with the number of processors allocated to
                each cluster
        """
        degeneracy = []
        for i, weight in enumerate(weights):
            degeneracy.append(int(weight*trajToDistribute))

        # divide remaining traj to distribute according to decimal part
        decimalPart = []
        decimalPart = [math.modf(weight*trajToDistribute)[0] for weight in weights]
        sortedDecimals = np.argsort(decimalPart)
        sortedDecimals = sortedDecimals[::-1]  # flip list

        leftProcessors = trajToDistribute-sum(degeneracy)
        for i in range(leftProcessors):
            degeneracy[sortedDecimals[i]] += 1

        return degeneracy

    def divideProportionalToArray(self, array, trajToDistribute):
        """
            Distribute the trajectories among the clusters according to the
            values of an array.

            :param Array: Weight of each cluster
            :type Array: np.Array
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int

            :returns: list -- List with the number of processors allocated to
                each cluster
        """
        if isinstance(array, list):
            array = np.array(array)
        weights = array/sum(array)
        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def divideInverselyProportionalToArray(self, array, trajToDistribute):
        """
            Distribute the trajectories among the clusters inversely proportional
            to the values of an array.

            :param Array: Weight of each cluster
            :type Array: np.Array
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int

            :returns: list -- List with the number of processors allocated to
                each cluster
        """
        if isinstance(array, list):
            array = np.array(array)
        weights = 1./array

        # Handle Nan cases
        weights[weights == np.inf] = 0

        # Handle all Nan cases
        if weights.any():
            weights /= sum(weights)
        else:
            weights[:] = 1./weights.shape[0]

        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def getMetrics(self, clusters):
        """
            Get the metric of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`

            :returns: np.Array -- Array containing the metric of the clusters
        """
        metrics = np.zeros(len(clusters))
        for i, cluster in enumerate(clusters):
            metrics[i] = cluster.getMetric()
        return metrics


class DensitySpawningCalculator(SpawningCalculator):
    """
        Subclass of Spawning calculator that ensures the definition of a density calculator.
    """

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        SpawningCalculator.__init__(self)
        self.type = "BaseDensityClass"  # change for abstract attribute
        self.densityCalculator = densityCalculator

    def calculateDensities(self, clusters):
        """
            Calculate the densities of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`

            :returns: np.Array -- Array containing the density of the clusters
        """
        densities = np.zeros(len(clusters))
        for i, cluster in enumerate(clusters):
            contacts = cluster.getContacts()
            cluster.density = self.densityCalculator.calculate(contacts, cluster.contactThreshold)
            densities[i] = cluster.density
        return densities


class IndependentRunsCalculator(SpawningCalculator):

    def __init__(self):
        SpawningCalculator.__init__(self)
        self.calculator = SameWeightDegeneracyCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.independent

    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        return self.calculator.calculate(clusters, trajToDistribute, spawningParams, currentEpoch)

    def log(self):
        """
            Log spawning information
        """
        pass

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration):
        """
            Write last trajectory structure as initial one for the next iteration

            :param outputPathConstants: Output constants that depend on the path
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param degeneracyOfRepresentatives: List with the degeneracy of
                each cluster (number of processors that will start from that state)
            :type degeneracyOfRepresentatives: list
            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`
            :param iteration: Number of epoch
            :type iteration: int

            :returns: int, list -- number of processors, list with the
                snapshot from which the trajectories will start in the next iteration
        """
        procMapping = []
        trajWildcard = os.path.join(outputPathConstants.epochOutputPathTempletized, constants.trajectoryBasename)
        trajectories = glob.glob(trajWildcard % (iteration-1))
        for num, trajectory in enumerate(trajectories):
            snapshots = utilities.getSnapshots(trajectory)
            lastSnapshot = snapshots[-1]
            nSnapshots = len(snapshots)
            del snapshots

            numTraj = int(os.path.splitext(trajectory.rsplit("_", 1)[-1])[0])
            outputFilename = outputPathConstants.tmpInitialStructuresTemplate % (iteration, num)
            procMapping.append((iteration-1, numTraj, nSnapshots-1))

            with open(outputFilename, 'w') as f:
                f.write(lastSnapshot)

        return len(trajectories), procMapping


class SameWeightDegeneracyCalculator(SpawningCalculator):

    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight

    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        numClusters = len(clusters)
        trajToDistribute = min(trajToDistribute, numClusters)
        samples = random.sample(range(numClusters), trajToDistribute)
        degeneracy = [0] * len(clusters)
        for sample in samples:
            degeneracy[sample] = 1
        return degeneracy

    def log(self):
        """
            Log spawning information
        """
        pass


class InverselyProportionalToPopulationCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional

    def log(self):
        """
            Log spawning information
        """
        pass

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        sizes = getSizes(clusters)
        densities = self.calculateDensities(clusters)

        if densities.any():
            weights = sizes/densities
        else:
            weights = sizes

        argweights = weights.argsort()
        weights_trimmed = np.zeros(len(sizes)) + 1e6
        weights_trimmed[argweights[:trajToDistribute]] = weights[argweights[:trajToDistribute]]
        return self.divideInverselyProportionalToArray(weights_trimmed, trajToDistribute)


class EpsilonDegeneracyCalculator(DensitySpawningCalculator):
    """
        It uses epsilon * numTraj trajectories proportional to their energy and the rest inversely proportional to each cluster's population
        We only consider the nclusters with best metric
    """

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator(), condition=blockNames.SpawningParams.minValue):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator(densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None
        self.condition = condition

    def log(self):
        """
            Log spawning information
        """
        if self.degeneracyTotal is not None:
            print("[SpawningLog] Total: %s" % str(self.degeneracyTotal))
        if self.degeneracyInverselyProportional is not None:
            print("[SpawningLog] Inversely prop: %s" % str(self.degeneracyInverselyProportional))
        if self.degeneracyMetricProportional is not None:
            print("[SpawningLog] Metric prop:    %s" % str(self.degeneracyMetricProportional))

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        trajToMetricProportional = int(spawningParams.epsilon * trajToDistribute)
        trajToInverselyProportional = trajToDistribute - trajToMetricProportional

        self.degeneracyInverselyProportional = self.inverselyProportionalCalculator.calculate(clusters, trajToInverselyProportional, spawningParams)
        self.degeneracyMetricProportional = self.divideProcessorsMetricProportional(clusters, trajToMetricProportional, spawningParams)

        self.degeneracyTotal = np.array(self.degeneracyInverselyProportional) + np.array(self.degeneracyMetricProportional)
        return self.degeneracyTotal.tolist()

    def divideProcessorsMetricProportional(self, clusters, trajToDistribute, spawningParams):
        """
            Distribute the trajectories among the clusters according to their
            metric.

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param weights: Weight of each cluster
            :type weights: np.Array
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int

            :returns: list -- List with the number of processors allocated to
                each cluster
        """

        metrics = self.getMetrics(clusters)

        if isinstance(metrics, list):
            metrics = np.array(metrics)

        # Shift so that differences become larger.
        # Also, we can now merge positive & negative values
        # Alternatives: Boltzmann weights
        if spawningParams.condition == blockNames.SpawningParams.minValue:
            shiftValue = np.max(metrics)
        else:
            shiftValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, shiftValue)
        bestClusters = shiftedMetrics.argsort()

        if spawningParams.condition == blockNames.SpawningParams.minValue:
            shiftedMetrics[bestClusters[spawningParams.nclusters:]] = 0  # only consider best ones
        else:
            shiftedMetrics[bestClusters[:-spawningParams.nclusters]] = 0  # only consider best ones

        metricWeights = spawningParams.metricWeights
        if metricWeights == blockNames.SpawningParams.linear:

            # all shiftedMetrics <= 0, sum(shiftedMetrics) < 0 => weights >= 0
            if abs(shiftedMetrics.sum()) < 1e-8:
                weights = np.ones(len(metrics))/len(metrics)
            else:
                weights = (1.*shiftedMetrics)/sum(shiftedMetrics)

        elif metricWeights == blockNames.SpawningParams.boltzmann:
            T = spawningParams.temperature
            kbT = 0.001987*T
            if abs(shiftedMetrics.sum()) < 1e-8:
                weights = np.ones(len(metrics))/len(metrics)
            else:
                weights = np.exp(-shiftedMetrics/kbT)
                weights /= sum(weights)
        else:
            raise ValueError("No appropiate value for the metricWeights "
                             "was found, please specify a correct value. The "
                             "default value of the metrics weighting is linear")

        return self.divideTrajAccordingToWeights(weights, trajToDistribute)


class VariableEpsilonDegeneracyCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.epsilonDegeneracyCalculator = EpsilonDegeneracyCalculator(densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.variableEpsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None
        self.maxContacts = None
        # print variable epsilon information
        epsilon_file = open("epsilon_values.txt", "w")
        epsilon_file.write("Iteration\tEpsilon\n")
        epsilon_file.close()

    def linearVariation(self, spawningParams, currentEpoch):
        """
            Calculate linear variation of epsilon with respect to the iteraetion

            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
        """
        if currentEpoch == 0:
            spawningParams.epsilon = spawningParams.minEpsilon
            return

        middleWindow = int(spawningParams.period/2)
        leftWindow = int(spawningParams.maxEpsilonWindow/2)
        rightWindow = leftWindow+middleWindow
        if currentEpoch == spawningParams.period-1:
            # Avoid negative epsilon for minEpsilon = 0
            return
        rateEpsilonVariation = [(spawningParams.maxEpsilon-spawningParams.minEpsilon)/(middleWindow-leftWindow), (spawningParams.maxEpsilon-spawningParams.minEpsilon)/(spawningParams.period-rightWindow-1)]
        spawningParams.epsilon += return_sign(currentEpoch, leftWindow,
                                              middleWindow, rightWindow) * rateEpsilonVariation[currentEpoch > middleWindow]

    def contactsVariation(self, clusters, spawningParams):
        """
            Calculate the variation of epsilon according to the contacts ratio

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
        """
        if self.maxContacts is None:
            self.maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        if spawningParams.epsilon < spawningParams.maxEpsilon:
            spawningParams.epsilon += calculateContactsVar(maxContacts-self.maxContacts, spawningParams.maxEpsilon)
        self.maxContacts = maxContacts

    def calculateEpsilonValue(self, spawningParams, currentEpoch, clusters):
        """
            Calculate variation of epsilon according to the selected parameters

            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
        """
        if spawningParams.varEpsilonType == blockNames.VariableEpsilonTypes.linearVariation:
            if currentEpoch is None or spawningParams.variationWindow < currentEpoch:
                spawningParams.epsilon = spawningParams.minEpsilon
                return
            self.linearVariation(spawningParams,
                                 (currentEpoch % spawningParams.period))
        elif spawningParams.varEpsilonType == blockNames.VariableEpsilonTypes.contactsVariation:
            self.contactsVariation(clusters, spawningParams)
        else:
            sys.exit("Unknown epsilon variation type! Choices are: " +
                     str(spawningTypes.EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY.values()))

    def logVariableEpsilon(self, epsilon, epoch):
        """
            Log spawning information
        """
        with open("epsilon_values.txt", "a") as epsilon_file:
            epsilon_file.write("%d\t%f\n" % (epoch, epsilon))

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        self.calculateEpsilonValue(spawningParams, currentEpoch, clusters)
        self.logVariableEpsilon(spawningParams.epsilon, currentEpoch)
        return self.epsilonDegeneracyCalculator.calculate(clusters, trajToDistribute, spawningParams, currentEpoch)


class SimulatedAnnealingCalculator(SpawningCalculator):

    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.simulatedAnnealing

    def log(self):
        """
            Log spawning information
        """
        pass

    def computeTemperature(self, params, epoch):
        T = params.temperature - params.decrement*epoch
        if T < 300:
            return 300
        else:
            return T

    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        metrics = self.getMetrics(clusters)

        minimumValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, minimumValue)

        T = self.computeTemperature(spawningParams, currentEpoch)
        kbT = 0.001987*T
        weights = np.exp(-shiftedMetrics/kbT)
        weights /= sum(weights)

        return self.divideProportionalToArray(weights, trajToDistribute)


class FASTDegeneracyCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.FAST
        self.densityCalculator = densityCalculator

    def normaliseArray(self, array):
        maxValue = float(np.max(array))
        minValue = float(np.min(array))

        if maxValue-minValue > 1e-5:
            normalisedArray = np.subtract(maxValue, array)/(maxValue - minValue)
        else:
            # all elements are equal
            n = len(array)
            normalisedArray = np.ones(n)/float(n)
        return normalisedArray

    def calculateNormalisedSizes(self, clusters):
        sizes = getSizes(clusters)

        densities = self.calculateDensities(clusters)
        weightedSizes = sizes/densities

        return self.normaliseArray(weightedSizes)

    def calculateNormalisedMetrics(self, clusters):
        metrics = self.getMetrics(clusters)
        return self.normaliseArray(metrics)

    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        normalisedSizes = self.calculateNormalisedSizes(clusters)
        normalisedMetrics = self.calculateNormalisedMetrics(clusters)

        weight = normalisedSizes + 1.*normalisedMetrics

        return self.divideProportionalToArray(weight, trajToDivide)

    def log(self):
        """
            Log spawning information
        """
        pass


class UCBCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.UCB
        self.prevMetrics = np.array([0.0])
        self.averages = []
        self.beta = 1.0
        self.averageMetric = 0
        self.epoch = np.array([0.0])

    def log(self):
        """
            Log spawning information
        """
        pass

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        self.epoch += 1
        sizes = np.array(getSizes(clusters))
        densities = self.calculateDensities(clusters)
        metrics = np.array(self.getMetrics(clusters))
        # self.averageMetric = np.mean(metrics)
        maximumValue = np.max(metrics)
        shiftedMetrics = np.subtract(metrics, maximumValue)
        if abs(shiftedMetrics.sum()) < 1e-8:
            weights = np.ones(len(metrics))/len(metrics)
        else:
            weights = (1.*shiftedMetrics)/min(shiftedMetrics)
        # for i, metric in enumerate(metrics):
        #     if i < len(self.prevMetrics):
        #         self.prevMetrics[i] = abs(metric-self.prevMetrics[i])/abs(metric)
        #         self.averages[i] += (self.prevMetrics[i]-self.averages[i])/sizes[i]
        #     else:
        #         self.prevMetrics.append(-(metric-self.averageMetric)/abs(metric))
        #         self.averages.append(-(metric-self.averageMetric)/abs(metric))
        l = self.prevMetrics.size
        n = weights.size
        try:
            self.prevMetrics = np.pad(self.prevMetrics, (0, n-l), str('constant'), constant_values=(0.0))
            self.epoch = np.pad(self.epoch, (0, n-l), str('constant'), constant_values=(1.0))
        except AttributeError:
            # Numpy version in life is too old to use pad function
            prevMetrics = np.zeros_like(weights)
            epochs = np.ones_like(weights)
            prevMetrics[:l] = self.prevMetrics
            self.prevMetrics = prevMetrics
            epochs[:l] = self.epoch
            self.epoch = epochs
        # avg[:l] = self.prevMetrics[:l]
        self.prevMetrics += (weights-self.prevMetrics)/self.epoch
        argweights = self.prevMetrics.argsort()
        weights_trimmed = np.zeros(len(sizes))
        weights_trimmed[argweights[-trajToDistribute:]] = self.prevMetrics[argweights[-trajToDistribute:]]
        # values = weights_trimmed+spawningParams.alpha*np.sqrt((1/sizes))
        # values = self.beta*weights_trimmed**2+spawningParams.alpha*(1/sizes**2)
        values = self.beta*weights_trimmed**2+spawningParams.alpha*(1/sizes)
        # values = self.beta*weights_trimmed**2+(spawningParams.alpha/((np.log2(currentEpoch+2))**(1/4.0)))*(1/sizes)
        # minVal = np.min(values)
        # if minVal < 0:
        #     # if there is a negative value shift all the values so that the min
        #     # value is zero
        #    values += abs(minVal)
        if densities.any():
            weights = values*densities
        else:
            weights = values
        argweights = weights.argsort()
        weights_trimmed = np.zeros(len(sizes))
        weights_trimmed[argweights[-trajToDistribute:]] = weights[argweights[-trajToDistribute:]]
        return self.divideProportionalToArray(weights_trimmed, trajToDistribute)


class REAPCalculator(DensitySpawningCalculator):
    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        """
            Spawning following the Reinforcement learning based Adaptive
            samPling (REAP) (Shamsi et al., arXiv, Oct 2017), where the reward
            given by the exploration on several reaction coordinates
            is maximized
        """
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.REAP
        self.weights = None
        self.metricInd = None
        self.rewards = None
        self.degeneracy = None
        # constraints so the weights have values between 0 and 1
        self.cons = ({'type': 'eq', 'fun': lambda x: np.array(x.sum()-1)})
        self.bounds = None

    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: list -- List containing the degeneracy of the clusters
        """
        population = []
        metrics = []
        if self.metricInd is None:
            if spawningParams.metricInd == -1:
                self.metricInd = list(range(3, clusters[0].metrics.size))
            else:
                self.metricInd = spawningParams.metricInd
            self.bounds = [(0, 1)]*len(self.metricInd)

        # Gather population and metrics data for all clusters
        for cluster in clusters:
            population.append(float(cluster.elements))
            metrics.append([cluster.metrics[i] for i in self.metricInd])
        self.degeneracy = np.zeros_like(population)
        metrics = np.array(metrics).T
        meanRew = np.mean(metrics, axis=1)
        stdRew = np.std(metrics, axis=1)
        # Filter top least populated clusters
        population = np.array(population)
        densities = self.calculateDensities(clusters)
        if densities.any():
            population /= densities
        argweights = np.argsort(population)
        metrics = metrics[:, argweights[:trajToDivide]]
        # Shift and scale all metrics to have mean 0 and std 1, so that the
        # weight of each metric is not affected by its magnitude (i.e. binding
        # energy ~ 10**2 while SASA <= 1)
        rewProv = np.abs(metrics-meanRew[:, np.newaxis])/stdRew[:, np.newaxis]

        if self.weights is None:
            self.weights = np.ones(len(self.metricInd))/len(self.metricInd)
        else:
            optimResult = optim.minimize(reward, self.weights, args=(rewProv,),
                                         method="SLSQP", constraints=self.cons,
                                         bounds=self.bounds)
            self.weights = optimResult.x
        self.rewards = (self.weights[:, np.newaxis]*rewProv).sum(axis=0)
        self.degeneracy[argweights[:trajToDivide]] = self.divideProportionalToArray(self.rewards, trajToDivide)
        return self.degeneracy.tolist()

    def log(self):
        """
            Log spawning information
        """
        if self.degeneracy is not None:
            print("[SpawningLog] Total: %s" % str(self.degeneracy))
        print("Metric indices")
        print(self.metricInd)
        print("Spawning weights")
        print(self.weights)


class NullSpawningCalculator(SpawningCalculator):
    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.null

    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        """
            Calculate the degeneracy of the clusters. In this particular class
            no spawning is performed, so this function just returns None

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :returns: None
        """
        return None
