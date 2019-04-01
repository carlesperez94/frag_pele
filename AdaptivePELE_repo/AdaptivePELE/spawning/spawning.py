from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import math
import sys
import numpy as np
import random
import scipy.optimize as optim
import os
import glob
from abc import abstractmethod
from scipy.linalg import lu, solve
from AdaptivePELE.constants import blockNames
from AdaptivePELE.constants import constants
from AdaptivePELE.utilities import utilities
from AdaptivePELE.spawning import spawningTypes
from AdaptivePELE.spawning import densitycalculator
try:
    # Check if the basestring type if available, this will fail in python3
    basestring
except NameError:
    basestring = str
PYEMMA = True
try:
    import pyemma.msm as msm
    from AdaptivePELE.freeEnergies import computeDeltaG as computedG
except ImportError:
    PYEMMA = False
MATPLOTLIB = True
try:
    import matplotlib.pyplot as plt
except ImportError:
    MATPLOTLIB = False
if MATPLOTLIB:
    try:
        # Minotauro version is too old
        plt.style.use("ggplot")
    except AttributeError:
        pass


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
        spawningParams = SpawningParams()
        spawningParams.buildSpawningParameters(spawningBlock)

        spawningCalculatorBuilder = SpawningBuilder()
        spawningCalculator = spawningCalculatorBuilder.buildSpawningCalculator(spawningBlock, spawningParams)

        return spawningCalculator


class SpawningBuilder:

    def buildSpawningCalculator(self, spawningBlock, spawningParams):
        """
            Build the selected spawning calculator object

            :param spawningBlock: Block of the control file with the spawning
                parameters
            :type spawningBlock: dict
            :param spawningParams: Object containing the parameters of the spawning
            :type spawningParams: :py:class:`.SpawningParams`

            :returns: :py:class:`.SpawningCalculator` -- SpawningCalculator
                object
        """
        densityBuilder = densitycalculator.DensityCalculatorBuilder()
        densityCalculator = densityBuilder.build(spawningBlock)

        spawningTypeString = spawningBlock[blockNames.StringSpawningTypes.type]
        if spawningTypeString == blockNames.StringSpawningTypes.sameWeight:
            spawningCalculator = SameWeightDegeneracyCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.independent:
            spawningCalculator = IndependentRunsCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.independentMetric:
            spawningCalculator = IndependentMetricCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.inverselyProportional:
            spawningCalculator = InverselyProportionalToPopulationCalculator(spawningParams, densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.epsilon:
            spawningCalculator = EpsilonDegeneracyCalculator(spawningParams, densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.fast:
            spawningCalculator = FASTDegeneracyCalculator(spawningParams, densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.variableEpsilon:
            spawningCalculator = VariableEpsilonDegeneracyCalculator(spawningParams, densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.UCB:
            spawningCalculator = UCBCalculator(spawningParams, densityCalculator)
        elif spawningTypeString == blockNames.StringSpawningTypes.REAP:
            spawningCalculator = REAPCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.null:
            spawningCalculator = NullSpawningCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.ProbabilityMSMCalculator:
            spawningCalculator = ProbabilityMSMCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.MetastabilityMSMCalculator:
            spawningCalculator = MetastabilityMSMCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.UncertaintyMSMCalculator:
            spawningCalculator = UncertaintyMSMCalculator(spawningParams)
        elif spawningTypeString == blockNames.StringSpawningTypes.IndependentMSMCalculator:
            spawningCalculator = IndependentMSMCalculator(spawningParams)
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
        self.lagtime = None
        self.minPos = None
        self.sasaColumn = None
        self.filter_value = None
        self.filter_col = None
        self.filterByMetric = None

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

        if spawningType not in spawningTypes.MSMSpawning:
            # reportFilename is now mandatory for all spawning not related to
            # MSM
            self.reportFilename = spawningParamsBlock[blockNames.SpawningParams.report_filename]
            # this paramaters are optional for many methods, some of the
            # assignments are redundant and might be later overriden, but since
            # there are many spawning methods we trade in a bit of efficient for
            # sanity
            self.condition = spawningParamsBlock.get(blockNames.SpawningParams.condition,
                                                     blockNames.SpawningParams.minValue)
            self.filterByMetric = spawningParamsBlock.get(blockNames.SpawningParams.filterByMetric, False)
            self.filter_col = spawningParamsBlock.get(blockNames.SpawningParams.filter_col)
            if self.filter_col is None and self.filterByMetric:
                raise utilities.RequiredParameterMissingException("Column not specified for cluster filtering")
            if self.filter_col is not None:
                self.filter_col -= 1
            self.filter_value = spawningParamsBlock.get(blockNames.SpawningParams.filter_value)
            if self.filter_value is None and self.filterByMetric:
                raise utilities.RequiredParameterMissingException("Filtering value not specified for cluster filtering")
            self.reportCol = spawningParamsBlock.get(blockNames.SpawningParams.report_col)
            if self.reportCol is not None:
                self.reportCol -= 1

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

        if spawningType == blockNames.StringSpawningTypes.independentMetric:
            # Start counting the columns by 1
            self.reportCol = spawningParamsBlock[blockNames.SpawningParams.report_col]-1
            self.condition = spawningParamsBlock.get(blockNames.SpawningParams.condition,
                                                     blockNames.SpawningParams.minValue)
        if spawningType in spawningTypes.MSMSpawning:
            self.lagtime = spawningParamsBlock[blockNames.SpawningParams.lagtime]
            self.condition = spawningParamsBlock.get(blockNames.SpawningParams.condition,
                                                     blockNames.SpawningParams.minValue)
            self.minPos = spawningParamsBlock.get(blockNames.SpawningParams.minPos)
            self.sasaColumn = spawningParamsBlock.get(blockNames.SpawningParams.SASA_column)
            if self.sasaColumn is not None:
                # start counting from one in the report file
                self.sasaColumn -= 1


class SpawningCalculator:
    """
        The purpose of this abstract class is to contain the behaviour of the different strategies for the spawning.
        Spawning is the way in which we split the different explorers at the begining of each epoch.
    """

    def __init__(self):
        self.type = "BaseClass"  # change for abstract attribute

    @abstractmethod
    def calculate(self, clusters, trajToDivide, currentEpoch=None, outputPathConstants=None):
        pass

    @abstractmethod
    def log(self):
        """
            Log spawning information
        """
        pass

    @abstractmethod
    def createPlots(self, outputPathConstants, currentEpoch, clusters):
        """
            Create the plots to do a quick analysis of the MSM and dG calculation

            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
        """
        pass

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration, topologies=None):
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
            :param topology_file: Topology file for non-pdb trajectories
            :type topology_file: str
            :param topology: Topology like structure to write PDB from xtc files
            :type topology: list

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

    def shouldWriteStructures(self):
        return True


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

    def __init__(self, parameters):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.independent
        self.parameters = parameters

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration, topologies=None):
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
            :param topology_file: Topology file for non-pdb trajectories
            :type topology_file: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`

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
            if isinstance(lastSnapshot, basestring):
                with open(outputFilename, 'w') as f:
                    f.write(lastSnapshot)
            else:
                utilities.write_mdtraj_object_PDB(lastSnapshot, outputFilename, topologies.getTopology(iteration-1, numTraj))

        return len(trajectories), procMapping


class IndependentMetricCalculator(SpawningCalculator):

    def __init__(self, parameters):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.independentMetric
        self.parameters = parameters

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration, topologies=None):
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
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`

            :returns: int, list -- number of processors, list with the
                snapshot from which the trajectories will start in the next iteration
        """
        procMapping = []
        trajWildcard = os.path.join(outputPathConstants.epochOutputPathTempletized, constants.trajectoryBasename)
        trajectories = glob.glob(trajWildcard % (iteration-1))
        for num in range(len(trajectories)):
            reportFilename = os.path.join(outputPathConstants.epochOutputPathTempletized % (iteration-1),
                                          "%s_%d" % (self.parameters.reportFilename, num+1))
            metric_array = np.genfromtxt(reportFilename, missing_values="--", filling_values=0)
            if len(metric_array.shape) < 2:
                metric_array = metric_array[np.newaxis, :]
            trajectory = glob.glob("%s_%d.*" % (trajWildcard % (iteration-1), num+1))
            assert len(trajectory) == 1, "Too many trajectories found in IndependentMetricCalculator"
            trajectory = trajectory[0]
            if self.parameters.condition == blockNames.SpawningParams.minValue:
                snapshot_ind = np.argmin(metric_array[:, self.parameters.reportCol])
            else:
                snapshot_ind = np.argmax(metric_array[:, self.parameters.reportCol])
            snapshots = utilities.getSnapshots(trajectory)
            snapshot = snapshots[snapshot_ind]
            del snapshots

            numTraj = int(os.path.splitext(trajectory.rsplit("_", 1)[-1])[0])
            outputFilename = outputPathConstants.tmpInitialStructuresTemplate % (iteration, num)
            procMapping.append((iteration-1, numTraj, snapshot_ind))
            if isinstance(snapshot, basestring):
                with open(outputFilename, 'w') as f:
                    f.write(snapshot)
            else:
                utilities.write_mdtraj_object_PDB(snapshot, outputFilename, topologies.getTopology(iteration-1, numTraj))

        return len(trajectories), procMapping


class SameWeightDegeneracyCalculator(SpawningCalculator):

    def __init__(self, parameters):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight
        self.parameters = parameters

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

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

    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional
        self.parameters = parameters

    def log(self):
        """
            Log spawning information
        """
        pass

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

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

    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator(parameters, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None
        self.parameters = parameters

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

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        trajToMetricProportional = int(self.parameters.epsilon * trajToDistribute)
        trajToInverselyProportional = trajToDistribute - trajToMetricProportional

        self.degeneracyInverselyProportional = self.inverselyProportionalCalculator.calculate(clusters, trajToInverselyProportional)
        self.degeneracyMetricProportional = self.divideProcessorsMetricProportional(clusters, trajToMetricProportional)

        self.degeneracyTotal = np.array(self.degeneracyInverselyProportional) + np.array(self.degeneracyMetricProportional)
        return self.degeneracyTotal.tolist()

    def divideProcessorsMetricProportional(self, clusters, trajToDistribute):
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
        if self.parameters.condition == blockNames.SpawningParams.minValue:
            shiftValue = np.max(metrics)
        else:
            shiftValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, shiftValue)
        bestClusters = shiftedMetrics.argsort()

        if self.parameters.condition == blockNames.SpawningParams.minValue:
            shiftedMetrics[bestClusters[self.parameters.nclusters:]] = 0  # only consider best ones
        else:
            shiftedMetrics[bestClusters[:-self.parameters.nclusters]] = 0  # only consider best ones

        metricWeights = self.parameters.metricWeights
        if metricWeights == blockNames.SpawningParams.linear:

            # all shiftedMetrics <= 0, sum(shiftedMetrics) < 0 => weights >= 0
            if abs(shiftedMetrics.sum()) < 1e-8:
                weights = np.ones(len(metrics))/len(metrics)
            else:
                weights = (1.*shiftedMetrics)/sum(shiftedMetrics)

        elif metricWeights == blockNames.SpawningParams.boltzmann:
            T = self.parameters.temperature
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

    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.epsilonDegeneracyCalculator = EpsilonDegeneracyCalculator(parameters, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.variableEpsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None
        self.maxContacts = None
        self.parameters = parameters
        # print variable epsilon information
        epsilon_file = open("epsilon_values.txt", "w")
        epsilon_file.write("Iteration\tEpsilon\n")
        epsilon_file.close()

    def linearVariation(self, currentEpoch):
        """
            Calculate linear variation of epsilon with respect to the iteraetion

            :param currentEpoch: Current iteration number
            :type currentEpoch: int
        """
        if currentEpoch == 0:
            self.parameters.epsilon = self.parameters.minEpsilon
            return

        middleWindow = int(self.parameters.period/2)
        leftWindow = int(self.parameters.maxEpsilonWindow/2)
        rightWindow = leftWindow+middleWindow
        if currentEpoch == self.parameters.period-1:
            # Avoid negative epsilon for minEpsilon = 0
            return
        rateEpsilonVariation = [(self.parameters.maxEpsilon-self.parameters.minEpsilon)/(middleWindow-leftWindow), (self.parameters.maxEpsilon-self.parameters.minEpsilon)/(self.parameters.period-rightWindow-1)]
        self.parameters.epsilon += return_sign(currentEpoch, leftWindow,
                                               middleWindow, rightWindow) * rateEpsilonVariation[currentEpoch > middleWindow]

    def contactsVariation(self, clusters):
        """
            Calculate the variation of epsilon according to the contacts ratio

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
        """
        if self.maxContacts is None:
            self.maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        if self.parameters.epsilon < self.parameters.maxEpsilon:
            self.parameters.epsilon += calculateContactsVar(maxContacts-self.maxContacts, self.parameters.maxEpsilon)
        self.maxContacts = maxContacts

    def calculateEpsilonValue(self, currentEpoch, clusters):
        """
            Calculate variation of epsilon according to the selected parameters

            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
        """
        if self.parameters.varEpsilonType == blockNames.VariableEpsilonTypes.linearVariation:
            if currentEpoch is None or self.parameters.variationWindow < currentEpoch:
                self.parameters.epsilon = self.parameters.minEpsilon
                return
            self.linearVariation(currentEpoch % self.parameters.period)
        elif self.parameters.varEpsilonType == blockNames.VariableEpsilonTypes.contactsVariation:
            self.contactsVariation(clusters)
        else:
            sys.exit("Unknown epsilon variation type! Choices are: " +
                     str(spawningTypes.EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY.values()))

    def logVariableEpsilon(self, epsilon, epoch):
        """
            Log spawning information
        """
        with open("epsilon_values.txt", "a") as epsilon_file:
            epsilon_file.write("%d\t%f\n" % (epoch, epsilon))

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        self.calculateEpsilonValue(currentEpoch, clusters)
        self.logVariableEpsilon(self.parameters.epsilon, currentEpoch)
        return self.epsilonDegeneracyCalculator.calculate(clusters, trajToDistribute, currentEpoch)


class SimulatedAnnealingCalculator(SpawningCalculator):

    def __init__(self, parameters):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.simulatedAnnealing
        self.parameters = parameters

    def log(self):
        """
            Log spawning information
        """
        pass

    def computeTemperature(self, epoch):
        T = self.parameters.temperature - self.parameters.decrement*epoch
        if T < 300:
            return 300
        else:
            return T

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        metrics = self.getMetrics(clusters)

        minimumValue = np.min(metrics)
        shiftedMetrics = np.subtract(metrics, minimumValue)

        T = self.computeTemperature(currentEpoch)
        kbT = 0.001987*T
        weights = np.exp(-shiftedMetrics/kbT)
        weights /= sum(weights)

        return self.divideProportionalToArray(weights, trajToDistribute)


class FASTDegeneracyCalculator(DensitySpawningCalculator):

    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.FAST
        self.densityCalculator = densityCalculator
        self.parameters = parameters

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

    def calculate(self, clusters, trajToDivide, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

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

    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.UCB
        self.prevMetrics = np.array([0.0])
        self.averages = []
        self.beta = 1.0
        self.averageMetric = 0
        self.epoch = np.array([0.0])
        self.parameters = parameters

    def log(self):
        """
            Log spawning information
        """
        pass

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

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
        # values = weights_trimmed+self.parameters.alpha*np.sqrt((1/sizes))
        # values = self.beta*weights_trimmed**2+self.parameters.alpha*(1/sizes**2)
        values = self.beta*weights_trimmed**2+self.parameters.alpha*(1/sizes)
        # values = self.beta*weights_trimmed**2+(self.parameters.alpha/((np.log2(currentEpoch+2))**(1/4.0)))*(1/sizes)
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
    def __init__(self, parameters, densityCalculator=densitycalculator.NullDensityCalculator()):
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
        self.parameters = parameters

    def calculate(self, clusters, trajToDivide, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        population = []
        metrics = []
        if self.metricInd is None:
            if self.parameters.metricInd == -1:
                self.metricInd = list(range(3, clusters[0].metrics.size))
            else:
                self.metricInd = self.parameters.metricInd
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
    def __init__(self, parameters):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.null
        self.parameters = parameters

    def calculate(self, clusters, trajToDivide, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters. In this particular class
            no spawning is performed, so this function just returns None

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: None
        """
        return None

    def shouldWriteStructures(self):
        return False


class MSMCalculator(SpawningCalculator):

    def __init__(self, parameters):
        if not PYEMMA:
            raise utilities.UnsatisfiedDependencyException("Pyemma module is necessary ot use MSM based spawning")
        SpawningCalculator.__init__(self)
        self.type = "BaseClass"  # change for abstract attribute
        self.parameters = parameters
        self.MSM = None

    def estimateMSM(self, dtrajs, outputPathConstants, currentEpoch):
        """
            Estimate and MSM using PyEMMA

            :param dtrajs: Discretized trajectories to estimate the Markov model
            :type dtrajs: np.ndarray
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

            :return: object -- Object containing the estimated MSM
        """
        self.MSM = msm.estimate_markov_model(dtrajs, self.parameters.lagtime)
        if outputPathConstants is not None and currentEpoch is not None:
            utilities.writeObject(outputPathConstants.MSMObjectEpoch % currentEpoch, self.MSM)

    def calculatedG(self, clusters, outputPathConstants, currentEpoch):
        """
            Calculate the ligand-binding dG from the estimated MSM

            :param clusters: Cluster centers
            :type clusters: np.ndarray
            :param trajPath: Path to the trajectories to analyse
            :type trajPath: str
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int

        """
        # need clusters for this step
        pi, clusters = computedG.ensure_connectivity(self.MSM, clusters)
        d = 0.75
        originalFilenames = glob.glob(os.path.join(outputPathConstants.allTrajsPath, "*traj*.dat"))
        originalCoordinates = computedG.gather_coordinates(originalFilenames)
        bins = computedG.create_box(clusters, originalCoordinates, d)
        microstateVolume = computedG.calculate_microstate_volumes_new(clusters, originalCoordinates, bins, d)
        gpmf, string = computedG.calculate_pmf(microstateVolume, pi)
        print("bound    Delta G     Delta W     Binding Volume:     Binding Volume contribution")
        print(string)

        pmf_xyzg = np.hstack((clusters, np.expand_dims(gpmf, axis=1)))
        np.savetxt(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "pmf_xyzg.dat"), pmf_xyzg)
        distance = None
        if self.parameters.minPos is not None:
            distance = np.linalg.norm(clusters-self.parameters.minPos, axis=1)
        return pi, gpmf, distance

    def getSASAvalues(self, clusters, outputPathConstants):
        """
            Get the SASA value for each cluster from the report file

            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
        """
        sasa = []
        for cl in clusters:
            epoch, traj, snapshot = cl.trajPosition
            report_filename = glob.glob(os.path.join(outputPathConstants.epochOutputPathTempletized % epoch, "*report*_%d" % traj))[0]
            report_values = utilities.loadtxtfile(report_filename)
            sasa.append(report_values[snapshot, self.parameters.sasaColumn])
        return sasa

    def createPlots(self, outputPathConstants, currentEpoch, clustering):
        """
            Create the plots to do a quick analysis of the MSM and dG calculation

            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param clustering: Existing clusters
            :type clusters: :py:class:`.Clustering`
        """
        if not MATPLOTLIB:
            raise utilities.UnsatisfiedDependencyException("Matplotlib installation not found!!")
        if self.parameters.minPos is None and self.parameters.sasaColumn is None:
            print("Plots can't be created since neither a reference minimum position nor SASA have been provided!!!")
            return
        clusters = clustering.pyemma_clustering.clusterCenters
        prob, G, distance = self.calculatedG(clusters, outputPathConstants, currentEpoch)
        sasa_values = None
        if self.parameters.sasaColumn is not None:
            sasa_values = self.getSASAvalues(clustering.clusters, outputPathConstants)
        if sasa_values is not None and distance is not None:
            f, axarr = plt.subplots(1, 2)
            axarr[0].scatter(distance, prob)
            axarr[0].set_xlabel("Distance to minimum")
            axarr[0].set_ylabel("Stationary distribution")
            axarr[1].scatter(sasa_values, prob)
            axarr[1].set_xlabel("SASA")
            f.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "eigenvector.png"))
            f, axarr = plt.subplots(1, 2)
            axarr[0].scatter(distance, G)
            axarr[0].set_xlabel("Distance to minimum")
            axarr[0].set_ylabel("PMF")
            axarr[1].scatter(sasa_values, G)
            axarr[1].set_xlabel("SASA")
            f.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "PMF.png"))
        elif distance is not None:
            plt.figure()
            plt.scatter(distance, prob)
            plt.xlabel("Distance to minimum")
            plt.ylabel("Stationary distribution")
            plt.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "eigenvector.png"))
            plt.figure()
            plt.scatter(distance, G)
            plt.xlabel("Distance to minimum")
            plt.ylabel("PMF")
            plt.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "PMF.png"))
        elif sasa_values is not None:
            plt.figure()
            plt.scatter(sasa_values, prob)
            plt.xlabel("SASA")
            plt.ylabel("Stationary distribution")
            plt.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "eigenvector.png"))
            plt.figure()
            plt.scatter(sasa_values, G)
            plt.xlabel("SASA")
            plt.ylabel("PMF")
            plt.savefig(os.path.join(outputPathConstants.epochOutputPathTempletized % currentEpoch, "PMF.png"))


class ProbabilityMSMCalculator(MSMCalculator):
    def __init__(self, parameters):
        MSMCalculator.__init__(self, parameters)
        self.type = spawningTypes.SPAWNING_TYPES.ProbabilityMSMCalculator
        self.parameters = parameters

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        # estimate MSM from clustering object
        self.estimateMSM(clusters.dtrajs, outputPathConstants, currentEpoch)
        nclusters = self.MSM.nstates_full
        # distribute seeds using the MSM
        probabilities = np.zeros(nclusters)
        for i, index in enumerate(self.MSM.active_set):
            probabilities[index] = self.MSM.stationary_distribution[i]
        if self.parameters.condition == blockNames.SpawningParams.minValue:
            sortedProbs = np.argsort(probabilities)
            probabilities = 1 - probabilities
        else:
            sortedProbs = np.argsort(probabilities)[::-1]
        neutral = 0.0
        probabilities[sortedProbs[trajToDistribute:]] = neutral
        if abs(probabilities.sum()) < 1e-8:
            probabilities = np.ones(nclusters)/nclusters
        else:
            probabilities /= sum(probabilities)
        return self.divideTrajAccordingToWeights(probabilities, trajToDistribute)


class MetastabilityMSMCalculator(MSMCalculator):
    def __init__(self, parameters):
        MSMCalculator.__init__(self, parameters)
        self.type = spawningTypes.SPAWNING_TYPES.MetastabilityMSMCalculator
        self.parameters = parameters

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        # estimate MSM from clustering object
        self.estimateMSM(clusters.dtrajs, outputPathConstants, currentEpoch)
        nclusters = self.MSM.nstates_full
        # distribute seeds using the MSM
        counts = self.MSM.count_matrix_full
        metastability = counts.diagonal()/counts.sum()

        if self.parameters.condition == blockNames.SpawningParams.minValue:
            sortedProbs = np.argsort(metastability)
            metastability = 1 - metastability
        else:
            sortedProbs = np.argsort(metastability)[::-1]
        metastability[sortedProbs[trajToDistribute:]] = 0.0
        if abs(metastability.sum()) < 1e-8:
            metastability = np.ones(nclusters)/nclusters
        else:
            metastability /= sum(metastability)
        return self.divideTrajAccordingToWeights(metastability, trajToDistribute)


class UncertaintyMSMCalculator(MSMCalculator):
    def __init__(self, parameters):
        MSMCalculator.__init__(self, parameters)
        self.type = spawningTypes.SPAWNING_TYPES.UncertaintyMSMCalculator
        self.parameters = parameters
        self.MSM = None

    def calculate_q(self, counts, nclusters):
        alpha = 1/float(nclusters)
        U_counts = counts + alpha
        w = U_counts.sum(axis=1)
        P = U_counts / w[:, np.newaxis]
        eigvalues, _ = np.linalg.eig(P)
        eigvalues.sort()
        eigvalues = eigvalues[::-1]
        ek = np.zeros(nclusters)
        ek[nclusters-1] = 1.0
        A = P - np.real(eigvalues[1])*np.eye(nclusters)
        perm, L, U = lu(A.T)
        x = solve(L.T, ek)
        xa = solve(U[:-1, :-1], -U[:-1, -1])
        xa = np.array(xa.tolist()+[1.0])
        norm_factor = xa.dot(perm.dot(x))
        si = np.outer(xa, perm.dot(x))/norm_factor
        q = []
        for i in range(nclusters):
            q.append(si[i].dot((np.diag(P[i])-np.outer(P[i], P[i])).dot(si[i])))
        return np.array(q), w

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        # estimate MSM from clustering object
        self.estimateMSM(clusters.dtrajs, outputPathConstants, currentEpoch)
        nclusters = self.MSM.nstates_full
        # distribute seeds using the MSM
        counts = self.MSM.count_matrix_full
        q, w = self.calculate_q(counts, nclusters)
        score = (q/(w+1))-(q/(w+1+trajToDistribute))
        score /= score.sum()
        sortedProbs = np.argsort(score)[::-1]
        score[sortedProbs[trajToDistribute:]] = 0.0
        if abs(score.sum()) < 1e-8:
            score = np.ones(nclusters)/nclusters
        else:
            score /= sum(score)
        return self.divideTrajAccordingToWeights(score, trajToDistribute)


class IndependentMSMCalculator(MSMCalculator):
    def __init__(self, parameters):
        MSMCalculator.__init__(self, parameters)
        self.type = spawningTypes.SPAWNING_TYPES.IndependentMSMCalculator
        self.parameters = parameters
        self.IndependentCalculator = IndependentRunsCalculator(parameters)

    def calculate(self, clusters, trajToDistribute, currentEpoch=None, outputPathConstants=None):
        """
            Calculate the degeneracy of the clusters

            :param clusters: Existing clusters
            :type clusters: :py:class:`.Clusters`
            :param trajToDistribute: Number of processors to distribute
            :type trajToDistribute: int
            :param currentEpoch: Current iteration number
            :type currentEpoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`

            :returns: list -- List containing the degeneracy of the clusters
        """
        # estimate MSM from clustering object
        self.estimateMSM(clusters.dtrajs, outputPathConstants, currentEpoch)

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration, topologies=None):
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
            :param topology_file: Topology file for non-pdb trajectories
            :type topology_file: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`

            :returns: int, list -- number of processors, list with the
                snapshot from which the trajectories will start in the next iteration
        """
        return self.IndependentCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clustering, iteration, topologies)
