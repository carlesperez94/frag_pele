import math
import sys
import numpy as np
import random
from AdaptivePELE.constants import blockNames
import spawningTypes
import densitycalculator
import os
import glob
from AdaptivePELE.constants import constants
from AdaptivePELE.utilities import utilities


def return_sign(i, m, n, r):
    """ Helper function, creates a three-piece step function"""
    if i <= n-m:
        return 1
    elif i <= r:
        return 0
    else:
        return -1


def calculateContactsVar(deltaR, epsMax):
    """ """
    if deltaR < 0.1:
        return 0
    elif deltaR > 1.0:
        return epsMax * 0.09
    else:
        return epsMax * 0.09 * deltaR


class SpawningAlgorithmBuilder:

    def build(self, spawningBlock):
        spawningCalculatorBuilder = SpawningBuilder()
        spawningCalculator = spawningCalculatorBuilder.buildSpawningCalculator(spawningBlock)

        spawningParams = SpawningParams()
        spawningParams.buildSpawningParameters(spawningBlock)

        return spawningCalculator, spawningParams


class SpawningBuilder:

    def buildSpawningCalculator(self, spawningBlock):

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
        self.nclusters = None #number of clusters to consider in epsilon

    def buildSpawningParameters(self, spawningBlock):
        spawningParamsBlock = spawningBlock[blockNames.SpawningParams.params]
        spawningType = spawningBlock[blockNames.StringSpawningTypes.type]

        # Params general to all
        # Params specific to epsilon related spawning
        if spawningType == blockNames.StringSpawningTypes.epsilon or \
                spawningType == blockNames.StringSpawningTypes.variableEpsilon:
            self.epsilon = spawningParamsBlock[blockNames.SpawningParams.epsilon]
            self.nclusters = spawningParamsBlock.get(blockNames.SpawningParams.nclusters, 5)
            self.metricWeights = spawningParamsBlock.get(blockNames.SpawningParams.metricWeights,
                                                         blockNames.SpawningParams.linear)

        if spawningType == blockNames.StringSpawningTypes.epsilon or \
                spawningType == blockNames.StringSpawningTypes.variableEpsilon or\
                spawningType == blockNames.StringSpawningTypes.fast or \
                spawningType == blockNames.StringSpawningTypes.simulatedAnnealing or \
                spawningType == blockNames.StringSpawningTypes.UCB:
            self.temperature = spawningParamsBlock.get(blockNames.SpawningParams.temperature,1000)
            self.reportFilename = spawningParamsBlock[blockNames.SpawningParams.report_filename]
            self.reportCol = spawningParamsBlock[blockNames.SpawningParams.report_col]

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


from abc import ABCMeta, abstractmethod
class SpawningCalculator:
    """
        The purpose of this abstract class is to contain the behaviour of the different strategies for the spawning.
        Spawning is the way in which we split the different explorers at the begining of each epoch.
    """

    def __init__(self):
        self.type = "BaseClass" #change for abstract attribute
        pass

    @abstractmethod
    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        pass

    @abstractmethod
    def log(self):
        pass

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration):
        """ Write initial structures for the next iteriation

            outputPathConstants [In] Output constants that depend on the path
            degeneracyOfRepresentatives [In] Array with the degeneracy of each
            cluster (i.e the number of processors that will be assigned to it)
            clustering [In] clustering object
            iteration [In] Number of epoch
        """
        tmpInitialStructuresTemplate = outputPathConstants.tmpInitialStructuresTemplate
        counts = 0
        procMapping = []
        for i, cluster in enumerate(clustering.clusters.clusters):
            for j in range(int(degeneracyOfRepresentatives[i])):
                outputFilename = tmpInitialStructuresTemplate % (iteration, counts)
                print 'Writing to ', outputFilename, 'cluster', i
                # cluster.writePDB(outputFilename)
                procMapping.append(cluster.writeSpawningStructure(outputFilename))

                counts += 1

        print "counts & cluster centers", counts, np.where(np.array(degeneracyOfRepresentatives) > 0)[0].size
        return counts, procMapping

    def divideTrajAccordingToWeights(self, weights, trajToDistribute):
        """
            weights must be normalized (i.e. sum(weights) = 1)
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
            Divides "trajToDistribute" proportionally to the values in "array"
        """
        if isinstance(array, list): array = np.array(array)
        weights = array/sum(array)
        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def divideInverselyProportionalToArray(self, array, trajToDistribute):
        """
            Divides "trajToDistribute" inversely proportional to the values in "array"
        """
        if isinstance(array, list): array = np.array(array)
        weights = 1./array

        #Handle Nan cases
        weights[weights == np.inf] = 0

        #Handle all Nan cases
        if weights.any():
            weights /= sum(weights)
        else:
            weights[:] = 1./weights.shape[0]

        return self.divideTrajAccordingToWeights(weights, trajToDistribute)

    def getMetrics(self, clusters):
        metrics = np.zeros(len(clusters))
        for i, cluster in enumerate(clusters):
            metrics[i] = cluster.getMetric()
        return metrics

    def getSizes(self, clusters):
        sizes = np.zeros(len(clusters))
        for i, cluster in enumerate(clusters):
            sizes[i] = cluster.elements
        return sizes

from abc import ABCMeta, abstractmethod
class DensitySpawningCalculator(SpawningCalculator):
    """
        Subclass of Spawning calculator, that ensures the definition of a density calculator.
    """

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        SpawningCalculator.__init__(self)
        self.type = "BaseDensityClass"  # change for abstract attribute
        self.densityCalculator = densityCalculator

    def calculateDensities(self, clusters):
        densities = np.zeros(len(clusters))
        for i, cluster in enumerate(clusters):
            contacts = cluster.getContacts()
            cluster.density = self.densityCalculator.calculate(contacts, cluster.contactThreshold)
            densities[i] = cluster.density
        return densities

class IndependentRunsCalculator(SpawningCalculator):

    def __init__(self):
        self.calculator = SameWeightDegeneracyCalculator()
        self.type = spawningTypes.SPAWNING_TYPES.independent

    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch=None):
        return self.calculator.calculate(clusters, trajToDistribute, spawningParams, currentEpoch)

    def log(self):
        pass

    def writeSpawningInitialStructures(self, outputPathConstants, degeneracyOfRepresentatives, clustering, iteration):
        """ Write last trajectory structure as initial one for the next iteriation

            outputPathConstants [In] Output constants that depend on the path
            degeneracyOfRepresentatives [In] Array with the degeneracy of each
            cluster (i.e the number of processors that will be assigned to it)
            clustering [In] clustering object
            iteration [In] Number of epoch
        """
        trajWildcard = os.path.join(outputPathConstants.epochOutputPathTempletized, constants.trajectoryBasename)
        trajectories = glob.glob(trajWildcard%(iteration-1))
        for trajectory in trajectories:
            lastSnapshot = utilities.getSnapshots(trajectory)[-1]

            num = int(trajectory.split("_")[-1][:-4])%len(trajectories)  #to start with 0
            outputFilename = outputPathConstants.tmpInitialStructuresTemplate % (iteration, num)

            f = open(outputFilename, 'w')
            f.write(lastSnapshot)
            f.close()

        return len(trajectories)


class SameWeightDegeneracyCalculator(SpawningCalculator):

    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.sameWeight

    # TODO: Don't back on pele, so that density calculator can be used
    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch=None):
        """
            We back on PELE to split equal trajectories
        """
        numClusters = len(clusters)
        trajToDistribute = min(trajToDistribute, numClusters)
        samples = random.sample(range(numClusters), trajToDistribute)
        degeneracy = [0] * len(clusters)
        for sample in samples:
            degeneracy[sample] = 1
        return degeneracy

    def log(self):
        pass


class InverselyProportionalToPopulationCalculator(DensitySpawningCalculator):

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.inverselyProportional

    def log(self):
        pass

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        sizes = self.getSizes(clusters)
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

    def __init__(self, densityCalculator=densitycalculator.NullDensityCalculator()):
        DensitySpawningCalculator.__init__(self, densityCalculator)
        self.inverselyProportionalCalculator = InverselyProportionalToPopulationCalculator(densityCalculator)
        self.type = spawningTypes.SPAWNING_TYPES.epsilon
        self.degeneracyInverselyProportional = None
        self.degeneracyMetricProportional = None
        self.degeneracyTotal = None

    # TODO add possibility for different pipes
    def log(self):
        if self.degeneracyTotal is not None:
            print "[SpawningLog] Total: %s" % str(self.degeneracyTotal)
        if self.degeneracyInverselyProportional is not None:
            print "[SpawningLog] Inversely prop: %s" % str(self.degeneracyInverselyProportional)
        if self.degeneracyMetricProportional is not None:
            print "[SpawningLog] Metric prop:    %s" % str(self.degeneracyMetricProportional)

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        trajToMetricProportional = int(spawningParams.epsilon * trajToDistribute)
        trajToInverselyProportional = trajToDistribute - trajToMetricProportional

        self.degeneracyInverselyProportional = self.inverselyProportionalCalculator.calculate(clusters, trajToInverselyProportional, spawningParams)
        self.degeneracyMetricProportional = self.divideProcessorsMetricProportional(clusters, trajToMetricProportional, spawningParams)

        self.degeneracyTotal = np.array(self.degeneracyInverselyProportional) + np.array(self.degeneracyMetricProportional)
        return self.degeneracyTotal

    def divideProcessorsMetricProportional(self, clusters, trajToDistribute, spawningParams):

        metrics = self.getMetrics(clusters)

        if isinstance(metrics, list): metrics = np.array(metrics)

        """
            Shift so that differences become larger.
            Also, we can now merge positive & negative values
            Alternatives: Boltzmann weights
        """
        maximumValue = np.max(metrics)
        shiftedMetrics = np.subtract(metrics, maximumValue)
        bestClusters = shiftedMetrics.argsort()

        shiftedMetrics[bestClusters[spawningParams.nclusters:]] = 0 #only consider best ones

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
        if self.maxContacts is None:
            self.maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        maxContacts = reduce(max, [cluster.contacts for cluster in clusters])
        if spawningParams.epsilon < spawningParams.maxEpsilon:
            spawningParams.epsilon += calculateContactsVar(maxContacts-self.maxContacts, spawningParams.maxEpsilon)
        self.maxContacts = maxContacts

    def calculateEpsilonValue(self, spawningParams, currentEpoch, clusters):
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
        epsilon_file = open("epsilon_values.txt", "a")
        epsilon_file.write("%d\t%f\n" % (epoch, epsilon))
        epsilon_file.close()

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        self.calculateEpsilonValue(spawningParams, currentEpoch, clusters)
        self.logVariableEpsilon(spawningParams.epsilon, currentEpoch)
        return self.epsilonDegeneracyCalculator.calculate(clusters, trajToDistribute, spawningParams, currentEpoch)


class SimulatedAnnealingCalculator(SpawningCalculator):

    def __init__(self):
        SpawningCalculator.__init__(self)
        self.type = spawningTypes.SPAWNING_TYPES.simulatedAnnealing

    def log(self):
        pass

    def computeTemperature(self, params, epoch):
        T = params.temperature - params.decrement*epoch
        if T < 300:
            return 300
        else:
            return T

    def calculate(self, clusters, trajToDistribute, spawningParams,
                  currentEpoch):

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
        sizes = self.getSizes(clusters)

        densities = self.calculateDensities(clusters)
        weightedSizes = sizes/densities

        return self.normaliseArray(weightedSizes)

    def calculateNormalisedMetrics(self, clusters):
        metrics = self.getMetrics(clusters)
        return self.normaliseArray(metrics)

    def calculate(self, clusters, trajToDivide, spawningParams, currentEpoch=None):
        normalisedSizes = self.calculateNormalisedSizes(clusters)
        normalisedMetrics = self.calculateNormalisedMetrics(clusters)

        weight = normalisedSizes + 1.*normalisedMetrics

        return self.divideProportionalToArray(weight, trajToDivide)

    def log(self):
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
        pass

    def calculate(self, clusters, trajToDistribute, spawningParams, currentEpoch=None):
        self.epoch += 1
        sizes = np.array(self.getSizes(clusters))
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
            self.prevMetrics = np.pad(self.prevMetrics, (0, n-l), 'constant', constant_values=(0.0))
            self.epoch = np.pad(self.epoch, (0, n-l), 'constant', constant_values=(1.0))
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
