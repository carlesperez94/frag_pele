from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import numpy as np
from AdaptivePELE.constants import blockNames
from AdaptivePELE.clustering import thresholdcalculatortypes
from abc import abstractmethod


class ThresholdCalculatorBuilder():
    def build(self, clusteringBlock):
        """
            Bulid the selecte thresholdCaulcualtor object

            :param clusteringBlock: Parameters block corresponding to the threshold calculator
            :type clusteringBlock: dict
            :returns: :py:class:`.ThresholdCalculator` -- thresholdCalculator object selected
        """
        try:
            thresholdCalculatorBlock = clusteringBlock[blockNames.ClusteringTypes.thresholdCalculator]
        except KeyError:
            # Default value if no threshold calculator block was defined
            return ThresholdCalculatorHeaviside()

        try:
            typeParam = thresholdCalculatorBlock[blockNames.ThresholdCalculator.type]
        except KeyError:
            sys.exit("Threshold calculator must have a type")

        if typeParam == blockNames.ThresholdCalculator.constant:
            try:
                paramsBlock = thresholdCalculatorBlock[blockNames.ThresholdCalculator.params]
                value = paramsBlock[blockNames.ThresholdCalculatorParams.value]
                return ThresholdCalculatorConstant(value)
            except KeyError:
                print("Using default parameters for constant threshold calculator")
                return ThresholdCalculatorConstant()
        elif typeParam == blockNames.ThresholdCalculator.heaviside:
            try:
                paramsBlock = thresholdCalculatorBlock[blockNames.ThresholdCalculator.params]
                values = paramsBlock[blockNames.ThresholdCalculatorParams.values]
                conditions = paramsBlock[blockNames.ThresholdCalculatorParams.conditions]
                return ThresholdCalculatorHeaviside(conditions, values)
            except KeyError:
                print("Using default parameters for Heaviside threshold calculator")
                return ThresholdCalculatorHeaviside()
        else:
            sys.exit("Unknown threshold calculator type! Choices are: " + str(thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPE_TO_STRING_DICTIONARY.values()))


class ThresholdCalculator():
    def __init__(self):
        self.type = "BaseClass"  # change for abstract attribute

    @abstractmethod
    def calculate(self, contacts):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    def __ne__(self, other):
        return not self.__eq__(other)


class ThresholdCalculatorConstant(ThresholdCalculator):
    def __init__(self, value=2):
        ThresholdCalculator.__init__(self)
        self.type = thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPES.constant
        self.value = value

    def calculate(self, contacts):
        """
            Calculate the threshold value of a cluster. In this case it is constant,
            the contacts ratio is only passed for compatibility purposes

            :param contacts: Contact ratio
            :type contacts: float
            :returns: float -- threshold value of the cluster
        """
        return self.value

    def __eq__(self, other):
        return self.type == other.type and self.value == other.value

    def getMaxThreshold(self):
        """
            Method that returns the maximum treshold possible, required for new
            distance-ordered clustering(in early development)

            :returns: float -- Maximum threshold possible
        """
        return self.value


class ThresholdCalculatorHeaviside(ThresholdCalculator):
    def __init__(self, conditions=None, values=None):
        ThresholdCalculator.__init__(self)
        self.type = thresholdcalculatortypes.THRESHOLD_CALCULATOR_TYPES.heaviside
        if conditions is None:
            conditions = [1.0, 0.75, 0.5]
        if values is None:
            values = [2, 3, 4, 5.0]

        if len(values) != len(conditions) and len(values) != len(conditions) + 1:
            raise ValueError('The number of values must be equal or one more, than the number of conditions')

        self.conditions = conditions
        self.values = values

    def calculate(self, contacts):
        """
            Calculate the threshold value of a cluster according to the contacts ratio
            and the selected conditions and values

            :param contacts: Contact ratio
            :type contacts: float
            :returns: float -- threshold value of the cluster
        """
        for i in range(len(self.conditions)):
            # change, so that whole condition is in array
            if contacts > self.conditions[i]:
                return self.values[i]
        # the way it's built, it makes more sense to return this value, but, should check that len(value) = len(conditions) + 1 in order to return the "else" value
        return self.values[-1]

    def getMaxThreshold(self):
        """
            Method that returns the maximum treshold possible, required for new
            distance-ordered clustering(in early development)

            :returns: float -- Maximum threshold possible
        """
        return max(self.values)

    def __eq__(self, other):
        return self.type == other.type and\
            np.allclose(self.conditions, other.conditions) and\
            np.allclose(self.values, other.values)
