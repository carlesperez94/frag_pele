from AdaptivePELE.constants import blockNames
import densitycalculatortypes
import sys


class DensityCalculatorBuilder():
    def build(self, spawningBlock):
        try:
            densityBlock = spawningBlock[blockNames.SpawningParams.density]
        except KeyError:
            print "Using null density calculator (no preference for any cluster)"
            return NullDensityCalculator()

        try:
            type = densityBlock[blockNames.DensityCalculator.type]
        except KeyError:
            sys.exit("Density calculator must have a type")

        if type == blockNames.DensityCalculator.null or type == blockNames.DensityCalculator.constant:
            print "Using constant density"
            return NullDensityCalculator()
        elif type == blockNames.DensityCalculator.heaviside:
            try:
                paramsBlock = densityBlock[blockNames.DensityCalculator.params]
                values = paramsBlock[blockNames.DensityCalculatorParams.values]
                conditions = paramsBlock[blockNames.DensityCalculatorParams.conditions]
                return DensityCalculatorHeaviside(conditions, values)
            except KeyError:
                print "Using default parameters for Heaviside density calculator"
                return DensityCalculatorHeaviside()
        elif type == blockNames.DensityCalculator.continuous:
            return ContinuousDensityCalculator()
        else:
            sys.exit("Unknown density calculator type! Choices are: " + str(densitycalculatortypes.DENSITY_CALCULATOR_TYPE_TO_STRING_DICTIONARY.values()))

from abc import ABCMeta, abstractmethod
class DensityCalculator():
    def __init__(self):
        self.type = "BaseClass"

    @abstractmethod
    def calculate(self, contacts, contactThreshold):
        pass


class DensityCalculatorHeaviside(DensityCalculator):
    # Mostly duplicated code with threshold calculator
    def __init__(self, conditions=[], values=[1.]):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside

        if len(values) != len(conditions) and len(values) != len(conditions) + 1:
            raise ValueError('The number of values must be equal or one more, than the number of conditions')

        self.conditions = conditions
        self.values = values

    def calculate(self, contacts, contactThreshold):
        for i in range(len(self.conditions)):
            # change, so that whole condition is in array
            if contacts > self.conditions[i]:
                return self.values[i]
        # the way it's built, it makes more sense to return this value, but,
        # should check that len(value) = len(conditions) + 1 in order to
        # return the "else" value
        return self.values[-1]


class NullDensityCalculator(DensityCalculator):
    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null

    def calculate(self, contacts, contactThreshold):
        return 1.


class ContinuousDensityCalculator(DensityCalculator):
    limit = {8: 1.0, 6: 0.4, 4: 0.1, 10: 4.0}
    slope = {8: 4.0, 6: 10.0, 4: 40.0, 10: 2.0}

    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.continuous

    def calculate(self, contacts, contactThreshold):
        # if contacts > 2.0:
        #     return 32
        # else:
        #     return 8.76571*contacts**2-2.44857*contacts+0.28829
        # if contacts > self.limit[contactThreshold]:
        #     return 8
        # else:
        #     return 64.0/(-self.slope[contactThreshold]*contacts+6)**3
        if contacts > 1.0:
            return 8
        else:
            return 64.0/(-4*contacts+6)**3
