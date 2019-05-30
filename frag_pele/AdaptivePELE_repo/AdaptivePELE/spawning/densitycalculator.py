from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.constants import blockNames
from abc import abstractmethod
from AdaptivePELE.spawning import densitycalculatortypes
import sys


def continousDensity(contacts):
    if contacts > 1.0:
        return 8.0
    else:
        return 64.0/(-4*contacts+6)**3


class DensityCalculatorBuilder():
    def build(self, spawningBlock):
        """
            Build the DensityCalculator object according to the parameters selcted

            :param spawningBlock: Spawning block of the control file
            :type spawningBlock: dict

            :returns: :py:class:`.DensityCalculator` -- DensityCalculator object selected
        """
        try:
            densityBlock = spawningBlock[blockNames.SpawningParams.density]
        except KeyError:
            print("Using null density calculator (no preference for any cluster)")
            return NullDensityCalculator()

        try:
            density_type = densityBlock[blockNames.DensityCalculator.type]
        except KeyError:
            sys.exit("Density calculator must have a type")

        if density_type == blockNames.DensityCalculator.null or density_type == blockNames.DensityCalculator.constant:
            print("Using constant density")
            return NullDensityCalculator()
        elif density_type == blockNames.DensityCalculator.heaviside:
            try:
                paramsBlock = densityBlock[blockNames.DensityCalculator.params]
                values = paramsBlock[blockNames.DensityCalculatorParams.values]
                conditions = paramsBlock[blockNames.DensityCalculatorParams.conditions]
                return DensityCalculatorHeaviside(conditions, values)
            except KeyError:
                print("Using default parameters for Heaviside density calculator")
                return DensityCalculatorHeaviside()
        elif density_type == blockNames.DensityCalculator.continuous:
            print("Using continuous density calculator")
            return ContinuousDensityCalculator()
        elif density_type == blockNames.DensityCalculator.exitContinuous:
            print("Using inverse continuous density calculator")
            return ExitContinousDensityCalculator()
        else:
            sys.exit("Unknown density calculator type! Choices are: " + str(densitycalculatortypes.DENSITY_CALCULATOR_TYPE_TO_STRING_DICTIONARY.values()))


class DensityCalculator():
    def __init__(self):
        self.type = "BaseClass"

    @abstractmethod
    def calculate(self, contacts, contactThreshold):
        pass


class DensityCalculatorHeaviside(DensityCalculator):
    # Mostly duplicated code with threshold calculator
    def __init__(self, conditions=None, values=None):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.heaviside
        if conditions is None:
            conditions = []
        if values is None:
            values = [1.]

        if len(values) != len(conditions) and len(values) != len(conditions) + 1:
            raise ValueError('The number of values must be equal or one more, than the number of conditions')

        self.conditions = conditions
        self.values = values

    def calculate(self, contacts, contactThreshold):
        """
            Calcuate the density value according to the contact ratio

            :param contacts: Contacts ratio
            :type contacts: float
            :param contactThreshold: Deprecated parameter
            :type contactThreshold: float

            :returns: float -- Density value for the value of the contacts ratio
        """
        for i in range(len(self.conditions)):
            # change, so that whole condition is in array
            if contacts > self.conditions[i]:
                # ensure that the value is a float
                return float(self.values[i])
        # the way it's built, it makes more sense to return this value, but,
        # should check that len(value) = len(conditions) + 1 in order to
        # return the "else" value
        return self.values[-1]


class NullDensityCalculator(DensityCalculator):
    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.null

    def calculate(self, contacts, contactThreshold):
        """
            Calcuate the density value according to the contact ratio, in this
            case is always one

            :param contacts: Contacts ratio
            :type contacts: float
            :param contactThreshold: Deprecated parameter
            :type contactThreshold: float

            :returns: float -- Density value for the value of the contacts ratio
        """
        return 1.


class ContinuousDensityCalculator(DensityCalculator):
    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.continuous

    def calculate(self, contacts, contactThreshold):
        """
            Calcuate the density value according to the contact ratio

            :param contacts: Contacts ratio
            :type contacts: float
            :param contactThreshold: Deprecated parameter
            :type contactThreshold: float

            :returns: float -- Density value for the value of the contacts ratio
        """
        return continousDensity(contacts)


class ExitContinousDensityCalculator(DensityCalculator):
    def __init__(self):
        DensityCalculator.__init__(self)
        self.type = densitycalculatortypes.DENSITY_CALCULATOR_TYPES.exitContinous

    def calculate(self, contacts, contactThreshold):
        # alternative possibility
        # if contacts > 1.0:
        #     return 1/4.0
        # elif contacts < 0.24:
        #     return 4.0
        # else:
        #     return (-4*contacts+6)**3/32.0

        return 1.0/continousDensity(contacts)
