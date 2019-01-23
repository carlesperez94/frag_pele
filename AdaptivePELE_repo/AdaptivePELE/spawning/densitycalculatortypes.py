from AdaptivePELE.constants import blockNames


class DENSITY_CALCULATOR_TYPES:
    null, heaviside, continuous, exitContinous = list(range(4))

DENSITY_CALCULATOR_TYPE_TO_STRING_DICTIONARY = {
    DENSITY_CALCULATOR_TYPES.heaviside: blockNames.DensityCalculator.heaviside,
    DENSITY_CALCULATOR_TYPES.null: blockNames.DensityCalculator.null,
    DENSITY_CALCULATOR_TYPES.continuous: blockNames.DensityCalculator.continuous,
    DENSITY_CALCULATOR_TYPES.exitContinous: blockNames.DensityCalculator.exitContinuous
}
