from AdaptivePELE.constants import blockNames

class THRESHOLD_CALCULATOR_TYPES:
    heaviside, constant = list(range(2))

THRESHOLD_CALCULATOR_TYPE_TO_STRING_DICTIONARY = {
    THRESHOLD_CALCULATOR_TYPES.heaviside : blockNames.ThresholdCalculator.heaviside,
    THRESHOLD_CALCULATOR_TYPES.constant : blockNames.ThresholdCalculator.constant
}
