from AdaptivePELE.constants import blockNames


class SIMULATION_TYPE:
    PELE, MD, TEST = list(range(3))

SIMULATION_TYPE_TO_STRING_DICTIONARY = {
    SIMULATION_TYPE.PELE: blockNames.SimulationType.pele,
    SIMULATION_TYPE.MD: blockNames.SimulationType.md,
    SIMULATION_TYPE.TEST: blockNames.SimulationType.test
}


class EXITCONDITION_TYPE:
    METRIC, METRICMULTIPLETRAJS, CLUSTERING = list(range(3))

EXITCONDITION_TYPE_TO_STRING_DICTIONARY = {
    EXITCONDITION_TYPE.METRIC: blockNames.ExitConditionType.metric,
    EXITCONDITION_TYPE.METRICMULTIPLETRAJS: blockNames.ExitConditionType.metricMultipleTrajs,
    EXITCONDITION_TYPE.CLUSTERING: blockNames.ExitConditionType.clustering
}
