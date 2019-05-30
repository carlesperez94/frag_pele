try:
    # Check if the basestring type if available, this will fail in python3
    basestring
except NameError:
    basestring = str


class ControlFileParams:
    generalParams = "GeneralParams"
    spawningBlockname = "SpawningParams"
    simulationBlockname = "SimulationParams"
    clusteringBlockname = "clusteringTypes"


class GeneralParams:
    mandatory = {
        "restart": "bool",
        "outputPath": "basestring",
        "initialStructures": "list"
    }
    params = {
        "restart": "bool",
        "outputPath": "basestring",
        "initialStructures": "list",
        "debug": "bool",
        "writeAllClusteringStructures": "bool",
        "nativeStructure": "basestring"
    }


class SpawningParams:
    params = {
        "epsilon": "numbers.Real",
        "T": "numbers.Real",
        "reportFilename": "basestring",
        "metricColumnInReport": "numbers.Real",
        "varEpsilonType": "basestring",
        "maxEpsilon": "numbers.Real",
        "minEpsilon": "numbers.Real",
        "variationWindow": "numbers.Real",
        "maxEpsilonWindow": "numbers.Real",
        "period": "numbers.Real",
        "alpha": "numbers.Real",
        "metricWeights": "basestring",
        "metricsInd": "list",
        "condition": "basestring",
        "n": "numbers.Real"
    }
    types = {
        "sameWeight": {
            "reportFilename": "basestring"
        },
        "independent": {
            "reportFilename": "basestring"
        },
        "inverselyProportional": {
            "reportFilename": "basestring"
        },
        "null": {
            "reportFilename": "basestring"
        },
        "epsilon": {
            "epsilon": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
        },
        "FAST": {
            "epsilon": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
        },
        "variableEpsilon": {
            "epsilon": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
            "varEpsilonType": "basestring",
            "maxEpsilon": "numbers.Real"
        },
        "UCB": {
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real"
        },
        "REAP": {
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real"
        }
    }
    density = {
        "types": {
            "heaviside": "basestring",
            "null": "basestring",
            "constant": "basestring",
            "exitContinuous": "basestring",
            "continuous": "basestring"
        },
        "params": {
            "heaviside": "basestring",
            "null": "basestring",
            "constant": "basestring",
            "values": "list",
            "conditions": "list",
            "exitContinuous": "basestring",
            "continuous": "basestring"
        }
    }


class SimulationParams:
    types = {
        "pele": {
            "processors": "numbers.Real",
            "controlFile": "basestring",
            "seed": "numbers.Real",
            "peleSteps": "numbers.Real",
            "iterations": "numbers.Real"
            },
        "test": {
            "destination": "basestring",
            "origin": "basestring",
            "processors": "numbers.Real",
            "seed": "numbers.Real",
            "peleSteps": "numbers.Real",
            "iterations": "numbers.Real"
            },
        "md": {}}
    params = {
        "executable": "basestring",
        "data": "basestring",
        "documents": "basestring",
        "destination": "basestring",
        "origin": "basestring",
        "processors": "numbers.Real",
        "controlFile": "basestring",
        "seed": "numbers.Real",
        "peleSteps": "numbers.Real",
        "iterations": "numbers.Real",
        "modeMovingBox": "basestring",
        "boxCenter": "list",
        "boxRadius": "numbers.Real",
        "runEquilibration": "bool",
        "equilibrationMode": "basestring",
        "equilibrationLength": "numbers.Real",
        "numberEquilibrationStructures": "numbers.Real",
        "useSrun": "bool",
        "exitCondition": "dict"
    }
    exitCondition = {
        "types": {
            "metric": "basestring",
            "clustering": "basestring",
            "metricMultipleTrajectories": "basestring"
        },
        "params": {
            "metricCol": "int",
            "exitValue": "numbers.Real",
            "condition": "basestring",
            "numTrajs": "numbers.Real"
        }
    }


class clusteringTypes:
    types = {
        "rmsd": {
        },
        "contactMap": {
            "similarityEvaluator": "basestring"
        },
        "lastSnapshot": {
        }
    }
    params = {
        "rmsd": "basestring",
        "contactMap": "basestring",
        "contactThresholdDistance": "numbers.Real",
        "lastSnapshot": "basestring",
        "ligandResname": "basestring",
        "ligandResnum": "numbers.Real",
        "ligandChain": "basestring",
        "similarityEvaluator": "basestring",
        "symmetries": "list",
        "alternativeStructure": "bool",
        "nclusters": "numbers.Real"
    }
    thresholdCalculator = {
        "types": {
            "heaviside": "basestring",
            "constant": "basestring"
        },
        "params": {
            "conditions": "list",
            "values": "list",
            "value": "numbers.Real",
            "heaviside": "basestring",
            "constant": "basestring"
        }
    }
