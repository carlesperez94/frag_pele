class ControlFileParams:
    generalParams = "GeneralParams"
    spawningBlockname = "SpawningParams"
    simulationBlockname = "SimulationParams"
    clusteringBlockname = "clusteringTypes"


class GeneralParams:
    mandatory = {
        "restart": "bool",
        "outputPath": "basestring",
        "initialStructures": "list",
        "debug": "bool",
        "writeAllClusteringStructures": "bool"
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
        "n": "numbers.Real"
    }
    types = {
        "sameWeight": {},
        "independent": {},
        "inverselyProportional": {},
        "epsilon": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
        },
        "FAST": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
        },
        "variableEpsilon": {
            "epsilon": "numbers.Real",
            "T": "numbers.Real",
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real",
            "varEpsilonType": "basestring",
            "maxEpsilon": "numbers.Real"
        },
        "UCB": {
            "reportFilename": "basestring",
            "metricColumnInReport": "numbers.Real"
        }
    }
    density = {
        "types": {
            "heaviside": "basestring",
            "null": "basestring",
            "constant": "basestring",
            "continuous": "basestring"
        },
        "params": {
            "heaviside": "basestring",
            "null": "basestring",
            "constant": "basestring",
            "values": "list",
            "conditions": "list",
            "continuous": "basestring"
        }
    }


class SimulationParams:
    types = {"pele": {
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
        "exitCondition": "dict"
    }
    exitCondition = {
        "types": {
            "metric": "basestring"
        },
        "params": {
            "metricCol": "int",
            "exitValue": "numbers.Real"
        }
    }


class clusteringTypes:
    types = {
        "rmsd": {
            "contactThresholdDistance": "numbers.Real",
        },
        "contactMap": {
            "contactThresholdDistance": "numbers.Real",
            "similarityEvaluator": "basestring"
        },
        "lastSnapshot": {
            "contactThresholdDistance": "numbers.Real",
        }
    }
    params = {
        "rmsd": "basestring",
        "contactMap": "basestring",
        "lastSnapshot": "basestring",
        "ligandResname": "basestring",
        "ligandResnum": "numbers.Real",
        "ligandChain": "basestring",
        "contactThresholdDistance": "numbers.Real",
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
