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
        "n": "numbers.Real",
        "lagtime": "numbers.Real",
        "minPos": "list",
        "SASA_column": "int",
        "filterByMetric": "bool",
        "filter_value": "numbers.Real",
        "filter_col": "int"
    }
    types = {
        "sameWeight": {
            "reportFilename": "basestring"
        },
        "independent": {
            "reportFilename": "basestring"
        },
        "independentMetric": {
            "metricColumnInReport": "numbers.Real",
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
        },
        "ProbabilityMSM": {
            "lagtime": "numbers.Real"
        },
        "MetastabilityMSM": {
            "lagtime": "numbers.Real"
        },
        "UncertaintyMSM": {
            "lagtime": "numbers.Real"
        },
        "IndependentMSM": {
            "lagtime": "numbers.Real"
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
        "md": {
            "processors": "numbers.Real",
            "seed": "numbers.Real",
            "productionLength": "numbers.Real",
            "iterations": "numbers.Real",
            "numReplicas": "numbers.Real"
        }}
    params = {
        "executable": "basestring",
        "data": "basestring",
        "documents": "basestring",
        "destination": "basestring",
        "origin": "basestring",
        "time": "numbers.Real",
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
        "srunParameters": "basestring",
        "mpiParameters": "basestring",
        "exitCondition": "dict",
        "trajectoryName": "basestring",
        "ligandCharge": "numbers.Real",
        "ligandName": "basestring",
        "nonBondedCutoff": "numbers.Real",
        "timeStep": "numbers.Real",
        "Temperature": "numbers.Real",
        "runningPlatform": "basestring",
        "minimizationIterations": "numbers.Real",
        "reporterFrequency": "numbers.Real",
        "productionLength": "numbers.Real",
        "WaterBoxSize": "numbers.Real",
        "forcefield": "basestring",
        "trajectoriesPerReplica": "numbers.Real",
        "equilibrationLengthNVT": "numbers.Real",
        "equilibrationLengthNPT": "numbers.Real",
        "devicesPerTrajectory": "int",
        "constraintsMinimization": "numbers.Real",
        "constraintsNVT": "numbers.Real",
        "constraintsNPT": "numbers.Real",
        "customparamspath": "basestring",
        "numReplicas": "numbers.Real",
        "maxDevicesPerReplica": "numbers.Real",
        "format": "basestring",
        "constraints": "list",
        "boxType": "basestring",
        "cylinderBases": "list"
    }
    exitCondition = {
        "types": {
            "metric": "basestring",
            "clustering": "basestring",
            "metricMultipleTrajectories": "basestring"
        },
        "params": {
            "metricCol": "numbers.Real",
            "exitValue": "numbers.Real",
            "condition": "basestring",
            "numTrajs": "numbers.Real"
        }
    }


class clusteringTypes:
    types = {
        "rmsd": {},
        "contactMap": {
            "similarityEvaluator": "basestring",
            "ligandResname": "basestring"
        },
        "lastSnapshot": {
            "ligandResname": "basestring"
        },
        "null": {},
        "MSM": {
            "ligandResname": "basestring",
            "nclusters": "numbers.Real"
        }
    }
    params = {
        "rmsd": "basestring",
        "contactMap": "basestring",
        "lastSnapshot": "basestring",
        "null": "basestring",
        "contactThresholdDistance": "numbers.Real",
        "ligandResname": "basestring",
        "ligandResnum": "numbers.Real",
        "ligandChain": "basestring",
        "similarityEvaluator": "basestring",
        "symmetries": "list",
        "alternativeStructure": "bool",
        "nclusters": "numbers.Real",
        "tica": "bool",
        "atom_Ids": "list",
        "writeCA": "bool",
        "sidechains": "bool",
        "tica_lagtime": "numbers.Real",
        "tica_nICs": "numbers.Real",
        "tica_kinetic_map": "bool",
        "tica_commute_map": "bool"
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
