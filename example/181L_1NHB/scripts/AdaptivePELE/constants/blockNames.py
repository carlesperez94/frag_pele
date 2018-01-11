class ClusteringTypes:
    type = "type"
    params = "params"
    rmsd = "rmsd"
    contactMap = "contactMap"
    lastSnapshot = "lastSnapshot"
    thresholdCalculator = "thresholdCalculator"
    ligandResname = "ligandResname"
    ligandResnum = "ligandResnum"
    ligandChain = "ligandChain"
    alternativeStructure = "alternativeStructure"
    contactThresholdDistance = "contactThresholdDistance"
    nclusters = "nclusters"
    similarityEvaluator = "similarityEvaluator"
    differenceDistance = "differenceDistance"
    Jaccard = "Jaccard"
    correlation = "correlation"
    symmetries = "symmetries"


class ThresholdCalculator:
    type = "type"
    params = "params"
    heaviside = "heaviside"
    constant = "constant"


class ThresholdCalculatorParams:
    conditions = "conditions"
    values = "values"
    value = "value"


class DensityCalculator:
    type = "type"
    params = "params"
    heaviside = "heaviside"
    null = "null"
    constant = "constant"
    continuous = "continuous"


class DensityCalculatorParams:
    conditions = "conditions"
    values = "values"


class StringSpawningTypes:
    type = "type"
    independent = "independent"
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    fast = "FAST"
    simulatedAnnealing = "simulatedAnnealing"
    # New parameters for variable epsilon(experimental)
    variableEpsilon = "variableEpsilon"
    UCB = "UCB"


class SpawningParams:
    params = "params"
    epsilon = "epsilon"
    temperature = "T"
    threshold = "threshold"
    report_filename = "reportFilename"
    report_col = "metricColumnInReport"
    # New parameters for variable epsilon(experimental)
    varEpsilonType = "varEpsilonType"
    maxEpsilon = "maxEpsilon"
    minEpsilon = "minEpsilon"
    variationWindow = "variationWindow"  # Last epoch of variable epsilon,if
    # current epoch > than variation Window, set epsilon to minEpsilon
    maxEpsilonWindow = "maxEpsilonWindow"
    period = "period"  # Only useful for periodic epsilon modes
    density = "density"
    metricWeights = "metricWeights"
    linear = "linear"
    boltzmann = "boltzmann"
    alpha = "alpha"
    nclusters = "n"


class SpawningDensity:
    values = "values"
    conditions = "conditions"


class VariableEpsilonTypes:
    linearVariation = "linearVariation"
    contactsVariation = "contactsVariation"


class SimulationType:
    type = "type"
    pele = "pele"
    md = "md"
    test = "test"


class SimulationParams:
    params = "params"
    processors = "processors"
    executable = "executable"
    templetizedControlFile = "controlFile"
    dataFolder = "data"
    documentsFolder = "documents"
    destination = "destination"
    origin = "origin"
    seed = "seed"
    peleSteps = "peleSteps"
    iterations = "iterations"
    exitCondition = "exitCondition"
    metricCol = "metricCol"
    exitValue = "exitValue"
    trajectories = "trajectories"


class ExitConditionType:
    type = "type"
    metric = "metric"
    clustering = "clustering"


class ControlFileParams:
    generalParams = "generalParams"
    spawningBlockname = "spawning"
    simulationBlockname = "simulation"
    clusteringBlockname = "clustering"


class GeneralParams:
    restart = "restart"
    outputPath = "outputPath"
    initialStructures = "initialStructures"
    debug = "debug"
    writeAllClustering = "writeAllClusteringStructures"
    nativeStructure = "nativeStructure"
