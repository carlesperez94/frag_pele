class ClusteringTypes:
    type = "type"
    params = "params"
    rmsd = "rmsd"
    contactMap = "contactMap"
    lastSnapshot = "lastSnapshot"
    null = "null"
    MSMClustering = "MSM"
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
    tica = "tica"
    atom_Ids = "atom_Ids"
    writeCA = "writeCA"
    sidechains = "sidechains"
    tica_lagtime = "tica_lagtime"
    tica_nICs = "tica_nICs"
    tica_kinetic_map = "tica_kinetic_map"
    tica_commute_map = "tica_commute_map"


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
    exitContinuous = "exitContinuous"


class DensityCalculatorParams:
    conditions = "conditions"
    values = "values"


class StringSpawningTypes:
    type = "type"
    independent = "independent"
    independentMetric = "independentMetric"
    sameWeight = "sameWeight"
    inverselyProportional = "inverselyProportional"
    epsilon = "epsilon"
    fast = "FAST"
    simulatedAnnealing = "simulatedAnnealing"
    variableEpsilon = "variableEpsilon"
    UCB = "UCB"
    REAP = "REAP"
    null = "null"
    ProbabilityMSMCalculator = "ProbabilityMSM"
    MetastabilityMSMCalculator = "MetastabilityMSM"
    UncertaintyMSMCalculator = "UncertaintyMSM"
    IndependentMSMCalculator = "IndependentMSM"


class SpawningParams:
    params = "params"
    epsilon = "epsilon"
    temperature = "T"
    threshold = "threshold"
    report_filename = "reportFilename"
    report_col = "metricColumnInReport"
    minValue = "min"
    maxValue = "max"
    condition = "condition"
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
    metricsInd = "metricsInd"
    lagtime = "lagtime"
    minPos = "minPos"
    SASA_column = "SASA_column"
    filterByMetric = "filterByMetric"
    filter_value = "filter_value"
    filter_col = "filter_col"


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
    time = "time"
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
    modeMovingBox = "modeMovingBox"
    modeMovingBoxBinding = "binding"
    modeMovingBoxUnBinding = "unbinding"
    equilibrationMode = "equilibrationMode"
    equilibrationLastSnapshot = "equilibrationLastSnapshot"
    equilibrationSelect = "equilibrationSelect"
    equilibrationCluster = "equilibrationCluster"
    numberEquilibrationStructures = "numberEquilibrationStructures"
    boxCenter = "boxCenter"
    boxRadius = "boxRadius"
    runEquilibration = "runEquilibration"
    condition = "condition"
    numTrajs = "numberTrajectories"
    equilibrationLength = "equilibrationLength"
    trajectoryName = "trajectoryName"
    srun = "useSrun"
    srunParameters = "srunParameters"
    mpiParameters = "mpiParameters"
    # params for MD
    ligandCharge = "ligandCharge"
    nonBondedCutoff = "nonBondedCutoff"
    Temperature = "temperature"
    runningPlatform = "runningPlatform"
    minimizationIterations = "minimizationIterations"
    repoterfreq = "reporterFrequency"
    productionLength = "productionLength"
    waterBoxSize = "WaterBoxSize"
    trajsPerReplica = "trajectoriesPerReplica"
    numReplicas = "numReplicas"
    timeStep = "timeStep"
    equilibrationLengthNVT = "equilibrationLengthNVT"
    equilibrationLengthNPT = "equilibrationLengthNPT"
    constraintsMin = "constraintsMinimization"
    constraintsNVT = "constraintsNVT"
    constraintsNPT = "constraintsNPT"
    devicesPerTrajectory = "devicesPerTrajectory"
    forcefield = "forcefield"
    customparamspath = "customparamspath"
    maxDevicesPerReplica = "maxDevicesPerReplica"
    format = "format"
    ligandName = "ligandName"
    constraints = "constraints"
    boxType = "boxType"
    cylinder = "cylinder"
    cylinderBases = "cylinderBases"
    sphere = "sphere"


class ExitConditionType:
    type = "type"
    metric = "metric"
    clustering = "clustering"
    metricMultipleTrajs = "metricMultipleTrajectories"


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
