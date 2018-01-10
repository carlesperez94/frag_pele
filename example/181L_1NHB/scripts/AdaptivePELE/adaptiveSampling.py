import numpy as np
import sys
import shutil
import os
import json
import time
import glob
import argparse
#import networkx as nx
import atexit
from AdaptivePELE.constants import blockNames, constants
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities
from AdaptivePELE.validator import controlFileValidator
from AdaptivePELE.spawning import spawning, spawningTypes
from AdaptivePELE.simulation import simulationrunner, simulationTypes
from AdaptivePELE.clustering import clustering, clusteringTypes


def parseArgs():
    parser = argparse.ArgumentParser(description="Perform several iterations of"
                                     " simulations using adaptive sampling to "
                                     "distribute the processors in order "
                                     "to optimize sampling")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    return args


def expandInitialStructuresWildcard(initialStructuresWildcard):
    """
        Returns the initial structures after expanding the initial structures wildcard

        :param initialStructureWildcard: Wildcard that matches the initial structures
        :type initialStructureWildcard: str

        :return totalInitialStructures: The expanded initial structures
        :rtype: list of str
    """
    totalInitialStructures = []
    for initialStructureWildcard in initialStructuresWildcard:
        expandedStructures = glob.glob(initialStructureWildcard)
        totalInitialStructures.extend(expandedStructures)
    return totalInitialStructures


def checkSymmetryDict(clusteringBlock, initialStructures, resname):
    """
        Check if the symmetries dictionary is valid for the ligand

        :param clusteringBlock: JSON block with the clustering-related parameters
        :type clusteringBlock: json
        :param initialStructures: List with initial structures
        :type initialStructures: list
        :param resname: Residue name of the ligand in the system pdb
        :type resname: str

        :raise AssertionError: If atoms are not found in the structure
     """
    symmetries = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.symmetries, {})
    for structure in initialStructures:
        PDB = atomset.PDB()
        PDB.initialise(str(structure), resname=resname)
        utilities.assertSymmetriesDict(symmetries, PDB)


def fixReportsSymmetry(outputPath, resname, nativeStructure, symmetries):
    """
        Adds a new column in the report file with the RMSD that takes into account symmetries.
        New reports are stored in the fixedReport_i where i is the number of the report

        :param outputPath: Path where trajectories are found
        :type outputPath: str
        :param resname: Residue name of the ligand in the pdb
        :type resname: str
        :param nativeStructure: Path to the native structure pdb
        :type nativeStructure: str
        :param symmetries: Dictionary containg the symmetries of the ligand
        :type symmetries: dict

        :raise IndexError: If original report file not found in output folder
    """
    outputFilename = "fixedReport_%d" #move to constants?
    trajName = "*traj*.pdb" #move to constants?
    reportName = "*report_%d" #move to constants?
    trajs = glob.glob(os.path.join(outputPath, trajName))
    nativePDB = atomset.PDB()
    nativePDB.initialise(str(nativeStructure), resname=resname)
    for traj in trajs:
        trajNum = utilities.getTrajNum(traj)
        rmsd = list(utilities.getRMSD(traj, nativePDB, resname, symmetries))
        try:
            reportFilename = glob.glob(os.path.join(outputPath, reportName) % trajNum)[0]
        except IndexError:
            raise IndexError("File %s not found in folder %s" % (reportName % trajNum, outputPath))
        rmsd.insert(0, "\tCorrected RMSD")
        with open(reportFilename, "r") as f:
            report = f.readlines()
        outfile = open(os.path.join(outputPath, outputFilename % trajNum), "w")
        for line, value in zip(report, rmsd):
            outfile.write(line.rstrip("\n")+str(value)+"\n")
        outfile.close()


def copyInitialStructures(initialStructures, tmpInitialStructuresTemplate, iteration):
    """
        Copies the initial structures from a certain iteration

        :param initialStructures: Name of the initial structures to copy
        :type initialStructures: list of str
        :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
        :type tmpInitialStructuresTemplate: str
        :param iteration: Epoch number
        :type iteration: int
    """

    for i, name in enumerate(initialStructures):
        shutil.copyfile(name, tmpInitialStructuresTemplate % (iteration, i))


def createMultipleComplexesFilenames(numberOfSnapshots, tmpInitialStructuresTemplate, iteration):
    """
        Creates the string to substitute the complexes in the PELE control file

        :param numberOfSnapshots: Number of complexes to write
        :type numberOfSnapshots: int
        :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
        :type tmpInitialStructuresTemplate: str
        :param iteration: Epoch number
        :type iteration: int

        :returns: jsonString to be substituted in PELE control file
        :rtype: str
    """
    jsonString = "\n"
    for i in range(numberOfSnapshots-1):
        jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, i)) + ",\n"
    jsonString += constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1))
    return jsonString


def generateTrajectorySelectionString(epoch, epochOutputPathTempletized):
    """
        Generates the template for the name of the trajectories in a given epoch
        :param Poch: Epoch number
        :type epoch: int
        :param epochOutputPathTempletized: Templetized path where the trajectories of any epoch are stored
        :type epochOutputPathTempletized: str
    """
    return "[\"" + os.path.join(epochOutputPathTempletized % epoch, constants.trajectoryBasename) + "\"]"


def findFirstRun(outputPath, clusteringOutputObject):
    """
        Find the last epoch that was properly simulated and clusterized and
        and return the first epoch to run in case of restart

        :param outputPath: Simulation output path
        :type outputPath: str
        :param clusteringOutputObject: Templetized name of the clustering object
        :type clusteringOutputObject: str

        :return: Current epoch
        :rtype: int 
    """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    objectsFound = []
    for epoch in epochFolders:
        if os.path.exists(clusteringOutputObject % epoch):
            objectsFound.append(epoch)
        if objectsFound and epoch < (objectsFound[0]-5):
            break
    while objectsFound:
        epoch = objectsFound.pop(0)
        if checkIntegrityClusteringObject(clusteringOutputObject % epoch):
            return epoch + 1
    return 0


def checkIntegrityClusteringObject(objectPath):
    """
        Test whether the found clustering object to reload is a valid object

        :param objectPath: Clustering object path
        :type objectPath: str

        :returns: True if the found clustering object is valid
        :rtype: bool
    """
    try:
        utilities.readClusteringObject(objectPath)
        return True
    except EOFError:
        return False


def loadParams(jsonParams):
    """
        Read the control file in JSON format and extract the blocks of simulation,
        general parameters, clustering and spawning

        :param jsonParams: Control file in JSON format from where the parameters will be read
        :type jsonParams: json str
    """
    jsonFile = open(jsonParams, 'r').read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
        parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]


def saveInitialControlFile(jsonParams, originalControlFile):
    """ Save the control file jsonParams in originalControlFile

        :param jsonParams: Input control file in JSON format
        :type jsonParams: str
        :param originalControlFile: Path where to save the control file
        :type originalControlFile: str
    """
    file = open(originalControlFile, 'w')
    jsonFile = open(jsonParams, 'r').read()
    file.write(jsonFile)


def needToRecluster(oldClusteringMethod, newClusteringMethod):
    """ 
        Check if the parameters have changed in a restart and we need to redo
        the clustering. In particular: type of clustering, theshold calculator
        or distance

        :param oldClusteringMethod: Clustering in a previous simulation before the restart
        :type oldClusteringMethod: :py:class:`.Clustering`
        :param newClusteringMethod: Clustering in the restarted simulation
        :type newClusteringMethod: :py:class:`.Clustering`

        :returns: If clustering needs to be redone
        :rtype: bool
    """

    # Check 1: change of type
    if oldClusteringMethod.type != newClusteringMethod.type:
        return True

    # Check 2: Change of thresholdCalculator and thresholdDistance
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.rmsd or\
       oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMap:
        return oldClusteringMethod.thresholdCalculator != newClusteringMethod.thresholdCalculator or\
                abs(oldClusteringMethod.contactThresholdDistance - newClusteringMethod.contactThresholdDistance) > 1e-7

    #TODO: add similarityEvaluator check for contactMap clustering


def clusterEpochTrajs(clusteringMethod, epoch, epochOutputPathTempletized):
    """
        Cluster the trajecotories of a given epoch

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering` 
        :param epoch: Number of the epoch to cluster
        :type epoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
"""

    snapshotsJSONSelectionString = generateTrajectorySelectionString(epoch, epochOutputPathTempletized)
    paths = eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No trajectories to cluster! Matching path:%s"%paths[-1])
    clusteringMethod.cluster(paths)


def clusterPreviousEpochs(clusteringMethod, finalEpoch, epochOutputPathTempletized, simulationRunner):
    """
        Cluster all previous epochs using the clusteringMethod object

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering`
        :param finalEpoch: Last epoch to cluster (not included)
        :type finalEpoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
        :param simulationRunner: Simulation runner object
"""
    for i in range(finalEpoch):
        simulationRunner.readMappingFromDisk(epochOutputPathTempletized % i)
        clusterEpochTrajs(clusteringMethod, i, epochOutputPathTempletized)


def getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner):
    """
        It reads the previous clustering method, and, if there are changes,
        it reclusters the previous trajectories. Returns the clustering object to use

        :param firstRun: New epoch to run
        :type firstRun: int
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params, to know what reportFile and column to read
        :type spawningParams: :py:class:`.SpawningParams`

        :returns clusteringMethod: The clustering method to use in the adaptive sampling simulation
        :rtype: :py:class:`.Clustering`
    """

    lastClusteringEpoch = firstRun - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % (lastClusteringEpoch)
    oldClusteringMethod = utilities.readClusteringObject(clusteringObjectPath)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)

    if needToRecluster(oldClusteringMethod, clusteringMethod):
        print "Reclustering!"
        startTime = time.time()
        clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)
        endTime = time.time()
        print "Reclustering took %s sec" % (endTime - startTime)
    else:
        clusteringMethod = oldClusteringMethod
        clusteringMethod.setCol(spawningParams.reportCol)

    return clusteringMethod


def buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock,
                                                         spawningParams, spawningCalculator, simulationRunner):
    """
        It reads the previous clustering method, and if there are changes (clustering method or related to thresholds),
        reclusters the previous trajectories. Returns the clustering object to use,
        and the initial structure filenames as strings

        :param firstRun: New epoch to run
        :type firstRun: int
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: str
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params
        :type spawningParams: :py:class:`.SpawningParams`
        :param spawningCalculator: Spawning calculator object
        :type spawningCalculator: :py:class:`.SpawningCalculator`
        :param simulationRunner: :py:class:`.SimulationRunner` Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`

        :returns: (clusteringMethod, initialStructuresAsString): 

        **clusteringMethod** (:py:class:`.Clustering`) The clustering method to use in the adaptive sampling simulation
        **initialStructuresAsString** (str) The initial structures filenames
    """

    clusteringMethod = getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner)

    if not hasattr(clusteringMethod, "conformationNetwork"):
        clusteringMethod.epoch = firstRun

    degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, firstRun)
    spawningCalculator.log()
    print "Degeneracy", degeneracyOfRepresentatives
    seedingPoints, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, firstRun)
    initialStructuresAsString = createMultipleComplexesFilenames(seedingPoints, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    simulationRunner.updateMappingProcessors(procMapping)

    return clusteringMethod, initialStructuresAsString


def buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, outputPath, controlFile, outputPathConstants, clusteringBlock, spawningParams, initialStructures):
    """
        Build the clustering object and copies initial structures from control file.
        Returns the clustering object to use and the initial structures filenames as string

        :param debug: In debug, it will not remove the simulations
        :type debug: bool
        :param outputPath: Simulation output path
        :type outputPath: str
        :param controlFile: Adaptive sampling control file
        :type controlFile: str
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
        :param clusteringBlock: Contains the new clustering block
        :type clusteringBlock: json
        :param spawningParams: Spawning params
        :type spawningParams: :py:class:`.SpawningParams`
        :param initialStructures: Control file initial structures
        :type initialStructures: list

        **clusteringMethod** (:py:class:`.Clustering`) The clustering method to use in the adaptive sampling simulation
        **initialStructuresAsString** (str) The initial structures filenames
    """
    if not debug:
        shutil.rmtree(outputPath)
    utilities.makeFolder(outputPath)
    saveInitialControlFile(controlFile, outputPathConstants.originalControlFile)

    firstRun = 0
    copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    initialStructuresAsString = createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)
    initialClusters = []
    return clusteringMethod, initialStructuresAsString, initialClusters


def preparePeleControlFile(epoch, outputPathConstants, simulationRunner, peleControlFileDictionary):
    """
        Substitute the parameters in the PELE control file specified with the
        provided in the control file

        :param epoch: Epoch number
        :type epoch: int
        :param outputPathConstants: Object that has as attributes constant related to the outputPath that will be used to create the working control file
        :type outputPathConstants: :py:class:`.OutputPathConstants`
        :param simulationRunner: Simulation runner object
        :type simulationRunner: `.SimulationRunner`
        :param peleControlFileDictionary: Dictonary containing the values of the parameters to substitute in the control file
        :type peleControlFileDictionary: dict
    """
    outputDir = outputPathConstants.epochOutputPathTempletized % epoch
    utilities.makeFolder(outputDir)
    peleControlFileDictionary["OUTPUT_PATH"] = outputDir
    peleControlFileDictionary["SEED"] = simulationRunner.parameters.seed + epoch*simulationRunner.parameters.processors
    simulationRunner.makeWorkingControlFile(outputPathConstants.tmpControlFilename % epoch, peleControlFileDictionary)


def main(jsonParams):
    """
        Main body of the adaptive sampling program.

        :params jsonParams: A string with the name of the control file to use
        :type jsonParams: str
    """

    controlFileValidator.validate(jsonParams)
    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    spawningCalculator, spawningParams = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    restart = generalParams[blockNames.GeneralParams.restart]
    debug = generalParams[blockNames.GeneralParams.debug]
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructuresWildcard = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams[blockNames.GeneralParams.writeAllClustering]
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = str(clusteringBlock[blockNames.ClusteringTypes.params][blockNames.ClusteringTypes.ligandResname])

    print "================================"
    print "            PARAMS              "
    print "================================"
    print "Restarting simulations", restart
    print "Debug:", debug

    print "Iterations: %d, Mpi processors: %d, Pele steps: %d" % (simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps)

    print "SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

    print "SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type]
    if simulationRunner.hasExitCondition():
        print "Exit condition:", simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY[simulationRunner.parameters.exitCondition.type]
    print "Clustering method:", clusteringBlock[blockNames.ClusteringTypes.type]

    print "Output path: ", outputPath
    print "Initial Structures: ", initialStructuresWildcard
    print "================================\n\n"

    initialStructures = expandInitialStructuresWildcard(initialStructuresWildcard)
    checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    if not debug: atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

    utilities.makeFolder(outputPath)
    utilities.makeFolder(outputPathConstants.tmpFolder)
    saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)

    startFromScratch = False
    if restart:
        firstRun = findFirstRun(outputPath, outputPathConstants.clusteringOutputObject)
        if firstRun == 0:
            startFromScratch = True
        else:
            clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock, spawningParams, spawningCalculator, simulationRunner)

    if startFromScratch or not restart:
        firstRun = 0  # if restart false, but there were previous simulations
        clusteringMethod, initialStructuresAsString, initialClusters = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, outputPath, jsonParams, outputPathConstants, clusteringBlock, spawningParams, initialStructures)

    peleControlFileDictionary = {"COMPLEXES": initialStructuresAsString, "PELE_STEPS": simulationRunner.parameters.peleSteps}

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print "Iteration", i

        print "Preparing control file..."
        preparePeleControlFile(i, outputPathConstants, simulationRunner, peleControlFileDictionary)

        print "Production run..."
        if not debug:
            startTime = time.time()
            simulationRunner.runSimulation(outputPathConstants.tmpControlFilename % i)
            endTime = time.time()
            print "PELE %s sec" % (endTime - startTime)

        simulationRunner.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i)

        print "Clustering..."
        startTime = time.time()
        clusterEpochTrajs(clusteringMethod, i, outputPathConstants.epochOutputPathTempletized)
        endTime = time.time()
        print "Clustering ligand: %s sec" % (endTime - startTime)

        degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, i)
        spawningCalculator.log()
        print "Degeneracy", degeneracyOfRepresentatives


        clusteringMethod.writeOutput(outputPathConstants.clusteringOutputDir % i,
                                     degeneracyOfRepresentatives,
                                     outputPathConstants.clusteringOutputObject % i, writeAll)

        if i > 0:
            # Remove old clustering object, since we already have a newer one
            try:
                os.remove(outputPathConstants.clusteringOutputObject % (i-1))
            except OSError:
                # In case of restart
                pass

        # Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            numberOfSeedingPoints, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, i+1)
            simulationRunner.updateMappingProcessors(procMapping)
            initialStructuresAsString = createMultipleComplexesFilenames(numberOfSeedingPoints, outputPathConstants.tmpInitialStructuresTemplate, i+1)
            peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

        if clusteringMethod.symmetries and nativeStructure:
            fixReportsSymmetry(outputPathConstants.epochOutputPathTempletized % i, resname,
                               nativeStructure, clusteringMethod.symmetries)

        # check exit condition, if defined
        if simulationRunner.hasExitCondition() and simulationRunner.checkExitCondition(clusteringMethod):
            print "Simulation exit condition met at iteration %d" % i
            break

if __name__ == '__main__':
    args = parseArgs()
    main(args.controlFile)
