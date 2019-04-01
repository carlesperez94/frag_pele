from __future__ import absolute_import, division, print_function, unicode_literals
import os
import sys
import ast
import json
import time
import glob
import atexit
import shutil
import signal
import errno
import argparse
import numpy as np
from builtins import range
from AdaptivePELE.constants import blockNames, constants
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities
from AdaptivePELE.utilities.synchronization import ProcessesManager
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
    arg = parser.parse_args()
    return arg


class InitialStructuresError(Exception):
    __module__ = Exception.__module__


def cleanProcessesFiles(folder):
    """
        Clean the processes files from a previous simulation
        :param folder: Folder where the files are stored
        :type folder: str
    """
    processes = glob.glob(os.path.join(folder, "*.proc"))
    for process_file in processes:
        try:
            os.remove(process_file)
        except OSError:
            pass


def printRunInfo(restart, debug, simulationRunner, spawningCalculator, clusteringBlock, outputPath, initialStructuresWildcard):
    """
        Print a summary of the run paramaters

        :param restart: Flag on whether to continue a previous simulation
        :type restart: bool
        :param debug: Flag to mark whether simulation is run in debug mode
        :type debug: bool
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param spawningCalculator: Spawning calculator object
        :type spawningCalculator: :py:class:`.SpawningCalculator`
        :param clusteringBlock: Block of the control file with clustering options
        :type clusteringBlock: dict
        :param outputPath: Path where to write the simulation output
        :type outputPath: str
        :param initialStructuresWildcard: Wildcard expression to find the initial structures
        :type initialStructuresWildcard: str

    """
    print("================================")
    print("            PARAMS              ")
    print("================================")
    print("Restarting simulations", restart)
    print("Debug:", debug)
    print("Iterations: %d, Mpi processors: %d, Pele steps: %d" % (simulationRunner.parameters.iterations, simulationRunner.parameters.processors, simulationRunner.parameters.peleSteps))
    print("SpawningType:", spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type])
    print("SimulationType:", simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY[simulationRunner.type])
    if simulationRunner.hasExitCondition():
        print("Exit condition:", simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY[simulationRunner.parameters.exitCondition.type])
    print("Clustering method:", clusteringBlock[blockNames.ClusteringTypes.type])
    print("Output path: ", outputPath)
    print("Initial Structures: ", initialStructuresWildcard)
    print("================================\n\n")
    sys.stdout.flush()


def cleanPreviousSimulation(output_path, allTrajs):
    """
        Clean the uneeded data from a previous simulation

        :param output_path: Path where the data is stored
        :type output_path: str
        :param allTrajs: Path where the discretized trajectories for MSM are stored
        :type allTrajs: str
    """
    equilibration_folders = glob.glob(os.path.join(output_path, "equilibration*"))
    for folder in equilibration_folders:
        try:
            shutil.rmtree(folder)
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                raise
            # If another process deleted the folder between the glob and the
            # actual removing an OSError is raised
    epochs = utilities.get_epoch_folders(output_path)
    for epoch in epochs:
        try:
            shutil.rmtree(os.path.join(output_path, epoch))
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                raise
    try:
        shutil.rmtree(allTrajs)
    except OSError as exc:
        # this folder may not exist, in which case we just carry on
        pass


def createMappingForFirstEpoch(initialStructures, topologies, processors):
    """
        Create the topology mapping for the first iteration

        :param initialStructures: List of the initial structures for the first iteration
        :type initialStructures: list
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param processors: Number of trajectories
        :type processors: int

    """
    topologyMapping = list(range(1, len(initialStructures)))+[0]
    topologyMapping = topologyMapping*int(np.ceil(processors/len(initialStructures)))
    topologies.topologyMap[0] = topologyMapping[:processors]


def writeTopologyFiles(topologies, destination):
    """
        Write the topology files to the desired destination

        :param topologies: List of topology files
        :type topologies: list
        :param destination: Path where to copy the toplogy files
        :type destination: str
    """
    destination = os.path.join(destination, "topology_%d.pdb")
    for i, topology in enumerate(topologies):
        shutil.copy(topology, destination % i)


def checkMetricExitConditionMultipleTrajsinRestart(firstRun, outputFolder, simulationRunner):
    """
        Check the previous simulation data when restarting a simulation with a
        multiple trajectory metric exit condition

        :param firstRun: First epoch to be run in restart simulation
        :type firstRun: int
        :param outputFolder: Folder of the previous simulation data
        :type outputFolder: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
    """
    if simulationRunner.hasExitCondition():
        if simulationRunner.parameters.exitCondition.type != blockNames.ExitConditionType.metricMultipleTrajs:
            return
        for i in range(firstRun):
            simulationRunner.parameters.exitCondition.checkExitCondition(outputFolder % i)


def mergeFilteredClustersAccordingToBox(degeneracy, clustersFiltering):
    """
        Merge the (possibly) partial degeneracy to obtain a complete list,
        the partial list comes from the clusters excluded by the moving box
        :param degeneracy: Degeneracy of the clusters
        :param degeneracy: list
        :param clustersFiltering: List of bools indicating whether a cluster was filtered or not
        :param clustersFiltering: list

        :returns list: -- complete list of cluster degeneracies
    """
    assert len(degeneracy) == sum(clustersFiltering)
    newDegeneracy = []
    for filtered in clustersFiltering:
        if filtered:
            newDegeneracy.append(degeneracy.pop(0))
        else:
            newDegeneracy.append(0)
    return np.array(newDegeneracy)


def expandInitialStructuresWildcard(initialStructuresWildcard):
    """
        Returns the initial structures after expanding the initial structures wildcard

        :param initialStructureWildcard: Wildcard that matches the initial structures
        :type initialStructureWildcard: str

        :return: list of str -- The expanded initial structures
    """
    totalInitialStructures = []
    for initialStructureWildcard in initialStructuresWildcard:
        expandedStructures = glob.glob(initialStructureWildcard)
        totalInitialStructures.extend(map(os.path.abspath, expandedStructures))
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
        PDB.initialise(structure, resname=resname)
        utilities.assertSymmetriesDict(symmetries, PDB)


def fixReportsSymmetry(outputPath, resname, nativeStructure, symmetries, topologies):
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
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`

        :raise IndexError: If original report file not found in output folder
    """
    outputFilename = "fixedReport_%d"  # move to constants?
    trajName = "*traj*"  # move to constants?
    reportName = "*report_%d"  # move to constants?
    epoch = int(os.path.basename(outputPath))
    trajs = glob.glob(os.path.join(outputPath, trajName))
    nativePDB = atomset.PDB()
    nativePDB.initialise(str(nativeStructure), resname=resname)
    for traj in trajs:
        trajNum = utilities.getTrajNum(traj)
        rmsd = list(utilities.getRMSD(traj, nativePDB, resname, symmetries, topology=topologies.getTopology(epoch, trajNum)))
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


def generateTrajectorySelectionString(epoch, epochOutputPathTempletized):
    """
        Generates the template for the name of the trajectories in a given epoch

        :param epoch: Epoch number
        :type epoch: int
        :param epochOutputPathTempletized: Templetized path where the trajectories of any epoch are stored
        :type epochOutputPathTempletized: str

        :returns: str -- Template for the name of the trajectories in a given
            epoch
    """
    return "[\"" + os.path.join(epochOutputPathTempletized % epoch, constants.trajectoryBasename) + "\"]"


def findFirstRun(outputPath, clusteringOutputObject, simulationRunner, restart):
    """
        Find the last epoch that was properly simulated and clusterized and
        and return the first epoch to run in case of restart

        :param outputPath: Simulation output path
        :type outputPath: str
        :param clusteringOutputObject: Templetized name of the clustering object
        :type clusteringOutputObject: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param restart: Whether to restart a previous simulation
        :type restart: bool

        :return: int -- Current epoch
    """

    folderWithSimulationData = outputPath
    allFolders = os.listdir(folderWithSimulationData)
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort(reverse=True)

    objectsFound = []
    for epoch in epochFolders:
        if simulationRunner.checkSimulationInterrupted(epoch, outputPath, restart):
            # this should only happen in MD simulations, where checkpoints are
            # periodically written in case the adaptive run dies at
            # mid-simulation, be able to use the already computed trajectories
            return epoch
        if os.path.exists(clusteringOutputObject % epoch):
            objectsFound.append(epoch)
        if objectsFound and epoch < (objectsFound[0]-5):
            break
    while objectsFound:
        epoch = objectsFound.pop(0)
        if checkIntegrityClusteringObject(clusteringOutputObject % epoch):
            return epoch + 1
    return None


def checkIntegrityClusteringObject(objectPath):
    """
        Test whether the found clustering object to reload is a valid object

        :param objectPath: Clustering object path
        :type objectPath: str

        :returns: bool -- True if the found clustering object is valid
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
    with open(jsonParams, 'r') as f:
        jsonFile = f.read()
    parsedJSON = json.loads(jsonFile)

    return parsedJSON[blockNames.ControlFileParams.generalParams], parsedJSON[blockNames.ControlFileParams.spawningBlockname],\
        parsedJSON[blockNames.ControlFileParams.simulationBlockname], parsedJSON[blockNames.ControlFileParams.clusteringBlockname]


def saveInitialControlFile(jsonParams, originalControlFile):
    """
        Save the adaptive control file jsonParams in originalControlFile

        :param jsonParams: Input control file in JSON format
        :type jsonParams: str
        :param originalControlFile: Path where to save the control file
        :type originalControlFile: str
    """
    with open(originalControlFile, 'w') as f:
        with open(jsonParams, 'r') as fr:
            jsonFile = fr.read()
        f.write(jsonFile)


def needToRecluster(oldClusteringMethod, newClusteringMethod):
    """
        Check if the parameters have changed in a restart and we need to redo
        the clustering. In particular: type of clustering, theshold calculator
        or distance

        :param oldClusteringMethod: Clustering in a previous simulation before the restart
        :type oldClusteringMethod: :py:class:`.Clustering`
        :param newClusteringMethod: Clustering in the restarted simulation
        :type newClusteringMethod: :py:class:`.Clustering`

        :returns: bool -- If clustering needs to be redone
    """

    # Check 1: change of type
    if oldClusteringMethod.type != newClusteringMethod.type:
        return True

    # Check 2: Change of thresholdCalculator and thresholdDistance
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.rmsd or\
       oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMap:
        return oldClusteringMethod.thresholdCalculator != newClusteringMethod.thresholdCalculator or\
                abs(oldClusteringMethod.contactThresholdDistance - newClusteringMethod.contactThresholdDistance) > 1e-7

    # Check 3: Change of similarity Evaluator in contactMap clustering
    if oldClusteringMethod.type == clusteringTypes.CLUSTERING_TYPES.contactMap:
        return oldClusteringMethod.similarityEvaluator.typeEvaluator != newClusteringMethod.similarityEvaluator.typeEvaluator


def clusterEpochTrajs(clusteringMethod, epoch, epochOutputPathTempletized, topologies, outputPathConstants=None):
    """
        Cluster the trajecotories of a given epoch

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering`
        :param epoch: Number of the epoch to cluster
        :type epoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
"""

    snapshotsJSONSelectionString = generateTrajectorySelectionString(epoch, epochOutputPathTempletized)
    paths = ast.literal_eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No trajectories to cluster! Matching path:%s" % paths[-1])
    clusteringMethod.cluster(paths, topology=topologies, epoch=epoch, outputPathConstants=outputPathConstants)


def clusterPreviousEpochs(clusteringMethod, finalEpoch, epochOutputPathTempletized, simulationRunner, topologies, outputPathConstants=None):
    """
        Cluster all previous epochs using the clusteringMethod object

        :param clusteringMethod: Clustering object
        :type clusteringMethod: :py:class:`.Clustering`
        :param finalEpoch: Last epoch to cluster (not included)
        :type finalEpoch: int
        :param epochOutputPathTempletized: Path where to find the trajectories
        :type epochOutputPathTempletized: str
        :param simulationRunner: Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param topologies: Topology object containing the set of topologies needed for the simulatioies
        :type topologies: :py:class:`.Topology`
        :param outputPathConstants: Contains outputPath-related constants
        :type outputPathConstants: :py:class:`.OutputPathConstants`
"""
    for i in range(finalEpoch):
        simulationRunner.readMappingFromDisk(epochOutputPathTempletized % i)
        topologies.readMappingFromDisk(epochOutputPathTempletized % i, i)
        clusterEpochTrajs(clusteringMethod, i, epochOutputPathTempletized, topologies, outputPathConstants)


def getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner, topologies, processManager):
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
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param processManager: Object to synchronize the possibly multiple processes
        :type processManager: :py:class:`.ProcessesManager`

        :returns: :py:class:`.Clustering` -- The clustering method to use in the
            adaptive sampling simulation
    """
    if not processManager.isMaster():
        for ij in range(firstRun):
            topologies.readMappingFromDisk(outputPathConstants.epochOutputPathTempletized % ij, ij)
        return
    lastClusteringEpoch = firstRun - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % (lastClusteringEpoch)
    oldClusteringMethod = utilities.readClusteringObject(clusteringObjectPath)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)

    clusteringMethod.setProcessors(simulationRunner.getWorkingProcessors())
    if needToRecluster(oldClusteringMethod, clusteringMethod):
        utilities.print_unbuffered("Reclustering!")
        startTime = time.time()
        clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner, topologies, outputPathConstants.allTrajsPath)
        endTime = time.time()
        utilities.print_unbuffered("Reclustering took %s sec" % (endTime - startTime))
    else:
        clusteringMethod = oldClusteringMethod
        clusteringMethod.setCol(spawningParams.reportCol)
        for ij in range(firstRun):
            topologies.readMappingFromDisk(outputPathConstants.epochOutputPathTempletized % ij, ij)

    return clusteringMethod


def buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock,
                                                         spawningParams, spawningCalculator, simulationRunner, topologies, processManager):
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
        :param topologies: Topology object containing the set of topologies needed for the simulation
        :type topologies: :py:class:`.Topology`
        :param processManager: Object to synchronize the possibly multiple processes
        :type processManager: :py:class:`.ProcessesManager`

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """
    processorManagerFilename = "procMapping.txt"
    clusteringMethod = getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner, topologies, processManager)
    if processManager.isMaster():
        degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.getClusterListForSpawning(), simulationRunner.getWorkingProcessors(), firstRun)
        spawningCalculator.log()
        _, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, firstRun, topologies=topologies)
        utilities.writeProcessorMappingToDisk(outputPathConstants.tmpFolder, processorManagerFilename, procMapping)
    else:
        clusteringMethod = None
    processManager.barrier()
    if not processManager.isMaster():
        procMapping = utilities.readProcessorMappingFromDisk(outputPathConstants.tmpFolder, processorManagerFilename)
    # for compatibility with old data
    procMapping = [element if element is not None else (0, 0, 0) for element in procMapping]
    topologies.mapEpochTopologies(firstRun, procMapping)
    simulationRunner.updateMappingProcessors(procMapping)
    processManager.barrier()
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(simulationRunner.getWorkingProcessors(), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    return clusteringMethod, initialStructuresAsString


def buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, controlFile, outputPathConstants, clusteringBlock, spawningParams, initialStructures, simulationRunner, processManager):
    """
        Build the clustering object and copies initial structures from control file.
        Returns the clustering object to use and the initial structures filenames as string

        :param debug: In debug, it will not remove the simulations
        :type debug: bool
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
        :param simulationRunner: :py:class:`.SimulationRunner` Simulation runner object
        :type simulationRunner: :py:class:`.SimulationRunner`
        :param processManager: Object to synchronize the possibly multiple processes
        :type processManager: :py:class:`.ProcessesManager`

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """
    firstRun = 0
    if processManager.isMaster():
        saveInitialControlFile(controlFile, outputPathConstants.originalControlFile)

        copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    processManager.barrier()
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

    if processManager.isMaster():
        clusteringBuilder = clustering.ClusteringBuilder()
        clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                             spawningParams.reportFilename,
                                                             spawningParams.reportCol)
    else:
        clusteringMethod = None
    return clusteringMethod, initialStructuresAsString


def main(jsonParams, clusteringHook=None):
    """
        Main body of the adaptive sampling program.

        :param jsonParams: A string with the name of the control file to use
        :type jsonParams: str
    """

    controlFileValidator.validate(jsonParams)
    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    spawningCalculator = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    restart = generalParams.get(blockNames.GeneralParams.restart, True)
    debug = generalParams.get(blockNames.GeneralParams.debug, False)
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructuresWildcard = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams.get(blockNames.GeneralParams.writeAllClustering, False)
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = clusteringBlock[blockNames.ClusteringTypes.params].get(blockNames.ClusteringTypes.ligandResname)

    initialStructures = expandInitialStructuresWildcard(initialStructuresWildcard)
    if not initialStructures:
        raise InitialStructuresError("No initial structures found!!!")

    if len(initialStructures) > simulationRunner.getWorkingProcessors():
        raise InitialStructuresError("Error: More initial structures than Working Processors found!!!")

    if resname is not None:
        checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    if not debug:
        atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

    simulationRunner.unifyReportNames(spawningCalculator.parameters.reportFilename)
    utilities.makeFolder(outputPath)
    utilities.makeFolder(outputPathConstants.tmpFolder)
    utilities.makeFolder(outputPathConstants.topologies)
    processManager = ProcessesManager(outputPath, simulationRunner.getNumReplicas())
    firstRun = findFirstRun(outputPath, outputPathConstants.clusteringOutputObject, simulationRunner, restart)
    if processManager.isMaster():
        printRunInfo(restart, debug, simulationRunner, spawningCalculator, clusteringBlock, outputPath, initialStructuresWildcard)
        saveInitialControlFile(jsonParams, outputPathConstants.originalControlFile)
    processManager.barrier()
    # once the replicas are properly syncronized there is no need for the
    # process files, and erasing them allows us to restart simulations
    cleanProcessesFiles(processManager.syncFolder)

    topologies = utilities.Topology(outputPathConstants.topologies)
    if restart and firstRun is not None:
        topology_files = glob.glob(os.path.join(outputPathConstants.topologies, "topology*.pdb"))
        topologies.setTopologies(topology_files)
        if firstRun == 0:
            createMappingForFirstEpoch(initialStructures, topologies, simulationRunner.getWorkingProcessors())
            clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, jsonParams, outputPathConstants, clusteringBlock, spawningCalculator.parameters, initialStructures, simulationRunner, processManager)
        else:
            clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInRestart(firstRun, outputPathConstants, clusteringBlock, spawningCalculator.parameters, spawningCalculator, simulationRunner, topologies, processManager)
        if processManager.isMaster():
            checkMetricExitConditionMultipleTrajsinRestart(firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)
        processManager.barrier()

    if firstRun is None or not restart:
        topologies.setTopologies(initialStructures)
        if processManager.isMaster():
            if not debug:
                cleanPreviousSimulation(outputPath, outputPathConstants.allTrajsPath)
            writeTopologyFiles(initialStructures, outputPathConstants.topologies)
        processManager.barrier()
        firstRun = 0  # if restart false, but there were previous simulations

        if simulationRunner.parameters.runEquilibration:
            initialStructures = simulationRunner.equilibrate(initialStructures, outputPathConstants, spawningCalculator.parameters.reportFilename, outputPath, resname, processManager, topologies)
            # write the equilibration structures for each replica
            processManager.writeEquilibrationStructures(outputPathConstants.tmpFolder, initialStructures)
            if processManager.isMaster() and simulationRunner.parameters.constraints:
                # write the new constraints for synchronization
                utilities.writeNewConstraints(outputPathConstants.tmpFolder, "new_constraints.txt", simulationRunner.parameters.constraints)
            processManager.barrier()

            if not processManager.isMaster() and simulationRunner.parameters.constraints:
                simulationRunner.parameters.constraints = utilities.readConstraints(outputPathConstants.tmpFolder, "new_constraints.txt")
            # read all the equilibration structures
            initialStructures = processManager.readEquilibrationStructures(outputPathConstants.tmpFolder)
            topologies.setTopologies(initialStructures, cleanFiles=processManager.isMaster())
            writeTopologyFiles(initialStructures, outputPathConstants.topologies)
        createMappingForFirstEpoch(initialStructures, topologies, simulationRunner.getWorkingProcessors())

        clusteringMethod, initialStructuresAsString = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, jsonParams, outputPathConstants, clusteringBlock, spawningCalculator.parameters, initialStructures, simulationRunner, processManager)

    if processManager.isMaster():
        clusteringMethod.setProcessors(simulationRunner.getWorkingProcessors())
    if simulationRunner.parameters.modeMovingBox is not None and simulationRunner.parameters.boxCenter is None:
        simulationRunner.parameters.boxCenter = simulationRunner.selectInitialBoxCenter(initialStructuresAsString, resname)
    for i in range(firstRun, simulationRunner.parameters.iterations):
        if processManager.isMaster():
            utilities.print_unbuffered("Iteration", i)
            outputDir = outputPathConstants.epochOutputPathTempletized % i
            utilities.makeFolder(outputDir)

            simulationRunner.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i)
            topologies.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i, i)

        processManager.barrier()
        if processManager.isMaster():
            utilities.print_unbuffered("Production run...")
        if not debug:
            simulationRunner.runSimulation(i, outputPathConstants, initialStructuresAsString, topologies, spawningCalculator.parameters.reportFilename, processManager)
        processManager.barrier()

        if processManager.isMaster():
            utilities.print_unbuffered("Clustering...")
            startTime = time.time()
            clusterEpochTrajs(clusteringMethod, i, outputPathConstants.epochOutputPathTempletized, topologies, outputPathConstants)
            endTime = time.time()
            utilities.print_unbuffered("Clustering ligand: %s sec" % (endTime - startTime))

            if clusteringHook is not None:
                clusteringHook(clusteringMethod, outputPathConstants, simulationRunner, i+1)
            clustersList = clusteringMethod.getClusterListForSpawning()
            clustersFiltered = [True for _ in clusteringMethod]

        if simulationRunner.parameters.modeMovingBox is not None:
            simulationRunner.getNextIterationBox(outputPathConstants.epochOutputPathTempletized % i, resname, topologies, i)
            if processManager.isMaster():
                clustersList, clustersFiltered = clusteringMethod.filterClustersAccordingToBox(simulationRunner.parameters)

        if processManager.isMaster():
            if spawningCalculator.parameters.filterByMetric:
                clustersList, clustersFiltered = clusteringMethod.filterClustersAccordingToMetric(clustersFiltered, spawningCalculator.parameters.filter_value, spawningCalculator.parameters.condition, spawningCalculator.parameters.filter_col)

            degeneracyOfRepresentatives = spawningCalculator.calculate(clustersList, simulationRunner.getWorkingProcessors(), i)
            spawningCalculator.log()
            # this method only does works with MSM-based spwaning methods,
            # creating a plot of the stationary distribution and the PMF, for
            # the rest of methods it does nothing
            spawningCalculator.createPlots(outputPathConstants, i, clusteringMethod)

            if degeneracyOfRepresentatives is not None:
                if simulationRunner.parameters.modeMovingBox is not None or spawningCalculator.parameters.filterByMetric:
                    degeneracyOfRepresentatives = mergeFilteredClustersAccordingToBox(degeneracyOfRepresentatives, clustersFiltered)
                utilities.print_unbuffered("Degeneracy", degeneracyOfRepresentatives)
                assert len(degeneracyOfRepresentatives) == len(clusteringMethod)
            else:
                # When using null or independent spawning the calculate method returns None
                assert spawningCalculator.type in spawningTypes.SPAWNING_NO_DEGENERACY_TYPES, "calculate returned None with spawning type %s" % spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

            clusteringMethod.writeOutput(outputPathConstants.clusteringOutputDir % i,
                                         degeneracyOfRepresentatives,
                                         outputPathConstants.clusteringOutputObject % i, writeAll)
            simulationRunner.cleanCheckpointFiles(outputPathConstants.epochOutputPathTempletized % i)

            if i > 0:
                # Remove old clustering object, since we already have a newer one
                try:
                    os.remove(outputPathConstants.clusteringOutputObject % (i-1))
                except OSError:
                    # In case of restart
                    pass

        # Prepare for next pele iteration
        if i != simulationRunner.parameters.iterations-1:
            # Differentiate between null spawning and the rest of spawning
            # methods
            if spawningCalculator.shouldWriteStructures():
                if processManager.isMaster():
                    _, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants,
                                                                                       degeneracyOfRepresentatives,
                                                                                       clusteringMethod,
                                                                                       i+1,
                                                                                       topologies=topologies)
                    utilities.writeProcessorMappingToDisk(outputPathConstants.tmpFolder, "processMapping.txt", procMapping)
                processManager.barrier()
                if not processManager.isMaster():
                    procMapping = utilities.readProcessorMappingFromDisk(outputPathConstants.tmpFolder, "processMapping.txt")
                simulationRunner.updateMappingProcessors(procMapping)
                topologies.mapEpochTopologies(i+1, procMapping)
                initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(simulationRunner.getWorkingProcessors(),
                                                                                              outputPathConstants.tmpInitialStructuresTemplate,
                                                                                              i+1)

        if processManager.isMaster():
            topologies.writeTopologyObject()
            if clusteringMethod.symmetries and nativeStructure:
                fixReportsSymmetry(outputPathConstants.epochOutputPathTempletized % i, resname,
                                   nativeStructure, clusteringMethod.symmetries, topologies)

            # check exit condition, if defined
            if simulationRunner.hasExitCondition():
                if simulationRunner.checkExitCondition(clusteringMethod, outputPathConstants.epochOutputPathTempletized % i):
                    utilities.print_unbuffered("Simulation exit condition met at iteration %d, stopping" % i)
                    # send a signal to all possible adaptivePELE copies to stop
                    for pid in processManager.lockInfo:
                        if pid != processManager.pid:
                            os.kill(pid, signal.SIGTERM)
                    break
                else:
                    utilities.print_unbuffered("Simulation exit condition not met at iteration %d, continuing..." % i)
        processManager.barrier()

if __name__ == '__main__':
    args = parseArgs()
    main(args.controlFile)
