from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import shutil
import os
import json
import time
import glob
import argparse
import atexit
import ast
import numpy as np
from builtins import range
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
    arg = parser.parse_args()
    return arg


class EmptyInitialStructuresError(Exception):
    __module__ = Exception.__module__


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


def filterClustersAccordingToBox(simulationRunnerParams, clusteringObject):
    """
        Filter the clusters to select only the ones whose representative
        structures will fit into the selected box

        :param simulationRunnerParams: :py:class:`.SimulationParameters` Simulation parameters object
        :type simulationRunnerParams: :py:class:`.SimulationParameters`
        :param clusteringObject: Clustering object
        :type clusteringObject: :py:class:`.Clustering`

        :returns list, list: -- list of the filtered clusters, list of bools flagging wether the cluster is selected or not

    """
    box_center = ast.literal_eval(simulationRunnerParams.boxCenter)
    box_radius = simulationRunnerParams.boxRadius
    clustersFiltered = []
    clustersSelected = []
    for cluster in clusteringObject.clusters.clusters:
        if utilities.distanceCOM(box_center, cluster.pdb.getCOM()) < (box_radius-1):
            clustersFiltered.append(cluster)
            clustersSelected.append(True)
        else:
            clustersSelected.append(False)
    return clustersFiltered, clustersSelected


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
        PDB.initialise(structure, resname=resname)
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
    outputFilename = "fixedReport_%d"  # move to constants?
    trajName = "*traj*.pdb"  # move to constants?
    reportName = "*report_%d"  # move to constants?
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


def generateTrajectorySelectionString(epoch, epochOutputPathTempletized):
    """
        Generates the template for the name of the trajectories in a given epoch

        :param Poch: Epoch number
        :type epoch: int
        :param epochOutputPathTempletized: Templetized path where the trajectories of any epoch are stored
        :type epochOutputPathTempletized: str

        :returns: str -- Template for the name of the trajectories in a given
            epoch
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

        :return: int -- Current epoch
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

        :returns: bool -- True if the found clustering object is valid
    """
    try:
        utilities.readClusteringObject(objectPath)
        return True
    except EOFError:
        return False


def __unicodeToStr(data):
    # convert dict
    if isinstance(data, dict):
        return {__unicodeToStr(key): __unicodeToStr(value) for key, value in data.items()}
    # convert list
    if isinstance(data, list):
        return [__unicodeToStr(val) for val in data]
    # convert unicode to str
    if isinstance(data, unicode):
        return data.encode('utf-8')

    return data


def loadParams(jsonParams):
    """
        Read the control file in JSON format and extract the blocks of simulation,
        general parameters, clustering and spawning

        :param jsonParams: Control file in JSON format from where the parameters will be read
        :type jsonParams: json str
    """
    jsonFile = open(jsonParams, 'r').read()
    # parsedJSON = json.loads(jsonFile, object_hook=__unicodeToStr)
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
    paths = ast.literal_eval(snapshotsJSONSelectionString)
    if len(glob.glob(paths[-1])) == 0:
        sys.exit("No trajectories to cluster! Matching path:%s" % paths[-1])
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
        :type simulationRunner: :py:class:`.SimulationRunner`
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

        :returns: :py:class:`.Clustering` -- The clustering method to use in the
            adaptive sampling simulation
    """

    lastClusteringEpoch = firstRun - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % (lastClusteringEpoch)
    oldClusteringMethod = utilities.readClusteringObject(clusteringObjectPath)

    clusteringBuilder = clustering.ClusteringBuilder()
    clusteringMethod = clusteringBuilder.buildClustering(clusteringBlock,
                                                         spawningParams.reportFilename,
                                                         spawningParams.reportCol)

    if needToRecluster(oldClusteringMethod, clusteringMethod):
        print("Reclustering!")
        startTime = time.time()
        clusterPreviousEpochs(clusteringMethod, firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)
        endTime = time.time()
        print("Reclustering took %s sec" % (endTime - startTime))
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

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """

    clusteringMethod = getWorkingClusteringObjectAndReclusterIfNecessary(firstRun, outputPathConstants, clusteringBlock, spawningParams, simulationRunner)

    degeneracyOfRepresentatives = spawningCalculator.calculate(clusteringMethod.clusters.clusters, simulationRunner.parameters.processors-1, spawningParams, firstRun)
    spawningCalculator.log()
    print("Degeneracy", degeneracyOfRepresentatives)
    seedingPoints, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, firstRun)
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(seedingPoints, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    simulationRunner.updateMappingProcessors(procMapping)

    return clusteringMethod, initialStructuresAsString


def buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, controlFile, outputPathConstants, clusteringBlock, spawningParams, initialStructures, simulationRunner):
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

        :returns: :py:class:`.Clustering`, str -- The clustering method to use in the adaptive sampling simulation and the initial structures filenames
    """
    saveInitialControlFile(controlFile, outputPathConstants.originalControlFile)

    firstRun = 0
    copyInitialStructures(initialStructures, outputPathConstants.tmpInitialStructuresTemplate, firstRun)
    initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(len(initialStructures), outputPathConstants.tmpInitialStructuresTemplate, firstRun)

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
    if simulationRunner.parameters.boxCenter is not None:
        peleControlFileDictionary["BOX_RADIUS"] = simulationRunner.parameters.boxRadius
        peleControlFileDictionary["BOX_CENTER"] = simulationRunner.parameters.boxCenter
    simulationRunner.makeWorkingControlFile(outputPathConstants.tmpControlFilename % epoch, peleControlFileDictionary)


def main(jsonParams, clusteringHook=None):
    """
        Main body of the adaptive sampling program.

        :param jsonParams: A string with the name of the control file to use
        :type jsonParams: str
    """

    controlFileValidator.validate(jsonParams)
    generalParams, spawningBlock, simulationrunnerBlock, clusteringBlock = loadParams(jsonParams)

    spawningAlgorithmBuilder = spawning.SpawningAlgorithmBuilder()
    spawningCalculator, spawningParams = spawningAlgorithmBuilder.build(spawningBlock)

    runnerbuilder = simulationrunner.RunnerBuilder()
    simulationRunner = runnerbuilder.build(simulationrunnerBlock)

    restart = generalParams.get(blockNames.GeneralParams.restart, True)
    debug = generalParams.get(blockNames.GeneralParams.debug, False)
    outputPath = generalParams[blockNames.GeneralParams.outputPath]
    initialStructuresWildcard = generalParams[blockNames.GeneralParams.initialStructures]
    writeAll = generalParams.get(blockNames.GeneralParams.writeAllClustering, False)
    nativeStructure = generalParams.get(blockNames.GeneralParams.nativeStructure, '')
    resname = str(clusteringBlock[blockNames.ClusteringTypes.params][blockNames.ClusteringTypes.ligandResname])

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

    initialStructures = expandInitialStructuresWildcard(initialStructuresWildcard)
    if not initialStructures:
        raise EmptyInitialStructuresError("No initial structures found!!!")
    checkSymmetryDict(clusteringBlock, initialStructures, resname)

    outputPathConstants = constants.OutputPathConstants(outputPath)

    if not debug:
        atexit.register(utilities.cleanup, outputPathConstants.tmpFolder)

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
            checkMetricExitConditionMultipleTrajsinRestart(firstRun, outputPathConstants.epochOutputPathTempletized, simulationRunner)

    if startFromScratch or not restart:
        firstRun = 0  # if restart false, but there were previous simulations
        # These two lines were previously inside the
        # buildNewClusteringAndWriteInitialStructuresInNewSimulation function,
        # however, the equilibration procedure has to be called in between
        if not debug:
            shutil.rmtree(outputPath)
        utilities.makeFolder(outputPath)

        if simulationRunner.parameters.runEquilibration:
            initialStructures = simulationRunner.equilibrate(initialStructures, outputPathConstants, spawningParams.reportFilename, outputPath, resname)
        clusteringMethod, initialStructuresAsString, _ = buildNewClusteringAndWriteInitialStructuresInNewSimulation(debug, jsonParams, outputPathConstants, clusteringBlock, spawningParams, initialStructures, simulationRunner)
    peleControlFileDictionary = {"COMPLEXES": initialStructuresAsString, "PELE_STEPS": simulationRunner.parameters.peleSteps,
                                 "BOX_RADIUS": simulationRunner.parameters.boxRadius}
    if simulationRunner.parameters.modeMovingBox is not None and simulationRunner.parameters.boxCenter is None:
        simulationRunner.parameters.boxCenter = simulationRunner.selectInitialBoxCenter(initialStructuresAsString, resname)

    for i in range(firstRun, simulationRunner.parameters.iterations):
        print("Iteration", i)

        print("Preparing control file...")
        preparePeleControlFile(i, outputPathConstants, simulationRunner, peleControlFileDictionary)

        print("Production run...")
        if not debug:
            startTime = time.time()
            simulationRunner.runSimulation(outputPathConstants.tmpControlFilename % i)
            endTime = time.time()
            print("PELE %s sec" % (endTime - startTime))

        simulationRunner.writeMappingToDisk(outputPathConstants.epochOutputPathTempletized % i)

        print("Clustering...")
        startTime = time.time()
        clusterEpochTrajs(clusteringMethod, i, outputPathConstants.epochOutputPathTempletized)
        endTime = time.time()
        print("Clustering ligand: %s sec" % (endTime - startTime))

        if clusteringHook is not None:
            clusteringHook(clusteringMethod, outputPathConstants, simulationRunner, i+1)

        if simulationRunner.parameters.modeMovingBox is not None:
            simulationRunner.getNextIterationBox(clusteringMethod, outputPathConstants.epochOutputPathTempletized % i, resname)
            clustersList, clustersFiltered = filterClustersAccordingToBox(simulationRunner.parameters, clusteringMethod)
        else:
            clustersList = clusteringMethod.clusters.clusters

        degeneracyOfRepresentatives = spawningCalculator.calculate(clustersList, simulationRunner.parameters.processors-1, spawningParams, i)
        spawningCalculator.log()

        if degeneracyOfRepresentatives is not None:
            if simulationRunner.parameters.modeMovingBox is not None:
                degeneracyOfRepresentatives = mergeFilteredClustersAccordingToBox(degeneracyOfRepresentatives, clustersFiltered)
            print("Degeneracy", degeneracyOfRepresentatives)
            assert len(degeneracyOfRepresentatives) == len(clusteringMethod.clusters.clusters)
        else:
            # When using null spawning the calculate method returns None
            assert spawningCalculator.type == spawningTypes.SPAWNING_TYPES.null, "calculate returned None with spawning type %s" % spawningTypes.SPAWNING_TYPE_TO_STRING_DICTIONARY[spawningCalculator.type]

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
            if degeneracyOfRepresentatives is not None:
                numberOfSeedingPoints, procMapping = spawningCalculator.writeSpawningInitialStructures(outputPathConstants, degeneracyOfRepresentatives, clusteringMethod, i+1)
                simulationRunner.updateMappingProcessors(procMapping)
                initialStructuresAsString = simulationRunner.createMultipleComplexesFilenames(numberOfSeedingPoints, outputPathConstants.tmpInitialStructuresTemplate, i+1)
                peleControlFileDictionary["COMPLEXES"] = initialStructuresAsString

        if clusteringMethod.symmetries and nativeStructure:
            fixReportsSymmetry(outputPathConstants.epochOutputPathTempletized % i, resname,
                               nativeStructure, clusteringMethod.symmetries)

        # check exit condition, if defined
        if simulationRunner.hasExitCondition():
            if simulationRunner.checkExitCondition(clusteringMethod, outputPathConstants.epochOutputPathTempletized % i):
                print("Simulation exit condition met at iteration %d, stopping" % i)
                break
            else:
                print("Simulation exit condition not met at iteration %d, continuing..." % i)

if __name__ == '__main__':
    args = parseArgs()
    main(args.controlFile)
