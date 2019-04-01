from __future__ import absolute_import, division, print_function, unicode_literals
import time
import os
import json
import shutil
import string
import sys
import numpy as np
import glob
import itertools
import multiprocessing as mp
from builtins import range
from AdaptivePELE.constants import constants, blockNames
from AdaptivePELE.simulation import simulationTypes
from AdaptivePELE.atomset import atomset, RMSDCalculator
from AdaptivePELE.utilities import utilities, PDBLoader
SKLEARN = True
OPENMM = True
try:
    from sklearn.cluster import KMeans
except ImportError:
    SKLEARN = False

try:
    import subprocess32 as subprocess
except ImportError:
    import subprocess

try:
    basestring
except NameError:
    basestring = str
try:
    from AdaptivePELE.simulation import openmm_simulations as sim
except ImportError:
    OPENMM = False
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


class SimulationParameters:
    def __init__(self):
        self.processors = 0
        self.executable = ""
        self.templetizedControlFile = ""
        self.dataFolder = ""
        self.documentsFolder = ""
        self.iterations = 0
        self.peleSteps = 0
        self.seed = 0
        self.exitCondition = None
        self.time = None
        self.boxCenter = None
        self.boxRadius = 20
        self.modeMovingBox = None
        self.runEquilibration = False
        self.equilibrationLength = None
        self.destination = None
        self.origin = None
        self.equilibrationMode = ""
        self.SASAforBox = None
        self.columnSASA = None
        self.reportName = None
        self.trajectoryName = None
        self.srun = False
        self.srunParameters = None
        self.mpiParameters = None
        self.numberEquilibrationStructures = 10
        self.reportName = ""
        # parameters needed for MD simulations and their defaults
        self.timeStep = 2
        self.ligandCharge = 0
        self.nonBondedCutoff = 8
        self.Temperature = 300
        self.runningPlatform = "CPU"
        self.minimizationIterations = 2000
        self.reporterFreq = None
        self.energyReport = True
        self.productionLength = 0
        self.ligandName = None
        self.waterBoxSize = 8
        self.trajsPerReplica = None
        self.numReplicas = 1
        self.equilibrationLengthNVT = 200000
        self.equilibrationLengthNPT = 500000
        self.devicesPerTrajectory = 1
        self.constraintsMin = 5
        self.constraintsNVT = 5
        self.constraintsNPT = 0.5
        self.maxDevicesPerReplica = None
        self.forcefield = "ff99SB"
        self.customparamspath = None
        self.format = None
        self.constraints = None
        self.boxType = None
        self.cylinderBases = None


class SimulationRunner:
    def __init__(self, parameters):
        self.parameters = parameters
        self.processorsToClusterMapping = []

    def runSimulation(self, epoch, outputPathConstants, initialStructuresAsString, topologies, reportFileName, processManager):
        pass

    def getWorkingProcessors(self):
        """
            Return the number of working processors, i.e. number of trajectories
        """
        return self.parameters.processors

    def getNumReplicas(self):
        """
            Return the number of replicas, only useful for MD simulations
        """
        return self.parameters.numReplicas

    def hasExitCondition(self):
        """
            Check if an exit condition has been set

            :returns: bool -- True if an exit condition is set
        """
        return self.parameters.exitCondition is not None

    def checkSimulationInterrupted(self, epoch, outputpath, restart):
        """
            Check wether the simulation was interrupted before finishing

            :param epoch: Epoch number
            :type epoch: int
            :param outputpath: Simulation output path
            :type outputpath: str
            :param restart: Whether to restart a previous simulation
            :type restart: bool

            :returns: bool -- True if the simulations where interrupted
        """
        # for Pele and Test simulation there is no proper way to check so return
        # False and rely on the clustering for such check
        return False

    def cleanCheckpointFiles(self, epoch):
        """
            Clean the restart files generated if the simulation was interrupted
            before finishing

            :param epoch: Epoch number
            :type epoch: int
        """
        # for Pele and Test simulation there is no proper way to restart, so
        # just pass
        pass

    def checkExitCondition(self, clustering, outputFolder):
        """
            Check if the exit condition has been met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- True if the exit condition is met
        """
        if self.parameters.exitCondition.type == simulationTypes.EXITCONDITION_TYPE.METRICMULTIPLETRAJS:
            return self.parameters.exitCondition.checkExitCondition(outputFolder)
        else:
            return self.parameters.exitCondition.checkExitCondition(clustering)

    def makeWorkingControlFile(self, workingControlFilename, dictionary, inputTemplate=None):
        """
            Substitute the values in the templetized control file

            :param workingControlFilename: Name of the template control file
            :type workingControlFilename: str
            :param dictionary: Dictonary containing the parameters to substitute
                in the control file
            :param inputFileTemplate: Template control file
            :type inputFileTemplate: str
            :type dictionary: dict
        """
        if inputTemplate is None:
            with open(self.parameters.templetizedControlFile, "r") as inputFile:
                inputFileContent = inputFile.read()
        else:
            inputFileContent = inputTemplate

        inputFileTemplate = string.Template(inputFileContent)
        outputFileContent = inputFileTemplate.substitute(dictionary)
        with open(workingControlFilename, "w") as outputFile:
            outputFileContent = outputFileContent.replace("'", '"')
            outputFile.write(outputFileContent)

    def updateMappingProcessors(self, mapping):
        """
            Update the value of the processorsToClusterMapping, a list with the
            snapshot from which the trajectories will start in the next iteration

            :param mapping: List with the snapshot from which the trajectories
                will start in the next iteration
            :type mapping: list
        """
        self.processorsToClusterMapping = mapping[1:]+[mapping[0]]

    def writeMappingToDisk(self, epochDir):
        """
            Write the processorsToClusterMapping to disk

            :param epochDir: Name of the folder where to write the
                processorsToClusterMapping
            :type epochDir: str
        """
        if len(self.processorsToClusterMapping) == 0:
            return
        utilities.writeProcessorMappingToDisk(epochDir, "processorMapping.txt", self.processorsToClusterMapping)

    def readMappingFromDisk(self, epochDir):
        """
            Read the processorsToClusterMapping from disk

            :param epochDir: Name of the folder where to write the
                processorsToClusterMapping
            :type epochDir: str
        """
        self.processorsToClusterMapping = utilities.readProcessorMappingFromDisk(epochDir, "processorMapping.txt")

    def setZeroMapping(self):
        """
            Set the processorsToClusterMapping to zero
        """
        self.processorsToClusterMapping = [0 for _ in range(1, self.parameters.processors)]

    def createMultipleComplexesFilenames(self, numberOfSnapshots, tmpInitialStructuresTemplate, iteration, equilibration=False):
        """
            Creates the string to substitute the complexes in the PELE control file

            :param numberOfSnapshots: Number of complexes to write
            :type numberOfSnapshots: int
            :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
            :type tmpInitialStructuresTemplate: str
            :param iteration: Epoch number
            :type iteration: int
            :param equilibration: Flag to mark wether the complexes are part of an
                equilibration run
            :type equilibration: bool

            :returns: str -- jsonString to be substituted in PELE control file
        """
        pass

    def prepareControlFile(self, epoch, outputPathConstants, peleControlFileDictionary):
        """
            Substitute the parameters in the PELE control file specified with the
            provided in the control file

            :param epoch: Epoch number
            :type epoch: int
            :param outputPathConstants: Object that has as attributes constant related to the outputPath that will be used to create the working control file
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param peleControlFileDictionary: Dictonary containing the values of the parameters to substitute in the control file
            :type peleControlFileDictionary: dict
        """
        outputDir = outputPathConstants.epochOutputPathTempletized % epoch
        peleControlFileDictionary["OUTPUT_PATH"] = outputDir
        peleControlFileDictionary["SEED"] = self.parameters.seed + epoch * self.parameters.processors
        if self.parameters.boxCenter is not None:
            peleControlFileDictionary["BOX_RADIUS"] = self.parameters.boxRadius
            peleControlFileDictionary["BOX_CENTER"] = self.parameters.boxCenter
        self.makeWorkingControlFile(outputPathConstants.tmpControlFilename % epoch, peleControlFileDictionary)

    def unifyReportNames(self, spawningReportName):
        """
            Ensure that the reportName in the simulation parameters is the same
            as the one provided in the spawning parameters

            :param spawningReportName: Name of the report file provided in the spawning parameters
            :type spawningReportName: str
        """
        if spawningReportName is None:
            return
        baseReportName = self.parameters.reportName.split("_%d")
        if baseReportName[0] != spawningReportName:
            baseReportName[0] = spawningReportName
            self.parameters.reportName = "_%d".join(baseReportName)


class PeleSimulation(SimulationRunner):
    def __init__(self, parameters):
        SimulationRunner.__init__(self, parameters)
        self.type = simulationTypes.SIMULATION_TYPE.PELE

    def createSymbolicLinks(self):
        """
            Create symbolic links to Data and Documents folder if they don't exist
        """
        if not os.path.islink("Data"):
            os.system("ln -s " + self.parameters.dataFolder + " Data")
        if not os.path.islink("Documents"):
            os.system("ln -s " + self.parameters.documentsFolder + " Documents")

    def getWorkingProcessors(self):
        """
            Return the number of working processors, i.e. number of trajectories
        """
        return self.parameters.processors-1

    def getNextIterationBox(self, outputFolder, resname, topologies=None, epoch=None):
        """
            Select the box for the next epoch, currently selecting the COM of
            the cluster with max SASA

            :param outputFolder: Folder to the trajectories
            :type outputFolder: str
            :param resname: Name of the ligand in the pdb
            :type resname: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`
            :param epoch: Epoch of the trajectories to analyse
            :type epoch: int

            :returns str: -- string to be substitued in PELE control file
        """
        metrics = utilities.getMetricsFromReportsInEpoch(self.parameters.reportName, outputFolder, self.parameters.processors)
        if self.parameters.modeMovingBox.lower() == blockNames.SimulationParams.modeMovingBoxBinding:
            SASAcluster = np.argmin(metrics[:, self.parameters.columnSASA])
            if metrics[SASAcluster, self.parameters.columnSASA] < self.parameters.SASAforBox:
                self.parameters.SASAforBox = metrics[SASAcluster, self.parameters.columnSASA]
            else:
                return
        elif self.parameters.modeMovingBox.lower() == blockNames.SimulationParams.modeMovingBoxUnBinding:
            SASAcluster = np.argmax(metrics[:, self.parameters.columnSASA])
            if metrics[SASAcluster, self.parameters.columnSASA] > self.parameters.SASAforBox:
                self.parameters.SASAforBox = metrics[SASAcluster, self.parameters.columnSASA]
            else:
                return
        else:
            raise ValueError("%s should be either binding or unbinding, but %s is provided!!!" % (blockNames.SimulationParams.modeMovingBox, self.parameters.modeMovingBox))
        # If this lines are reached then a new extreme SASA value was
        # identified and we proceed to extract the corresponding center of mass
        trajNum = int(metrics[SASAcluster, -2])
        snapshotNum = int(metrics[SASAcluster, -1])
        snapshot = utilities.getSnapshots(os.path.join(outputFolder, self.parameters.trajectoryName % trajNum))[snapshotNum]
        snapshotPDB = atomset.PDB()
        snapshotPDB.initialise(snapshot, resname=resname, topology=topologies.getTopology(epoch, trajNum))
        self.parameters.boxCenter = str(snapshotPDB.getCOM())
        return

    def selectInitialBoxCenter(self, initialStructuresAsString, resname):
        """
            Select the coordinates of the first box, currently as the center of
            mass of the first initial structure provided

            :param initialStructuresAsString: String containing the files of the initial structures
            :type initialStructuresAsString: str
            :param resname: Residue name of the ligand in the system pdb
            :type resname: str

            :returns str: -- string to be substitued in PELE control file
        """
        # This is a dictionary because it's prepared to be subtitued in the PELE
        # control files (i.e. JSON format)
        try:
            initialStructuresDict = json.loads(initialStructuresAsString.split(",")[0])
            initialStruct = str(initialStructuresDict['files'][0]['path'])
        except ValueError:
            # If a json valid string is not passed, interpret the input string
            # as a path to a file
            initialStruct = initialStructuresAsString

        # Note: this lines were not changed when adding support for xtc
        # trajectories because it was assumed that initial structures would
        # still be pdbs
        PDBinitial = atomset.PDB()
        PDBinitial.initialise(initialStruct, resname=resname)
        return repr(PDBinitial.getCOM())

    def runEquilibrationPELE(self, runningControlFile):
        """
        Run a short PELE equilibration simulation

        :param runningControlFile: Path of the control file to run
        :type runningControlFile: str
        """

        self.createSymbolicLinks()
        if self.parameters.srun:
            toRun = ["srun", "-n", str(self.parameters.processors)] + self.parameters.srunParameters + [self.parameters.executable, runningControlFile]
        else:
            toRun = ["mpirun", "-np", str(self.parameters.processors)] + self.parameters.mpiParameters + [self.parameters.executable, runningControlFile]

        utilities.print_unbuffered(" ".join(toRun))
        startTime = time.time()
        proc = subprocess.Popen(toRun, shell=False, universal_newlines=True)
        (out, err) = proc.communicate()
        if out:
            print(out)
        if err:
            print(err)

        endTime = time.time()
        utilities.print_unbuffered("PELE took %.2f sec" % (endTime - startTime))

    def runSimulation(self, epoch, outputPathConstants, initialStructuresAsString, topologies, reportFileName, processManager):
        """
            Run a short PELE simulation

            :param epoch: number of the epoch
            :type epoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`
            :param reportFileName: Name of the report file
            :type reportFileName: str
            :param processManager: Object to synchronize the possibly multiple processes
            :type processManager: :py:class:`.ProcessesManager`
        """
        trajName = "".join(self.parameters.trajectoryName.split("_%d"))
        ControlFileDictionary = {"COMPLEXES": initialStructuresAsString,
                                 "PELE_STEPS": self.parameters.peleSteps,
                                 "BOX_RADIUS": self.parameters.boxRadius,
                                 "REPORT_NAME": reportFileName,
                                 "TRAJECTORY_NAME": trajName}
        self.prepareControlFile(epoch, outputPathConstants, ControlFileDictionary)
        self.createSymbolicLinks()
        runningControlFile = outputPathConstants.tmpControlFilename % epoch
        if self.parameters.srun:
            toRun = ["srun", "-n", str(self.parameters.processors)] + self.parameters.srunParameters + [self.parameters.executable, runningControlFile]
        else:
            toRun = ["mpirun", "-np", str(self.parameters.processors)] + self.parameters.mpiParameters + [self.parameters.executable, runningControlFile]
        utilities.print_unbuffered(" ".join(toRun))
        startTime = time.time()
        if self.parameters.time:
            try:
                proc = subprocess.Popen(toRun, shell=False, universal_newlines=True)
                (out, err) = proc.communicate(timeout=self.parameters.time)
            except subprocess.TimeoutExpired:
                utilities.print_unbuffered("Simulation has reached the established time limit, exiting now!!")
                out = err = None
                proc.kill()
        else:
            proc = subprocess.Popen(toRun, shell=False, universal_newlines=True)
            (out, err) = proc.communicate()
        if out:
            print(out)
        if err:
            print(err)

        endTime = time.time()
        utilities.print_unbuffered("PELE took %.2f sec" % (endTime - startTime))

    def getEquilibrationControlFile(self, peleControlFileDict):
        """
            Filter unnecessary parameters and return a minimal PELE control
            file for equilibration runs

            :param peleControlFileDict: Dictionary with pele control file options
            :type peleControlFileDict: dict

            :returns: dict -- Dictionary with pele control file options
        """
        # Set small rotations and translations
        peleControlFileDict["commands"][0]["Perturbation"]["parameters"]["translationRange"] = 0.5
        peleControlFileDict["commands"][0]["Perturbation"]["parameters"]["rotationScalingFactor"] = 0.05
        # Remove dynamical changes in control file
        peleControlFileDict["commands"][0]["PeleTasks"][0].pop("exitConditions", None)
        peleControlFileDict["commands"][0]["PeleTasks"][0].pop("parametersChanges", None)
        # Set box_radius to 2
        peleControlFileDict["commands"][0]["Perturbation"]["Box"]["fixedCenter"] = "$BOX_CENTER"
        peleControlFileDict["commands"][0]["Perturbation"]["Box"]["radius"] = 2
        # Ensure random tags exists in metrics
        metricsBlock = peleControlFileDict["commands"][0]["PeleTasks"][0]["metrics"]
        nMetrics = len(metricsBlock)
        randomIndexes = sorted([i for i in range(nMetrics) if metricsBlock[i]['type'] == "random"])
        # Delete random numbers from metrics
        for index in randomIndexes[::-1]:
            del metricsBlock[index]
        nMetrics = len(metricsBlock)
        # Add new random number to metrics
        metricsBlock.extend([{'tag': 'rand', 'type': 'random'}])
        # Add equilibration dynamical changes
        changes = {
            "doThesechanges": {
                "Perturbation::parameters": {"rotationScalingFactor": 0.050}
            },
            "ifAnyIsTrue": ["rand >= .5"],
            "otherwise": {
                "Perturbation::parameters": {"rotationScalingFactor": 0.15}
            }
            }
        peleControlFileDict["commands"][0]["PeleTasks"][0]["parametersChanges"] = [changes]

        return peleControlFileDict

    def getMetricColumns(self, JSONdict):
        """
            Extract the column of a similarity distance (RMSD or distance) from the pele control file

            :param JSONdict: Dictionary containing a parsed PELE control file
            :type JSONdict: dict

            :returns: int -- Column index of the similarity metric
        """
        hasRMSD = False
        RMSDCol = None
        distanceCol = None
        hasDistance = False
        for i, metricBlock in enumerate(JSONdict["commands"][0]["PeleTasks"][0]['metrics']):
            if 'rmsd' in metricBlock['type'].lower():
                hasRMSD = True
                RMSDCol = i+4
            elif 'distance' in metricBlock['type'].lower():
                hasDistance = True
                distanceCol = i+4
        if hasRMSD:
            return RMSDCol
        elif hasDistance:
            return distanceCol
        else:
            return None

    def calculateEquilibrationLength(self):
        """
            Calculate the number of steps for the equilibration lenght accoding
            to the available processors
        """
        # Total steps is an approximate number of total steps to produce
        totalSteps = 1000
        # Take at least 5 steps
        stepsPerProc = max(int(totalSteps/float(self.parameters.processors)), 5)
        # but no more than 50
        return min(stepsPerProc, 50)

    def equilibrate(self, initialStructures, outputPathConstants, reportFilename, outputPath, resname, processManager, topologies=None):
        """
            Run short simulation to equilibrate the system. It will run one
            such simulation for every initial structure and select appropiate
            structures to start the simulation

            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: list of str
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param reportBaseFilename: Name of the file that contains the metrics of the snapshots to cluster
            :type reportBaseFilename: str
            :param outputPath: Path where trajectories are found
            :type outputPath: str
            :param resname: Residue name of the ligand in the system pdb
            :type resname: str
            :param processManager: Object to synchronize the possibly multiple processes
            :type processManager: :py:class:`.ProcessesManager`
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`

            :returns: list --  List with initial structures
        """
        if resname is None:
            raise utilities.RequiredParameterMissingException("Resname not specified in clustering block!!!")
        newInitialStructures = []
        newStructure = []
        if self.parameters.equilibrationLength is None:
            self.parameters.equilibrationLength = self.calculateEquilibrationLength()
        trajName = "".join(self.parameters.trajectoryName.split("_%d"))
        equilibrationPeleDict = {"PELE_STEPS": self.parameters.equilibrationLength, "SEED": self.parameters.seed}
        peleControlFileDict, templateNames = utilities.getPELEControlFileDict(self.parameters.templetizedControlFile)
        peleControlFileDict = self.getEquilibrationControlFile(peleControlFileDict)
        similarityColumn = self.getMetricColumns(peleControlFileDict)

        for i, structure in enumerate(initialStructures):
            equilibrationOutput = os.path.join(outputPath, "equilibration_%d" % (i+1))
            equilibrationControlFile = outputPathConstants.tmpControlFilenameEqulibration % (i+1)
            utilities.makeFolder(equilibrationOutput)
            shutil.copyfile(structure, outputPathConstants.tmpInitialStructuresEquilibrationTemplate % (i+1))
            initialStructureString = self.createMultipleComplexesFilenames(1, outputPathConstants.tmpInitialStructuresEquilibrationTemplate, i+1, equilibration=True)
            equilibrationPeleDict["OUTPUT_PATH"] = equilibrationOutput
            equilibrationPeleDict["COMPLEXES"] = initialStructureString
            equilibrationPeleDict["BOX_CENTER"] = self.selectInitialBoxCenter(structure, resname)
            equilibrationPeleDict["BOX_RADIUS"] = 2
            equilibrationPeleDict["REPORT_NAME"] = reportFilename
            equilibrationPeleDict["TRAJECTORY_NAME"] = trajName
            for name in ["BOX_CENTER", "BOX_RADIUS"]:
                # If the template PELE control file is not templetized with the
                # box information, include it manually
                if name not in templateNames:
                    templateNames[name] = '"$%s"' % name
            utilities.print_unbuffered("Running equilibration for initial structure number %d" % (i+1))
            peleControlString = json.dumps(peleControlFileDict, indent=4)
            for key, value in templateNames.items():
                # Remove double quote around template keys, so that PELE
                # understands the options
                peleControlString = peleControlString.replace(value, "$%s" % key)

            self.makeWorkingControlFile(equilibrationControlFile, equilibrationPeleDict, peleControlString)
            self.runEquilibrationPELE(equilibrationControlFile)
            # Extract report, trajnames, metrics columns from pele control file
            reportNames = os.path.join(equilibrationOutput, self.parameters.reportName)
            trajNames = os.path.join(equilibrationOutput, self.parameters.trajectoryName)
            if len(initialStructures) == 1 and self.parameters.equilibrationMode == blockNames.SimulationParams.equilibrationLastSnapshot:
                newStructure.extend(self.selectEquilibrationLastSnapshot(self.parameters.processors, trajNames, topology=topologies.topologies[i]))
            elif self.parameters.equilibrationMode == blockNames.SimulationParams.equilibrationSelect:
                newStructure.extend(self.selectEquilibratedStructure(self.parameters.processors, similarityColumn, resname, trajNames, reportNames, topology=topologies.topologies[i]))
            elif self.parameters.equilibrationMode == blockNames.SimulationParams.equilibrationCluster:
                newStructure.extend(self.clusterEquilibrationStructures(resname, trajNames, reportNames, topology=topologies.topologies[i]))

        if len(newStructure) > self.getWorkingProcessors():
            # if for some reason the number of selected structures exceeds
            # the number of available processors, randomly sample the initial
            # structures
            newStructure = np.random.choice(newStructure, self.getWorkingProcessors(), replace=False)
        for j, struct in enumerate(newStructure):
            newStructurePath = os.path.join(equilibrationOutput, 'equilibration_struc_%d_%d.pdb' % (i+1, j+1))
            with open(newStructurePath, "w") as fw:
                fw.write(struct)
            newInitialStructures.append(newStructurePath)
        return newInitialStructures

    def clusterEquilibrationStructures(self, resname, trajWildcard, reportWildcard, topology=None):
        """
            Cluster the equilibration run

            :param resname: Name of the ligand in the pdb
            :type resname: str
            :param trajWildcard: Templetized path to trajectory files"
            :type trajWildcard: str
            :param reportWildcard: Templetized path to report files"
            :type reportWildcard: str
            :param topology: Topology for non-pdb trajectories
            :type topology: list

            :returns: list -- List with the pdb snapshots of the representatives structures
        """
        if not SKLEARN:
            raise utilities.UnsatisfiedDependencyException("No installation of scikit-learn found. Please, install scikit-learn or select a different equilibrationMode.")
        energyColumn = 3
        # detect number of trajectories available
        nTrajs = len(glob.glob(trajWildcard.rsplit("_", 1)[0]+"*"))+1
        data = []
        for i in range(1, nTrajs):
            report = utilities.loadtxtfile(reportWildcard % i)
            snapshots = utilities.getSnapshots(trajWildcard % i)
            for nSnap, (line, snapshot) in enumerate(zip(report, snapshots)):
                conformation = atomset.PDB()
                conformation.initialise(snapshot, resname=resname, topology=topology)
                com = conformation.getCOM()
                data.append([line[energyColumn], i, nSnap]+com)
        data = np.array(data)
        data = data[data[:, 0].argsort()]
        nPoints = max(self.parameters.numberEquilibrationStructures, data.shape[0]//4)
        data = data[:nPoints]
        n_clusters = min(self.parameters.numberEquilibrationStructures, data.shape[0])
        kmeans = KMeans(n_clusters=n_clusters).fit(data[:, 3:])
        utilities.print_unbuffered("Clustered equilibration output into %d clusters!" % n_clusters)
        clustersInfo = {x: {"structure": None, "minDist": 1e6} for x in range(self.parameters.numberEquilibrationStructures)}
        for conf, cluster in zip(data, kmeans.labels_):
            dist = np.linalg.norm(kmeans.cluster_centers_[cluster]-conf[3:])
            if dist < clustersInfo[cluster]["minDist"]:
                clustersInfo[cluster]["minDist"] = dist
                clustersInfo[cluster]["structure"] = tuple(conf[1:3].astype(int))
        initialStructures = []
        for cl in range(self.parameters.numberEquilibrationStructures):
            if clustersInfo[cl]["structure"] is None:
                # If a cluster has no structure assigned, skip it
                continue
            traj, snap = clustersInfo[cl]["structure"]
            conf = utilities.getSnapshots(trajWildcard % traj)[snap]
            if not isinstance(conf, basestring):
                conf = utilities.get_mdtraj_object_PDBstring(conf, topology)
            initialStructures.append(conf)
        return initialStructures

    def selectEquilibrationLastSnapshot(self, nTrajs, trajWildcard, topology=None):
        """
            Select the last snapshot of each trajectory as a  representative
            initial structure from the equilibration run

            :param nTrajs: Number of trajectories
            :type nTrajs: int
            :param trajWildcard: Templetized path to trajectory files"
            :type trajWildcard: str
            :param topology: Topology for non-pdb trajectories
            :type topology: list

            :returns: list -- List with the pdb snapshots of the representatives structures
        """
        initialStructures = []
        for ij in range(1, nTrajs):
            conf = utilities.getSnapshots(trajWildcard % ij)[-1]
            if not isinstance(conf, basestring):
                conf = utilities.get_mdtraj_object_PDBstring(conf, topology)
            initialStructures.append(conf)

        return initialStructures

    def selectEquilibratedStructure(self, nTrajs, similarityColumn, resname, trajWildcard, reportWildcard, topology=None):
        """
            Select a representative initial structure from the equilibration
            run

            :param nTrajs: Number of trajectories
            :type nTrajs: int
            :param similarityColumn: Column number of the similarity metric
                (RMSD or distance)
            :type similarityColumn: int
            :param resname: Name of the ligand in the pdb
            :type resname: str
            :param trajWildcard: Templetized path to trajectory files"
            :type trajWildcard: str
            :param reportWildcard: Templetized path to report files"
            :type reportWildcard: str
            :param topology: Topology for non-pdb trajectories
            :type topology: list

            :returns: list -- List with the pdb snapshots of the representatives structures
        """
        energyColumn = 3
        values = []
        indices = []
        rowIndex = 0
        if similarityColumn is None:
            RMSDCalc = RMSDCalculator.RMSDCalculator()
            initial = atomset.PDB()
        else:
            cols = sorted([energyColumn, similarityColumn])

        for i in range(1, nTrajs):
            indices.append(rowIndex)
            report = utilities.loadtxtfile(reportWildcard % i)
            if similarityColumn is None:
                snapshots = utilities.getSnapshots(trajWildcard % i)
                report_values = []
                if i == 1:
                    initial.initialise(snapshots.pop(0), resname=resname, topology=topology)
                    report_values.append([0, report[0, energyColumn]])
                    report = report[1:, :]
                for j, snap in enumerate(snapshots):
                    pdbConformation = atomset.PDB()
                    pdbConformation.initialise(snap, resname=resname, topology=topology)
                    report_values.append([RMSDCalc.computeRMSD(initial, pdbConformation), report[j, energyColumn]])
            else:
                report_values = report[:, cols]
            values.extend(report_values)
            rowIndex += len(report_values)

        values = np.array(values)
        if energyColumn > similarityColumn or similarityColumn is None:
            similarityColumn, energyColumn = list(range(2))
        else:
            energyColumn, similarityColumn = list(range(2))
        maxEnergy = values.max(axis=0)[energyColumn]
        # Substract the max value so all values will be negative (avoid sign problems)
        values[:, energyColumn] -= maxEnergy
        minEnergy = values.min(axis=0)[energyColumn]
        # Normalize to the minimum value
        values[:, energyColumn] /= minEnergy
        freq, bins = np.histogram(values[:, similarityColumn])
        freq = np.divide(freq, np.float(np.sum(freq)))
        # Get center of the bins
        centers = (bins[:-1]+bins[1:])/2.0
        for row in values:
            row[similarityColumn] = freq[np.argmin(np.abs(row[similarityColumn] - centers))]
        distance = values.sum(axis=1)
        indexSelected = np.argmax(distance)
        # Get trajectory and snapshot number from the index selected
        trajNum = len(indices)-1
        for i, index in enumerate(indices):
            if index == indexSelected:
                trajNum = i
                break
            elif indexSelected < index:
                trajNum = i-1
                break
        snapshotNum = indexSelected-indices[trajNum]
        # trajectories are 1-indexed
        trajNum += 1
        # return as list for compatibility with selectEquilibrationLastSnapshot
        conf = utilities.getSnapshots(trajWildcard % trajNum)[snapshotNum]
        if not isinstance(conf, basestring):
            conf = utilities.get_mdtraj_object_PDBstring(conf, topology)
        return [conf]

    def createMultipleComplexesFilenames(self, numberOfSnapshots, tmpInitialStructuresTemplate, iteration, equilibration=False):
        """
            Creates the string to substitute the complexes in the PELE control file

            :param numberOfSnapshots: Number of complexes to write
            :type numberOfSnapshots: int
            :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
            :type tmpInitialStructuresTemplate: str
            :param iteration: Epoch number
            :type iteration: int
            :param equilibration: Flag to mark whether the complexes are part of an
                equilibration run
            :type equilibration: bool

            :returns: str -- jsonString to be substituted in PELE control file
        """
        jsonString = ["\n"]
        for i in range(numberOfSnapshots-1):
            if equilibration:
                jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration)) + ",\n")
            else:
                jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, i)) + ",\n")
        if equilibration:
            jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration)))
        else:
            jsonString.append(constants.inputFileTemplate % (tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1)))
        return "".join(jsonString)


class MDSimulation(SimulationRunner):
    def __init__(self, parameters):
        SimulationRunner.__init__(self, parameters)
        self.type = simulationTypes.SIMULATION_TYPE.MD
        self.antechamberTemplate = constants.AmberTemplates.antechamberTemplate
        self.parmchkTemplate = constants.AmberTemplates.parmchk2Template
        self.tleapTemplate = constants.AmberTemplates.tleapTemplate
        self.prmtopFiles = []
        self.restart = False

        if not OPENMM:
            raise utilities.UnsatisfiedDependencyException("No installation of OpenMM found. Please, install OpenMM to run MD simulations.")
        if not self.checkAmbertools():
            raise utilities.UnsatisfiedDependencyException("No installation of AmberTools found. Please, install AmberTools to run MD simulations.")

    def checkAmbertools(self):
        try:
            subprocess.Popen(['antechamber'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return True
        except OSError:
            # If antechamber is not defined (i.e. ambertools not available)
            # popen raises an OSError
            return False

    def getWorkingProcessors(self):
        """
            Return the number of working processors, i.e. number of trajectories
        """
        return self.parameters.processors

    def equilibrate(self, initialStructures, outputPathConstants, reportFilename, outputPath, resname, processManager, topologies=None):
        """
            Run short simulation to equilibrate the system. It will run one
            such simulation for every initial structure
            (some of the arguments are not needed, such as (reportFilename, outputPath and topology) but they are kept
            to follow the same structure has the Super class definition)

            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: list of str
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param reportBaseFilename: Name of the file that contains the metrics of the snapshots to cluster
            :type reportBaseFilename: str
            :param outputPath: Path where trajectories are found
            :type outputPath: str
            :param resname: Residue name of the ligand in the system pdb
            :type resname: str
            :param processManager: Object to synchronize the possibly multiple processes
            :type processManager: :py:class:`.ProcessesManager`
            :param topology: Topology object
            :type topology: :py:class:

            :returns: list -- List with initial structures
        """
        if self.parameters.trajsPerReplica*processManager.id >= len(initialStructures):
            # Only need to launch as many simulations as initial structures
            # synchronize the replicas that will not run equilibration with the
            # replicas that will do, i.e with the synchronize in the middle of
            # this method
            processManager.barrier()
            return []
        initialStructures = processManager.getStructureListPerReplica(initialStructures, self.parameters.trajsPerReplica)
        # the new initialStructures list contains tuples in the form (i,
        # structure) where i is the index of structure in the original list
        newInitialStructures = []
        solvatedStrcutures = []
        equilibrationFiles = []
        equilibrationOutput = outputPathConstants.equilibrationDir
        utilities.makeFolder(equilibrationOutput)
        # AmberTools generates intermediate files in the current directory, change to the tmp folder
        workingdirectory = os.getcwd()
        os.chdir(outputPathConstants.tmpFolder)
        temporalFolder = os.getcwd()
        ligandPDB = self.extractLigand(initialStructures[0][1], resname, "", processManager.id)
        ligandmol2 = "%s.mol2" % resname
        ligandfrcmod = "%s.frcmod" % resname
        Tleapdict = {"RESNAME": resname, "MOL2": ligandmol2, "FRCMOD": ligandfrcmod, "DUM": ""}
        antechamberDict = {"LIGAND": ligandPDB, "OUTPUT": ligandmol2, "CHARGE": self.parameters.ligandCharge}
        parmchkDict = {"MOL2": ligandmol2, "OUTPUT": ligandfrcmod}
        if processManager.isMaster() and not self.parameters.customparamspath and resname is not None:
            self.prepareLigand(antechamberDict, parmchkDict)
        amber_file_path = ""
        # Change the Mol2 and Frcmod path to the new user defined path
        if resname is None:
            self.tleapTemplate = constants.AmberTemplates.tleapTemplatenoLigand
        if self.parameters.customparamspath and resname is not None:
            if not os.path.exists(self.parameters.customparamspath):
                # As the working directory has changed, eval if the given path is a global or relative one.
                self.parameters.customparamspath = os.path.join(workingdirectory, self.parameters.customparamspath)
                if not os.path.exists(self.parameters.customparamspath):
                    raise FileNotFoundError("No custom parameters found in the given path")
            Tleapdict["MOL2"] = os.path.join(self.parameters.customparamspath, Tleapdict["MOL2"])
            Tleapdict["FRCMOD"] = os.path.join(self.parameters.customparamspath, Tleapdict["FRCMOD"])
            amber_file_path = self.parameters.customparamspath
        if self.parameters.boxCenter or self.parameters.cylinderBases:
            with open(os.path.join(amber_file_path, "%s.prep" % constants.AmberTemplates.DUM_res), "w") as fw:
                fw.write(constants.AmberTemplates.DUM_prep)
            with open(os.path.join(amber_file_path, "%s.frcmod" % constants.AmberTemplates.DUM_res), "w") as fw:
                fw.write(constants.AmberTemplates.DUM_frcmod)

        processManager.barrier()
        if self.parameters.constraints is not None:
            prev_constraints = None
        new_constraints = None
        for i, structure in initialStructures:
            TleapControlFile = "tleap_equilibration_%d.in" % i
            pdb = PDBLoader.PDBManager(structure, resname)
            constraints_map = pdb.preparePDBforMD(constraints=self.parameters.constraints, boxCenter=self.parameters.boxCenter, cylinderBases=self.parameters.cylinderBases)
            if constraints_map is not None:
                new_constraints = updateConstraints(self.parameters.constraints, constraints_map)
                if prev_constraints is None:
                    prev_constraints = new_constraints
                assert new_constraints == prev_constraints
                prev_constraints = new_constraints
            structure = pdb.writeAll(outputpath=temporalFolder, outputname="initial_%d.pdb" % i)
            prmtop = os.path.join(workingdirectory, outputPathConstants.topologies, "system_%d.prmtop" % i)
            inpcrd = os.path.join(workingdirectory, equilibrationOutput, "system_%d.inpcrd" % i)
            finalPDB = os.path.join(workingdirectory, equilibrationOutput, "system_%d.pdb" % i)
            Tleapdict["FORCEFIELD"] = constants.AmberTemplates.forcefields[self.parameters.forcefield]

            Tleapdict["BOXSIZE"] = pdb.compute_water_box(waterBoxSize=self.parameters.waterBoxSize,
                                                         boxCenter=self.parameters.boxCenter,
                                                         boxRadius=self.parameters.boxRadius)
            Tleapdict["COMPLEX"] = structure
            Tleapdict["PRMTOP"] = prmtop
            Tleapdict["INPCRD"] = inpcrd
            Tleapdict["SOLVATED_PDB"] = finalPDB
            Tleapdict["BONDS"] = pdb.getDisulphideBondsforTleapTemplate()
            Tleapdict["MODIFIED_RES"] = pdb.getModifiedResiduesTleapTemplate()
            if self.parameters.boxCenter or self.parameters.cylinderBases:
                Tleapdict["DUM"] = "loadamberprep %s.prep\nloadamberparams %s.frcmod\n" % (constants.AmberTemplates.DUM_res, constants.AmberTemplates.DUM_res)
            self.makeWorkingControlFile(TleapControlFile, Tleapdict, self.tleapTemplate)
            self.runTleap(TleapControlFile)
            shutil.copy("leap.log", os.path.join(workingdirectory, equilibrationOutput, "leap_%d.log" % i))
            solvatedStrcutures.append(finalPDB)
            if not os.path.isfile(inpcrd):
                raise FileNotFoundError("Error While running Tleap, check %s/leap_%d.log for more information." %
                                        (os.path.join(workingdirectory, equilibrationOutput), i))
            self.prmtopFiles.append(prmtop)
            equilibrationFiles.append((prmtop, inpcrd))
        if new_constraints is not None:
            self.parameters.constraints = new_constraints
        assert len(equilibrationFiles) == len(initialStructures), "Equilibration files and initial structures don't match"
        assert len(equilibrationFiles) <= self.parameters.trajsPerReplica, "Too many equilibration structures per replica"
        os.chdir(workingdirectory)
        pool = mp.Pool(len(equilibrationFiles))
        workers = []
        startTime = time.time()
        utilities.print_unbuffered("Equilibrating System")
        for i, equilibrationFilePair in enumerate(equilibrationFiles):
            reportName = os.path.join(equilibrationOutput, "equilibrated_system_%d.pdb" % (i+processManager.id*self.parameters.trajsPerReplica))
            workers.append(pool.apply_async(sim.runEquilibration, args=(equilibrationFilePair, reportName, self.parameters, i)))

        for worker in workers:
            newInitialStructures.append(worker.get())
        endTime = time.time()
        utilities.print_unbuffered("Equilibration took %.2f sec" % (endTime - startTime))
        return newInitialStructures

    def runTleap(self, TleapControlFile):
        """
        Method that runs the Tleap software form Ambertools

        :param TleapControlFile: Path to the Tleap.in file
        :type TleapControlFile: str

        """
        tleapCommand = "tleap -f %s" % TleapControlFile
        print("System Preparation")
        startTime = time.time()
        proc = subprocess.Popen(tleapCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
        (out, err) = proc.communicate()
        print(out)
        if err:
            print("Error Found: %s" % err)
            raise utilities.UnsatisfiedDependencyException("Error Runing Tleap. Please check your installation of Ambertools.")
        endTime = time.time()
        print("System preparation took %.2f sec" % (endTime - startTime))

    def extractLigand(self, PDBtoOpen, resname, outputpath, id_replica):
        """
            Extracts the ligand from a given PDB

            :param PDBtoOpen: string with the pdb to prepare
            :type PDBtoOpen: str
            :param resname: string with the code of the ligand
            :type resname: str
            :param outputPath: Path where the pdb is written
            :type outputPath: str
            :param id_replica: Id of the current replica
            :type id_replica: int

            :returns: str -- string with the ligand pdb
        """
        ligandpdb = os.path.join(outputpath, "raw_ligand.pdb")
        if resname is None:
            return ""
        if id_replica:
            return ligandpdb
        with open(ligandpdb, "w") as out:
            with open(PDBtoOpen, "r") as inp:
                for line in inp:
                    if resname in line:
                        out.write(line)
        return ligandpdb

    def prepareLigand(self, antechamberDict, parmchkDict):
        """
        Runs antechamber and parmchk2 to obtain the mol2 and frcmod of the ligand

        :param antechamberDict: Dictonary containing the parameters to substitute in the antechamber command
        :type antechamberDict: dict
        :param parmchkDict: Dictonary containing the parameters to substitute in the parmchk2 command
        :type parmchkDict: dict
        """
        antechamberCommand = string.Template(self.antechamberTemplate)
        antechamberCommand = antechamberCommand.substitute(antechamberDict)
        parmchkCommand = string.Template(self.parmchkTemplate)
        parmchkCommand = parmchkCommand.substitute(parmchkDict)
        print(antechamberCommand)
        startTime = time.time()
        proc = subprocess.Popen(antechamberCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
        (out, err) = proc.communicate()
        if out:
            print(out)
        if err:
            print("Error Found: %s" % err)
            raise utilities.UnsatisfiedDependencyException("Error Runing Antechamber. Please check your installation of Ambertools.")
        print(parmchkCommand)
        proc = subprocess.Popen(parmchkCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
        (out, err) = proc.communicate()
        if out:
            print(out)
        if err:
            print("Error Found: %s" % err)
            raise utilities.UnsatisfiedDependencyException("Error Runing Parmchk2. Please check your installation of Ambertools.")
        endTime = time.time()
        print("Ligand preparation took %.2f sec" % (endTime - startTime))

    def runSimulation(self, epoch, outputPathConstants, initialStructuresAsString, topologies, reportFileName, processManager):
        """
            Run a MD simulation using OpenMM

            :param epoch: number of the epoch
            :type epoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`
            :param reportFileName: Name of the report file
            :type reportFileName: str
            :param processManager: Object to synchronize the possibly multiple processes
            :type processManager: :py:class:`.ProcessesManager`
        """
        outputDir = outputPathConstants.epochOutputPathTempletized % epoch
        structures_to_run = initialStructuresAsString.split(":")
        if self.restart:
            if epoch == 0:
                # if the epoch is 0 the original equilibrated pdb files are taken as intial structures
                equilibrated_structures = glob.glob(os.path.join(outputPathConstants.topologies, "top*pdb"))
                structures_to_run = sorted(equilibrated_structures, key=utilities.getTrajNum)
            checkpoints = glob.glob(os.path.join(outputDir, "checkpoint*.chk"))
            checkpoints = sorted(checkpoints, key=utilities.getTrajNum)
        # always read the prmtop files from disk to serve as communication
        # between diffrent processses
        prmtops = glob.glob(os.path.join(outputPathConstants.topologies, "*prmtop"))
        # sort the prmtops according to the original topology order
        self.prmtopFiles = sorted(prmtops, key=utilities.getPrmtopNum)
        # To follow the same order as PELE (important for processor mapping)
        structures_to_run = structures_to_run[1:]+[structures_to_run[0]]
        structures_to_run = [structure for i, structure in zip(range(self.parameters.processors), itertools.cycle(structures_to_run))]
        structures_to_run = processManager.getStructureListPerReplica(structures_to_run, self.parameters.trajsPerReplica)
        startingFilesPairs = [(self.prmtopFiles[topologies.getTopologyIndex(epoch, structure[0]+1)], structure[1]) for structure in structures_to_run]
        utilities.print_unbuffered("Starting OpenMM Production Run of %d steps..." % self.parameters.productionLength)
        startTime = time.time()
        pool = mp.Pool(self.parameters.trajsPerReplica)
        workers = []
        seed = self.parameters.seed + epoch * self.parameters.processors
        for i, startingFiles in enumerate(startingFilesPairs):
            checkpoint = None
            if self.restart:
                checkpoint = checkpoints[i + processManager.id * self.parameters.trajsPerReplica]
            workerNumber = i
            workers.append(pool.apply_async(sim.runProductionSimulation, args=(startingFiles, workerNumber, outputDir, seed, self.parameters, reportFileName, checkpoint, self.parameters.ligandName, processManager.id, self.parameters.trajsPerReplica, self.restart)))
        for worker in workers:
            worker.get()
        endTime = time.time()
        self.restart = False
        utilities.print_unbuffered("OpenMM took %.2f sec" % (endTime - startTime))

    def unifyReportNames(self, spawningReportName):
        """
            Ensure that the reportName in the simulation parameters is the same
            as the one provided in the spawning parameters

            :param spawningReportName: Name of the report file provided in the spawning parameters
            :type spawningReportName: str
        """
        pass

    def createMultipleComplexesFilenames(self, numberOfSnapshots, tmpInitialStructuresTemplate, iteration, equilibration=False):
        """
            Creates the string to substitute the complexes in the PELE control file

            :param numberOfSnapshots: Number of complexes to write
            :type numberOfSnapshots: int
            :param tmpInitialStructuresTemplate: Template with the name of the initial strutctures
            :type tmpInitialStructuresTemplate: str
            :param iteration: Epoch number
            :type iteration: int
            :param equilibration: Flag to mark wether the complexes are part of an
                equilibration run
            :type equilibration: bool

            :returns: str with the files to be used
        """
        initialStructures = []
        for i in range(numberOfSnapshots-1):
            initialStructures.append(tmpInitialStructuresTemplate % (iteration, i)+":")
        initialStructures.append((tmpInitialStructuresTemplate % (iteration, numberOfSnapshots-1)))
        return "".join(initialStructures)

    def checkSimulationInterrupted(self, epoch, outputpath, restart):
        """
            Check wether the simulation was interrupted before finishing

            :param epoch: Epoch number
            :type epoch: int
            :param outputpath: Simulation output path
            :type outputpath: str
            :param restart: Whether to restart a previous simulation
            :type restart: bool

            :returns: bool -- True if the simulations where interrupted
        """
        # to be implemented depending on implementation details
        simulationpath = os.path.join(outputpath, str(epoch), "checkpoint*")
        if glob.glob(simulationpath):
            # if the simulation does not have the restart option set we don't
            # want to restart
            self.restart = True and restart
            return True
        else:
            return False

    def cleanCheckpointFiles(self, checkpointDir):
        """
            Clean the restart files generated if the simulation was interrupted
            before finishing

            :param checkpointDir: Directory with the checkpoints
            :type checkpointDir: str
        """
        # to be implemented depending on implementation details
        checkpointsToRemove = glob.glob(os.path.join(checkpointDir, "checkpoint*.chk"))
        for files in checkpointsToRemove:
            os.remove(files)


class TestSimulation(SimulationRunner):
    """
        Class used for testing
    """
    def __init__(self, parameters):
        SimulationRunner.__init__(self, parameters)
        self.type = simulationTypes.SIMULATION_TYPE.TEST
        self.copied = False
        self.parameters = parameters

    def getWorkingProcessors(self):
        """
            Return the number of working processors, i.e. number of trajectories
        """
        return self.parameters.processors-1

    def runSimulation(self, epoch, outputPathConstants, initialStructuresAsString, topologies, reportFileName, processManager):
        """
            Copy file to test the rest of the AdaptivePELE procedure

            :param epoch: number of the epoch
            :type epoch: int
            :param outputPathConstants: Contains outputPath-related constants
            :type outputPathConstants: :py:class:`.OutputPathConstants`
            :param initialStructures: Name of the initial structures to copy
            :type initialStructures: str
            :param topologies: Topology object containing the set of topologies needed for the simulation
            :type topologies: :py:class:`.Topology`
            :param reportFileName: Name of the report file
            :type reportFileName: str
            :param processManager: Object to synchronize the possibly multiple processes
            :type processManager: :py:class:`.ProcessesManager`
        """
        ControlFileDictionary = {"COMPLEXES": initialStructuresAsString,
                                 "PELE_STEPS": self.parameters.peleSteps,
                                 "BOX_RADIUS": self.parameters.boxRadius}
        self.prepareControlFile(epoch, outputPathConstants, ControlFileDictionary)
        if not self.copied:
            tmp_sync = os.path.join(outputPathConstants.tmpFolder, os.path.split(processManager.syncFolder)[1])
            shutil.copytree(processManager.syncFolder, tmp_sync)
            if os.path.exists(self.parameters.destination):
                shutil.rmtree(self.parameters.destination)
            shutil.copytree(self.parameters.origin, self.parameters.destination)
            shutil.copytree(tmp_sync, processManager.syncFolder)
            self.copied = True

    def makeWorkingControlFile(self, workingControlFilename, dictionary, inputTemplate=None):
        pass


class ClusteringExitCondition:
    def __init__(self, ntrajs):
        self.clusterNum = 0
        self.ntrajs = ntrajs
        self.type = simulationTypes.EXITCONDITION_TYPE.CLUSTERING

    def checkExitCondition(self, clustering):
        """
            Iterate over all unchecked cluster and check if the exit condtion
            is met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- Returns True if the exit condition has been met
        """
        newClusterNum = clustering.getNumberClusters()
        clusterDiff = newClusterNum - self.clusterNum
        self.clusterNum = newClusterNum
        return clusterDiff < 0.1*self.ntrajs


class MetricExitCondition:
    def __init__(self, metricCol, metricValue, condition):
        self.metricCol = metricCol
        self.metricValue = metricValue
        self.type = simulationTypes.EXITCONDITION_TYPE.METRIC
        if condition == ">":
            self.condition = lambda x, y: x > y
        else:
            self.condition = lambda x, y: x < y

    def checkExitCondition(self, clustering):
        """
            Iterate over all clusters and check if the exit condtion
            is met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- Returns True if the exit condition has been met
        """
        for cluster in clustering.clusters.clusters:
            metric = cluster.getMetricFromColumn(self.metricCol)

            if metric is not None and self.condition(metric, self.metricValue):
                return True
        return False


class MetricMultipleTrajsExitCondition:
    def __init__(self, metricCol, metricValue, condition, reportWildCard, numTrajs, nProcessors):
        self.metricCol = metricCol
        self.metricValue = metricValue
        self.type = simulationTypes.EXITCONDITION_TYPE.METRICMULTIPLETRAJS
        self.numTrajs = numTrajs
        self.nProcessors = nProcessors
        self.trajsFound = 0
        if condition == ">":
            self.condition = lambda x, y: x.max() > y
        else:
            self.condition = lambda x, y: x.min() < y
        self.report = reportWildCard

    def checkExitCondition(self, outputFolder):
        """
            Iterate over all reports and check if the exit condtion
            is met

            :param clustering: Clustering object
            :type clustering: :py:class:`.Clustering`

            :returns: bool -- Returns True if the exit condition has been met
        """
        for j in range(1, self.nProcessors):
            report = utilities.loadtxtfile(os.path.join(outputFolder, self.report % j))
            if self.condition(report[:, self.metricCol], self.metricValue):
                self.trajsFound += 1
        return self.trajsFound >= self.numTrajs


class RunnerBuilder:
    def build(self, simulationRunnerBlock):
        """
            Build the selected  SimulationRunner object

            :param simulationRunnerBlock: Block of the control file
                corresponding to the simulation step
            :type simulationRunnerBlock: dict

            :returns: :py:class:`.SimulationRunner` -- SimulationRunner object
                selected
        """
        simulationType = simulationRunnerBlock[blockNames.SimulationType.type]
        paramsBlock = simulationRunnerBlock[blockNames.SimulationParams.params]
        params = SimulationParameters()
        if simulationType == blockNames.SimulationType.pele:
            params.processors = paramsBlock[blockNames.SimulationParams.processors]
            params.time = paramsBlock.get(blockNames.SimulationParams.time, None)
            params.dataFolder = paramsBlock.get(blockNames.SimulationParams.dataFolder, constants.DATA_FOLDER)
            params.documentsFolder = paramsBlock.get(blockNames.SimulationParams.documentsFolder, constants.DOCUMENTS_FOLDER)
            params.executable = paramsBlock.get(blockNames.SimulationParams.executable, constants.PELE_EXECUTABLE)
            if params.dataFolder is None or params.documentsFolder is None or params.executable is None:
                raise utilities.ImproperParameterValueException("PELE parameters not defined! Please ensure that you have defined the path to the PELE executable, the Data and Documents paths")
            params.templetizedControlFile = paramsBlock[blockNames.SimulationParams.templetizedControlFile]
            params.iterations = paramsBlock[blockNames.SimulationParams.iterations]
            params.peleSteps = paramsBlock[blockNames.SimulationParams.peleSteps]
            params.seed = paramsBlock[blockNames.SimulationParams.seed]
            params.trajectoryName = paramsBlock.get(blockNames.SimulationParams.trajectoryName)
            peleDict, _ = utilities.getPELEControlFileDict(params.templetizedControlFile)
            params.reportName, trajectoryName = utilities.getReportAndTrajectoryWildcard(peleDict)
            if params.trajectoryName is None:
                params.trajectoryName = trajectoryName
            else:
                params.trajectoryName = "_%d".join(os.path.splitext(params.trajectoryName))
            params.modeMovingBox = paramsBlock.get(blockNames.SimulationParams.modeMovingBox)
            if params.modeMovingBox is not None:
                if params.modeMovingBox.lower() == blockNames.SimulationParams.modeMovingBoxBinding:
                    params.SASAforBox = 1.0
                elif params.modeMovingBox.lower() == blockNames.SimulationParams.modeMovingBoxUnBinding:
                    params.SASAforBox = 0.0
                params.columnSASA = utilities.getSASAcolumnFromControlFile(peleDict)
            params.boxCenter = paramsBlock.get(blockNames.SimulationParams.boxCenter)
            params.boxRadius = paramsBlock.get(blockNames.SimulationParams.boxRadius, 20)
            params.runEquilibration = paramsBlock.get(blockNames.SimulationParams.runEquilibration, False)
            params.equilibrationMode = paramsBlock.get(blockNames.SimulationParams.equilibrationMode, blockNames.SimulationParams.equilibrationSelect)
            params.equilibrationLength = paramsBlock.get(blockNames.SimulationParams.equilibrationLength)
            params.numberEquilibrationStructures = paramsBlock.get(blockNames.SimulationParams.numberEquilibrationStructures, 10)
            params.srun = paramsBlock.get(blockNames.SimulationParams.srun, False)
            params.trajsPerReplica = params.processors
            params.numReplicas = 1
            params.srunParameters = paramsBlock.get(blockNames.SimulationParams.srunParameters, None)
            if params.srunParameters is not None:
                params.srunParameters = params.srunParameters.strip().split()
            else:
                params.srunParameters = []
            params.mpiParameters = paramsBlock.get(blockNames.SimulationParams.mpiParameters, None)
            if params.mpiParameters is not None:
                params.mpiParameters = params.mpiParameters.strip().split()
            else:
                params.mpiParameters = []
            exitConditionBlock = paramsBlock.get(blockNames.SimulationParams.exitCondition, None)
            if exitConditionBlock:
                exitConditionBuilder = ExitConditionBuilder()
                params.exitCondition = exitConditionBuilder.build(exitConditionBlock, params.templetizedControlFile, params.processors)

            return PeleSimulation(params)
        elif simulationType == blockNames.SimulationType.md:
            params.iterations = paramsBlock[blockNames.SimulationParams.iterations]
            params.processors = paramsBlock[blockNames.SimulationParams.processors]
            params.productionLength = paramsBlock[blockNames.SimulationParams.productionLength]
            params.seed = paramsBlock[blockNames.SimulationParams.seed]
            params.reporterFreq = paramsBlock[blockNames.SimulationParams.repoterfreq]
            params.numReplicas = paramsBlock[blockNames.SimulationParams.numReplicas]
            params.devicesPerTrajectory = paramsBlock.get(blockNames.SimulationParams.devicesPerTrajectory, 1)
            params.trajsPerReplica = int(params.processors/params.numReplicas)
            assert params.trajsPerReplica*params.numReplicas == params.processors, "Number of trajectories requested does not match the number of replicas"
            params.maxDevicesPerReplica = paramsBlock.get(blockNames.SimulationParams.maxDevicesPerReplica)
            params.runEquilibration = True
            params.equilibrationLengthNVT = paramsBlock.get(blockNames.SimulationParams.equilibrationLengthNVT, 200000)
            params.equilibrationLengthNPT = paramsBlock.get(blockNames.SimulationParams.equilibrationLengthNPT, 500000)
            params.format = paramsBlock.get(blockNames.SimulationParams.format, "xtc")
            if params.format not in constants.md_supported_formats:
                raise utilities.ImproperParameterValueException("Not supported %s format specified, supported formats are %s" % (params.format, constants.formats_md_string))
            params.timeStep = paramsBlock.get(blockNames.SimulationParams.timeStep, 2)
            params.boxRadius = paramsBlock.get(blockNames.SimulationParams.boxRadius, 20)
            params.boxCenter = paramsBlock.get(blockNames.SimulationParams.boxCenter)
            params.boxType = paramsBlock.get(blockNames.SimulationParams.boxType, blockNames.SimulationParams.sphere)
            params.cylinderBases = paramsBlock.get(blockNames.SimulationParams.cylinderBases)
            if params.boxType == blockNames.SimulationParams.cylinder and params.cylinderBases is None:
                raise utilities.RequiredParameterMissingException("To use a cylinder box you need to specify the basis points for the cylinder")
            if params.cylinderBases is not None:
                # round all coordinates to 3 decimals to fit into pdb format and
                # also ensure that they are numpy arrays
                params.cylinderBases = [np.round(el, decimals=3) for el in params.cylinderBases]
            if params.boxType is not None and params.boxType not in (blockNames.SimulationParams.sphere, blockNames.SimulationParams.cylinder):
                raise utilities.ImproperParameterValueException("Unknown %s box type, supported formats are %s" % (params.boxType, " ".join([blockNames.SimulationParams.sphere, blockNames.SimulationParams.cylinder])))
            params.ligandCharge = paramsBlock.get(blockNames.SimulationParams.ligandCharge, 1)
            params.waterBoxSize = paramsBlock.get(blockNames.SimulationParams.waterBoxSize, 8)
            params.forcefield = paramsBlock.get(blockNames.SimulationParams.forcefield, "ff99SB")
            params.nonBondedCutoff = paramsBlock.get(blockNames.SimulationParams.nonBondedCutoff, 8)
            params.Temperature = paramsBlock.get(blockNames.SimulationParams.Temperature, 300)
            params.runningPlatform = paramsBlock.get(blockNames.SimulationParams.runningPlatform, "CPU")
            params.minimizationIterations = paramsBlock.get(blockNames.SimulationParams.minimizationIterations, 2000)
            params.constraintsMin = paramsBlock.get(blockNames.SimulationParams.constraintsMin, 5)
            params.constraintsNVT = paramsBlock.get(blockNames.SimulationParams.constraintsNVT, 5)
            params.constraintsNPT = paramsBlock.get(blockNames.SimulationParams.constraintsNPT, 0.5)
            params.customparamspath = paramsBlock.get(blockNames.SimulationParams.customparamspath)
            params.ligandName = paramsBlock.get(blockNames.SimulationParams.ligandName)
            params.constraints = paramsBlock.get(blockNames.SimulationParams.constraints)
            if params.ligandName is None and params.boxCenter is not None:
                raise utilities.ImproperParameterValueException("Ligand name is necessary to establish the box")
            return MDSimulation(params)
        elif simulationType == blockNames.SimulationType.test:
            params.processors = paramsBlock[blockNames.SimulationParams.processors]
            params.destination = paramsBlock[blockNames.SimulationParams.destination]
            params.origin = paramsBlock[blockNames.SimulationParams.origin]
            params.iterations = paramsBlock[blockNames.SimulationParams.iterations]
            params.peleSteps = paramsBlock[blockNames.SimulationParams.peleSteps]
            params.seed = paramsBlock[blockNames.SimulationParams.seed]
            params.trajsPerReplica = params.processors
            params.numReplicas = 1
            return TestSimulation(params)
        else:
            sys.exit("Unknown simulation type! Choices are: " + str(simulationTypes.SIMULATION_TYPE_TO_STRING_DICTIONARY.values()))


class ExitConditionBuilder:
    def build(self, exitConditionBlock, templetizedControlFile, nProcessors):
        """
            Build the selected exit condition object

            :param exitConditionBlock: Block of the control file
                corresponding to the exit condition
            :type exitConditionBlock: dict

            :returns: :py:class:`.MetricExitCondition` -- MetricExitCondition object
                selected
        """
        exitConditionType = exitConditionBlock[blockNames.ExitConditionType.type]
        exitConditionParams = exitConditionBlock[blockNames.SimulationParams.params]
        if exitConditionType == blockNames.ExitConditionType.metric:
            # Start counting the columns by 1
            metricCol = exitConditionParams[blockNames.SimulationParams.metricCol]-1
            metricValue = exitConditionParams[blockNames.SimulationParams.exitValue]
            condition = exitConditionParams.get(blockNames.SimulationParams.condition, "<")
            if condition not in [">", "<"]:
                raise ValueError("In MetricExitCondition the parameter condition only accepts > or <, but %s was passed" % condition)
            return MetricExitCondition(metricCol, metricValue, condition)
        elif exitConditionType == blockNames.ExitConditionType.metricMultipleTrajs:
            # Start counting the columns by 1
            metricCol = exitConditionParams[blockNames.SimulationParams.metricCol]-1
            metricValue = exitConditionParams[blockNames.SimulationParams.exitValue]
            numTrajs = exitConditionParams[blockNames.SimulationParams.numTrajs]
            condition = exitConditionParams.get(blockNames.SimulationParams.condition, "<")
            if condition not in [">", "<"]:
                raise ValueError("In MetricMultipleTrajsExitCondition the parameter condition only accepts > or <, but %s was passed" % condition)
            peleControlFileDict, _ = utilities.getPELEControlFileDict(templetizedControlFile)
            reportWildCard, _ = utilities.getReportAndTrajectoryWildcard(peleControlFileDict)
            return MetricMultipleTrajsExitCondition(metricCol, metricValue, condition, reportWildCard, numTrajs, nProcessors)
        elif exitConditionType == blockNames.ExitConditionType.clustering:
            ntrajs = exitConditionParams[blockNames.SimulationParams.trajectories]
            return ClusteringExitCondition(ntrajs)
        else:
            sys.exit("Unknown exit condition type! Choices are: " + str(simulationTypes.EXITCONDITION_TYPE_TO_STRING_DICTIONARY.values()))


def updateConstraints(constraints_orig, constraints_map):
    new_const = []
    for constr in constraints_orig:
        atom1, atom2, dist = constr
        atom1 = atom1.split(":")
        atom2 = atom2.split(":")
        atom1[2] = int(atom1[2])
        atom2[2] = int(atom2[2])
        atom1[2] = constraints_map[tuple(atom1[1:])]
        atom2[2] = constraints_map[tuple(atom2[1:])]
        new_const.append([":".join([str(i) for i in atom1]), ":".join([str(i) for i in atom2]), str(dist)])
    return new_const
