import time
import os
import pickle
import argparse
import json
from AdaptivePELE.atomset import RMSDCalculator
import AdaptivePELE.atomset.atomset as atomset
from AdaptivePELE.clustering import clustering, clusteringTypes, thresholdcalculator
from AdaptivePELE.spawning import spawning, densitycalculator
from AdaptivePELE.analysis import analyseClustering
from AdaptivePELE.simulation import simulationrunner
import AdaptivePELE.adaptiveSampling as adaptiveSampling


def parseArgs():
    parser = argparse.ArgumentParser(description="Perform two runs of adaptive"
                                     "sampling ")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    jsonFile = open(args.controlFile, "r").read()
    parsedJSON = json.loads(jsonFile)
    return parsedJSON


class OrderedContactsClustering(clustering.Clustering):
    def __init__(self, thresholdCalculator, resname=None,
                 reportBaseFilename=None, columnOfReportFile=None,
                 contactThresholdDistance=8, symmetries=[]):
        clustering.Clustering.__init__(self, resname, reportBaseFilename,
                                       columnOfReportFile,
                                       contactThresholdDistance)
        self.type = clusteringTypes.CLUSTERING_TYPES.contacts
        self.thresholdCalculator = thresholdCalculator
        self.symmetries = symmetries
        self.initialPDB = None
        self.maxThreshold = self.thresholdCalculator.getMaxThreshold()
        self.distancesList = []
        self.RMSDCalculator = RMSDCalculator.RMSDCalculator(symmetries)

    def getOptimalMetric(self):
        optimalMetric = 100
        optimalMetricIndex = 0
        for i, cluster in enumerate(self.clusters.clusters):
            if cluster.getMetric() < optimalMetric:
                optimalMetric = cluster.getMetric()
                optimalMetricIndex = i
        return optimalMetricIndex

    def addSnapshotToCluster(self, snapshot, metrics=[], col=None):
        pdb = atomset.PDB()
        pdb.initialise(snapshot, resname=self.resname)
        if not self.initialPDB:
            self.initialPDB = pdb

        contacts = pdb.countContacts(self.resname,
                                     self.contactThresholdDistance)
        numberOfLigandAtoms = pdb.getNumberOfAtoms()
        contactsPerAtom = float(contacts)/numberOfLigandAtoms

        threshold = self.thresholdCalculator.calculate(contactsPerAtom)
        threshold2 = threshold**2
        threshold2 = 10

        initial_scd = atomset.computeSquaredCentroidDifference(self.initialPDB,
                                                               pdb)
        snapshotPosition, closerClusterInd, closerClusterRMSD = self.getCloserCluster(initial_scd, threshold2, pdb)
        try:
            clusterToAssign = self.clusters.clusters[closerClusterInd]
            if closerClusterRMSD < clusterToAssign.threshold:
                clusterToAssign.addElement(metrics)
                return
        except TypeError:
            # When the first snapshot is read the index of the closest cluster
            # center is set to None
            pass

        # if made it here, the snapshot was not added into any cluster

        cluster = clustering.Cluster(pdb, thresholdRadius=threshold,
                                     contacts=contactsPerAtom, metrics=metrics,
                                     metricCol=col, density=1)
        if snapshotPosition is not None:
            self.clusters.insertCluster(snapshotPosition, cluster)
            self.distancesList.insert(snapshotPosition, initial_scd)
        else:
            self.clusters.addCluster(cluster)
            self.distancesList.append(initial_scd)

    def getCloserCluster(self, initial_scd, threshold2, pdb):
        minimumCloseCluster = 0
        nclusters = len(self.clusters.clusters)
        snapshotPosition = None
        for distance_initial in self.distancesList:
            if (snapshotPosition is None and distance_initial > initial_scd):
                snapshotPosition = minimumCloseCluster
            if (abs(distance_initial-initial_scd) < threshold2):
                break
            minimumCloseCluster += 1
        closerClusterInd = None
        closerClusterRMSD = 100

        while (minimumCloseCluster < nclusters):
            cluster = self.clusters.clusters[minimumCloseCluster]
            distance = self.distancesList[minimumCloseCluster]
            if (snapshotPosition is None and distance > initial_scd):
                snapshotPosition = minimumCloseCluster
            if (abs(distance-initial_scd) > threshold2):
                break
            clusterRMSD = self.RMSDCalculator.computeRMSD(cluster.pdb, pdb)
            if (clusterRMSD < closerClusterRMSD):
                closerClusterRMSD = clusterRMSD
                closerClusterInd = minimumCloseCluster
            minimumCloseCluster += 1
        return snapshotPosition, closerClusterInd, closerClusterRMSD


def printNumberSnapshotsEpoch(paths_report, i):
    trajs = clustering.getAllTrajectories(paths_report)
    total_snapshots = 0
    for traj in trajs:
        for line in open(traj, "r"):
            total_snapshots += 1
        total_snapshots -= 1
    print "Total snapsthots for epoch %d: %d" % (i, total_snapshots)


def clusterSnapshotsEpoch(ClOrd, path, i, ntrajs, clusteringObject):
    startTimeOrd = time.time()
    ClOrd.cluster(path)
    endTimeOrd = time.time()
    print "Total time of clustering ordered contacts, epoch %d: %.6f" % (i, endTimeOrd-startTimeOrd)
    print "Number of clusters ordered contacts epoch %d: %d" % (i, len(ClOrd.clusters.clusters))
    # degeneraciesOrd = spawningObject.calculate(ClOrd.clusters.clusters, ntrajs, {})
    degeneraciesOrd = [0] * len(ClOrd.clusters.clusters)
    ClOrd.writeOutput("clsummary", degeneraciesOrd, clusteringObject, False)
    os.remove("clsummary/summary.txt")


def createPathway(initial_cluster, final_cluster, ClOrd, distanceThreshold):
    pathway = [initial_cluster]
    rowind = final_cluster
    while (rowind > 0):
        pathway.insert(1, rowind)
        clusterInit = ClOrd.clusters.clusters[rowind]
        cluster_distance = ClOrd.distancesList[rowind]
        minimumCloseCluster = rowind-1
        while (minimumCloseCluster > 0):
            distance_initial = ClOrd.distancesList[minimumCloseCluster]
            if (abs(distance_initial-cluster_distance) < distanceThreshold):
                break
            minimumCloseCluster -= 1

        minimumRMSD = 1000
        closerCluster = rowind-1
        while (minimumCloseCluster > 0):
            cluster = ClOrd.clusters.clusters[minimumCloseCluster]
            distance = ClOrd.distancesList[minimumCloseCluster]
            if (abs(distance-cluster_distance) > distanceThreshold):
                break
            clusterRMSD = ClOrd.RMSDCalculator.computeRMSD(cluster.pdb,
                                                           clusterInit.pdb)
            if (clusterRMSD < minimumRMSD):
                minimumRMSD = clusterRMSD
                closerCluster = minimumCloseCluster
            minimumCloseCluster -= 1
        rowind = closerCluster
    return pathway


def writePathwayTrajectory(ClOrd, pathway, filename):
    pathwayFile = open(filename, "w")
    pathwayFile.write("REMARK 000 File created using PELE++\n")
    pathwayFile.write("REMARK 000 Pathway trajectory created using findPathway program\n")
    pathwayFile.write("REMARK 000 List of cluster belonging to the pathway %s\n" % ' '.join(map(str, pathway)))
    for i, step_cluster in enumerate(pathway):
        cluster = ClOrd.clusters.clusters[step_cluster]
        pathwayFile.write("MODEL %d\n" % (i+1))
        pdbStr = cluster.pdb.pdb
        pdbList = pdbStr.split("\n")
        for line in pdbList:
            line = line.strip()
            # Avoid writing previous REMARK block
            if line.startswith("REMARK ") or line.startswith("MODEL "):
                continue
            elif line:
                pathwayFile.write(line+"\n")
        pathwayFile.write("ENDMDL\n")
    pathwayFile.close()


def plotMetricsDegeneracy(ClPath, resname, degeneracies):
    title = "Degeneracies along pathway"
    title2 = "Metrics along pathway"
    comCoord, metrics = analyseClustering.extractCOMMatrix(ClPath.clusters.clusters,
                                                           resname)[:2]
    analyseClustering.plotClusters(comCoord, degeneracies, title)
    analyseClustering.plotClusters(comCoord, metrics, title2)
    import matplotlib.pyplot as plt
    plt.show()


def clusterTrajectories(resname, trajFolder, clusteringObject,
                        clusteringThreshold, reportBaseFilename="report",
                        symmetries=[]):
    thresholdCalculatorBuilder = thresholdcalculator.ThresholdCalculatorBuilder()
    thresholdCalculator = thresholdCalculatorBuilder.build({
            "thresholdCalculator": {
                "type": "constant",
                "params": {
                    "value": clusteringThreshold
                }
            }
    })
    ClOrd = OrderedContactsClustering(thresholdCalculator, resname=resname,
                                      reportBaseFilename=reportBaseFilename,
                                      columnOfReportFile=4,
                                      symmetries=symmetries)

    if os.path.exists(clusteringObject):
        with open(clusteringObject, "r") as f:
            ClOrd = pickle.load(f)
    else:
        allFolders = os.listdir(trajFolder)
        Epochs = [epoch for epoch in allFolders if epoch.isdigit()]
        for i in range(len(Epochs)):
            path = [trajFolder+"/%d/traj*" % i]
            paths_report = [trajFolder+"/%d/report*" % i]
            printNumberSnapshotsEpoch(paths_report, i)
            clusterSnapshotsEpoch(ClOrd, path, i, 100, clusteringObject)
    return ClOrd, thresholdCalculator


def makeNewControlFile(degeneracies, ClPath, templetizedInitialName,
                       secondControlFileTemp, secondControlFile):
    initialConfList = []
    for j, nproc in enumerate(degeneracies):
        if nproc:
            for i in range(nproc):
                nameInitialStructure = templetizedInitialName % (j, i)
                ClPath.clusters.clusters[j].writePDB(nameInitialStructure)
                initialConfList.append(nameInitialStructure)
    controlFileDict = {"INITIAL": initialConfList}
    simulationParameters = simulationrunner.SimulationParameters()
    simulationParameters.templetizedControlFile = secondControlFileTemp
    simulationRunner = simulationrunner.SimulationRunner(simulationParameters)
    simulationRunner.makeWorkingControlFile(secondControlFile, controlFileDict)


def main(jsonBlock):
    # Parameters
    firstControlFile = jsonBlock["firstControlFile"]
    ntrajs = jsonBlock["ntrajs"]
    resname = jsonBlock["resname"]
    symmetries = jsonBlock["symmetries"]
    trajFolder = jsonBlock["trajFolder"]
    clusteringObject = jsonBlock["clusteringObject"]
    clusteringThreshold = jsonBlock["clusteringThreshold"]
    pathwayFilename = jsonBlock["pathwayFilename"]
    templetizedInitialName = jsonBlock["templetizedInitialName"].encode()
    secondControlFileTemp = jsonBlock["secondControlFileTemp"]
    secondControlFile = jsonBlock["secondControlFile"]
    distanceThreshold = jsonBlock["distanceThreshold"]

    if firstControlFile:
        # Run first adaptive
        adaptiveSampling.main(firstControlFile)

    # Cluster trajectories

    ClOrd, thresholdCalculator = clusterTrajectories(resname, trajFolder,
                                                     clusteringObject,
                                                     clusteringThreshold,
                                                     symmetries=symmetries)

    # use graph algorithm to establish a path
    initial_cluster = 0
    final_cluster = ClOrd.getOptimalMetric()
    pathway = createPathway(initial_cluster, final_cluster, ClOrd,
                            distanceThreshold)

    # write pathway into a single trajectory
    writePathwayTrajectory(ClOrd, pathway, pathwayFilename)

    # create clustering object with only the pathway clusters
    ClPath = OrderedContactsClustering(thresholdCalculator, resname=resname,
                                       reportBaseFilename="report",
                                       columnOfReportFile=4,
                                       symmetries=symmetries)
    ClPath.clusters.clusters = map(lambda x: ClOrd.clusters.clusters[x],
                                   pathway)

    # spawning along the trajectory
    spawningParams = spawning.SpawningParams()
    densityCalculatorBuilder = densitycalculator.DensityCalculatorBuilder()
    densityCalculator = densityCalculatorBuilder.build({})
    spawningPathway = spawning.InverselyProportionalToPopulationCalculator(densityCalculator)
    # Set a least 1 processors from the extrems of the path
    degeneracies = spawningPathway.calculate(ClPath.clusters.clusters,
                                              ntrajs-2, spawningParams)
    degeneracies[0] += 1
    degeneracies[-1] += 1
    print "degeneracies over pathway:"
    print degeneracies
    print ""
    plots = False
    if plots:
        plotMetricsDegeneracy(ClPath, resname, degeneracies)

    # Prepare second adaptive
    makeNewControlFile(degeneracies, ClPath, templetizedInitialName,
                       secondControlFileTemp, secondControlFile)
    adaptiveSampling.main(secondControlFile)

if __name__ == "__main__":
    jsonBlock = parseArgs()
    main(jsonBlock)
