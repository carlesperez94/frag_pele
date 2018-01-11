import argparse
from AdaptivePELE.spawning import spawning, densitycalculator
from AdaptivePELE.clustering import clustering
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import RMSDCalculator
import AdaptivePELE.atomset.atomset as atomset
from scipy.sparse.csgraph import shortest_path
import numpy as np


class PathwayError(Exception):
    pass


def parseArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('clusteringObj', type=str)
    parser.add_argument('native', type=str)
    parser.add_argument('ntrajs', type=int)
    parser.add_argument('threshold', type=float)
    parser.add_argument('pathwayFilename', type=str, default="pathway.pdb", nargs='?')
    args = parser.parse_args()
    return args


def getOptimalCluster(clusteringObj, native, RMSDCalc):
    optimalMetric = 1000000
    optimalIndex = 0
    nativePDB = atomset.PDB()
    nativePDB.initialise(native)
    for i, cluster in enumerate(clusteringObj.clusters.clusters):
        metric = RMSDCalc.computeRMSD(cluster.pdb, nativePDB)
        if metric < optimalMetric:
            optimalMetric = metric
            optimalIndex = i
    print optimalMetric
    return optimalIndex


def createNetworkMatrix(clusteringObj, threshold, RMSDCalc):
    n = clusteringObj.clusters.getNumberClusters()
    matrix = np.zeros((n, n))
    for i, cluster in enumerate(clusteringObj.clusters.clusters):
        for j in range(i, n):
            if i == j:
                continue
            cluster2 = clusteringObj.getCluster(j)
            rmsd =  RMSDCalc.computeRMSD(cluster.pdb, cluster2.pdb)
            if rmsd > threshold:
                matrix[i, j] = matrix[j, i] = np.inf
            else:
                matrix[i, j] = matrix[j, i] = RMSDCalc.computeRMSD(cluster.pdb, cluster2.pdb)
    return matrix


def obtainShortestPath(distanceMatrix):
    foo, predecessors = shortest_path(distanceMatrix, return_predecessors=True)
    return predecessors


def createPathway(initial_cluster, final_cluster, predecessors):
    pathway = [final_cluster]
    i, j = initial_cluster, final_cluster
    while (i != j) and j >= 0:
        j = predecessors[i, j]
        pathway.insert(0, j)
    if pathway[0] < 0:
        # The network defined by the distanceMatrix is not connected
        raise PathwayError("Unable to find a path since network is not "
                           "connected, use a bigger threshold")

    return pathway


def removeRemarksPDB(native, modelNum):
    pdbStr = ["MODEL %d\n" % modelNum]
    with open(native) as f:
        for line in f:
            if line.startswith("MODEL"):
                break
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                pdbStr.append(line)
        for line in f:
            pdbStr.append(line)
    pdbStr.append("ENDMDL\n")
    return ''.join(pdbStr)


def writePathwayTrajectory(ClOrd, pathway, filename, native):
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
    pathwayFile.write(removeRemarksPDB(native, i+2))
    pathwayFile.close()


def main(args):

    # Parameters
    clusteringObj = utilities.readClusteringObject(args.clusteringObj)
    native = args.native
    pathwayFilename = args.pathwayFilename
    ntrajs = args.ntrajs
    threshold = args.threshold
    RMSDCalc = RMSDCalculator.RMSDCalculator(clusteringObj.symmetries)

    # use graph algorithm to establish a path
    initial_cluster = 0
    final_cluster = getOptimalCluster(clusteringObj, native, RMSDCalc)
    distanceMatrix = createNetworkMatrix(clusteringObj, threshold, RMSDCalc)
    predecessors = obtainShortestPath(distanceMatrix)
    pathway = createPathway(initial_cluster, final_cluster, predecessors)
    print "Pathway clusters:"
    print pathway

    # write pathway into a single trajectory
    writePathwayTrajectory(clusteringObj, pathway, pathwayFilename, native)

    # create clustering object with only the pathway clusters
    ClPath = clustering.Clustering()
    ClPath.clusters.clusters = map(lambda x: clusteringObj.clusters.clusters[x],
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

if __name__ == "__main__":
    args = parseArgs()
    main(args)
