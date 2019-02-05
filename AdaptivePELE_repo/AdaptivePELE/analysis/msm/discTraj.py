import os
import glob
import atomset
import numpy as np
from utilities import utilities
from constants import constants
from atomset import atomset
import adaptiveSampling

def getLastClusteringObject():
    outputPath = "."
    outputPathConstants = constants.OutputPathConstants(outputPath)
    lastEpoch = adaptiveSampling.findFirstRun(outputPath, outputPathConstants.clusteringOutputObject) - 1
    clusteringObjectPath = outputPathConstants.clusteringOutputObject % lastEpoch
    clusteringObject = utilities.readClusteringObject(clusteringObjectPath)
    return clusteringObject

class Constants():
    def __init__(self):
        self.folder = "discretized"
        self.file = "traj_%s_%d.disctraj" #epoch_trajNum
        self.traj = "*traj*.pdb"

def discretizeAllTrajs(epoch, constants, clusteringObject):
    allTrajs = glob.glob(os.path.join(epoch, constants.traj))
    allDiscTrajs = []
    trajNums = []
    for traj in allTrajs:
        discTraj = discretizeTraj(traj, clusteringObject)
        allDiscTrajs.append(discTraj)
        trajNums.append(utilities.getTrajNum(traj))
    return trajNums, allDiscTrajs

def findBelongingCluster(snapshot, clusteringObject):
    """
        Returns the closest cluster in terms of rmsd
    """
    pdb = atomset.PDB()
    pdb.initialise(snapshot, resname=clusteringObject.resname)

    clusters = clusteringObject.clusters.clusters
    minRmsd = 1e10
    minRmsd2 = minRmsd*minRmsd
    minArg = None

    distances = []
    for i, cluster in enumerate(clusters):
        #centroid distance is a minimum bound for RMSD distance
        if atomset.computeSquaredCentroidDifference(cluster.pdb, pdb) > minRmsd2:
            continue

        rmsd = atomset.computeRMSD(cluster.pdb, pdb, clusteringObject.symmetries)
        if rmsd < minRmsd:
            minRmsd = rmsd
            minRmsd2 = rmsd*rmsd
            minArg = i
    return minArg


def discretizeTraj(traj, clusteringObject):
    snapshots = utilities.getSnapshots(traj)

    discTraj = np.zeros(len(snapshots))
    for i, snapshot in enumerate(snapshots):
        discTraj[i] = findBelongingCluster(snapshot, clusteringObject)
    return discTraj

def writeDiscTrajs(epoch, trajNums, discTrajs, constants):
    for trajNum, discTraj in zip(trajNums, discTrajs):
        output = os.path.join(constants.folder, constants.file%(epoch, trajNum))
        np.savetxt(output, discTraj, fmt="%d")

def discretize():
    constants = Constants()

    clusteringObject = getLastClusteringObject()

    allFolders = os.listdir(".")
    epochFolders = [int(epoch) for epoch in allFolders if epoch.isdigit()]
    epochFolders.sort()

    utilities.cleanup(constants.folder)
    utilities.makeFolder(constants.folder)

    for folder in epochFolders:
        print "Epoch", folder
        trajNums, discTrajs = discretizeAllTrajs(epoch, constants, clusteringObject)
        writeDiscTrajs(folder, trajNums, discTrajs, constants)


if __name__ == "__main__":
    discretize()


