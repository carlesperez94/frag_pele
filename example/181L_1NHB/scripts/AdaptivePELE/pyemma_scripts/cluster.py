import os
import numpy as np
import glob
import pyemma.coordinates as coor
from pyemma.coordinates.clustering import AssignCenters

class Cluster:
    def __init__(self, numClusters, trajectoryFolder, trajectoryBasename, stride=1, alwaysCluster=True):
        """
            alwaysCluster: clusterize regardless of wether discretized/clusterCenters.dat exists or not
        """

        self.discretizedFolder = "discretized"
        self.clusterCentersFile = os.path.join(self.discretizedFolder, "clusterCenters.dat")
        self.dTrajTemplateName = os.path.join(self.discretizedFolder, "%s.disctraj")
        self.clusteringFile = "clustering_object.pkl"
        self.stride = stride
        self.trajFilenames = []
        self.dtrajs = []
        self.alwaysCluster = alwaysCluster
        self.numClusters = numClusters

        self.clusterTrajectories(trajectoryFolder, trajectoryBasename)

    
    def cluster(self, trajectories):
        """ Cluster the trajectories into numClusters clusters using kmeans
        algorithm.
        Returns a KmeansClusteringObject
        """
        return coor.cluster_kmeans(data=trajectories, k=self.numClusters, max_iter=500, stride=self.stride)

    def assignNewTrajecories(self, trajs):
        assign = AssignCenters(self.clusterCentersFile)
        dTrajs = assign.assign(trajs)
        return dTrajs

    def clusterTrajectories(self, trajectoryFolder, trajectoryBasename):
        print "Loading trajectories..."
        self.x, self.trajFilenames = loadCOMFiles(trajectoryFolder, trajectoryBasename)

        # cluster & assign
        if self.alwaysCluster or not os.path.exists(self.clusterCentersFile):
            print "Clustering data..."
            cl = self.cluster(self.x)
            makeFolder(self.discretizedFolder)
            self.writeClusterCenters(cl, self.clusterCentersFile)
            print "Assigning data..."
            self.dtrajs = cl.dtrajs[:]
        else:
            print "Assigning data (clustering exists)..."
            self.dtrajs = self.assignNewTrajecories(self.x)

        print "Writing clustering data..."
        self.writeDtrajs(self.trajFilenames, self.dtrajs, self.dTrajTemplateName)

    def writeClusterCenters(self, cl, outputFilename):
        np.savetxt(outputFilename, cl.clustercenters, fmt="%.5f %.5f %.5f")

    def writeDtrajs(self, filenames, dtrajs, filenameTemplate="%s.disctraj"):
        for filename, dtraj in zip(filenames, dtrajs):
            fname = os.path.split(filename)[-1][:-4]
            dtrajfname = filenameTemplate%(fname)
            np.savetxt(dtrajfname, dtraj, fmt="%d")


"""
#Standalone functions
"""

def loadCOMFiles(trajectoryFolder, trajectory_basename):
    trajectoryBasename = os.path.join(trajectoryFolder, trajectory_basename)

    # load traj
    files = glob.glob(trajectoryBasename)
    x = len(files)*[0]
    for i, file in enumerate(files):
        currentX = np.loadtxt(file, usecols=(1, 2, 3))
        x[i] = currentX
    if not x:
        raise ValueError("Didn't find any trajectory files in the specified path!!!")
    return x, files

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
