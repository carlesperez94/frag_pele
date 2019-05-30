from __future__ import absolute_import, division, print_function, unicode_literals
import os
import numpy as np
import glob
import pyemma.coordinates as coor
from pyemma.coordinates.clustering import AssignCenters
import scipy


class Cluster:
    def __init__(self, numClusters, trajectoryFolder, trajectoryBasename, stride=1, alwaysCluster=True):
        """
            alwaysCluster: clusterize regardless of whether discretized/clusterCenters.dat exists or not
        """

        self.discretizedFolder = "discretized"
        self.clusterCentersFile = os.path.join(self.discretizedFolder, "clusterCenters.dat")
        self.clusterCenters = np.array([])
        self.dTrajTemplateName = os.path.join(self.discretizedFolder, "%s.disctraj")
        self.clusteringFile = "clustering_object.pkl"
        self.stride = stride
        self.trajFilenames = []
        self.dtrajs = []
        self.alwaysCluster = alwaysCluster
        self.numClusters = numClusters
        self.trajectoryFolder = trajectoryFolder
        self.trajectoryBasename = trajectoryBasename
        self.x = []

    def cluster(self, trajectories):
        """ Cluster the trajectories into numClusters clusters using kmeans
        algorithm.
        Returns a KmeansClusteringObject
        """
        return coor.cluster_kmeans(data=trajectories, k=self.numClusters, max_iter=500, stride=self.stride)

    def assignNewTrajectories(self, trajs):
        # wrap the clusterCentersFile argument in a str call to pass pyemma
        # assign check if isinstance of str
        assign = AssignCenters(str(self.clusterCentersFile))
        dTrajs = assign.assign(trajs)
        return dTrajs

    def clusterTrajectories(self):
        print("Loading trajectories...")
        self.x, self.trajFilenames = loadTrajFiles(self.trajectoryFolder, self.trajectoryBasename)

        # cluster & assign
        if self.alwaysCluster or not os.path.exists(self.clusterCentersFile):
            print("Clustering data...")
            cl = self.cluster(self.x)  # cl: pyemma's clusteringObject
            makeFolder(self.discretizedFolder)
            self.clusterCenters = cl.clustercenters
            self._writeClusterCenters(self.clusterCenters, self.clusterCentersFile)
            print("Assigning data...")
            self.dtrajs = cl.dtrajs[:]
        else:
            print("Assigning data (clustering exists)...")
            self.clusterCenters = np.loadtxt(self.clusterCentersFile)
            self.dtrajs = self.assignNewTrajectories(self.x)

        print("Writing clustering data...")
        self._writeDtrajs(self.trajFilenames, self.dtrajs, self.dTrajTemplateName)

    def eliminateLowPopulatedClusters(self, clusterCountsThreshold, tau=None):
        if self.dtrajs == []:
            print("Call clusterTrajectories() first!")
            return

        dtrajs = np.array(self.dtrajs).copy()
        try:
            dtrajs = np.concatenate(dtrajs[:, :-tau])
        except:
            dtrajs = np.concatenate(dtrajs)

        dummy = np.zeros(dtrajs.size)
        data = np.ones(dtrajs.size)
        # using sparse is fast and sucint
        counts = np.ravel(scipy.sparse.coo_matrix((data, (dtrajs, dummy)), shape=(self.numClusters, 1)).toarray())

        clustersToDelete = np.argwhere(counts < clusterCountsThreshold)
        if clustersToDelete.shape[0] > 0:
            print("Removing %d clusters due to a small number of counts (less than %d)" % (clustersToDelete.shape[0], clusterCountsThreshold))
            self.clusterCenters = np.delete(self.clusterCenters, clustersToDelete, axis=0)
            self._writeClusterCenters(self.clusterCenters, self.clusterCentersFile)
            print("Reassigning trajectories")
            self.dtrajs = self.assignNewTrajectories(self.x)

    def _writeClusterCenters(self, clusterCenters, outputFilename):
        np.savetxt(outputFilename, clusterCenters, fmt=b"%.5f")

    def _writeDtrajs(self, filenames, dtrajs, filenameTemplate="%s.disctraj"):
        for filename, dtraj in zip(filenames, dtrajs):
            fname = os.path.split(filename)[-1][:-4]
            dtrajfname = filenameTemplate % (fname)
            np.savetxt(dtrajfname, dtraj, fmt=b"%d")


# Standalone functions

def loadTrajFiles(trajectoryFolder, trajectory_basename):
    trajectoryBasename = os.path.join(trajectoryFolder, trajectory_basename)

    # load traj
    files = glob.glob(trajectoryBasename)
    x = len(files)*[0]
    for i, f in enumerate(files):
        currentX = np.loadtxt(f, ndmin=2)[:, 1:]
        x[i] = currentX
    if not x:
        raise ValueError("Didn't find any trajectory files in the specified path!!!")
    return x, files


def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
