import pyemma.coordinates as coor
import glob
import os
import numpy as np


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


def loadCoordinates(path, trajectories_basename='*traj*.pdb', topfile=''):
    """ Load the coordinates from the simulation into an object that can be
    used with the rest of PyEMMA tools.
    Returns an array with the trajectories"""

    feat = coor.featurizer(topfile)

    path_to_file = os.path.join(path, trajectories_basename)
    files = glob.glob(path_to_file)
    return coor.load(files, feat)

# More preprocessing staff?


def clusterTrajectories(trajectories, numClusters, stride=1):
    """ Cluster the trajectories into numClusters clusters using kmeans
    algorithm.
    Returns a KmeansClusteringObject
    """
    return coor.cluster_kmeans(data=trajectories, k=numClusters, max_iter=20, stride=stride)

def clusterRegularSpace(trajectories, dmin, stride=1):
    """
        Cluster the trajectories using Regular Space clustering, which is a
        modified version of Hartigan's leader algorithm
    """
    return coor.cluster_regspace(data=trajectories, dmin = dmin, stride=stride)


