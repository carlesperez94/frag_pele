import os
import numpy as np
import cPickle
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def writeClusterCenters(cl, outputFilename):
    np.savetxt(outputFilename, cl.clustercenters)

def writeDtrajs(filenames, dtrajs, filenameTemplate="%s.disctraj"):
    for filename, dtraj in zip(filenames, dtrajs):
        fname = os.path.split(filename)[-1][:-4]
        dtrajfname = filenameTemplate%(fname)
        np.savetxt(dtrajfname, dtraj, fmt="%d")

def makeFolder(outputDir):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

def saveMSM(MSM_object):
    """Save the MSM object to avoid having to run again
    the more computationally expensive part"""
    with open("MSM_object.pkl","w") as MSMfile:
        cPickle.dump(MSM_object, MSMfile, -1)

def saveClustering(clustering_object, clusteringFile):
    """Save the clustering object to avoid having to run again
    the more computationally expensive part"""
    with open(clusteringFile,"w") as clustering_file:
        cPickle.dump(clustering_object, clustering_file, -1)

def loadMSM(MSMFile):
    with open(MSMFile) as MSMfile:
        MSM_object = cPickle.load(MSMfile)
    return MSM_object

def loadClustering(clusteringFile):
    with open(clusteringFile,"r") as clustering_file:
        cl = cPickle.load(clustering_file)
    return cl

def plot3DScatter(matrix):
    fig = plt.figure()
    ax = Axes3D(fig)
    ccx = matrix[:, 0]
    ccy = matrix[:, 1]
    ccz = matrix[:, 2]
    ax.scatter(ccx, ccy, zs=ccz)
    # scatter1 = ax.scatter(ccx, ccy, zs=ccz, c=metrics)
    # fig.colorbar(scatter1)
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    return fig
