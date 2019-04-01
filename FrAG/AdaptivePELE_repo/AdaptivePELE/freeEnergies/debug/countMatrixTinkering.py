from AdaptivePELE.utilities import utilities
from AdaptivePELE.pyemma_scripts import runMarkovChainModel as run
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def getDW(dist, volume, node=None):
    kb = 0.0019872041
    T = 300
    gpmf = -kb*T*np.log(dist/volume)
    gpmf -= gpmf.min()
    if node is not None:
        np.savetxt("GPMF/gmpf_%d.dat" % node, gpmf)
    return gpmf.max()

def getPi(counts):
    trans = run.buildRevTransitionMatrix(counts)
    eiv, eic = run.getSortedEigen(trans)
    return run.getStationaryDistr(eic[:,0])

m1 = utilities.readClusteringObject("MSM_object_0.pkl.bckp")

vol1 = np.loadtxt("volumeOfClusters_0.dat")

c1 = m1.count_matrix_full
c1 += 1/100.

c1[11,94] *= 10

pi = getPi(c1)
np.savetxt("stat.dat", pi)
stat = m1.stationary_distribution
dw = getDW(pi, vol1)
print "DW", dw
