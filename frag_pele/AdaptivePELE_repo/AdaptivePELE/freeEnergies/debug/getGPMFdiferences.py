from AdaptivePELE.utilities import utilities
from AdaptivePELE.pyemma_scripts import runMarkovChainModel as run
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def getGPMF(dist, volume, node=None):
    kb = 0.0019872041
    T = 300
    gpmf = -kb*T*np.log(dist/volume)
    gpmf -= gpmf.min()
    if node is not None:
        np.savetxt("GPMF/gmpf_%d.dat" % node, gpmf)
    return gpmf.max()

filename = "differences_400_gpmf.dat"
if os.path.exists(filename):
    differences = np.loadtxt(filename)
    plt.plot(differences)
    plt.show()
    sys.exit()

m1 = utilities.readClusteringObject("MSM_object_0.pkl")
# m2 = utilities.readClusteringObject("/home/jgilaber/3PTB_free_energies/3ptb_PELE_short_steps/run2/400/MSM_0/MSM_object_0.pkl")

vol1 = np.loadtxt("volumeOfClusters_0.dat")
vol2 = np.loadtxt("/home/jgilaber/3PTB_free_energies/3ptb_PELE_short_steps/run2/400/MSM_0/volumeOfClusters_0.dat")

if len(m1.active_set) != 100:
    c1 = m1.count_matrix_full+1/float(100)
    # trans = run.buildRevTransitionMatrix(c1)
    # eiv, eic = run.getSortedEigen(trans)
    # pi = run.getStationaryDistr(eic[:,0])
else:
    # pi = np.loadtxt("stationaryDistribution.dat")
    c1 = m1.count_matrix_full

if len(m2.active_set) != 100:
    c2 = m2.count_matrix_full+1/float(100)
    trans2 = run.buildRevTransitionMatrix(c2)
    eiv, eic2 = run.getSortedEigen(trans2)
    pi2 = run.getStationaryDistr(eic2[:,0])
else:
    pi2 = np.loadtxt("/home/jgilaber/3PTB_free_energies/3ptb_PELE_short_steps/run2/400/MSM_0/stationaryDistribution.dat")
    c2 = m2.count_matrix_full

# gmpf1 = getGPMF(pi, vol1)
gmpf2 = getGPMF(pi2, vol2)
differences = []

# counts = c1.copy()
# for state in [8, 24, 49, 82, 94]:
#     counts[state,:] = c2[state,:]
#     counts[:,state] = c2[:,state]
#
# trans = run.buildRevTransitionMatrix(counts)
# eiv, eic = run.getSortedEigen(trans)
# dist = run.getStationaryDistr(eic[:,0])
# gmpf = getGPMF(dist, vol1)
# differences.append(gmpf-gmpf2)
# print differences
# acum = False
# if acum:
#     counts = c1.copy()
# for state in xrange(c1.shape[0]):
#     if not acum:
#         counts = c1.copy()
#     print state
#     counts[state,:] = c2[state,:]
#     counts[:,state] = c2[:,state]
#     trans = run.buildRevTransitionMatrix(counts)
#     eiv, eic = run.getSortedEigen(trans)
#     dist = run.getStationaryDistr(eic[:,0])
#     gmpf = getGPMF(dist, vol1, state)
#     differences.append(gmpf-gmpf2)
# np.savetxt(filename, differences)
# plt.plot(differences)
# plt.show()
