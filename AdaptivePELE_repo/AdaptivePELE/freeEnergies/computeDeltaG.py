"""
    Script that computes the absolute binding free energy of an MSM.
    It needs:
        1) Trajectory wildcard in order to compute cluster volume
        2) A discretized/clusterCenters.dat file with the cluster center
        3) A MSM_object.pkl obtained with pyemma in order to obtain the stationary distribution
        4) For the moment, it needs of a reweightingT, in order to do a histogram reweighting, but does not seem to work that well
"""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import numpy as np
import glob
import sys
import argparse
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from pyemma.coordinates.clustering import AssignCenters
from AdaptivePELE.freeEnergies import runMarkovChainModel as run
from AdaptivePELE.freeEnergies import utils
import itertools


def assignNewTrajectories(trajs, clusterCenters):
    assign = AssignCenters(clusterCenters)
    dTrajs = assign.assign(trajs)
    return dTrajs


def expandTrajs(trajList):
    d = 0.5
    combinations = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0],
                             [-1, 0, 0], [0, 0, -1]])
    return list(trajList+combinations*d)


def parseArgs():
    parser = argparse.ArgumentParser(description="Script that computes delta G")
    parser.add_argument('trajectory', type=str)
    parser.add_argument('-r', '--reweight', type=float, default=1000)
    args = parser.parse_args()
    return args.trajectory, args.reweight


def writePDB(pmf_xyzg, title="clusters.pdb"):
    # templateLine = "HETATM%s  H%sCLT L 502    %s%s%s  0.75%s           H\n"
    templateLine = "HETATM%s H%sCLT L 502    %s%s%s  0.75%s                H\n"

    content = ""
    for i, line in enumerate(pmf_xyzg):
        number = str(i).rjust(5)
        number3 = str(i).ljust(3)
        x = ("%.3f" % line[0]).rjust(8)
        y = ("%.3f" % line[1]).rjust(8)
        z = ("%.3f" % line[2]).rjust(8)
        g = ("%.3f" % line[-1]).rjust(8)

        content += templateLine % (number, number3, x, y, z, g)
    f = open(title, 'w')
    f.write(content)
    f.close()


def calcQ(T, gpmf):
    kb = 0.0019872041
    Q = np.sum(np.exp(-gpmf/(kb*T)))
    return Q


def reweightProbabilities(T, Torig, origProb):
    """
        Torig: temperature at which simulations were run
    """
    kb = 0.0019872041
    gpmf = -kb*Torig*np.log(origProb)
    gpmf[gpmf == -np.inf] = np.inf
    gpmf -= gpmf.min()  # It does not make any difference here

    print(gpmf[gpmf == np.inf])
    print(gpmf[gpmf == -np.inf])
    a = (T-Torig) / (kb * T * Torig)
    correction1 = np.exp(a * gpmf)
    beta = 1/kb/T
    betaOrig = 1/kb/Torig
    Qorig = np.sum(np.exp(-betaOrig*gpmf))
    Qnew = np.sum(np.exp(-beta*gpmf))
    correction2 = Qorig / Qnew

    return correction1 * correction2 * origProb


def loadMSM(MSMFile):
    with open(MSMFile, "rb") as MSMfile:
        MSM_object = cPickle.load(MSMfile)
    return MSM_object


def reestimate_transition_matrix(count_matrix):
    count_matrix += 1/float(count_matrix.shape[0])
    return utils.buildRevTransitionMatrix(count_matrix)


def ensure_connectivity(MSMObject, allClusters):
    if len(allClusters) == MSMObject.stationary_distribution.size:
        pi = MSMObject.stationary_distribution
        clusters = allClusters[MSMObject.connected_sets[0]]
    else:
        ######################
        # Reconstruct stationary distribution with pseudocounts to ensure
        # connectivity
        print("Adding pseudocounts to enforce connectivity")
        counts = MSMObject.count_matrix_full
        trans = reestimate_transition_matrix(counts)
        _, eic = run.getSortedEigen(trans)
        pi = run.getStationaryDistr(eic[:, 0])
        clusters = allClusters
    return pi, clusters


def gather_coordinates(originalFilenames):
    originalCoordinates = []
    for originalFilename in originalFilenames:
        trajOriginalCoordinates = list(np.loadtxt(originalFilename, ndmin=2)[:, 1:])
        if np.random.random() < 0.0:
            # Add artificial points nearby to improve volume estimation, set
            # randomly since its very slow
            sys.stderr.write("Introducing artificial neighbours\n")
            newCoords = map(expandTrajs, trajOriginalCoordinates)
            trajOriginalCoordinates.extend(list(itertools.chain.from_iterable(newCoords)))
        originalCoordinates.append(np.array(trajOriginalCoordinates))
    return originalCoordinates


def create_box(clusters, originalCoordinates, d):
    """
        Create the box discretization for the estimation of volumes
    """
    dimensions = clusters.shape[1]
    maxval = dimensions*[-np.inf]
    minval = dimensions*[np.inf]
    for coord in originalCoordinates:
        cmaxval = coord.max(axis=0)
        cminval = coord.min(axis=0)
        maxval = np.maximum(cmaxval, maxval)
        minval = np.minimum(cminval, minval)

    print("Maximum bounds", maxval, "Minimum bounds", minval)

    # Rounded floor and ceiling in intervals of "d" (e.g., floor of 1.73 with d = 0.5, will be 1.5 instead of 1.0, in order to optimize box creation.
    # An extra box is included in the ceiling, so that all the points are contained in the range given by arange
    bins = np.array([np.arange(np.floor(minval[i]) + d*int((minval[i] - np.floor(minval[i]))/d),
                               np.ceil(maxval[i]) + d*(int((maxval[i] - np.ceil(maxval[i]))/d) + 1),
                               d) for i in range(3)])
    return bins


def calculate_microstate_volumes(clusters, originalCoordinates, bins, d):
    """
        Estimate the clusters volumes using a cubic discretization of volumes
    """
    numberOfClusters = clusters.shape[0]
    print("Number of clusters", numberOfClusters)
    histogram = np.array([])
    histograms = []
    microstateVolume = np.zeros(numberOfClusters)

    # dtrajs = clusteringObject.assign(originalCoordinates)
    dtrajs = assignNewTrajectories(originalCoordinates, clusters)
    for i in range(numberOfClusters):
        allCoords = []
        for trajOriginalCoordinates, dtraj in zip(originalCoordinates, dtrajs):
            assert dtraj.shape[0] == trajOriginalCoordinates.shape[0]
            belongingFrames = np.argwhere(dtraj == i)
            trajCoords = trajOriginalCoordinates[belongingFrames, :3]
            trajCoords = trajCoords.flatten().tolist()

            allCoords.extend(trajCoords)

        allCoords = np.reshape(allCoords, (-1, 3))

        current_hist, _ = np.histogramdd(allCoords, bins=bins)
        histograms.append(current_hist)

        if histogram.size == 0:
            histogram = np.copy(current_hist)
        else:
            histogram += current_hist

    nRows, nCols, nDepth = histogram.shape
    pseudo = False
    for i in range(numberOfClusters):
        histogramCluster = histograms[i]
        if pseudo:
            # Add "pseudocounts" to try to fill the holes that lead to volume
            # underestimation compared to Matlab script for free energies
            histogramTotal = histogram.copy()
            for x, y, z in zip(*np.where(histogramCluster)):
                upBound = max(x-1, 0)
                lowBound = min(x+2, nRows)
                leftBound = max(0, y-1)
                rightBound = min(y+2, nCols)
                topBound = max(z-1, 0)
                botBound = min(z+2, nDepth)
                signsCluster = np.sign(histogramCluster[upBound:lowBound, leftBound:rightBound, topBound:botBound])
                # signs = np.sign(histogramTotal[upBound:lowBound, leftBound:rightBound, topBound:botBound])
                histogramCluster[upBound:lowBound, leftBound:rightBound, topBound:botBound] += (1-signsCluster)*d/8  # + signsCluster*d/2
                histogramTotal[upBound:lowBound, leftBound:rightBound, topBound:botBound] += (1-signsCluster)  # + signs * d/2
            histogramTotal = histogramTotal[histogramCluster > 0]
        else:
            histogramTotal = histogram[histogramCluster > 0]
        histogramCluster = histogramCluster[histogramCluster > 0]
        microstateVolume[i] = (histogramCluster/histogramTotal).sum() * d**3
    return microstateVolume


def calculate_microstate_volumes_new(clusters, originalCoordinates, bins, d):
    """
        Estimate the clusters volumes using a cubic discretization of volumes
    """
    numberOfClusters = clusters.shape[0]
    print("Number of clusters", numberOfClusters)

    allCoords = []
    # The coordinates array is built through lists in order to be able to
    # process trajectories of different lenght
    for coord in originalCoordinates:
        allCoords.extend(coord[:, :3].tolist())
    allCoords = np.array(allCoords)

    histogram, _ = np.histogramdd(allCoords, bins=bins)

    centers_trajs = []
    for i, matrix in enumerate(histogram):
        centers_x = []
        x_val = bins[0][i]
        for j, row in enumerate(matrix):
            y_val = bins[1][j]
            indices, = np.nonzero(row)
            if indices.size > 0:
                init = indices[0]
                end = indices[-1]
                centers_x.extend([[x_val, y_val, bins[2][zi]] for zi in range(init, end+1)])
        if len(centers_x):
            centers_trajs.append(np.array(centers_x))

    dtrajs = assignNewTrajectories(centers_trajs, clusters[:, :3])
    microstateVolume = np.zeros(numberOfClusters)
    for traj in dtrajs:
        counts, _ = np.histogram(traj, bins=np.arange(-0.5, numberOfClusters))
        microstateVolume += counts
    microstateVolume *= d**3
    # provisional return for debugging and testing script
    # return microstateVolume, centers_trajs, dtrajs
    return microstateVolume


def calculate_pmf(microstateVolume, pi):
    """
        Compute a potential of mean force given a stationary distribution
        (probabilities) and cluster volumes
    """
    kb = 0.0019872041
    T = 300
    beta = 1 / (kb * T)
    newDist = pi/microstateVolume
    newDist /= newDist[newDist != np.inf].sum()
    gpmf = -kb*T*np.log(newDist)
    print(gpmf[gpmf == -np.inf])
    print(gpmf[gpmf == np.inf])
    gpmf[gpmf == -np.inf] = np.inf  # to avoid contribution later
    gpmf -= gpmf.min()

    deltaW = -gpmf[gpmf != np.inf].max()
    print("bound    Delta G     Delta W     Binding Volume:     Binding Volume contribution")

    upperGpmfValues = np.arange(0, -deltaW, 0.5)

    # Initialize string variable in case loop is not accessed
    string = ""

    for upperGpmfValue in upperGpmfValues:
        bindingVolume = 0
        for g, volume in zip(gpmf, microstateVolume):
            if g <= upperGpmfValue:
                bindingVolume += np.exp(-beta * g) * volume
        deltaG = deltaW - kb*T*np.log(bindingVolume/1661)
        string = "%.1f\t%.3f\t%.3f\t%.3f\t%.3f" % (upperGpmfValue, deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))
        print(string)
    return gpmf, string


def main(trajWildcard, reweightingT=1000):
    allClusters = np.loadtxt("discretized/clusterCenters.dat")
    MSMObject = loadMSM('MSM_object.pkl')

    pi, clusters = ensure_connectivity(MSMObject, allClusters)
    # radius of the cube for volume determination
    d = 0.75

    originalFilenames = glob.glob(trajWildcard)
    # originalFilenames = glob.glob("rawData/"+trajWildcard)
    originalCoordinates = gather_coordinates(originalFilenames)

    bins = create_box(clusters, originalCoordinates, d)
    method = "new"
    if method == "new":
        print("Using new volume estimation with radius %.2f" % d)
        microstateVolume = calculate_microstate_volumes_new(clusters, originalCoordinates, bins, d)
    else:
        print("Using old volume estimation with radius %.2f" % d)
        microstateVolume = calculate_microstate_volumes(clusters, originalCoordinates, bins, d)
    np.savetxt("volumeOfClusters.dat", microstateVolume)

    gpmf, string = calculate_pmf(microstateVolume, pi)

    pmf_xyzg = np.hstack((clusters, np.expand_dims(gpmf, axis=1)))
    np.savetxt("pmf_xyzg.dat", pmf_xyzg)

    writePDB(pmf_xyzg)
    return string

if __name__ == "__main__":
    trajWildcard_name, reweight = parseArgs()
    main(trajWildcard_name, reweight)
