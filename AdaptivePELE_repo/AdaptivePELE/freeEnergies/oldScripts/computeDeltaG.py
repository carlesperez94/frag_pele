import numpy as np
import helper
import glob
import sys
import argparse
from scipy.ndimage import filters
from pyemma.coordinates.clustering import AssignCenters
import itertools

"""
    Script that computes the absolute binding free energy of an MSM.
    It needs:
        1) Trajectory wildcard in order to compute cluster volume
        2) A discretized/clusterCenters.dat file with the cluster center
        3) A MSM_object.pkl obtained with pyemma in order to obtain the stationary distribution
        4) For the moment, it needs of a reweightingT, in order to do a histogram reweighting, but does not seem to work that well
"""

def assignNewTrajecories(trajs, clusterCenters):
    assign = AssignCenters(clusterCenters)
    dTrajs = assign.assign(trajs)
    return dTrajs

def expandTrajs(trajList):
    d = 0.5
    combinations = np.array([[0,1,0], [1,0,0], [0,0,1], [0,-1,0], [-1,0,0], [0,0,-1]])
    return list(trajList+combinations*d)

def parseArgs():
    parser = argparse.ArgumentParser(description="Script that computes delta G")
    parser.add_argument('trajectory', type=str)
    parser.add_argument('-r', '--reweight', type=float, default=1000)
    args = parser.parse_args()
    return args.trajectory, args.reweight

def writePDB(pmf_xyzg, title="clusters.pdb"):
    templateLine = "HETATM%s  H1  CLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for i, line in enumerate(pmf_xyzg):
        number = str(i).rjust(5)
        x = ("%.3f"%line[0]).rjust(8)
        y = ("%.3f"%line[1]).rjust(8)
        z = ("%.3f"%line[2]).rjust(8)
        g = ("%.3f"%line[3]).rjust(8)

        content += templateLine%(number, x, y, z, g)
    f = open(title, 'w')
    f.write(content)
    f.close()

def calcQ(T, gpmf):
    Q = np.sum(np.exp(-gpmf/(kb*T)))
    return Q

def reweightProbabilities(T, Torig, origProb):
    """
        Torig: temperature at which simulations were run
    """
    kb = 0.0019872041
    gpmf = -kb*Torig*np.log(origProb)
    gpmf[gpmf == -np.inf] = np.inf
    gpmf -= gpmf.min() #It does not make any difference here

    print gpmf[gpmf == np.inf]
    print gpmf[gpmf == -np.inf]
    a = (T-Torig) / (kb * T * Torig)
    correction1 = np.exp(a * gpmf)
    beta = 1/kb/T
    betaOrig = 1/kb/Torig
    Qorig = np.sum(np.exp(-betaOrig*gpmf))
    Qnew = np.sum(np.exp(-beta*gpmf))
    correction2 = Qorig / Qnew

    return correction1 * correction2 * origProb


def main(trajWildcard, reweightingT=1000):
    #clusteringObject = helper.loadMSM('clustering_object.pkl')
    #allClusters = clusteringObject.clustercenters
    allClusters = np.loadtxt("discretized/clusterCenters.dat")

    MSMObject = helper.loadMSM('MSM_object.pkl')
    pi = MSMObject.stationary_distribution
    np.savetxt("stationaryDist_small.dat", pi)


    r = allClusters[MSMObject.connected_sets[0]]


    #filename = "output.txt"
    #data = np.loadtxt(filename)

    #r = data[:,0:3]
    #pi = data[:,3]

    d = 0.75

    originalFilenames = glob.glob(trajWildcard)

    originalCoordinates = []
    for i, originalFilename in enumerate(originalFilenames):
        trajOriginalCoordinates = list(np.loadtxt(originalFilename, usecols=(1,2,3)))
        if np.random.random() < 0.0:
            # Add artificial neighbours to improve volume estimation, set
            # randomly since its very slow
            sys.stderr.write("Introducing artificial neighbours\n")
            newCoords = map(expandTrajs, trajOriginalCoordinates)
            trajOriginalCoordinates.extend(list(itertools.chain.from_iterable(newCoords)))
        originalCoordinates.append(np.array(trajOriginalCoordinates))

    maxval = 3*[-np.inf]
    minval = 3*[np.inf]
    for coord in originalCoordinates:
        cmaxval = coord.max(axis=0)
        cminval = coord.min(axis=0)
        maxval = np.maximum(cmaxval, maxval)
        minval = np.minimum(cminval, minval)

    print "Maximum bounds", maxval, "Minimum bounds", minval

    #Rounded floor and ceiling in intervals of "d" (e.g., floor of 1.73 with d = 0.5, will be 1.5 instead of 1.0, in order to optimize box creation.
    #An extra box is included in the ceiling, so that all the points are contained in the range given by arange
    bins = np.array([np.arange(np.floor(minval[i]) + d*int((minval[i] - np.floor(minval[i]))/d),
                        np.ceil(maxval[i]) + d*(int((maxval[i] - np.ceil(maxval[i]))/d) + 1),
            d) for i in range(3)])

    """
    bins = [np.arange(np.floor(minval[i]),
                        np.ceil(maxval[i]),
            d) for i in range(3)]
    """

    numberOfClusters = r.shape[0]
    histogram = np.array([])
    histogramFreq = np.array([])
    histograms = []

    microstateVolume = np.zeros(numberOfClusters)

    print "Number of clusters", numberOfClusters


    #dtrajs = clusteringObject.assign(originalCoordinates)
    clusterCenters = r
    dtrajs = assignNewTrajecories(originalCoordinates, clusterCenters)
    for i in range(numberOfClusters):
        allCoords = []
        for j,(trajOriginalCoordinates,dtraj) in enumerate(zip(originalCoordinates, dtrajs)):
            assert dtraj.shape[0] == trajOriginalCoordinates.shape[0]
            belongingFrames = np.argwhere(dtraj==i)
            trajCoords = trajOriginalCoordinates[belongingFrames, :]
            trajCoords = trajCoords.flatten().tolist()

            """
            if allCoords.size == 0:
                allCoords = np.copy(trajCoords)
            else:
                allCoords = np.vstack((allCoords, trajCoords))
            """
            allCoords.extend(trajCoords)

        allCoords = np.reshape(allCoords, (-1,3))

        current_hist, edges = np.histogramdd(allCoords, bins=bins)
        histograms.append(current_hist)

        #filtered_hist = filters.gaussian_filter(current_hist, sigma=1)
        #if current_hist.sum() == 0:
        #    filtered_hist = np.zeros(current_hist.shape)
        #else:
        #    filtered_hist = current_hist/current_hist.sum()

        #nonZeroIndices = np.argwhere(filtered_hist > 0)
        #microstateVolume[i] = len(nonZeroIndices) * d**3

        if histogram.size == 0:
            #histogramFreq = pi[i]*filtered_hist
            histogram = np.copy(current_hist)
        else:
            #histogramFreq += pi[i]*filtered_hist
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
                signs = np.sign(histogramTotal[upBound:lowBound, leftBound:rightBound, topBound:botBound])
                histogramCluster[upBound:lowBound, leftBound:rightBound, topBound:botBound] += (1-signsCluster)*d/8  # + signsCluster*d/2
                histogramTotal[upBound:lowBound, leftBound:rightBound, topBound:botBound] += (1-signsCluster)  # + signs * d/2
            histogramTotal = histogramTotal[histogramCluster > 0]
        else:
            histogramTotal = histogram[histogramCluster > 0]
        histogramCluster = histogramCluster[histogramCluster > 0]
        microstateVolume[i] = (histogramCluster/histogramTotal).sum() * d**3


    np.savetxt("volumeOfClusters.dat", microstateVolume)

    # microstateVolume = np.loadtxt("volumeOfClusters.dat")

    #Torig = 1000
    #Tnew = reweightingT
    #newProb = reweightProbabilities(Tnew, Torig, pi)
    #print Torig, Tnew
    #print "normalization", np.sum(newProb)
    newProb = pi


    kb = 0.0019872041
    T = 300
    beta = 1 / (kb * T)
    gpmf = -kb*T*np.log(newProb/microstateVolume)
    print gpmf[gpmf == -np.inf]
    print gpmf[gpmf == np.inf]
    gpmf[gpmf == -np.inf] = np.inf #to avoid contribution later
    gpmf -= gpmf.min()

    deltaW = -gpmf[gpmf != np.inf].max()
    print "bound    Delta G     Delta W     Binding Volume:     Binding Volume contribution"

    upperGpmfValues = np.arange(0,-deltaW,0.5)

    for upperGpmfValue in upperGpmfValues:
        bindingVolume = 0
        for g, volume in zip(gpmf, microstateVolume):
            if g <= upperGpmfValue:
                bindingVolume += np.exp(-beta * g) * volume
        deltaG = deltaW - kb*T*np.log(bindingVolume/1661)
        string = "%.1f\t%.3f\t%.3f\t%.3f\t%.3f" % (upperGpmfValue, deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))
        print string


    pmf_xyzg = np.hstack((r, np.expand_dims(gpmf,axis=1)))
    np.savetxt("pmf_xyzg.dat", pmf_xyzg)

    writePDB(pmf_xyzg)
    return string


    #sys.exit()
    """
    """

    gpmf = -kb*T*np.log(histogram)
    gpmf -= gpmf.min()

    #inExplorationRange = np.argwhere(gpmf != np.inf)
    deltaW = -gpmf[gpmf != np.inf].max()

    for indices in np.argwhere(gpmf != np.inf):
        i = indices[0]
        j = indices[1]
        k = indices[2]
        #print bins[0][i], bins[1][j], bins[2][k], gpmf[i,j,k]

    upperGpmfValues = np.arange(0,-deltaW,0.5)
    bindingVolumes = []
    deltaGs = []
    print "Delta G     Delta W     Binding Volume:     Binding Volume contribution"
    for upperGpmfValue in upperGpmfValues:
        bindingVolume = 0
        for i,g in np.ndenumerate(gpmf[gpmf <= upperGpmfValue]):
            bindingVolume += np.exp(-beta*g)
        bindingVolume *= d**3

        deltaG = deltaW - kb*T*np.log(bindingVolume/1661)

        deltaGs.append(deltaG)
        bindingVolumes.append(bindingVolume)

        #np.savetxt("pmf_xyzg.dat", np.hstack((r, np.expand_dims(gpmf,axis=1))))

        string = "%.3f\t%.3f\t%.3f\t%.3f" % (deltaG, deltaW, bindingVolume, -kb*T*np.log(bindingVolume/1661))
        print string
    return string

    import matplotlib.pyplot as plt
    plt.figure(1)
    plt.plot(upperGpmfValues, deltaGs)
    plt.figure(2)
    plt.plot(upperGpmfValues, bindingVolumes)
    plt.show()

if __name__ == "__main__":
    trajWildcard, reweight = parseArgs()
    main(trajWildcard, reweight)
