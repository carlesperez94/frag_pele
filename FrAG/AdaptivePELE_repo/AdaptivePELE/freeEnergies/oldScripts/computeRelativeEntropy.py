import json
import os
import sys
import helper
import runMarkovChainModel as markov
import trajectories
from pyemma.coordinates.clustering import AssignCenters
from pyemma import msm as pyemmaMSM
import msm
import numpy as np
import matplotlib.pyplot as plt

from AdaptivePELE.pyemma_scripts import revTransitionMatrix #cython implementation

def readJSON(controlFile):
    with open(controlFile, "r") as f:
        paramsJSON = json.load(f)
    return paramsJSON

def readParams(control_file):
    params = readJSON(control_file)
    disctrajFolder = params["disctrajFolder"]
    trajectoryFolder = params["goldenTrajFolder"]
    trajectoryFolder2 = params["trajFolder"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    lagtime = params.get("lagtime", 0)
    sampleSize = params.get("sampleSize", None)
    numRuns = params.get("numRuns", 1)
    dTraj = params.get("dtraj", None)
    maxNtraj = params.get("maxntraj", None)
    minNtraj = params.get("minntraj", None)
    dlength = params.get("dlength", None)
    maxlength = params.get("maxlength", None)
    return disctrajFolder, trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtime, sampleSize, numRuns, dTraj, maxNtraj, minNtraj, dlength, maxlength

def getStationaryDistributionAndTransitionMatrix(dTrajs, nclusters, lagtime):
    counts = markov.estimateCountMatrix(dTrajs, nclusters, lagtime)
    counts += 1.0/counts.shape[0]

    realT = revTransitionMatrix.buildRevTransitionMatrix_fast(counts, iterations=1000)
    T = markov.buildTransitionMatrix(counts) #fastest
    #other (slower) options
    #transition = markov.buildRevTransitionMatrix_fast(counts, iterations=5)
    #transition = markov.buildRevTransitionMatrix(counts)

    eigenvals, eigenvec = markov.getSortedEigen(realT)
    pi = markov.getStationaryDistr(eigenvec[:,0])

    return pi, T


def makeRandomSampleOfNtrajs(X, ntrajs=None, length=None):
    if ntrajs:
        indices = np.array(np.random.choice(range(len(X)), ntrajs))
    else:
        indices = range(len(X))

    try:
        Xsample = map(lambda x: X[x][:length,:], indices)
    except:
        for index in indices:
            try:
                X[index][:length,:]
            except:
                import sys
                sys.exit("There is a problem with the trajectory!")

    return Xsample

def makeRandomSampleOfdtrajs(dtrajs, ntrajs=None, length=None):
    if ntrajs:
        indices = np.array(np.random.choice(range(len(dtrajs)), ntrajs))
    else:
        indices = range(len(dtrajs))

    try:
        Xsample = map(lambda x: np.array(dtrajs[x][:length]), indices)
    except:
        for index in indices:
            try:
                dtrajs[index][:length]
            except:
                import sys
                sys.exit("There is a problem with the trajectory!")

    return Xsample

def assignTrajectories(goldenMSMClusterCenters, X):
    assign = AssignCenters(goldenMSMClusterCenters)
    dTrajs = assign.assign(X)
    return dTrajs

def estimateTWithDiscTrajs(dTrajs, nclusters, lagtime):
    #own estimation
    counts = markov.estimateCountMatrix(dTrajs, nclusters, lagtime)
    counts += 1.0/counts.shape[0]
    #transition = revTransitionMatrix.buildRevTransitionMatrix_fast(counts, iterations=100)
    transition = markov.buildTransitionMatrix(counts)
    #other (slower) options
    #transition = markov.buildRevTransitionMatrix_fast(counts, iterations=5)
    #transition = markov.buildRevTransitionMatrix(counts)
    return transition

def plotIsocostLines(extent, allTrajLengths, numberOfTrajs, steps=10):
    minCost = allTrajLengths[0]*numberOfTrajs[0]
    maxCost = allTrajLengths[-1]*numberOfTrajs[-1]
    d = (maxCost - minCost) / steps
    for cost in np.arange(minCost, maxCost, d):
        x = np.arange(extent[0], extent[1], 1)
        y = cost / x
        plt.plot(x,y, color="black")


def readDTrajOrAssign(filename, goldenMSMClusterCenters, X):
    if os.path.exists(filename):
        print "Reading current disc trajs from file: ", filename
        discTrajs = np.load(filename).tolist()
    else:
        discTrajs = assignTrajectories(goldenMSMClusterCenters, X)
        np.save(filename, np.array(discTrajs))
    return discTrajs

def buildDTraj(trajectoryFolder, trajectoryBasename, disctrajFolder, filename="golden_dtraj.npy"):
    """
        Builds discretized trajectories.
        If it finds g_dtraj.npy, it builds them reading the file.
        Otherwise, it reads X, reads cluster centers, and assigns them using the Voronoi of EMMA.
    """
    clusterCenters = os.path.join(disctrajFolder, "discretized/clusterCenters.dat")
    try:
        goldenMSMClusterCenters = np.loadtxt(clusterCenters)
    except:
        try:
            clusterCenters = os.path.join(disctrajFolder, "clusterCenters.dat")
            goldenMSMClusterCenters = np.loadtxt(clusterCenters)
        except:
            sys.exit("Didn't find cluster centers")

    if os.path.exists(filename):
        print "Reading current disc trajs from file: ", filename
        dTrajs = np.load(filename).tolist()
    else: #TODO: refactor readDTrajOrAssign
        goldX, unused = trajectories.loadCOMFiles(trajectoryFolder, trajectoryBasename)
        dTrajs = readDTrajOrAssign(filename, goldenMSMClusterCenters, goldX)
    return dTrajs, goldenMSMClusterCenters

def main(controlFile):
    """
        Takes cluster centers file, builds dtrajs, and computes relative entropy
    """
    try:
        disctrajFolder, trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtime, sampleSize, numRuns, dtrajs, maxNtraj, minNtraj, dlength, maxlength = readParams(controlFile)
    except IOError:
        print "Removing discretized trajectory files"
        if controlFile == "rm":
            os.remove("golden_dtraj.npy")
            os.remove("c_dtraj.npy")

    #Deprecated, to avoid dependencies on pyemma
    #refTransition, refStationaryDist, lagtime = getTransitionMatrix(trajectoryFolder, trajectoryBasename, numClusters, lagtimes, itsOutput, numberOfITS, itsErrors, lagtime, stride)

    dTrajs, goldenMSMClusterCenters = buildDTraj(trajectoryFolder, trajectoryBasename, disctrajFolder, "golden_dtraj.npy")
    pi, T = getStationaryDistributionAndTransitionMatrix(dTrajs, goldenMSMClusterCenters.shape[0], lagtime)

    #np.random.seed(250793)

    seq = True
    entropies = []

    if seq:
        try:
            X,unused = trajectories.loadCOMFiles(trajectoryFolder2, trajectoryBasename)
            discTrajs = readDTrajOrAssign("c_dtraj.npy", goldenMSMClusterCenters, X)

            if dtrajs is None:
                dtrajs = 100
            if maxNtraj is None:
                maxNtraj = len(discTrajs)
            numberOfTrajs = range(minNtraj, maxNtraj, dtrajs)

            #dTrajs = 100
            # numberOfTrajs = range(50, sampleSize, 50)

            #only trying different traj lengths if sampleSize is defined in control file
            shortestTrajSize = min([len(i) for i in X])
            lowerLimit = 2*lagtime

            if maxlength is None:
                maxlength = shortestTrajSize
            if dlength is None:
                dlength = 100
            allTrajLengths = range(lowerLimit, maxlength, dlength)
        except TypeError:
            numberOfTrajs = [None]
            allTrajLengths = [None]
    else:
        # epochFolders = [int(i) for i in os.listdir(trajectoryFolder2) if i.isdigit()]
        dTrajs = 1
        epochFolders = range(3)
        epochFolders.sort()
        # import pdb as debug
        # debug.set_trace()
        lowerLimit = 200+50
        numberOfTrajs = epochFolders
        upperLimit = 401
        dTrajLengths = 50
        allTrajLengths = range(lowerLimit, upperLimit, dTrajLengths)


    for length in allTrajLengths:
        if length: print "Working with trajectories of length: %d" % length
        lengthEntropies = []
        if not seq:
            X = []
        for ntrajs in numberOfTrajs:

            if not seq:
                currentEpochFolder = os.path.join(trajectoryFolder2, str(ntrajs))
                currentX,unused = trajectories.loadCOMFiles(currentEpochFolder, trajectoryBasename)
                X.extend(currentX)

            if not ntrajs % 25:
                print "Starting loop for sample of %d trajs" % ntrajs
            relativeEntropy = 0
            for j in range(numRuns):
                if seq:
                    #Xsample = makeRandomSampleOfNtrajs(X, ntrajs, length)
                    discTrajsSample = makeRandomSampleOfdtrajs(discTrajs,  ntrajs, length)
                else:
                    discTrajsSample = map(lambda x: x[:length,:], discTrajs)

                transitionMatrix = estimateTWithDiscTrajs(discTrajsSample, goldenMSMClusterCenters.shape[0], lagtime)
                try:
                    s = markov.getRelativeEntropy(pi, T, transitionMatrix)
                    relativeEntropy += s
                except ValueError:
                    j -= 1
            lengthEntropies.append(relativeEntropy/float(numRuns))
            print length, ntrajs, relativeEntropy/float(numRuns)
        print lengthEntropies
        entropies.append(lengthEntropies)

    np.save("matrix_adaptive.npy", entropies)
    #entropies = np.load("matrix_adaptive.npy")
    entropies = np.log10(entropies)

    for i, length in enumerate(allTrajLengths):
        for j, ntrajs  in enumerate(numberOfTrajs):
            print length, ntrajs, entropies[i][j]
        print ""

    if ntrajs and length:
        plt.figure(1)
        if seq:
            extent = [numberOfTrajs[0] - dtrajs/2, numberOfTrajs[-1] + dtrajs/2,
                    allTrajLengths[0] - dlength/2, allTrajLengths[-1] + dlength/2]
            #plot isocost lines
            plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 9)
            plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
        else:
            plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
        import os
        cwd = os.getcwd()
        cwd = cwd.replace("/", "_")
        #plt.save(cwd + ".eps")
        #plt.show()
        #plt.imshow(entropies, interpolation="nearest", extent=[numberOfTrajs[0], numberOfTrajs[-1], 800, 801])
        plt.colorbar()
        plt.savefig(cwd + ".eps")
        plt.figure(2)
        if seq:
            plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 9)
            plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=extent)
        else:
            plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=extent)
        plt.colorbar()
        import os
        cwd = os.getcwd()
        cwd = cwd.replace("/", "_")
        plt.savefig(cwd + "2.eps")
        plt.show()
    elif ntrajs:
        print numberOfTrajs, entropies[0]
        plt.plot(numberOfTrajs, entropies[0])
        plt.show()
    else:
        print entropies

if __name__ == "__main__":
    controlFile = sys.argv[1]
    MSM_object = main(controlFile)
