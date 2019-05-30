from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import numpy as np
import shutil
import glob
import sys
import matplotlib.pyplot as plt
from AdaptivePELE.freeEnergies import checkDetailedBalance
from AdaptivePELE.freeEnergies import ownBuildMSM
from AdaptivePELE.freeEnergies import computeDeltaG
from AdaptivePELE.freeEnergies import cluster


class Parameters:
    def __init__(self, ntrajs, length, lagtime, nclusters, nruns, useAllTrajInFirstRun, computeDetailedBalance, trajWildcard, folderWithTraj, lagtimes=None, skipFirstSteps=0, clusterCountsThreshold=0, clusteringStride=1):
        # If ntrajs/length = None, all trajs/lengths will be used
        self.trajWildcard = trajWildcard
        self.folderWithTraj = folderWithTraj
        self.ntrajs = ntrajs
        self.length = length
        self.lagtime = lagtime
        self.nclusters = nclusters
        self.nruns = nruns
        self.useAllTrajInFirstRun = useAllTrajInFirstRun
        self.computeDetailedBalance = computeDetailedBalance
        if lagtimes is None:
            self.lagtimes = []
        else:
            self.lagtimes = lagtimes
        self.skipFirstSteps = skipFirstSteps
        self.clusterCountsThreshold = clusterCountsThreshold
        self.clusteringStride = clusteringStride


def __rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def __rmFiles(trajWildcard):
    allfiles = glob.glob(trajWildcard)
    for f in allfiles:
        __rm(f)


def __prepareWorkingControlFile(lagtime, clusters, trajectoryFolder, trajectoryBasename, workingControlFile, lagtimes, clusterCountsThreshold=0, clusteringStride=1):
    """
    #Unused alternative #1, need of a templetized control file
    simulationParameters = simulationrunner.SimulationParameters()
    simulationParameters.templetizedControlFile = controlFile
    sr = simulationrunner.SimulationRunner(simulationParameters)

    controlFileDictionary = {"lagtime": lagtime, "clusters": clusters}
    sr.makeWorkingControlFile(workingControlFile, controlFileDictionary)
    """

    workingFolder = os.path.split(trajectoryFolder)[0]  # note that we assume a workingFolder/origTrajs structure (typically origTrajs=rawData)
    try:
        string = "{\"trajectoryFolder\":\"%s\", \"trajectoryBasename\":\"%s\", \"numClusters\":%d, \"stride\":%d, \"lagtime\":%d, \"itsOutput\":\"its.png\", \"lagtimes\":%s, \"clusterCountsThreshold\":%d}" % (workingFolder, trajectoryBasename, clusters, clusteringStride, lagtime, lagtimes, clusterCountsThreshold)
    except TypeError:
        string = "{\"trajectoryFolder\":\"%s\", \"trajectoryBasename\":\"%s\", \"numClusters\":%d, \"stride\":%d, \"itsOutput\":\"its.png\", \"lagtimes\":%s, \"clusterCountsThreshold\":%d}" % (workingFolder, trajectoryBasename, clusters, clusteringStride, lagtimes, clusterCountsThreshold)
    with open(workingControlFile, 'w') as f:
        f.write(string)


def __constructMSM(workingControlFile):
    ownBuildMSM.main(workingControlFile)


def __computeDG(trajWildcard):
    deltaGLine = computeDeltaG.main(trajWildcard)
    return deltaGLine


def __getDstName(bootstrap, i, trajFile):
    # Equiv lambda
    # getDstName = lambda bootstrap, i, trajFile: "traj_.%d.dat"%i if bootstrap else os.path.split(trajFile)[-1]
    if bootstrap:
        return "traj_.%d.dat" % i
    else:
        return os.path.split(trajFile)[-1]


def copyWorkingTrajectories(fileWildcard, length=None, ntrajs=None, bootstrap=True, skipFirstSteps=0):
    """
        Function that copies trajectories that match "fileWildcard" into the current directory.
        It may copy a subset and a part of them (length)

        Warning! If not using bootstrap, it uses all the trajectories

        :param fileWildcard: Wildcard to match original files
        :type fileWildcard: str
        :param length: Trajectory length to consider, if None (default value), the full trajectory will be considered
        :type length: int
        :param ntrajs: Number of trajs to consider. If None (default value), a number equal to the total  will be considered
        :type ntrajs: int
        :param bootstrap: Bootstrap ntrajs from the original (default is True)
        :type bootstrap: bool
        :param skipFirstSteps: Skip first trajectory steps (default value is 0)
        :type skipFirstSteps: int

        :returns: list -- writenFiles, in order to ease a posterior cleanup

    """
    allFiles = glob.glob(fileWildcard)

    if bootstrap is False:
        trajFiles = allFiles
    else:
        if ntrajs is None:
            ntrajs = len(allFiles)
        trajFiles = np.random.choice(allFiles, ntrajs)

    writenFiles = []
    for i, trajFile in enumerate(trajFiles):
        dst = __getDstName(bootstrap, i, trajFile)
        writenFiles.append(dst)
        traj = np.loadtxt(trajFile)
        if length is None:
            length = -2  # so that later eveything is copied
        try:
            trimmedTraj = traj[skipFirstSteps:length+1, :]
            if len(trimmedTraj) > 0:
                np.savetxt(dst, trimmedTraj, fmt=b"%.4f", delimiter="\t")
        except:
            sys.exit("There is a problem with %s" % trajFile)
    print("Boostraping trajectories", ntrajs, len(trajFiles), len(set(trajFiles)))
    # trajFiles.sort()
    # print(trajFiles)
    return writenFiles


def __cleanupFiles(trajWildcard, cleanupClusterCenters=True):
    __rmFiles("clustering_object.pkl")
    __rmFiles("MSM_object.pkl")
    __rmFiles("discretized/traj_*")
    __rmFiles(trajWildcard)
    if cleanupClusterCenters:
        __rmFiles("discretized/clusterCenter*")


def __setVariablesForFirstIteration(useAllTrajInFirstRun, i, ntrajs):
    if useAllTrajInFirstRun and i == 0:
        print("Not using bootstrap in iteration %d" % i)
        bootstrap = False
        nWorkingTrajs = None  # Not necessary, just to make it explicit that all of them are used
    else:
        print("Using bootstrap in iteration %d" % i)
        bootstrap = True
        nWorkingTrajs = ntrajs

    return bootstrap, nWorkingTrajs


def __copyMSMDataFromRun(i):
    try:  # it may not exist
        shutil.copyfile("its.png", "its_%d.png" % i)
    except IOError:
        pass
    shutil.copyfile("discretized/clusterCenters.dat", "clusterCenters_%d.dat" % i)
    shutil.copyfile("volumeOfClusters.dat", "volumeOfClusters_%d.dat" % i)
    shutil.copyfile("clusters.pdb", "clusters_%d.pdb" % i)
    shutil.copyfile("pmf_xyzg.dat", "pmf_xyzg_%d.dat" % i)
    shutil.copyfile("MSM_object.pkl", "MSM_object_%d.pkl" % i)
    if i == 0:
        try:
            shutil.copyfile("db_frobenius.eps", "db_frobenius_%d.eps" % i)
            shutil.copyfile("db_abs_diff.eps", "db_abs_diff_%d.eps" % i)
            shutil.copyfile("db_flux.eps", "db_flux_%d.eps" % i)
        except IOError:
            pass


def __printList(l, label):
    print(label)
    print("=====")
    for el in l:
        print(el)


def __getMeanAndStdFromList(l, accessFunction=lambda x: x):
    values = [float(accessFunction(element)) for element in l]
    return np.mean(values), np.std(values)


def getRepresentativePDBs(filesWildcard, run):
    files = glob.glob(filesWildcard)
    trajs = [np.loadtxt(f)[:, 1:] for f in files]
    cl = cluster.Cluster(0, "", "")
    cl.clusterCenters = np.loadtxt(cl.clusterCentersFile)
    dtrajs = cl.assignNewTrajectories(trajs)
    numClusters = cl.clusterCenters.shape[0]
    centersInfo = {x: {"structure": None, "minDist": 1e6} for x in range(numClusters)}
    for i, traj in enumerate(trajs):
        traj_name = files[i]
        _, epochNum, trajNum = os.path.splitext(traj_name)[0].split("_", 2)
        for nSnap, snapshot in enumerate(traj):
            clusterInd = dtrajs[i][nSnap]
            dist = np.linalg.norm(cl.clusterCenters[clusterInd]-snapshot)
            if dist < centersInfo[clusterInd]['minDist']:
                centersInfo[clusterInd]["minDist"] = dist
                centersInfo[clusterInd]["structure"] = (epochNum, trajNum, str(nSnap))

    if not os.path.exists("representative_structures"):
        os.makedirs("representative_structures")
    with open("representative_structures/representative_structures_%d.dat" % run, "w") as fw:
        fw.write("Cluster\tEpoch\tTrajectory\tSnapshot\n")
        for clNum in range(numClusters):
            fw.write("%d\t" % clNum+"\t".join(centersInfo[clNum]["structure"])+"\n")


def estimateDG(parameters, cleanupClusterCentersAtStart=False):
    """
        Estimates the absolute binding free energy using the parameters in the Parameters object.

        It copies the trajectory files from "folderWithTraj" into the current folder makes a certain number of iterations

        Documentation needs to be expanded, but the code style aims to help readability
    """

    workingControlFile = "control_MSM.conf"
    origFilesWildcard = os.path.join(parameters.folderWithTraj, parameters.trajWildcard)
    origFilesNonRepeatedWildcard = os.path.join(parameters.folderWithTraj, "extractedCoordinates", parameters.trajWildcard)

    __prepareWorkingControlFile(parameters.lagtime, parameters.nclusters, parameters.folderWithTraj, parameters.trajWildcard, workingControlFile, parameters.lagtimes, parameters.clusterCountsThreshold, parameters.clusteringStride)

    deltaGs = []
    detailedBalance = []
    __cleanupFiles(parameters.trajWildcard, cleanupClusterCentersAtStart)

    for i in range(parameters.nruns):
        bootstrap, nWorkingTrajs = __setVariablesForFirstIteration(parameters.useAllTrajInFirstRun, i, parameters.ntrajs)

        copyWorkingTrajectories(origFilesWildcard, parameters.length, nWorkingTrajs, bootstrap, parameters.skipFirstSteps)

        __constructMSM(workingControlFile)

        deltaG = __computeDG(parameters.trajWildcard)
        deltaGs.append(deltaG)

        getRepresentativePDBs(origFilesNonRepeatedWildcard, i)

        if parameters.computeDetailedBalance:
            avgAsymmetricFlux = checkDetailedBalance.main(folder="discretized", countsThreshold=0, lagtime=parameters.lagtime, printFigs=False)
            detailedBalance.append(avgAsymmetricFlux)

        __copyMSMDataFromRun(i)

        __cleanupFiles(parameters.trajWildcard, True)

        # Close all open matplotlib windows, apparently the its plot opens a lot
        # of windows that are not closed, which consumes a lot of memory (not
        # sure how much exactly)
        plt.close("all")
    # PLOT RESULTS
    # FIX TO WORK WITH NONES
    # print("clusters: %d, ntrajs: %d, trajLength: %d, lagtime: # % d"%(parameters.nclusters, parameters.ntrajs, parameters.length, # parameters.lagtime))
    __printList(deltaGs, "dG")
    meanDG, stdDG = __getMeanAndStdFromList(deltaGs, lambda element: element.split()[1])
    print("dG = %f +- %f" % (meanDG, stdDG))
    __printList(detailedBalance, "Asymmetric fluxes (see D.Lecina PhD thesis for more info)")
    meanDB, stdDB = __getMeanAndStdFromList(detailedBalance)  # DB from detailed balance
    print("Asymmetric flux = %f +- %f" % (meanDB, stdDB))
    with open("results_summary.txt", "w") as fw:
        fw.write("Lagtime %d\n" % parameters.lagtime)
        fw.write("Number of clusters %d\n" % parameters.nclusters)
        fw.write("\n")
        fw.write("=====\n")
        fw.write("dG\n")
        fw.write("bound    Delta G     Delta W     Binding Volume:     Binding Volume contribution\n")
        for el in deltaGs:
            fw.write("%s\n" % el)
        fw.write("=====\n")
        fw.write("Asymmetric fluxes (see D.Lecina PhD thesis for more info)\n")
        for el in detailedBalance:
            fw.write("%.3f\n" % el)
        fw.write("\n")
        fw.write("dG = %f +- %f\n" % (meanDG, stdDG))
        fw.write("Asymmetric flux = %f +- %f\n" % (meanDB, stdDB))

    return meanDG, stdDG, meanDB, stdDB

if __name__ == "__main__":
    params = Parameters(ntrajs=None, length=None, lagtime=25, nclusters=100,
                        nruns=1, skipFirstSteps=0, useAllTrajInFirstRun=True,
                        computeDetailedBalance=True, trajWildcard="traj_*",
                        folderWithTraj="rawData", lagtimes=[1, 10, 25, 50],
                        clusterCountsThreshold=0, clusteringStride=1)
    estimateDG(params, cleanupClusterCentersAtStart=False)
