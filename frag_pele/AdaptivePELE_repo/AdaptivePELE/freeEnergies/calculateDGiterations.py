from __future__ import absolute_import, division, print_function, unicode_literals
from six import reraise as raise_
import os
import sys
import shutil
from AdaptivePELE.freeEnergies import estimateDG, prepareMSMFolders

def isfinished(folders):
    for folder in folders:
        if not os.path.exists(os.path.join(folder, "results_summary.txt")):
                return False
    return True


# list of tuples with format (lagtime, #clusters)
iterations = [(25, 100), (50, 100), (100, 100), (200, 100), (400, 100),
              (25, 200), (50, 200), (100, 200), (200, 200), (400, 200),
              (25, 400), (50, 400), (100, 400), (200, 400), (400, 400)]
runFolder = os.getcwd()
print("Running from " + runFolder)
for tau, k in iterations:
    destFolder = "%d/%dcl" % (tau, k)
    if not os.path.exists(destFolder):
        os.makedirs(destFolder)
    os.chdir(destFolder)
    folders_MSM = glob.glob("MSM_*")
    if isfinished(folders_MSM):
        print("Skipping run with lagtime %d, clusters %d" % (tau, k))
        os.chdir(runFolder)
        continue
    else:
        for folder in folders_MSM:
            shutil.rmtree(folder)
    prepareMSMFolders.main(trajsPath=runFolder)
    print("***************")
    print("Estimating dG value in folder" + os.getcwd())
    parameters = estimateDG.Parameters(ntrajs=None, length=None, lagtime=tau,
                                       nclusters=k, nruns=10, skipFirstSteps=0,
                                       useAllTrajInFirstRun=True,
                                       computeDetailedBalance=True,
                                       trajWildcard="traj_*",
                                       folderWithTraj="rawData",
                                       lagtimes=[1, 10, 25, 50, 100, 250, 500, 1000],
                                       clusterCountsThreshold=0)
    try:
        estimateDG.estimateDG(parameters, cleanupClusterCentersAtStart=True)
    except Exception as err:
        if "distribution contains entries smaller" in str(err):
            print("Caught exception in step with lag %d and k %d, moving to next iteration" % (tau, k))
            with open("error.txt", "w") as fe:
                fe.write("Caught exception in step with lag %d and k %d, moving to next iteration\n" % (tau, k))
        else:
            t, v, tb = sys.exc_info()
            raise_(t, v, tb)
    os.chdir(runFolder)
