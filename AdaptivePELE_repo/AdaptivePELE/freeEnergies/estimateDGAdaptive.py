from __future__ import absolute_import, division, print_function, unicode_literals
import os
import sys
import glob
import numpy as np
from six import reraise as raise_
from AdaptivePELE.freeEnergies import estimateDG


def main(trajsPerEpoch, lagtime, nclusters, clusteringStride=1, nruns=10, lagtimes=[1, 10, 25, 50, 100, 250, 400, 500, 600, 1000]):
    allFolders = np.array(glob.glob("MSM_*"))
    epochs = [int(folder[4:]) for folder in allFolders]
    args = np.argsort(epochs)
    sortedFolders = allFolders[args]
    origDir = os.getcwd()
    resultsFile = os.path.join(origDir, "results.txt")

    with open(resultsFile, "a") as f:
        f.write("#Epoch DG StdDG Db StdDb\n")
        f.write("#=======================\n")

    resultsEpoch = []
    initialEpoch = 0
    for i, folder in enumerate(sortedFolders[initialEpoch:]):
        epoch = i + initialEpoch
        print(epoch, folder)
        os.chdir(folder)
        parameters = estimateDG.Parameters(ntrajs=trajsPerEpoch*(epoch+1),
                                           length=None,
                                           lagtime=lagtime,
                                           lagtimes=lagtimes,
                                           nclusters=nclusters,
                                           nruns=nruns,
                                           useAllTrajInFirstRun=True,
                                           computeDetailedBalance=True,
                                           trajWildcard="traj_*",
                                           folderWithTraj="rawData",
                                           clusterCountsThreshold=0,
                                           clusteringStride=clusteringStride)
        try:
            dG, stdDg, db, stdDb = estimateDG.estimateDG(parameters, cleanupClusterCentersAtStart=True)
            print("FINAL RESULTS EPOCH %d: dG: %f +- %f, asymmetric fluxes: %f +- %f" % (epoch, dG, stdDg, db, stdDb))
            resultsEpoch.append([dG, stdDg, db, stdDb])
            with open(resultsFile, "a") as f:
                f.write("%d %.3f %.3f %.3f %.3f\n" % (epoch, dG, stdDg, db, stdDb))
        except Exception as err:
            resultsEpoch.append(["Estimation in this epoch crashed"])
            if "distribution contains entries smaller" in str(err):
                print("Caught exception in folder %s with lag %d and k %d, moving to next iteration" % (folder, lagtime, nclusters))
                with open("error.txt", "w") as fe:
                    fe.write("Caught exception in folder %s with lag %d and k %d, moving to next iteration\n" % (folder, lagtime, nclusters))
            else:
                raise_(*sys.exc_info())
        os.chdir("..")

    print("Results")
    print("epoch, DG, stdDg, DB, stdDb")
    print("=====")
    for i, results in enumerate(resultsEpoch):
        print(i, results)

if __name__ == "__main__":
    trajs_Epoch = 50
    lag_time = 50
    n_clusters = 100
    clustering_Stride = 10
    main(trajs_Epoch, lag_time, n_clusters, clusteringStride=clustering_Stride)
