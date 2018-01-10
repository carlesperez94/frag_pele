import os
import numpy as np
import estimateDG
import shutil
import matplotlib.pyplot as plt


def lengthVsNtrajs(nruns, lagtime, clusters, lengths, ntrajs, outputFilename, cache, skipFirstSnaphots):
    nLengths = len(lengths)
    nNtrajs = len(ntrajs)
    results = np.zeros((nLengths, nNtrajs))
    stdDev = np.zeros((nLengths, nNtrajs))
    db = np.zeros((nLengths, nNtrajs))
    stdDb = np.zeros((nLengths, nNtrajs))
    for i, length in enumerate(lengths):
        for j, ntraj in enumerate(ntrajs):
            if (length, ntraj) in cache:
                print "Loading cached computation for length:%d and ntrajs:%d"%(length, ntraj)
                results[i][j], stdDev[i][j], db[i][j], stdDb[i][j] = cache[(length, ntraj)]
                with open(outputFilename, 'a') as f:
                    f.write("%d %d %f %f %f %f\n" %(length, ntraj, results[i][j], stdDev[i][j], db[i][j], stdDb[i][j]))
                continue
            print "Computing for length:%d and ntrajs:%d"%(length, ntraj)
            parameters = estimateDG.Parameters(ntrajs=ntraj, length=length, lagtime=lagtime, nclusters=clusters, nruns=nruns,
                            useAllTrajInFirstRun=False, computeDetailedBalance=True, trajWildcard="traj_*",
                            folderWithTraj="rawData", skipFirstSteps=skipFirstSnaphots)

            origFilesWildcard = os.path.join(parameters.folderWithTraj, parameters.trajWildcard)
            copiedFiles = estimateDG.copyWorkingTrajectories(origFilesWildcard, length, ntraj, bootstrap=True)
            results[i][j], stdDev[i][j], db[i][j], stdDb[i][j] = estimateDG.estimateDG(parameters)
            with open(outputFilename, 'a') as f:
                f.write("%d %d %f %f %f %f\n" %(length, ntraj, results[i][j], stdDev[i][j], db[i][j], stdDb[i][j]))
    return results, stdDev, db, stdDb

def saveResultsFileBckp(outputFilename):
    i = 1
    bckpFilename = outputFilename+".%d.bckp"
    while os.path.isfile(bckpFilename%i): i += 1
    try:
        shutil.copy(outputFilename, bckpFilename%i)
    except IOError:
        pass

def main():
    lagtime = 25
    clusters = 100

    ilengths = 200
    flengths = 4000
    dlengths = 200
    lengths = range(ilengths, flengths, dlengths)
    itrajs = 32-1
    ftrajs = 256
    dtrajs = 32
    ntrajs = range(itrajs, ftrajs, dtrajs)
    nruns = 10
    skipFirstSnaphots = 0
    outputFilename = "results.txt"
    cache = {}
    if os.path.exists(outputFilename):
        with open(outputFilename) as fr:
            for line in fr:
                contents = line.rstrip().split()
                if not contents[0].isdigit():
                    continue
                l, traj, dg, stdDG, db, stdDB = map(float, contents)
                cache[(int(l), int(traj))] = (dg, stdDG, db, stdDB)

    saveResultsFileBckp(outputFilename)

    with open(outputFilename, 'w') as f:
        f.write("Computing DG, stdDG, DB and stdDB for different lengths and number of trajectories\n")
        f.write("Lengths: %s\n"% lengths)
        f.write("Ntrajs: %s\n"% ntrajs)
        f.write("Skipping first: %d snapshots of each trajectory\n"%skipFirstSnaphots)
    results,stdDev,db,stdDb = lengthVsNtrajs(nruns, lagtime, clusters, lengths, ntrajs, outputFilename, cache, skipFirstSnaphots)
    np.save("results.npy", results)
    np.save("stdDev.npy", stdDev)
    np.save("db.npy", db)
    np.save("stdDb.npy", stdDb)

    #results = np.load("results.npy")


    extent = [itrajs+1-dtrajs/2,ftrajs+dtrajs/2,ilengths-dlengths/2,flengths+dlengths/2] # +1 for aesthetical purposes
    plt.figure(1)
    plt.imshow(results, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.figure(2)
    plt.imshow(results, interpolation="bilinear", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    main()
