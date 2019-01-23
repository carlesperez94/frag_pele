import os
import numpy as np
import firstSnapshots
from pyemma_scripts import buildMSM, computeDeltaG
from simulation import simulationrunner
import shutil
import matplotlib.pyplot as plt

def rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def computeDG(lagtime, clusters):
    simulationParameters = simulationrunner.SimulationParameters()
    controlFile = "templetized_control_MSM.conf"
    simulationParameters.templetizedControlFile = controlFile
    sr = simulationrunner.SimulationRunner(simulationParameters)

    controlFileDictionary = {"lagtime": lagtime, "clusters": clusters}
    sr.makeWorkingControlFile("control_MSM.conf", controlFileDictionary)
    buildMSM.main("control_MSM.conf")
    deltaGLine = computeDeltaG.main("traj_*")
    return deltaGLine 

def estimateDGValue(nruns, lagtime, clusters, length, ntraj, col=1):
    deltaGs = []
    for i in range(nruns):
        firstSnapshots.main(length, ntraj)
        #rm("clustering_object.pkl") 
        rm("MSM_object.pkl") 
        deltaGLine = computeDG(lagtime, clusters)
        dG = float(deltaGLine.split()[1])
        deltaGs.append(dG)
    return np.mean(deltaGs), np.std(deltaGs)

def lengthVsNtrajs(nruns, lagtime, clusters, lengths, ntrajs, col=1):
    nLengths = len(lengths)
    nNtrajs = len(ntrajs)
    results = np.zeros((nLengths, nNtrajs))
    stdDev = np.zeros((nLengths, nNtrajs))
    for i, length in enumerate(lengths):
        for j, ntraj in enumerate(ntrajs):
            print "computing for length:%d and ntrajs:%d"%(length, ntraj)
            results[i][j], stdDev[i][j] = estimateDGValue(nruns, lagtime, clusters,\
                                                         length, ntraj, col)
    return results, stdDev

def main():
    """
        Makes a plot à la relative entropy, but with free energies
    """

    """
    lagtime = 200
    clusters = 100
    lengths = [400, 500, 600, 700, 800, 900, 1000]
    lengths = [500, 600, 700, 800, 900, 1000]
    col = 1 #delta G
    ntrajs = [32-1, 64-1, 128-1, 256-1, 512-1]
    ntrajs = range(64-1, 512, 64)
    nruns = 10
    results,stdDev = lengthVsNtrajs(nruns, lagtime, clusters, lengths, ntrajs, col)
    np.save("results.npy", results)
    np.save("stdDev.npy", stdDev)
    """
    results = np.load("results.npy")

    
    extent = [64-32,512+32,500-50,1000+50]
    plt.figure(1)
    plt.imshow(results, interpolation="nearest", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.figure(2)
    plt.imshow(results, interpolation="bilinear", origin="lower", aspect="auto", extent=extent)
    plt.colorbar()
    plt.show()
    

if __name__ == "__main__":
    main()
