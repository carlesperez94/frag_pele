import computeRelativeEntropy
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def plot(extent, entropies, allTrajLengths, numberOfTrajs, vmin, vmax):
    #plot isocost lines
    computeRelativeEntropy.plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 9)
    plt.imshow(entropies, interpolation="nearest", origin="lower", aspect="auto", extent=extent, vmin=vmin, vmax=vmax)
    cwd = os.getcwd()
    cwd = cwd.replace("/", "_")
    #plt.save(cwd + ".eps")
    #plt.show()
    #plt.imshow(entropies, interpolation="nearest", extent=[numberOfTrajs[0], numberOfTrajs[-1], 800, 801])
    plt.colorbar()
    plt.savefig(cwd + ".eps")
    plt.figure(2)
    computeRelativeEntropy.plotIsocostLines(extent, allTrajLengths, numberOfTrajs, 9)
    plt.imshow(entropies, interpolation="bilinear", origin="lower", aspect="auto",  extent=extent, vmin=vmin, vmax=vmax)
    plt.colorbar()
    cwd = os.getcwd()
    cwd = cwd.replace("/", "_")
    #plt.save(cwd + "2.eps")
    plt.savefig(cwd + "2.eps")

def main(controlFile):
    """
        Takes cluster centers file, builds dtrajs, and computes relative entropy
    """
    disctrajFolder, trajectoryFolder, trajectoryFolder2, trajectoryBasename, numClusters, lagtime, sampleSize, numRuns, dtrajs, maxNtraj, minNtraj, dlength, maxlength = computeRelativeEntropy.readParams(controlFile)

    """
    Params
    """
    if dtrajs is None:
        dtrajs = 100
    if maxNtraj is None:
        maxNtraj = len(discTrajs)
    if minNtraj is None:
        minNtraj = 100
    numberOfTrajs = range(minNtraj, maxNtraj, dtrajs)

    #dTrajs = 100
    # numberOfTrajs = range(50, sampleSize, 50)

    #only trying different traj lengths if sampleSize is defined in control file
    lowerLimit = 2*lagtime
    
    if maxlength is None:
        sys.exit("Set a maximum length in control file")
    if dlength is None:
        dlength = 100
    allTrajLengths = range(lowerLimit, maxlength, dlength)


    entropies = np.load("/Users/daniel/Desktop/1o3p_fullAdaptive/matrix_adaptive.npy")
    entropies = np.load("/Users/daniel/Desktop/az/PR_COORD/PR_cortiCOORD_3/matrix_adaptive.npy")
    #entropies = np.load("/Users/daniel/1f5k_1o3p_1sqa/1o3p_OBC/matrix_adaptive.npy")
    entropies = np.load("/Users/daniel/1f5k_1o3p_1sqa/1o3p_OBC_benchmark/matrix_adaptive.npy")
    entropies2 = np.load("/Users/daniel/Desktop/az/PR_COORD/PR_cortiCOORD_3/matrix_adaptive.npy")
    entropies = np.log10(entropies)
    entropies2 = np.log10(entropies2)
    vmin = min(np.min(entropies), np.min(entropies2))
    vmax = max(np.max(entropies), np.max(entropies2))

    for i, length in enumerate(allTrajLengths):
        for j, ntrajs  in enumerate(numberOfTrajs):
            print length, ntrajs, entropies[i][j]
        print ""

    extent = [numberOfTrajs[0] - dtrajs/2, numberOfTrajs[-1] + dtrajs/2,
                allTrajLengths[0] - dlength/2, allTrajLengths[-1] + dlength/2]
    plt.figure(1)
    plot(extent, entropies, allTrajLengths, numberOfTrajs, vmin=vmin, vmax=vmax)
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1])
