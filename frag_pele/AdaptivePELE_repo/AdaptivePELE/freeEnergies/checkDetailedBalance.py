from __future__ import absolute_import, division, print_function, unicode_literals
import glob
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg

FOLDER = "discretized"
CLUSTER_CENTERS = "clusterCenters.dat"
TRAJECTORY_MATCHING_PATTERN = "*.disctraj"


def parseArguments():
    desc = "Program that analyses a column of data, printing different statistical values and a histogram if desired."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", default='discretized', help="Folder with cluster contents")
    parser.add_argument("-t", "--threshold", default=0, type=float, help="If Cij < threshold and Cji < threshold, Cij = Cji = 0")
    parser.add_argument("-l", "--lagtime", default=1, type=int, help="Lagtime")
    parser.add_argument("-p", "--printFigs", action="store_true")
    args = parser.parse_args()

    return args.f, args.threshold, args.lagtime, args.printFigs


def getNumberOfClusters(folder_name):
    clusterCentersFilename = os.path.join(folder_name, CLUSTER_CENTERS)
    clusterCentersFile = open(clusterCentersFilename, 'r')
    clusterCenters = clusterCentersFile.readlines()
    clusterCentersFile.close()
    return len(clusterCenters)


def normalizePopulations(populations):
    return populations / populations.sum()


def normalizeTransitions(transitions):
    transitionsSum = transitions.sum(axis=1, dtype='float')
    transitionsSum = np.transpose(transitionsSum)
    normalizedPopulations = np.divide(transitions, transitionsSum[:, np.newaxis])
    return normalizedPopulations


def computeCountsAndCountMatrix(trajectories, numberOfClusters, lag_time=1):

    transitions = np.zeros((numberOfClusters, numberOfClusters))
    populations = np.zeros(numberOfClusters)

    for trajectoryFilename in trajectories:
        dtraj = np.loadtxt(trajectoryFilename, dtype=int, ndmin=1)
        # can be done much faster with sparse matrices (see runMarkovChain script)
        for i in range(len(dtraj) - lag_time):
            fromCluster = dtraj[i]
            toCluster = dtraj[i + lag_time]
            transitions[fromCluster][toCluster] += 1
        populations[dtraj] += 1

    return populations, transitions


def removeNoise(countMatrix, threshold):
    for i in range(len(countMatrix)):
        for j in range(i+1, len(countMatrix[i])):
            if countMatrix[i][j] < threshold or countMatrix[j][i] < threshold:
                countMatrix[i][j] = 0
                countMatrix[j][i] = 0
    return countMatrix


def computeCountsAndCountMatrixRemovingNoise(folder_name, countsThreshold, lag_time):
    """
        Computes the number of counts per cluster (#times it's been visited)
        and the count matrix, setting to 0 those counts below the threshold
    """
    trajectoryMatchingPattern = os.path.join(folder_name, TRAJECTORY_MATCHING_PATTERN)
    trajectories = glob.glob(trajectoryMatchingPattern)

    nclusters = getNumberOfClusters(folder_name)
    countsPerCluster, countMatrix = computeCountsAndCountMatrix(trajectories, nclusters, lag_time)

    countMatrixWithoutNoise = removeNoise(countMatrix, countsThreshold)
    countMatrixWithoutNoise += 1./nclusters  # added pseudo counts to avoid nans in relative entropy calc

    return countsPerCluster, countMatrixWithoutNoise


def computePopulationsAndTransitionProbabilities(folder, countsThreshold, lagtime):
    counts, countMatrix = computeCountsAndCountMatrixRemovingNoise(folder, countsThreshold, lagtime)

    # countMatrix = np.array(transitions)
    # sparseCountMatrix = sparse.csr_matrix(countMatrix)
    # sparseInString = sparseCountMatrix.__str__()
    # outputFile = open("countMatrix.dat", "w")
    # outputFile.write(sparseInString)
    # outputFile.close()

    populations = normalizePopulations(counts)
    transitions = normalizeTransitions(countMatrix)

    # for i in range(len(transitions)):
    #     for j in range(i+1, len(transitions[i])):
    #         if transitions[i][j] < probabilityThershold and transitions[j][i] < probabilityThershold:
    #             transitions[i][j] = 0
    #             transitions[j][i] = 0

    return populations, transitions


def plotMatrix(figureNumber, titleString, matrix, cmap=''):
    fig = plt.figure(figureNumber)
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    plt.title(titleString, fontsize=20)
    if cmap != '':
        cax = plt.pcolor(matrix, alpha=0.9, cmap=cmap)
    else:
        cax = plt.pcolor(matrix, alpha=0.9)

    # add colorbar if does not exist
    # if exists, remove previous, and add new one
    if len(plt.gcf().axes) > 1:
        # if so, then the last axes must be the colorbar.
        # we get its extent
        pts = plt.gcf().axes[-1].get_position().get_points()
        # and its label
        label = plt.gcf().axes[-1].get_ylabel()
        # and then remove the axes
        plt.gcf().axes[-1].remove()
        # then we draw a new axes a the extents of the old one
        cax = plt.gcf().add_axes([pts[0][0], pts[0][1], pts[1][0]-pts[0][0], pts[1][1]-pts[0][1]])
        # and add a colorbar to it
        cbar = plt.colorbar(cax=cax)
        cbar.ax.set_ylabel(label)
        # unfortunately the aspect is different between the initial call to colorbar
        #   without cax argument. Try to reset it (but still it's somehow different)
        cbar.ax.set_aspect(20)
    else:
        plt.colorbar(cax)


def barPlot(figureNumber, titleString, array):
    fig = plt.figure(figureNumber)
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    plt.title(titleString)
    ax.bar(np.arange(len(array)), array, width=1, alpha=0.4, color='b')


def getRelativeEntropy(goldenStationary, goldenT, T):
    return np.dot(goldenStationary, goldenT*np.log(goldenT/T)).sum()


def main(folder, countsThreshold, lagtime, printFigs=False):
    populations, transitions = computePopulationsAndTransitionProbabilities(folder, countsThreshold, lagtime)

    if printFigs:
        cmap = matplotlib.cm.coolwarm
        cmap.set_bad('w', 1.)
        matplotlib.rc('font', family='Times New Roman')

    if printFigs:
        # transitionsWithoutDiagonal = np.copy(transitions)
        # for i in range(len(transitions)):
        #     transitionsWithoutDiagonal[i][i] = 0
        # plotMatrix(6, r'$P_{ij}, P_{ii} = 0 \forall i$', transitionsWithoutDiagonal,cmap)
        pass

    # p_i * P_ij
    detailedBalanceComponents = np.copy(transitions)
    for i in range(len(detailedBalanceComponents)):
        detailedBalanceComponents[i] = np.multiply(populations[i], transitions[i])

    if printFigs:
        # barPlot(1, 'Population', populations)

        # plotMatrix(2, r'$P_{ij}$', transitions, cmap)
        plotMatrix(3, r'$\pi_i P_{ij}$', detailedBalanceComponents, cmap)
        # plt.savefig("db_flux.eps")

    detailedBalanceComponentsAbsoluteDifference = np.absolute(detailedBalanceComponents - detailedBalanceComponents.T) / 2.  # factor 2 to avoid the metric to go from 0 to 2, but from 0 to 1
    detailedBalanceComponentsAverage = np.multiply(detailedBalanceComponents + detailedBalanceComponents.T, 0.5)
    frobeniusAvg = linalg.norm(detailedBalanceComponentsAbsoluteDifference) / linalg.norm(detailedBalanceComponentsAverage)
    print("|semidiff| / |average|", frobeniusAvg)

    np.seterr(divide='ignore', invalid='ignore')
    detailedbalanceComponentsRelativeDifference = np.divide(detailedBalanceComponentsAbsoluteDifference, detailedBalanceComponentsAverage)
    maskedDetailedbalanceComponentsRelativeDifference = np.ma.array(detailedbalanceComponentsRelativeDifference, mask=np.isnan(detailedbalanceComponentsRelativeDifference))

    if printFigs:
        plotMatrix(4, r'$|F_{ij} - F_{ji}|$', detailedBalanceComponentsAbsoluteDifference, cmap)
        plt.savefig("db_abs_diff.eps")
        plotMatrix(5, r'$M(F)$', maskedDetailedbalanceComponentsRelativeDifference, cmap)
        plt.savefig("db_frobenius.eps")

    return frobeniusAvg


if __name__ == '__main__':
    folder_path, counts_Threshold, lagTime, print_Figs = parseArguments()
    main(folder_path, counts_Threshold, lagTime, print_Figs)
