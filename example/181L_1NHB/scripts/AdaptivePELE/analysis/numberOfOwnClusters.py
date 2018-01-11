import socket
machine = socket.gethostname()
import matplotlib
if machine == "bsccv03":
   matplotlib.use('wxagg')
import matplotlib.pyplot as plt
import numpy as np
import os
import collections
import argparse


def printHelp():
    desc = "Program that prints the number of clusters throughout an adaptive sampling simulation. "\
            "It must be run in the root folder. "
    parser = argparse.ArgumentParser(description=desc)
    args = parser.parse_args()
    return args


def getClusteringSummaryContent(file):
    if os.path.isfile(file):
        summaryContent = np.genfromtxt(file)
        if summaryContent.ndim > 1:
            return summaryContent
        elif summaryContent.ndim == 1:
            # file has only one line
            return np.array([summaryContent])
        else:
            # file existed but was empty
            return []
    else:
        return []


def getTotalNumberOfClustersPerEpoch(templetizedClusteringSummaryFile, folder):
    allFolders = os.listdir(folder)
    numberOfEpochs = len([epoch for epoch in allFolders if epoch.isdigit() and os.path.isfile(templetizedClusteringSummaryFile%int(epoch))])

    totalNumberOfClustersPerEpoch = []
    for epoch in range(numberOfEpochs):
        clusteringSummary = getClusteringSummaryContent(templetizedClusteringSummaryFile % epoch)

        if clusteringSummary != []:
            totalNumberOfClustersPerEpoch.append(len(clusteringSummary))

    return totalNumberOfClustersPerEpoch

def findDifferentClustersInEpoch(column, summaryFile):
    """
        Returns a dictionary with the set of different elements in column and the number of elements in this epoch
    """
    clusteringSummary = getClusteringSummaryContent(summaryFile)

    epochDictionary = {}
    if clusteringSummary != []:
        for line in clusteringSummary:
            value = line[column]
            if not value in epochDictionary:
                epochDictionary[value] = len(np.argwhere(clusteringSummary[:,column] == value))
    return epochDictionary

def findDifferentClustersForAllEpochs(column, templetizedClusteringSummaryFile, numberOfEpochs):
    """
        Returns a list with dictionaries for all epochs. The dictionary has the set of different values (according to column) and their number
    """
    clustersPerEpoch = []
    for epoch in range(numberOfEpochs):
        summaryFile = templetizedClusteringSummaryFile%epoch
        epochDictionary = findDifferentClustersInEpoch(column, summaryFile)

        clustersPerEpoch.append(epochDictionary)
    return clustersPerEpoch

def getAllDifferentValues(clustersPerEpoch):
    allValues = set()
    for epochSummary in clustersPerEpoch:
        for value, numClusters in epochSummary.iteritems():
            allValues.update([value])
    return allValues

def buildClustersPerValue(clustersPerEpoch, numberOfEpochs):
    """
        Returns dictionary with lists for the different values. The length of the list is equal to the number of "clustering/summary.txt" files found
    """
    clustersPerValue = collections.defaultdict(list)

    allValues = getAllDifferentValues(clustersPerEpoch)

    for epochSummary in clustersPerEpoch:
        foundValues = set()
        for value, numClusters in epochSummary.iteritems():
            clustersPerValue[value].append(numClusters)
            foundValues.update([value])

        for value in allValues - foundValues:
            clustersPerValue[value].append(0)

    return clustersPerValue

def getNumberOfClustersPerEpochForGivenColumn(column, templetizedClusteringSummaryFile, folder):
    allFolders = os.listdir(folder)
    numberOfEpochs = len([epoch for epoch in allFolders if epoch.isdigit() and os.path.isfile(templetizedClusteringSummaryFile%int(epoch))])

    clustersPerEpoch = findDifferentClustersForAllEpochs(column, templetizedClusteringSummaryFile, numberOfEpochs)

    return buildClustersPerValue(clustersPerEpoch, numberOfEpochs)


def plotClustersPerValue(clustersPerValue):
    values = clustersPerValue.keys()
    sortedValues = np.sort(values)
    for value in sortedValues:
        plt.plot(clustersPerValue[value], label=str(value))


def plotContactsHistogram(folder, templetizedClusteringSummaryFile):
    allFolders = os.listdir(folder)
    lastEpoch = len([epoch for epoch in allFolders if epoch.isdigit() and os.path.isfile(templetizedClusteringSummaryFile%int(epoch))]) - 1
    lastSummary = templetizedClusteringSummaryFile%lastEpoch
    contactsColumn = 3
    allContacts = np.loadtxt(lastSummary, usecols=(contactsColumn,), ndmin=1)
    plt.hist(allContacts)

def main():
    printHelp()

    #Params
    clusteringFileDensityColumn = 5
    clusteringFileThresholdColumn = 4
    clusteringFolder = "clustering"
    summaryFile = "summary.txt"
    folder = "."
    # end params

    clusteringSummaryFile = os.path.join(clusteringFolder, summaryFile)
    templetizedClusteringSummaryFile = os.path.join("%d", clusteringSummaryFile)


    totalNumberOfClustersPerEpoch = getTotalNumberOfClustersPerEpoch(templetizedClusteringSummaryFile, folder)
    clustersPerDensityValue = getNumberOfClustersPerEpochForGivenColumn(clusteringFileDensityColumn, templetizedClusteringSummaryFile, folder)
    clustersPerThresholdValue = getNumberOfClustersPerEpochForGivenColumn(clusteringFileThresholdColumn, templetizedClusteringSummaryFile, folder)

    plt.figure(1)
    plt.plot(totalNumberOfClustersPerEpoch, label="All clusters")

    plotClustersPerValue(clustersPerDensityValue)
    plt.legend(loc=2)
    # plt.title("n=64, different thresholds, variable density")
    # plt.savefig("../3ptb_4_64_numberOfClusters_density_corner.png")

    plt.figure(2)
    plt.plot(totalNumberOfClustersPerEpoch, label="All clusters")
    plotClustersPerValue(clustersPerThresholdValue)
    plt.legend(loc=2)
    # plt.title("n=64, different thresholds, variable density")
    # plt.savefig("set2.png")

    plt.figure(3)
    plotContactsHistogram(folder, templetizedClusteringSummaryFile)

    plt.show()

if __name__ == "__main__":
    main()
