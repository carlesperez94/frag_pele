from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import os
import glob


templetizedEpochFolder = "%d"


def findTrajectoryFirstBindingEvent(metrics, thresholdValue, unBinding=False):
    if unBinding:
        firstBindingEvent = np.argmax(metrics[:, 1] >= thresholdValue)
        if metrics[firstBindingEvent, 1] >= thresholdValue:
            return metrics[firstBindingEvent, 0]
        else:
            return None
    else:
        # argmax returns "0" if no argument matches
        firstBindingEvent = np.argmax(metrics[:, 1] <= thresholdValue)
        if metrics[firstBindingEvent, 1] <= thresholdValue:
            return metrics[firstBindingEvent, 0]
        else:
            return None


def findEpochFirstBindingEvent(thresholdValue, columnInReport, reportWildcard="*report_*", unBinding=False):
    foundBindingEvent = False
    minNumberOfSteps = 1e10

    reportFiles = glob.glob(reportWildcard)
    for reportFile in reportFiles:
        metrics = np.loadtxt(reportFile, usecols=(1, columnInReport), ndmin=2)
        if metrics.shape[0] == 0:
            continue
        firstBindingEvent = findTrajectoryFirstBindingEvent(metrics, thresholdValue, unBinding=unBinding)
        if firstBindingEvent is not None:
            foundBindingEvent = True
            if firstBindingEvent < minNumberOfSteps:
                minNumberOfSteps = firstBindingEvent
    if foundBindingEvent:
        return minNumberOfSteps
    else:
        return None


def getAllSortedEpochs():
    allFolders = os.listdir(".")
    epochFolders = sorted([int(epoch) for epoch in allFolders if epoch.isdigit()])
    return epochFolders


def findFirstBindingEvent(stepsPerEpoch, columnInReport, thresholdValue, unBinding=False):
    epochFolders = getAllSortedEpochs()

    for epoch in epochFolders:
        os.chdir(str(epoch))
        epochFirstBindingEvent = findEpochFirstBindingEvent(thresholdValue, columnInReport, unBinding=unBinding)

        if epochFirstBindingEvent is not None:
            return epochFirstBindingEvent + epoch * stepsPerEpoch

        os.chdir("..")

    return None


def main(stepsPerEpoch=10, columnInReport=4, thresholdValue=2):
    return findFirstBindingEvent(stepsPerEpoch, columnInReport, thresholdValue)


if __name__ == "__main__":
    stepsForBE = main()
    print(stepsForBE)
