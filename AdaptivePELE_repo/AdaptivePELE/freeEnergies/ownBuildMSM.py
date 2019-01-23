import argparse
import json
import sys
from AdaptivePELE.freeEnergies import cluster
from AdaptivePELE.freeEnergies import estimate


def parseArgs():
    parser = argparse.ArgumentParser(description="Build a MSM using PyEMMA "
                                     "package")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    return args


def readParams(control_file):
    try:
        with open(control_file, "r") as f:
            params = json.load(f)
    except:
        sys.exit("Sth wrong with the control file: %s" % control_file)

    trajectoryFolder = params["trajectoryFolder"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    stride = params.get("stride", 1)
    lagtimes = params.get("lagtimes", [])
    lagtime = params.get("lagtime", None)
    numPCCA = params.get("numPCCA", None)
    itsOutput = params.get("itsOutput", None)
    numberOfITS = params.get("numberOfITS", -1)
    errors = params.get("errors", False)
    mlags = params.get("mlags", None)
    clusterCountsThreshold = params.get("clusterCountsThreshold", 0)

    return trajectoryFolder, trajectoryBasename, numClusters, stride, lagtimes, numPCCA, itsOutput, numberOfITS, errors, mlags, lagtime, clusterCountsThreshold


def main(control_file):

    # parameters
    trajectoryFolder, trajectoryBasename, numClusters, stride, lagtimes, _, _, numberOfITS, _, _, lagtime, clusterCountsThreshold = readParams(control_file)

    # program
    clusteringObject = cluster.Cluster(numClusters, trajectoryFolder, trajectoryBasename, alwaysCluster=False, stride=stride)
    clusteringObject.clusterTrajectories()
    clusteringObject.eliminateLowPopulatedClusters(clusterCountsThreshold)
    calculateMSM = estimate.MSM(error=False, dtrajs=clusteringObject.dtrajs)
    calculateMSM.estimate(lagtime=lagtime, lagtimes=lagtimes, numberOfITS=numberOfITS)


if __name__ == "__main__":
    arguments = parseArgs()
    main(arguments.controlFile)
