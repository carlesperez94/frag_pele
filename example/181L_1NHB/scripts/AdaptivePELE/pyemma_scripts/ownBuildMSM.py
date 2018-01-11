import argparse
import cluster
import estimate
import json

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
        import sys
        sys.exit("Sth wrong with the control file: %s" % control_file)

    trajectoryFolder = params["trajectoryFolder"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    lagtimes = params.get("lagtimes", [])
    lagtime = params.get("lagtime", None)
    numPCCA = params.get("numPCCA", None)
    itsOutput = params.get("itsOutput", None)
    numberOfITS = params.get("numberOfITS", -1)
    errors = params.get("errors", False)
    mlags = params.get("mlags", None)

    return trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, errors, mlags, lagtime

def main(control_file):

    # parameters
    trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, errors, mlags, lagtime = readParams(control_file)

    # program
    clusteringObject = cluster.Cluster(numClusters, trajectoryFolder, trajectoryBasename, alwaysCluster=False)
    calculateMSM = estimate.MSM(error=False, dtrajs=clusteringObject.dtrajs)
    calculateMSM.estimate(lagtime=lagtime, lagtimes=lagtimes, numberOfITS=numberOfITS)



if __name__ == "__main__":
    args = parseArgs()
    main(args.controlFile)
