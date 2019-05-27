import MSMblocks
import numpy as np
import argparse


def parseArgs():
    parser = argparse.ArgumentParser(description="Build a MSM using PyEMMA "
                                     "package")
    parser.add_argument('controlFile', type=str)
    args = parser.parse_args()
    return args


def readParams(control_file):
    params = MSMblocks.readParams(control_file)
    trajectoryFolder = params["trajectoryFolder"]
    trajectoryBasename = params["trajectoryBasename"]
    numClusters = params["numClusters"]
    lagtimes = params["lagtimes"]
    lagtime = params.get("lagtime", None)
    numPCCA = params["numPCCA"]
    itsOutput = params["itsOutput"]
    numberOfITS = params["numberOfITS"]
    itsErrors = params["itsErrors"]
    error_estimationCK = params["error_estimationCK"]
    state_labels = params["state_labels"]
    mlags = params["mlags"]
    stride = params.get("stride", 1)
    if state_labels is None:
        state_labels = 'auto'  # json returns string as
    # unicode, and this breaks some code in pyemma
    outfile_fluxTPT = params["outfile_fluxTPT"]
    return trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, itsErrors, error_estimationCK, state_labels, outfile_fluxTPT, mlags, lagtime, stride


def main(control_file):

    # parameters
    trajectoryFolder, trajectoryBasename, numClusters, lagtimes, numPCCA, itsOutput, numberOfITS, itsErrors, error_estimationCK, state_labels, outfile_fluxTPT, mlags, lagtime, stride = readParams(control_file)

    # program
    prepareMSM = MSMblocks.PrepareMSM(numClusters, trajectoryFolder, trajectoryBasename, stride=stride)
    cl = prepareMSM.getClusteringObject()
    calculateMSM = MSMblocks.MSM(cl, lagtimes, numPCCA, itsOutput, numberOfITS,
                                 itsErrors, error_estimationCK, mlags, lagtime, prepareMSM.dtrajs)
    calculateMSM.estimate()
    MSM_object = calculateMSM.getMSM_object()
    #calculateMSM.writeClustersForWMD()
    # TPTinstance = MSMblocks.TPT(MSM_object, cl, outfile_fluxTPT, state_labels)
    # TPT_Object = TPTinstance.getTPTObject()
    # coarseTPT_Object = TPTinstance.getCoarseTPTObject()



if __name__ == "__main__":
    args = parseArgs()
    main(args.controlFile)
