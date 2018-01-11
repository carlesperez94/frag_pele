import AdaptivePELE.atomset.atomset as atomset
import os
import glob
import numpy as np
from AdaptivePELE.utilities import utilities
import argparse
import json

def extendReportWithRmsd(reportFile, rmsds):
    newShape = reportFile.shape
    fixedReport = np.zeros((newShape[0], newShape[1]+1))
    fixedReport[:,:-1] = reportFile
    fixedReport[:,-1] = rmsds
    return fixedReport

def parseArguments():
    desc = "Program that fixes RMSD symmetries of a PELE report file."\
            "Control file is a JSON file that contains \"resname\", \"native\", "\
            "symmetries, and, optionally, the column to substitute in report. "\
            "Example of content:"\
            "{"\
            "\"resname\" : \"K5Y\","\
            "\"native\" : \"native.pdb\","\
            "\"symmetries\" : {[\"4122:C12:K5Y\":\"4123:C13:K5Y\", \"4120:C10:K5Y\":\"4127:C17:K5Y\"]},"\
            "\"column\" = 5"\
            "}"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("controlFile", type=str, help="Control File name")
    args = parser.parse_args()

    return  args.controlFile

def readControlFile(controlFile):
    jsonFile = open(controlFile, 'r').read()
    parsedJSON = json.loads(jsonFile)
    resname = parsedJSON["resname"]
    nativeFilename = parsedJSON["native"]
    symmetries = parsedJSON["symmetries"]
    rmsdColInReport = parsedJSON.get("column")
    if not rmsdColInReport:
        #append to the end
        rmsdColInReport = -1

    return resname, nativeFilename, symmetries, rmsdColInReport

def main(controlFile):
    #Constants
    folder = "."
    outputFilename = "fixedReport_%d"
    trajName = "*traj*.pdb"
    reportName = "*report_%d"
    #end constants

    resname, nativeFilename, symmetries, rmsdColInReport = readControlFile(controlFile)


    nativePDB = atomset.PDB()
    nativePDB.initialise(nativeFilename, resname=resname)

    allFolders = os.listdir(folder)
    epochs = [epoch for epoch in allFolders if epoch.isdigit()]

    for epoch in epochs:
        print "Epoch", epoch
        os.chdir(epoch)
        allTrajs = glob.glob(trajName)

        for traj in allTrajs:
            rmsds = utilities.getRMSD(traj, nativePDB, resname, symmetries)
            trajNum = utilities.getTrajNum(traj)
            try:
                reportFilename = glob.glob(reportName%trajNum)[0]
            except IndexError:
                raise IndexError("File %s not found in folder %s"%(reportName%trajNum, epoch))

            reportFile = np.loadtxt(reportFilename, ndmin=2)

            if rmsdColInReport > 0 and rmsdColInReport < reportFile.shape[1]:
                reportFile[:,rmsdColInReport] = rmsds
                fixedReport = reportFile
            else:
                fixedReport = extendReportWithRmsd(reportFile, rmsds)

            #print fixedReport
            np.savetxt(outputFilename%trajNum, fixedReport, fmt='%.4f')

        os.chdir("..")

if __name__ == "__main__":
    controlFile = parseArguments()
    main(controlFile)
