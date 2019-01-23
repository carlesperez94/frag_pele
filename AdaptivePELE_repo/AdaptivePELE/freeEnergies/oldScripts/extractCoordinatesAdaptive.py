#! /usr/bin/env python
import glob 
import sys
import os
import argparse

def loadAllAtomsInPdb(filename, atoms):
    snapshotNumber = 0

    allCoordinates = []
    coordinates = []

    file = open(filename,'r')
    for line in file:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            line = line.split()
            if line[1] in atoms:
                coordinates.append(line[6])
                coordinates.append(line[7])
                coordinates.append(line[8])

        if line[:6] == "ENDMDL":
            allCoordinates.append([str(snapshotNumber)] + coordinates)
            coordinates = []
            snapshotNumber += 1

    file.close()

    return allCoordinates

def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('.', 0, last)

    return filename[first+1:last]

def getOutputFilename(directory, filename, baseOutputFilename):
    filenumber = extractFilenumber(filename)
    return directory + baseOutputFilename + "_" + filenumber + ".dat"

def writeToFile(allCoordinates, outputFilename):

    with open(outputFilename, 'w') as f:
        for line in allCoordinates:
            for i in range(len(line) - 1):
                f.write(line[i] + ' ')
            f.write(line[-1] + '\n')


def parseArguments():
    desc = "Extracts coordinates in <currentFolder>/extractedCoordinates/coord*"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("folderTraj",
                        help="Folder where the epochs trajs are stored")
    parser.add_argument("-f", default="traj_*", nargs='+', help="Traj Files wildcard")
    parser.add_argument("-s", nargs='+', help="List of atoms")
    args = parser.parse_args()

    return  args.f, args.s, args.folderTraj




def main():
    originalPDBfileWildcard, atoms, folderTraj = parseArguments()

    allFolders = os.listdir(folderTraj)
    Epochs = [epoch for epoch in allFolders if epoch.isdigit()]
    for folder in Epochs:
        currentEpochFolder = os.path.join(folderTraj, folder)
        directory = 'extractedCoordinates/'
        directory = os.path.join(currentEpochFolder, directory)
        if not os.path.exists(directory):
                os.makedirs(directory)

        baseOutputFilename = 'coord'

        print 'Atoms to extract coordinates from ', atoms

        originalPDBfile = glob.glob(os.path.join(currentEpochFolder,originalPDBfileWildcard)) 
        for filename in originalPDBfile:
            print filename
            allCoordinates = loadAllAtomsInPdb(filename, atoms)
            outputFilename = getOutputFilename(directory, filename, baseOutputFilename)
            writeToFile(allCoordinates, outputFilename)




if __name__ == "__main__":
    main()

