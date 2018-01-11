import os
import argparse
import glob
# import minimumGCluster.atomset as atomset
from AdaptivePELE.atomset import atomset
import re

DIRECTORY = '%s/extractedCoordinates'
BASEOUTPUTFILENAME = 'coord'


def loadAllAtomsInPdb(filename):
    fileContent = open(filename).read()
    fileContent = fileContent.split('ENDMDL')
    return fileContent


def parseArguments():
    desc = "Extracts coordinates in <currentFolder>/extractedCoordinates/coord*.\
            It either extracts the resname COM coordinates or those of an atomId, depending on the input"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-folderTraj", default = ".",
                        help="Folder where original trajs are found and stored")
    parser.add_argument("-atomId", default="", help="serial:atomName:resname, e.g. 2048:C1:AIN")
    parser.add_argument("-resname", default="", help="Ligand resname")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()
    return args.folderTraj, args.atomId, args.resname


def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('_')
    number = re.sub("[^0-9]", "", filename[first+1:last])

    return number


def getOutputFilename(directory, filename, baseOutputFilename):
    filenumber = extractFilenumber(filename)
    return os.path.join(directory, baseOutputFilename+"_"+filenumber+".dat")


def getPDBCOM(allCoordinates, resname):
    COMs = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        COMs.append(pdb.extractCOM())
    return COMs

def getAtomCoord(allCoordinates, resname, atomId):
    coords = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=resname, heavyAtoms=True)
        coords.append(pdb.getAtom(atomId).getAtomCoords())
    return coords

def writeToFile(COMs, outputFilename):
    with open(outputFilename, 'w') as f:
        for i, line in enumerate(COMs):
            f.write(str(i) + ' ')
            for i in range(len(line) - 1):
                f.write(str(line[i]) + ' ')
            f.write(str(line[-1]) + '\n')


def main():
    folderTrajs, atomId, resname = parseArguments()

    if atomId == "" and resname == "":
        sys.exit("Both resname or atomId should be provided")
    elif resname == "":
        resname = atomId.split(":")[-1] #the atom Id last element is the resname


    allFolders = os.listdir(folderTrajs)
    Epochs = [epoch for epoch in allFolders if epoch.isdigit()]
    for folder in Epochs:
        pathFolder = os.path.join(folderTrajs, folder)
        if not os.path.exists(DIRECTORY % pathFolder):
            os.makedirs(DIRECTORY % pathFolder)

        originalPDBfiles = glob.glob(pathFolder+'/*traj*.pdb')
        for filename in originalPDBfiles:
            print filename
            allCoordinates = loadAllAtomsInPdb(filename)
            # because of the way it's split, the last element is empty
            if atomId != "":
                coords = getAtomCoord(allCoordinates[:-1], resname, atomId)
            else:
                coords = getPDBCOM(allCoordinates[:-1], resname)

            outputFilename = getOutputFilename(DIRECTORY, filename,
                                               BASEOUTPUTFILENAME)
            writeToFile(coords, outputFilename % pathFolder)


if __name__ == "__main__":
    main()
