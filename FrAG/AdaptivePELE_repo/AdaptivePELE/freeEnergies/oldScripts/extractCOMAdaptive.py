import os
import argparse
import glob
# import minimumGCluster.atomset as atomset
import AdaptivePELE.atomset.atomset as atomset

DIRECTORY = '%s/extractedCoordinates'
BASEOUTPUTFILENAME = 'coord'


def loadAllAtomsInPdb(filename):
    fileContent = open(filename).read()
    fileContent = fileContent.split('ENDMDL')
    return fileContent


def parseArguments():
    desc = "Extracts coordinates in <currentFolder>/extractedCoordinates/coord*"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("folderTraj",
                        help="Folder where the epochs trajs are stored")
    parser.add_argument("resname", help="Name of the ligand in the pdb")
    # parser.add_argument("-f", nargs='+', help="Files to get coordinates")
    args = parser.parse_args()
    return args.folderTraj, args.resname


def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('_')

    return filename[first+1:last]


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


def writeToFile(COMs, outputFilename):
    with open(outputFilename, 'w') as f:
        for i, line in enumerate(COMs):
            f.write(str(i) + ' ')
            for i in range(len(line) - 1):
                f.write(str(line[i]) + ' ')
            f.write(str(line[-1]) + '\n')


def main():
    folderTrajs, resname = parseArguments()

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
            COMs = getPDBCOM(allCoordinates[:-1], resname)

            outputFilename = getOutputFilename(DIRECTORY, filename,
                                               BASEOUTPUTFILENAME)
            writeToFile(COMs, outputFilename % pathFolder)


if __name__ == "__main__":
    main()
