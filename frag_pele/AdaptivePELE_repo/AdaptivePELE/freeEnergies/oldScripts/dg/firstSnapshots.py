import glob
import os
import argparse
import numpy as np

def parseArguments():
    desc = "Copies first n steps of trajectories in data folder to current folder"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-s", "--steps", type=int, help="Number of steps")
    #Not ready for bootstrap yet, need to copy files to non-clashing destination
    parser.add_argument("-b", nargs='?', type=int, const=-1, metavar='bootstrap', help="Bootstrap trajectories from data folder")
    args = parser.parse_args()

    return  args.steps, args.b


def main(steps, bootstrap):
    dataFolder = "rawData"
    fileWildcard = "traj_*.dat"
    fullFileWildcard = os.path.join(dataFolder, fileWildcard)

    allFiles = glob.glob(fullFileWildcard)

    if not bootstrap is None:
        if bootstrap == -1:
            numberOfFiles = len(allFiles)
        else:
            numberOfFiles = bootstrap#len(allFiles)
        trajFiles = np.random.choice(allFiles, numberOfFiles)
        outputFilename = "traj_.%d.dat"
    else:
        trajFiles = allFiles

    for i,trajFile in enumerate(trajFiles):
        traj = np.loadtxt(trajFile)
        if bootstrap:
            outputPath = outputFilename%i
        else:
            outputPath = os.path.split(trajFile)[-1]
        try:
            np.savetxt(outputPath, traj[:steps+1,:], fmt="%d\t%.4f\t%.4f\t%.4f")
        except:
            import sys
            sys.exit("There is a problem with %s"%trajFile)

if __name__ == "__main__":
    steps, bootstrap = parseArguments()
    main(steps, bootstrap)
