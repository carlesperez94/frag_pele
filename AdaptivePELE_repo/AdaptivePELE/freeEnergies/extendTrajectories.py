from __future__ import absolute_import, division, print_function, unicode_literals
import numpy as np
import glob
import os
import re
import sys
import argparse
import ast


def parseArguments():
    desc = "Program that extends trajectories.\n\
            Trajectories are joined with those from which the spawning cluster was discovered\n\
            Two options are available:\n\
                *) full: Trajectories are traced back to epoch 0\n\
                *) prev: Trajectories are extended with up to the last 'lagtime' snapshots of a previous trajectory\n"

    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--type", choices=['full', 'prev'], default='prev', help="Traj reconstruction type")
    parser.add_argument("-l", "--lagtime", default=0, type=int, help="Lagtime to be used (if proceeds)")
    parser.add_argument("--outputDir", default="allTrajs_reconstructed", help="Output directory")
    parser.add_argument("--inputDir", default="allTrajs", help="Directory with trajectories")
    args = parser.parse_args()

    return args.type, args.lagtime, args.outputDir, args.inputDir


def tryToOpenMapping(mapFilename):
    try:
        with open(mapFilename) as f:
            mapping = f.read().split(":")
            opened = True
    except:
        mapping = None
        opened = False

    return mapping, opened


def sameCoords(coords1, coords2):
    threshold = 1e-4
    diff = np.abs(coords1 - coords2)
    return (diff < threshold).all()


def checkFirstMatchingSnapshot(traj, snapshot, coords):
    for i in range(snapshot, traj.shape[0]):
        if sameCoords(coords[1:], traj[i, 1:]):  # coords[0] is the snapshot num, which is not important and is discarded in MSM
            firstMatchingSnapshot = i
            return firstMatchingSnapshot
    raise IndexError


def findSnapshotAndOpenTraj(trajName, lastSnapshot, coords, firstSnapshot=0):
    if lastSnapshot is None or coords is None:
        return np.loadtxt(trajName)[firstSnapshot:]
    else:
        traj = np.loadtxt(trajName)
        try:
            snapshot = checkFirstMatchingSnapshot(traj, lastSnapshot, coords)
            return traj[firstSnapshot:snapshot]
        except IndexError:
            sys.exit("Did not find matching traj in trajName: %s; coords:%s, from snapshot:%d" % (trajName, coords, lastSnapshot))


def reconstructFullTrajectory(mapping, thisTrajMap, trajNameTempletized, coords):
    """
        thisTrajMap contains the exact point at which a cluster was discovered
        Note that the number of snapshot corresponds to the accepted steps and
        not absolute steps. There are different ways to overcome the limitation.
        The fastest is looking at the report file. A slower way is looking at
        the exact coordinates. It is slower, but the main advantage is that we
        do not need any extra file.
    """
    (epoch, num, snapshot) = thisTrajMap

    try:
        thisTraj = findSnapshotAndOpenTraj(trajNameTempletized % (epoch, num), snapshot, coords)
    except:  # this is due to an error in adaptiveSampling. Once the bug is found, please remove the except block
        epoch += 1
        thisTraj = findSnapshotAndOpenTraj(trajNameTempletized % (epoch, num), snapshot, coords)

    if epoch == 0:
        return thisTraj
    else:
        prevTrajMap = ast.literal_eval(mapping[epoch][num-1])
        return np.vstack((reconstructFullTrajectory(mapping, prevTrajMap, trajNameTempletized, thisTraj[0]), thisTraj))


def addUncountedSnapshots(mapping, thisTrajMap, trajNameTempletized, coords, lagtime):
    """
        This function adds all possible previous uncounted snapshots
        (i.e. those in the last lagtime snapshots) to the current traj

        thisTrajMap contains the exact point at which a cluster was discovered
        Note that the number of snapshot corresponds to the accepted steps and
        not absolute steps. There are different ways to overcome the limitation.
        The fastest is looking at the report file. A slower way is looking at
        the exact coordinates. It is slower, but the main advantage is that we
        do not need any extra file.
    """

    (epoch, num, snapshot) = thisTrajMap
    thisTraj = findSnapshotAndOpenTraj(trajNameTempletized % (epoch, num), snapshot, None)

    if epoch == 0:
        return thisTraj

    prevTrajMap = ast.literal_eval(mapping[epoch][num-1])
    (epoch, num, snapshot) = prevTrajMap
    try:
        # only consider the last "lagtime" snapshots
        # if the initial point was found before the last lagtime snapshots, then: prevTraj = []
        prevTraj = findSnapshotAndOpenTraj(trajNameTempletized % (epoch, num), snapshot, thisTraj[0], firstSnapshot=-lagtime)
    except:
        epoch += 1
        prevTraj = findSnapshotAndOpenTraj(trajNameTempletized % (epoch, num), snapshot, thisTraj[0], firstSnapshot=-lagtime)

    return np.vstack((prevTraj, thisTraj))


def main():
    choice, lagtime, outputDir, inputDir = parseArguments()

    mapFilename = "%d/processorMapping.txt"

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    trajWildcard = "traj_%d_*.dat"  # %d for the epoch
    trajName = "traj_%d_%d.dat"  # %d for the epoch and number

    trajNameTempletized = os.path.join(inputDir, trajName)

    allFolders = os.listdir(".")
    epochFolders = [int(re.sub("MSM_", "", epoch)) for epoch in allFolders if epoch.startswith("MSM")]

    numberOfEpochs = max(epochFolders)
    mappings = []
    for epoch in range(0, numberOfEpochs):
        epochMapping, _ = tryToOpenMapping(mapFilename % epoch)
        mappings.append(epochMapping)

    newSizes = []
    for epoch in range(0, numberOfEpochs):
        allFiles = glob.glob(os.path.join(inputDir, trajWildcard % epoch))
        for source in allFiles:
            print(source)
            num = int(source.split("_")[-1][:-4])
            if choice == "full":
                fullTraj = reconstructFullTrajectory(mappings, (epoch, num, None), trajNameTempletized, None)
            elif choice == "prev":
                fullTraj = addUncountedSnapshots(mappings, (epoch, num, None), trajNameTempletized, None, lagtime)

            newSizes.append(fullTraj.shape[0])

            fname = os.path.split(source)[-1]
            dst = os.path.join(outputDir, fname)
            np.savetxt(dst, fullTraj)

    newSizes = np.array(newSizes)
    avgNewSize = np.average(newSizes)
    print("")
    print("Avg new size: %.2f +/- %.2f" % (avgNewSize, np.std(newSizes)))
    try:
        origSize = np.loadtxt(allFiles[0]).shape[0]
        print("Assuming orig trajectories of %d steps" % origSize)
        print("New trajectories are {0:.2f}% larger".format((avgNewSize/origSize - 1)*100))
    except:
        pass

if __name__ == "__main__":
    main()
