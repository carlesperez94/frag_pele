from __future__ import absolute_import, division, print_function, unicode_literals
import glob
import os
import shutil


class Constants():
    def __init__(self, trajs_path=None):
        self.trajFolder = "allTrajs"
        if trajs_path is not None:
            self.trajFolder = os.path.join(trajs_path, self.trajFolder)
        self.origTrajFiles = os.path.join(self.trajFolder, "traj_*")
        self.trajFileEpoch = os.path.join(self.trajFolder, "traj_%d_*")
        self.nonRepeatedTrajEpoch = os.path.join(self.trajFolder, "extractedCoordinates", "traj_%d_*")
        self.msmFolder = "MSM_%d"
        self.rawDataFolder = os.path.join(self.msmFolder, "rawData")
        self.nonRepeatedCoordsFolder = os.path.join(self.rawDataFolder, "extractedCoordinates")
        self.templetizedControlFileMSM = "templetized_control_MSM.conf"


def extractEpoch(f):
    # Split the filename blablb_0_1.dat into [balblb, 0, 1.dat]
    return f.rsplit("_", 2)[1]


def getAllDifferentEpochs(origTrajFiles):
    trajFiles = glob.glob(origTrajFiles)
    epochs = set([])
    for f in trajFiles:
        epoch = extractEpoch(f)
        epochs.add(int(epoch))
    epochs = sorted(epochs)
    return epochs


def makeMSMFolders(epochs, msmFolder):
    for epoch in epochs:
        if not os.path.exists(msmFolder % epoch):
            os.makedirs(msmFolder % epoch)


def makeRawDataFolders(epochs, rawDataFolder, nonRepeatedCoordsFolder):
    """
        This folder contains symbolic links to the corresponding trajectories
    """
    for epoch in epochs:
        folder = rawDataFolder % epoch
        nonRepeatedFolder = nonRepeatedCoordsFolder % epoch
        if not os.path.exists(folder):
            os.makedirs(folder)
        if not os.path.exists(nonRepeatedFolder):
            os.makedirs(nonRepeatedFolder)


def makeSymbolicLinksForFiles(filelist, destination):
    for src in filelist:
        [srcFolder, filename] = os.path.split(src)
        # srcFolder = os.path.abspath(srcFolder)
        # Switch the symbolic links to use a relative path, so if the
        # folders containing the data are moved they will not break
        srcFolder = os.path.relpath(srcFolder, start=destination)
        src = os.path.join(srcFolder, filename)
        dest = os.path.join(destination, filename)
        try:
            if not os.path.isfile(dest):
                os.symlink(src, dest)
        except OSError:
            pass


def makeSymbolicLinks(epochs, rawDataFolder, trajFileEpoch, trajNonRepeatedEpoch, nonRepeatedRawData):
    for epoch in epochs:
        destFolder = rawDataFolder % epoch
        destNonRepeatedFolder = nonRepeatedRawData % epoch
        for prevEpoch in range(epoch+1):
            sourcesPrevEpoch = glob.glob(trajFileEpoch % prevEpoch)
            makeSymbolicLinksForFiles(sourcesPrevEpoch, destFolder)
            nonRepeatedPrevEpoch = glob.glob(trajNonRepeatedEpoch % prevEpoch)
            makeSymbolicLinksForFiles(nonRepeatedPrevEpoch, destNonRepeatedFolder)


def copyMSMcontrolFile(epochs, msmFolder, templetizedControlFileMSM):
    scriptsFolder = os.path.dirname(os.path.realpath(__file__))
    scriptsFile = os.path.join(scriptsFolder, templetizedControlFileMSM)
    print(scriptsFile)
    for epoch in epochs:
        dst = os.path.join(msmFolder % epoch, templetizedControlFileMSM)
        shutil.copyfile(scriptsFile, dst)


def main(trajsPath=None):
    constants = Constants(trajsPath)

    epochs = getAllDifferentEpochs(constants.origTrajFiles)

    makeMSMFolders(epochs, constants.msmFolder)
    makeRawDataFolders(epochs, constants.rawDataFolder, constants.nonRepeatedCoordsFolder)
    makeSymbolicLinks(epochs, constants.rawDataFolder, constants.trajFileEpoch, constants.nonRepeatedTrajEpoch, constants.nonRepeatedCoordsFolder)
    # copyMSMcontrolFile(epochs, constants.msmFolder, constants.templetizedControlFileMSM)


if __name__ == "__main__":
    main()
