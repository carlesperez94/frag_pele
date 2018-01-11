import glob
import os
import shutil

class Constants():
    def __init__(self):
        self.trajFolder = "allTrajs"
        self.origTrajFiles = os.path.join(self.trajFolder, "traj_*")
        self.trajFileEpoch = os.path.join(self.trajFolder, "traj_%d_*")
        self.msmFolder = "MSM_%d"
        self.rawDataFolder = os.path.join(self.msmFolder, "rawData")
        self.templetizedControlFileMSM = "templetized_control_MSM.conf"

def extractEpoch(f):
    first = f.find("_")
    second = f.rfind("_")
    epoch = f[first+1:second]
    return epoch

def getAllDifferentEpochs(origTrajFiles):
    trajFiles = glob.glob(origTrajFiles)
    epochs = set([])
    for f in trajFiles:
        epoch = extractEpoch(f)
        epochs.add(int(epoch))
    epochs = list(epochs)
    epochs.sort()
    return epochs

def makeMSMFolders(epochs, msmFolder):
    for epoch in epochs:
        if not os.path.exists(msmFolder%epoch):
            os.makedirs(msmFolder%epoch)

def makeRawDataFolders(epochs, rawDataFolder):
    """
        This folder contains symbolic links to the corresponding trajectories
    """
    for epoch in epochs:
        folder = rawDataFolder%epoch
        if not os.path.exists(folder):
            os.makedirs(folder)

def makeSymbolicLinks(epochs, rawDataFolder, trajFileEpoch):
    for epoch in epochs:
        destFolder = rawDataFolder%epoch
        for prevEpoch in range(epoch+1):
            sourcesPrevEpoch = glob.glob(trajFileEpoch%prevEpoch)
            for src in sourcesPrevEpoch:
                [srcFolder, filename] = os.path.split(src)
                srcFolder = os.path.abspath(srcFolder)
                src = os.path.join(srcFolder, filename)
                dest = os.path.join(destFolder, filename)
                try:
                    if not os.path.isfile(dest): os.symlink(src, dest)
                except OSError:
                    pass

def copyMSMcontrolFile(epochs, msmFolder, templetizedControlFileMSM):
    scriptsFolder = os.path.dirname(os.path.realpath(__file__))
    scriptsFile = os.path.join(scriptsFolder, templetizedControlFileMSM)
    print scriptsFile
    for epoch in epochs:
        dst = os.path.join(msmFolder%epoch, templetizedControlFileMSM) 
        shutil.copyfile(scriptsFile, dst)

def main():
    constants = Constants()

    epochs = getAllDifferentEpochs(constants.origTrajFiles)

    makeMSMFolders(epochs, constants.msmFolder)
    makeRawDataFolders(epochs, constants.rawDataFolder)
    makeSymbolicLinks(epochs, constants.rawDataFolder, constants.trajFileEpoch)
    #copyMSMcontrolFile(epochs, constants.msmFolder, constants.templetizedControlFileMSM)


if __name__ == "__main__":
    main()
