import os
import shutil
import glob
import re

if not os.path.exists("allTrajs"):
    os.makedirs("allTrajs")

outputTrajName = "allTrajs/traj_%s_%s.dat"
origTrajectories = "%s/repeatedExtractedCoordinates/coord_"
allFolders = os.listdir(".")
Epochs = [epoch for epoch in allFolders if epoch.isdigit()]
for folder in Epochs:
    trajectories = glob.glob(origTrajectories%folder + "*")
    for inputTrajectory in trajectories:
        trajectoryNumber = re.sub('\.dat$', '', inputTrajectory)
        trajectoryNumber = re.sub(origTrajectories%folder, '', trajectoryNumber)
        shutil.copyfile(inputTrajectory, outputTrajName%(folder, trajectoryNumber))
