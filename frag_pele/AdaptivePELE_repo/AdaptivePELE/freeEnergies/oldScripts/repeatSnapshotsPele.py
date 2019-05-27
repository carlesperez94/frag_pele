import os
import numpy
import glob
import re
import argparse


def parseArguments():
    desc = "Extracts coordinates in <currentFolder>/extractedCoordinates/coord*"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("folderTraj",
                        help="Folder where the epochs trajs are stored")
    args = parser.parse_args()
    return args.folderTraj

folderTraj = parseArguments()
inputTrajectoryFolder = "extractedCoordinates/"
baseTrajectoryName = "coord_"
reportName = '*report_'

outputTrajectoryFolder = "repeatedExtractedCoordinates/"
if not os.path.exists(outputTrajectoryFolder):
    os.makedirs(outputTrajectoryFolder)

inputTrajectories = glob.glob(inputTrajectoryFolder + baseTrajectoryName + '*')

for inputTrajectory in inputTrajectories:
    trajectoryNumber = re.sub('\.dat$','',inputTrajectory)
    trajectoryNumber = re.sub(inputTrajectoryFolder + baseTrajectoryName,'',trajectoryNumber)
    try:
        reportFile = glob.glob(reportName + trajectoryNumber)[0]
    except:
        print "Couldn't find file that matches: ", reportName + trajectoryNumber
        continue


    with open(inputTrajectory) as f:
        trajectory = f.read().splitlines()

    steps = numpy.loadtxt(reportFile, dtype='int', comments='#', usecols=(1,))

    #Improvement: Whenever the initial step is added in the trajectory, we can add steps[0] times
    # the initial structure
    completeTrajectory = []
    counter=0
    """
    if len(trajectory) == 0: #one step traj, we could repeat #steps/epoch times
        snapshot = trajectory[i].split()
        snapshot[0] = str(counter)
        snapshot = ' '.join(snapshot)
        completeTrajectory.append(snapshot)
    """
    if len(trajectory) > 0:
        for i in range(len(trajectory) - 1):
            try:
                repeated = steps[i+1] - steps[i]
            except:
                #sys.exit("sth wrong in trajectory " + inputTrajectory)
                print "sth wrong in trajectory " + inputTrajectory
                #continue

                #Changed behavior, write until the end of the information in report file
                snapshot = trajectory[i].split()
                snapshot[0] = str(counter)
                snapshot = ' '.join(snapshot)
                completeTrajectory.append(snapshot)

                break
            for j in range(repeated):
                snapshot = trajectory[i].split()
                snapshot[0] = str(counter)
                snapshot = ' '.join(snapshot)
                completeTrajectory.append(snapshot)
                counter += 1

        #!!!!
        #WARNING!!! Add last snapshot DID NOT CHECK when report/traj don't match
        #!!!!
        snapshot = trajectory[-1].split()
        snapshot[0] = str(counter)
        snapshot = ' '.join(snapshot)
        completeTrajectory.append(snapshot)

        outputFilename = outputTrajectoryFolder + baseTrajectoryName + trajectoryNumber + '.dat'
        outputFile = open(outputFilename, 'w')
        for snapshot in completeTrajectory:
            outputFile.write("%s\n" % snapshot)
        outputFile.close()
