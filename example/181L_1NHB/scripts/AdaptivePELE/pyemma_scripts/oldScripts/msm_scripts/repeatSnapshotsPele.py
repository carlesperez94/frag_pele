import os
import numpy
import glob
import re
import sys

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

    #Whenever the initial step is added in the trajectory, we can add steps[0] times
    # the initial structure
    completeTrajectory = []
    counter=0
    for i in range(len(trajectory) - 1):
        try:
            repeated = steps[i+1] - steps[i]
        except:
            #sys.exit("sth wrong in trajectory " + inputTrajectory)
            print "sth wrong in trajectory " + inputTrajectory
            continue
        for j in range(repeated):
            snapshot = trajectory[i].split()
            snapshot[0] = str(counter)
            snapshot = ' '.join(snapshot)
            completeTrajectory.append(snapshot)
            counter += 1

    outputFilename = outputTrajectoryFolder + baseTrajectoryName + trajectoryNumber + '.dat'
    outputFile = open(outputFilename, 'w')
    for snapshot in completeTrajectory:
        outputFile.write("%s\n" % snapshot)
    outputFile.close()
