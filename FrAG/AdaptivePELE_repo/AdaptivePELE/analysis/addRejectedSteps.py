"""
    Adds rejected steps in report

    To be run in simulation root folder

    Needs refactoring
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob

# to be changed by user
reportFilename = "report_"
reportOutputFilename = "allReport_"
stepsPerEpoch = 4

allFolders = os.listdir(".")
epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]

for folder in epochFolders:
    print("Writing reports in folder", folder)
    os.chdir(folder)

    reports = glob.glob(reportFilename + "*")
    for report in reports:
        with open(report, "r") as f:
            newReportContent = []
            reportContent = f.readlines()

        newReportContent.append(reportContent[0])
        count = 0
        for i in range(1, len(reportContent)-1):  # the last one will be checked against steps per epoch
            splitContent = reportContent[i].split()
            thisAccStep = int(float(splitContent[1]))
            nextAccStep = int(float(reportContent[i+1].split()[1]))
            for j in range(nextAccStep-thisAccStep):
                splitContent[1] = str(count)
                newReportContent.append(" ".join(splitContent))
                count += 1

        # last one is always the steps per epoch (+1, since 0 is initial)
        splitContent = reportContent[-1].split()
        lastAccStep = int(float(splitContent[1]))
        lastTotalStep = stepsPerEpoch + 1
        for i in range(lastTotalStep - lastAccStep):
            splitContent[1] = str(count)
            newReportContent.append(" ".join(splitContent))
            count += 1

        # write file
        reportNum = report.split("_")[-1]
        with open(reportOutputFilename + reportNum, "w") as f:
            for line in newReportContent:
                f.write(line + "\n")

    os.chdir("..")
