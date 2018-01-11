import os
import argparse
import glob
import types



def parseArguments():
    desc = "Generates output for gnuplot\n"\
            "It MUST be run from the root epoch folder (i.e., where it can find the folders 0/, 1/, 2/, ... lastEpoch/"\
            "To be run for example like: \">python generateGnuplotFile.py | gnuplot -persist\""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("steps", type=int, default=25, help="Pele steps per run") 
    parser.add_argument("xcol", type=int, default=2, help="xcol") 
    parser.add_argument("ycol", type=int, default=4, help="ycol") 
    parser.add_argument("filename", type=str, default="report_", help="Report filename") 
    parser.add_argument("-be", action="store_true", help="Points")
    parser.add_argument("-rmsd", action="store_true", help="Lines")
    args = parser.parse_args()
    return args.steps, args.xcol, args.ycol, args.filename, args.be, args.rmsd

def generateNestedString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps=False, replotFirst=False):
    """
        TotalNumberOfSteps -> not only considering steps in current epoch, but steps with all previous epochs
    """
    allFolders = os.listdir('.')
    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
    numberOfEpochs=int(len(epochFolders))


    dictionary = {'reportName':reportName, 'col2':column2, 'numberOfEpochs':numberOfEpochs, 'withLines':''}
    
    #runs of epoch 0, assumed constant
    numberOfRunsPerEpoch=len( glob.glob( os.path.join(str(0), reportName+"*") ) )
    dictionary['runsPerEpoch'] = numberOfRunsPerEpoch


    if printWithLines:
        dictionary['withLines'] = "w l"

    if type(column1) == types.IntType:
        if totalNumberOfSteps:
            dictionary['col1'] = "($" + str(column1) + "+ (%d*j))"%stepsPerRun #adds steps per runs, so that it mathes the total number of steps
        else:
            dictionary['col1'] = str(column1)
    elif type(column1) == types.LambdaType:
        dictionary['col1'] = "(" + str(column1(i)) + ")"

    return gnuplotString % dictionary + "\n"


def generateForLoopString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps=False, replotFirst=False):
    """
        TotalNumberOfSteps -> not only considering steps in current epoch, but steps with all previous epochs
    """
    allFolders = os.listdir('.')
    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
    numberOfEpochs=int(len(epochFolders))


    dictionary = {'reportName':reportName, 'col2':column2, 'numberOfEpochs':numberOfEpochs, 'withLines':''}
    plottingString = ""
    for i in range(numberOfEpochs):
        dictionary['epoch'] = i

        numberOfRunsPerEpoch=len( glob.glob( os.path.join(str(i), reportName+"*") ) )
        dictionary['runsPerEpoch'] = numberOfRunsPerEpoch

        if printWithLines:
            dictionary['withLines'] = "w l"

        if type(column1) == types.IntType:
            if totalNumberOfSteps:
                dictionary['col1'] = "($" + str(column1) + "+" + str(stepsPerRun*i) + ")" #adds steps per runs, so that it mathes the total number of steps
            else:
                dictionary['col1'] = str(column1)
        elif type(column1) == types.LambdaType:
            dictionary['col1'] = "(" + str(column1(i)) + ")"

        if i != 0 or replotFirst:
            plottingString += "re"

        plottingString += gnuplotString % dictionary + "\n"
        #plottingString += "pause 0.25\n"

    print plottingString



def generatePrintString(stepsPerRun, xcol, ycol, reportName, kindOfPrint):

    if kindOfPrint == "PRINT_RMSD_STEPS":
        printWithLines = True
        totalNumberOfSteps=True
    elif kindOfPrint == "PRINT_BE_RMSD":
        printWithLines = False
        totalNumberOfSteps=False

    representativeReportName='clustering/reports'

    gnuplotString = "plot for [j=0:%(numberOfEpochs)d-1] for [i=1:%(runsPerEpoch)d] \"\".j.\"/%(reportName)s\".i u %(col1)s:%(col2)d lt 6 lc palette frac j/%(numberOfEpochs)d. notitle %(withLines)s"
    return generateNestedString(gnuplotString, reportName, xcol, ycol, stepsPerRun, printWithLines, totalNumberOfSteps, False)


if __name__ == "__main__":
    stepsPerRun, xcol, ycol, filename, be, rmsd = parseArguments()
    #VARIABLES TO SET WHEN PRINTING
    if be:
        kindOfPrint = "PRINT_BE_RMSD"
    elif rmsd:
        kindOfPrint = "PRINT_RMSD_STEPS"

    print generatePrintString(stepsPerRun, xcol, ycol, filename, kindOfPrint)
