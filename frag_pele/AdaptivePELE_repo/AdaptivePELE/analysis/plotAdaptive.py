from __future__ import absolute_import, division, print_function, unicode_literals
import os
import argparse
import glob


def parseArguments():
    """
        Parse command line arguments

        :returns: int, int, int, str, bool, bool -- Number of steps per epoch,
            column to plot in the X axis, column to plot in the Y axis, name of
            the files containing the simulation data, whether to plot the data
            as points, wether to plot the data as lines
    """
    desc = "Generates output for gnuplot\n"\
           "It MUST be run from the root epoch folder (i.e., where it can find the folders 0/, 1/, 2/, ... lastEpoch/"\
           "To be run for example like: \">python generateGnuplotFile.py | gnuplot -persist\""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("steps", type=int, default=4, help="Pele steps per run")
    parser.add_argument("xcol", type=int, default=2, help="xcol")
    parser.add_argument("ycol", type=int, default=4, help="ycol")
    parser.add_argument("filename", type=str, default="report_", help="Report filename")
    parser.add_argument("-points", action="store_true", help="Plot using points")
    parser.add_argument("-lines", action="store_true", help="Plot using lines")
    parser.add_argument("-zcol", type=int, default=None, help="Column to define color according to metric")
    parser.add_argument("-t", "--traj_range", type=str, default=None, help="Range of trajs to select, e.g to select trajs from 1 to 10, 1:10")

    args = parser.parse_args()
    return args.steps, args.xcol, args.ycol, args.filename, args.points, args.lines, args.zcol, args.traj_range


def generateNestedString(gnuplotString, reportName, column1, column2, stepsPerRun, printWithLines, totalNumberOfSteps=False, replotFirst=False, paletteModifier=None, trajs_range=None):
    """
        Generate a string to be passed to gnuplot

        :param gnuplotString: Template string for gnuplot
        :type gnuplotString: str
        :param reportName: Name of the files containing the simulation data
        :type reportName: str
        :param column1: Column to plot in the X axis
        :type column1: int
        :param column2: Column to plot in the Y axis
        :type column2: int
        :param stepsPerRun: Number of steps per epoch,
        :type stepsPerRun: int
        :param TotalNumberOfSteps: Not only considering steps in current epoch,
            but steps with all previous epochs
        :type TotalNumberOfSteps: bool
        :param replotFirst: Deprecated parameter
        :type replotFirst: bool
        :param paletteModifier: Wheter to use the epoch as color or a column
        :type paletteModifier: int

        :returns: str -- String to plot using gnuplot
    """
    allFolders = os.listdir('.')
    epochFolders = [epoch for epoch in allFolders if epoch.isdigit()]
    numberOfEpochs = int(len(epochFolders))

    dictionary = {'reportName': reportName, 'col2': column2, 'numberOfEpochs': numberOfEpochs, 'withLines': ''}

    if trajs_range is not None:
        start, end = map(int, traj_range.split(":"))
        dictionary['startTraj'] = start
        dictionary['runsPerEpoch'] = end
    else:
        dictionary['startTraj'] = 1
        # runs of epoch 0, assumed constant
        numberOfRunsPerEpoch = len(glob.glob(os.path.join(str(0), reportName+"*")))
        dictionary['runsPerEpoch'] = numberOfRunsPerEpoch

    if printWithLines:
        dictionary['withLines'] = "w l"

    if isinstance(column1, int):
        if totalNumberOfSteps:
            dictionary['col1'] = "($" + str(column1) + "+ (%d*j))" % stepsPerRun  # adds steps per runs, so that it mathes the total number of steps
        else:
            dictionary['col1'] = str(column1)

    return gnuplotString % dictionary + "\n"


def generatePrintString(stepsPerRun, xcol, ycol, reportName, kindOfPrint, paletteModifier, trajs_range):
    """
        Generate a template string to use with gnuplot

        :param stepsPerRun: Number of steps per epoch,
        :type stepsPerRun: int
        :param xcol: Column to plot in the X axis
        :type xcol: int
        :param ycol: Column to plot in the Y axis
        :type ycol: int
        :param reportName: Name of the files containing the simulation data
        :type reportName: str
        :param kindOfPrint:  Kind of lines to plot (solid or points)
        :type kindOfPrint: bool
        :param paletteModifier: Third column to specify color
        :type paletteModifier: int
        :trajs_range: Range of trajectories to plot
        :type trajs_range: str

        :returns: str -- String to plot using gnuplot
    """
    if kindOfPrint == "PRINT_RMSD_STEPS":
        printWithLines = True
        totalNumberOfSteps = True
    elif kindOfPrint == "PRINT_BE_RMSD":
        printWithLines = False
        totalNumberOfSteps = False
    if paletteModifier is None:
        stringPalette = "frac j/%(numberOfEpochs)d. "
        colorMetric = ""
    else:
        stringPalette = ""
        colorMetric = ":%d" % paletteModifier

    gnuplotString = "".join(["plot for [i=%(startTraj)d:%(runsPerEpoch)d] for [j=0:%(numberOfEpochs)d-1] \'\'.j.\'/%(reportName)s\'.i u %(col1)s:%(col2)d", colorMetric, " lt 6 lc palette ", stringPalette, "notitle %(withLines)s"])
    return generateNestedString(gnuplotString, reportName, xcol, ycol, stepsPerRun, printWithLines, totalNumberOfSteps, False, paletteModifier, trajs_range)


if __name__ == "__main__":
    steps_Run, Xcol, Ycol, filename, be, rmsd, colModifier, traj_range = parseArguments()
    # VARIABLES TO SET WHEN PRINTING
    if be:
        kind_Print = "PRINT_BE_RMSD"
    elif rmsd:
        kind_Print = "PRINT_RMSD_STEPS"

    printLine = generatePrintString(steps_Run, Xcol, Ycol, filename, kind_Print, colModifier, traj_range)
    print(printLine)
