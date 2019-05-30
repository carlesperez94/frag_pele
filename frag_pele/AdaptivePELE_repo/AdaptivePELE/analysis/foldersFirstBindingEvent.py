from __future__ import absolute_import, division, print_function, unicode_literals
import argparse
import os
from AdaptivePELE.analysis import findfirstbindingevent
from AdaptivePELE.analysis import analyse


def parseArguments():
    """
        Parse the command-line options

        :returns: list, int, float, int, bool -- List of folders, column with
            binding event related metric, threshold for a binding event to be
            considered, number of steps per epoch to be consdidered, wether the
            simulation to analyse is and adaptive or sequential simulation
    """
    desc = "Program that computes the first binding event for a series of adaptive sampling runs"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("column", type=int, help="Column with binding event related metric")
    parser.add_argument("threshold", type=float, help="Threshold for a binding event to be considered")
    parser.add_argument("stepsPerEpoch", type=int, help="StepsPerEpoch")
    parser.add_argument("-seq", action='store_true', help="Use a sequential run, instead of adaptive")
    parser.add_argument("-u", action='store_true', help="Look for unbinding event, instead of binding")
    parser.add_argument("folders", nargs='+', default=".", help="Folders with adaptive sampling runs")
    args = parser.parse_args()

    return args.folders, args.column, args.threshold, args.stepsPerEpoch, args.seq, args.u


def main(folders, column, threshold, stepsPerEpoch, sequential, unbinding):
    """
        Calculate first binding event statistics (mean, median, std)

        :param folders: List of folders
        :type folders: list
        :param column: Column with binding event related metric
        :type column: int
        :param threshold: Threshold for a binding event to be considered
        :type threshold: float
        :param stepsPerEpoch: Number of steps per epoch to be consdidered
        :type stepsPerEpoch: int
        :param sequential: Whether the simulation to analyse is and adaptive or
            sequential simulation
        :type sequential: bool
    """
    print_label = "binding"
    if unbinding:
        print_label = "unbinding"
    firstBE = []
    for folder in folders:
        cwd = os.getcwd()
        os.chdir(folder)
        if sequential:
            stepsFirstBindingEvent = findfirstbindingevent.findEpochFirstBindingEvent(threshold, column, unBinding=unbinding)
        else:
            stepsFirstBindingEvent = findfirstbindingevent.findFirstBindingEvent(stepsPerEpoch, column, threshold, unBinding=unbinding)

        if stepsFirstBindingEvent is None:
            print("Didn't see any first %s event in folder:" % print_label, folder)
        else:
            firstBE.append(stepsFirstBindingEvent)

        os.chdir(cwd)

    if len(firstBE) > 1:
        print(firstBE)
        analyse.analyseData(firstBE)
    elif len(firstBE) == 1:
        print(firstBE[0])
    else:
        print("No %s event found" % print_label)


if __name__ == "__main__":
    folders_name, col, thresholds, steps_epoch, seq, unbind = parseArguments()
    # We count columns starting by 1
    col -= 1
    main(folders_name, col, thresholds, steps_epoch, seq, unbind)
