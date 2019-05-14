from __future__ import absolute_import, division, print_function, unicode_literals
import math
import os
import numpy
import argparse


def parseArguments():
    desc = "Program that analyses a column of data, printing different statistical values and a histogram if desired."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("col", type=int, default=1, help="Column with the values (Integer values in [1, N]).")
    parser.add_argument("files", nargs='+', help="Files with the data to analyse.")
    parser.add_argument("-o", metavar="OUTLIERRANGE", type=float, default=1, help="Factor that widens the range for a value not to be considered as an outlier.")
    parser.add_argument("-p", nargs='?', type=int, const=10, metavar='bins', help="Plots histogram. Can take as argument the number of bins.")
    parser.add_argument("-e", action="store_true", default=False, help="If plotting, draw the error bars.")
    args = parser.parse_args()

    return args.col, args.files, args.p, args.o, args.e


def analyseData(data):
    numberOfElements = len(data)
    average = sum(data)
    average /= numberOfElements

    data.sort()
    minVal = data[0]
    maxVal = data[-1]
    median = data[numberOfElements//2]

    variance = 0
    for number in data:
        variance += math.pow(number - average, 2)

    variance = variance/(numberOfElements-1)

    standardDeviation = math.sqrt(variance)

    print("------------------------------")
    print("RESULTS")
    print("------------------------------")
    print("# of elements: ", numberOfElements)
    print("Average value: ", average)
    print("Std deviation: ", standardDeviation)
    print("Maximum vaule: ", maxVal)
    print("Minumum value: ", minVal)
    print("Median  value: ", median)
    print("------------------------------")


def returnListWithoutOutliers(data, outlierRange):
    """
        An outlier is defiend as a datapoint not in [Q1 - 1.5*IQR*outlierRange, Q3 + 1.5*IQR*outlierRange],
        where IQR is the interquartile range: Q3 - Q1
    """
    data.sort()

    dataPointsBefore = len(data)
    Q1 = data[dataPointsBefore//4]
    Q3 = data[3*dataPointsBefore//4]
    IQR = Q3 - Q1

    lowerFence = Q1 - 1.5 * IQR * outlierRange
    upperFence = Q3 + 1.5 * IQR * outlierRange

    filteredData = [i for i in data if i >= lowerFence and i <= upperFence]

    dataPointsAfter = len(filteredData)
    print('Removed ' + str(dataPointsBefore - dataPointsAfter) + ' outliers')

    return filteredData


def printFilenamesAbsolutePath(files):
    print("files:")
    for f in files:
        print(os.path.abspath(f))


def readDataFromFiles(files, column):
    data = []
    for f in files:
        tmpdata = numpy.loadtxt(f, unpack=True, usecols=[column])
        try:
            data.extend(tmpdata)
        except:
            data.append(tmpdata)
    return data


def main():
    column, files, bins, outlierRange, plotErrorBars = parseArguments()
    # We count columns starting with 1
    column -= 1

    data = readDataFromFiles(files, column)

    printFilenamesAbsolutePath(files)
    print("---------------------------")
    print("Before removing outliers:")

    analyseData(data)

    filteredData = returnListWithoutOutliers(data, outlierRange)

    print("---------------------------")
    print("After removing outliers:")

    analyseData(filteredData)

    if bins:
        import histogram
        plotErrorBars = plotErrorBars not in ['False', 'false', 0]
        histogram.plot_histogram(data, bins, plotErrorBars, 1, "Before removing outliers")
        histogram.plot_histogram(filteredData, bins, plotErrorBars, 2, "After removing outliers (%.2f * IQR)" % (outlierRange*1.5))
        # To keep plots "alive"
        raw_input()


if __name__ == '__main__':
    main()
