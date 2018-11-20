import sys
from AdaptivePELE.clustering import clustering, thresholdcalculator
from AdaptivePELE.spawning import spawning, densitycalculator
from AdaptivePELE.constants import constants
from AdaptivePELE.utilities import utilities
import pandas as pd
import glob
import os

def cluster_traject(resname, trajToDistribute, columnToChoose, distance_contact, clusterThreshold, path_to_cluster,
                    output_path, epsilon=0.5, report_basename="report", condition ="min", metricweights="linear",
                    nclusters=5):

    outputPathConst = constants.OutputPathConstants(output_path)
    outputPathConst.tmpFolder = output_path
    outputPathConst.buildTmpFolderConstants(outputPathConst.tmpFolder)
    utilities.makeFolder(outputPathConst.tmpFolder)

    thresholdCalc = thresholdcalculator.ThresholdCalculatorConstant(value=clusterThreshold)
    similarityEval = clustering.CMSimilarityEvaluator("Jaccard")
    clusteringObject = clustering.ContactMapAccumulativeClustering(thresholdCalc, similarityEval, resname=resname,
                                                                   reportBaseFilename=report_basename,
                                                                   columnOfReportFile=columnToChoose,
                                                                   contactThresholdDistance=distance_contact,
                                                                   altSelection=True)

    clusteringObject.cluster([path_to_cluster], ignoreFirstRow=True)
    spawning_params = spawning.SpawningParams()
    spawning_params.reportFilename = report_basename
    spawning_params.epsilon = epsilon
    spawning_params.nclusters = nclusters
    spawning_params.metricWeights = metricweights
    spawning_params.condition = condition

    density = densitycalculator.NullDensityCalculator()
    spawningObject = spawning.EpsilonDegeneracyCalculator(density)
    degeneracy = spawningObject.calculate(clusteringObject.clusters, trajToDistribute, spawning_params)
    spawningObject.log()

    spawningObject.writeSpawningInitialStructures(outputPathConst, degeneracy, clusteringObject, 0)


def get_column_num(path, header_column, report_basename="report"):
    reports = glob.glob(os.path.join(path, "*{}*".format(report_basename)))
    try:
        reports[0]
    except IndexError:
        raise IndexError("Not report file found. Check you are in adaptive's or Pele root folder")

    data = pd.read_csv(reports[0], sep='    ', engine='python')
    header_list = data.columns.values.tolist()
    column_number = header_list.index(header_column)
    return column_number

