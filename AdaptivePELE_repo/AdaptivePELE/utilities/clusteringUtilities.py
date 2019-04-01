from __future__ import absolute_import, division, print_function, unicode_literals
import os
from AdaptivePELE.utilities import utilities


def writeStructures(clusteringObject, listStructures, checker=lambda x: True, outputPath="cluster.pdb"):
    """
        Print all clusters in listStructures that meet the condition specified
        by the checker

        :param clusteringObject: Clustering object with clusters to print
        :type clusteringObject: :py:class:`.Clustering`
        :param checker: Lambda function with the checker that should evaluate to True for intersted structures
        :type checker: function
        :param outputPath: Output cluster pdb filename
        :type outputPath: str
    """
    clObject = utilities.readClusteringObject(clusteringObject)
    nameStructure = os.path.splitext(outputPath)
    outputName = nameStructure[0]+'_%d'+nameStructure[1]
    path = os.path.split(outputName)
    pathToWrite = path[1]
    if path[0]:
        utilities.makeFolder(path[0])
        pathToWrite = os.path.join(path[0], path[1])

    if listStructures is None or len(listStructures) == 0:  # If no listStructures, write all
        listStructures = range(len(clObject.clusters.clusters))

    for element in listStructures:
        cluster = clObject.clusters.clusters[element]
        if checker is None or checker(cluster):
            print("Writing", pathToWrite % element)
            cluster.writePDB(pathToWrite % element)
