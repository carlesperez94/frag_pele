import os
import pickle
import utilities


def writeStructures(clusteringObject, listStructures, checker=lambda x: True, outputPath="cluster.pdb"):
    """
        Function that prints all clusters in listStructures so that checker evaluates to true

        clusteringObject [In] Clustering object with clusters to print
        checker [In] Lambda function with the checker that should evaluate to True for intersted structures
        outputPath [In] Output cluster pdb filename
    """
    with open(clusteringObject, "rb") as f:
        clObject = pickle.load(f)
    nameStructure = os.path.splitext(outputPath)
    outputName = nameStructure[0]+'_%d'+nameStructure[1]
    path = os.path.split(outputName)
    pathToWrite = path[1]
    if path[0]:
        utilities.makeFolder(path[0])
        pathToWrite = os.path.join(path[0], path[1])

    if listStructures is None or len(listStructures) == 0: #If no listStructures, write all
        listStructures = range(len(clObject.clusters.clusters))

    output = ""
    for element in listStructures:
        cluster = clObject.clusters.clusters[element]
        if checker is None or checker(cluster):
            print "Writing", pathToWrite%element
            cluster.pdb.pdb += "\nENDMDL\n"
            output += cluster.pdb.pdb
            cluster.writePDB(pathToWrite % element)

