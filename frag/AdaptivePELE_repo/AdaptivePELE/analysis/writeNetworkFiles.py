import os
import sys
import argparse
from AdaptivePELE.utilities import utilities
import matplotlib.pyplot as plt
try:
    # This might fail for older versions of matplotlib (e.g in life cluster)
    plt.style.use("ggplot")
except:
    pass


def parseArguments():
    desc = "Write the information related to the conformation network to file\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("clusteringObject", type=str, help="Path to the clustering object")
    parser.add_argument("suffix", type=str, help="Suffix to append to file names")
    parser.add_argument("metricCol", type=int, help="Column of the metric of interest")
    parser.add_argument("-o", type=str, default=None, help="Output path where to write the files")
    parser.add_argument("-c", "--cond", type=str, default="min", help="Condition on the metric optimality, options are max or min")
    parser.add_argument("-b", "--bindEn", type=int, default=None, help="Column of the binding energy in the report file")
    args = parser.parse_args()
    return args.clusteringObject, args.suffix, args.metricCol, args.o, args.cond, args.bindEn


if __name__ == "__main__":
    clusteringObject, suffix, metricCol, outputPath, metricOptimization, bindingEnergy = parseArguments()
    if outputPath is not None:
        outputPath = os.path.join(outputPath, "")
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)
    else:
        outputPath = ""
    sys.stderr.write("Reading clustering object...\n")
    cl = utilities.readClusteringObject(clusteringObject)
    if cl.conformationNetwork is None:
        sys.exit("Clustering object loaded has no conformation network!!")
    conf = cl.conformationNetwork
    optimalCluster = cl.getOptimalMetric(metricCol, simulationType=metricOptimization)
    pathway = conf.createPathwayToCluster(optimalCluster)
    if not os.path.exists(outputPath+"conformationNetwork%s.edgelist" % suffix):
        sys.stderr.write("Writing conformation network...\n")
        conf.writeConformationNetwork(outputPath+"conformationNetwork%s.edgelist" % suffix)
    if not os.path.exists(outputPath+"FDT%s.edgelist" % suffix):
        sys.stderr.write("Writing FDT...\n")
        conf.writeFDT(outputPath+"FDT%s.edgelist" % suffix)
    if not os.path.exists(outputPath+"pathwayFDT%s.pdb" % suffix):
        sys.stderr.write("Writing pathway to optimal cluster...\n")
        # cl.writePathwayOptimalCluster(outputPath+"pathwayFDT%s.pdb" % suffix)
        cl.writePathwayTrajectory(pathway, outputPath+"pathwayFDT%s.pdb" % suffix)
    if not os.path.exists(outputPath+"nodesPopulation%s.txt" % suffix):
        sys.stderr.write("Writing nodes population...\n")
        cl.writeConformationNodePopulation(outputPath+"nodesPopulation%s.txt" % suffix)
    if not os.path.exists(outputPath+"nodesMetric%s.txt" % suffix):
        sys.stderr.write("Writing nodes metrics...\n")
        cl.writeClusterMetric(outputPath+"nodesMetric%s.txt" % suffix, metricCol)
    if bindingEnergy is not None:
        plt.figure()
        plt.plot(pathway, [cl.clusters.clusters[i].getMetricFromColumn(bindingEnergy) for i in pathway])
        plt.xlabel("Cluster number")
        plt.ylabel("Binding energy(kcal/mol)")
        plt.savefig(outputPath+"bindingEnergy_%s.png" % suffix)
    plt.figure()
    plt.plot(pathway, [cl.clusters.clusters[i].contacts for i in pathway])
    plt.xlabel("Cluster number")
    plt.ylabel("Contacts ratio")
    plt.savefig(outputPath+"contacts_%s.png" % suffix)
    plt.figure()
    plt.plot(pathway, [cl.clusters.clusters[i].getMetricFromColumn(3) for i in pathway])
    plt.xlabel("Cluster number")
    plt.ylabel("Energy(kcal/mol)")
    plt.savefig(outputPath+"totalEnergy_%s.png" % suffix)
    plt.show()
