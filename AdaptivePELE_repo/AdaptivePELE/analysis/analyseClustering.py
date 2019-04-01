from __future__ import absolute_import, division, print_function, unicode_literals
import pickle
import numpy as np
import matplotlib.pyplot as plt
from AdaptivePELE.atomset import atomset
from mpl_toolkits.mplot3d import Axes3D
import argparse


def extractCOMMatrix(clusters, resname, topology=None):
    """ Extract a matrix contaning the coordinates of the center of mass of
    the ligand for each cluster structure

        clusters [In] List of clusters
        resname [In] Residue name of the ligand in the pdb
    """
    n = len(clusters)
    cluster_matrix = np.zeros((n, 3))
    metrics = np.zeros(n)
    population = np.zeros(n)
    total_elements = 0
    contacts = np.zeros(n)
    if clusters[0].pdb.isFromPDBFile() and topology is None:
        raise ValueError("Need to pass a topology file to process non-pdb trajectories")
    for index, cluster in enumerate(clusters):
        metrics[index] = cluster.metrics[cluster.metricCol]
        contacts[index] = cluster.contacts
        ligandPDB = atomset.PDB()
        ligandPDB.initialise(cluster.pdb.get_pdb_string(), resname=resname, topology=topology)
        cluster_matrix[index, :] = ligandPDB.extractCOM()
        population[index] = cluster.elements
        total_elements += cluster.elements
    return cluster_matrix, metrics, total_elements, population, contacts


def plotClusters2D(cluster_matrix, metrics, title):
    """ Create all combination of xyz projections in 2D of the scatter plot
    of the center of mass of the ligand with a colormap given by a certain
    quantity (usually a metric or the clusters population)

        cluster_matrix [In] matrix contaning the coordinates of the center of
        mass of the ligand for each cluster structure
        metrics [In] Array with the quantity that will be used to create the
        colormap
        title [In] Title for the plot figure
    """
    ccx = cluster_matrix[:, 0]
    ccy = cluster_matrix[:, 1]
    ccz = cluster_matrix[:, 2]
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex='col',
                             sharey='row')
    fig.suptitle(title)
    scatter1 = axes[0][0].scatter(ccx, ccy, c=metrics)  # ,label="Set %d"%index)
    axes[0][1].scatter(ccz, ccy, c=metrics)  # ,label="Set %d"%index)
    axes[1][0].scatter(ccx, ccz, c=metrics)  # ,label="Set %d"%index)
    fig.colorbar(scatter1)

#    axes[1][0].legend(loc='center right', bbox_to_anchor=[1.8,0.5])
    axes[0][0].set_ylabel('y')
    axes[1][0].set_ylabel('z')
    axes[1][0].set_xlabel('x')
    axes[1][1].axis('off')
    axes[0][1].set_xticks(axes[1][1].get_xticks())
    axes[0][1].set_xticklabels(axes[1][1].get_xticklabels())
    axes[0][1].set_xlabel('z')
    return fig


def plotClusters(cluster_matrix, metrics, title):
    """ Create a 3D scatter plot of the center of mass of the ligand with a
    colormap given by a certain quantity
    (usually a metric or the clusters population)

        cluster_matrix [In] matrix contaning the coordinates of the center of
        mass of the ligand for each cluster structure
        metrics [In] Array with the quantity that will be used to create the
        colormap
        title [In] Title for the plot figure
    """
    fig = plt.figure()
    ax = Axes3D(fig)
    ccx = cluster_matrix[:, 0]
    ccy = cluster_matrix[:, 1]
    ccz = cluster_matrix[:, 2]
    fig.suptitle(title)
    scatter1 = ax.scatter(ccx, ccy, zs=ccz, c=metrics)
    fig.colorbar(scatter1)
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlabel('x')
    return fig


def extractInfo(inputFile):
    clusterInfo = np.loadtxt(inputFile)
    return clusterInfo[:, 1]


def plotClusteringData(pklObjectFilename, resname, titlemetric, titlepopulation,
                       titlecontacts, metricPlotFilename="",
                       populationPlotFilename="", contactsPlotFilename="",
                       metricFlag=False, populationFlag=False,
                       contactsFlag=False, inputFile=None, topology=None):

    with open(pklObjectFilename, "r") as f:
        clObject = pickle.load(f)

    comCoord, metrics, totalElements, population, contacts = extractCOMMatrix(clObject.clusters.clusters, resname, topology=topology)

    if inputFile is not None:
        clustersInfo = extractInfo(inputFile)
        plotClusters(comCoord, clustersInfo, "")

    if metricFlag:
        plot = plotClusters(comCoord, metrics, titlemetric)
        if metricPlotFilename:
            plot.savefig(metricPlotFilename)

    if populationFlag:
        plotContpop = plotClusters(comCoord, population, titlepopulation)
        if populationPlotFilename:
            plotContpop.savefig(populationPlotFilename)

    if contactsFlag:
        plotContcont = plotClusters(comCoord, contacts, titlecontacts)
        if contactsPlotFilename:
            plotContcont.savefig(contactsPlotFilename)

    print("Number of elements", totalElements)


def parseArguments():
    desc = "3D visualization of clustering"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("clusteringObject", type=str, help="Clustering object")
    parser.add_argument("resname", type=str, help="Resname in the pdb")
    parser.add_argument("-metrics", action="store_true", help="Wether to plot the metric of the clusters as color")
    parser.add_argument("-population", action="store_true", help="Wether to plot the population of the clusters as color")
    parser.add_argument("-contacts", action="store_true", help="Wether to plot the contacts of the clusters as color")
    parser.add_argument("-i", type=str, default=None, help="File with cluster-associated information")
    parser.add_argument("-top", type=str, default=None, help="PDB file with topology information, necessary if working with xtc files")
    args = parser.parse_args()

    return args.clusteringObject, args.resname, args.metrics, args.population, args.contacts, args.i

if __name__ == "__main__":
    pklObject_filename, lig_resname, metricsFlag, population_flag, contacts_flag, input_file = parseArguments()

    metricPlot_filename = ""  # "results/contactClusters.png"
    populationPlot_filename = ""  # "results/contactClusterspop.png"
    contactsPlot_filename = ""  # "results/contactClustersContacts.png"
    title_metric = "Metrics Contacts"
    title_population = "Population Contacts"
    title_contacts = "Number of contacts Contacts"
    topology = None
    if top is not None:
        topology = utilities.getTopologyFile(top)

    plotClusteringData(pklObject_filename, lig_resname, title_metric, title_population,
                       title_contacts, metricPlot_filename,
                       populationPlot_filename, contactsPlot_filename, metricsFlag,
                       population_flag, contacts_flag, input_file, topology=topology)
    plt.show()
