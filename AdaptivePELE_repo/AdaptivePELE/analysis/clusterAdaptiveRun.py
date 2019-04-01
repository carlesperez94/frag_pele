import matplotlib
matplotlib.use('Agg')
import sys
from functools import partial
import os
import multiprocessing as mp
import glob
from pylab import rcParams
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from AdaptivePELE.utilities import utilities
from AdaptivePELE.freeEnergies import cluster, extractCoords
from AdaptivePELE.analysis import splitTrajectory, simulationToCsv 

def parseArgs():
    parser = argparse.ArgumentParser(description="Script that reclusters the Adaptive clusters")
    parser.add_argument('nClusters', type=int, help="Number of clusters")
    parser.add_argument('crit1', type=int, help="Metric1 to calculate from clusters")
    parser.add_argument('crit2', type=int, help="Metric2 to calculate from clusters")
    parser.add_argument("ligand_resname", type=str, help="Name of the ligand in the PDB")
    parser.add_argument("-atomId", nargs="*", default="", help="Atoms to use for the coordinates of the conformation, if not specified use the center of mass")
    parser.add_argument('--o', type=str, help="Output folder", default="Cluster_analisis")
    parser.add_argument('--top', type=str, help="Topology file", default=None)
    parser.add_argument('--cpus', type=int, help="Cpus to use", default=1)
    parser.add_argument('--report', type=str, help="Report filenames i.e. run_report_", default="report_")
    parser.add_argument('--traj', type=str, help="Trajectory filenames i.e. run_trajectory_", default="trajectory_")
    parser.add_argument('--use_pdb', action="store_true", help="To use when having pdb files with .xtc extension")
    parser.add_argument('--png', action="store_true", help="Save plot in png format")
    parser.add_argument('--CA', action="store_true", help="Cluster by CA")
    parser.add_argument('--sidechains', action="store_true", help="Cluster by sidechain RMSD")
    parser.add_argument('--restart', action="store_true", help="Restart analysis from previous clusters")
    args = parser.parse_args()
    return args.nClusters, args.crit1, args.crit2, args.ligand_resname, args.atomId, args.o, args.top, args.cpus, args.report, args.traj, args.use_pdb, args.png, args.CA, args.sidechains, args.restart


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def write_snapshot(snap_num, trajectory, filename, topology=None, use_pdb=False):
    if not topology:
        snapshots = utilities.getSnapshots(trajectory, topology=topology, use_pdb=use_pdb)
        with open(filename, "w") as fw:
            fw.write(snapshots[snap_num])
    else:
        splitTrajectory.main("", [trajectory, ], topology, [snap_num+1,],template=filename, use_pdb=use_pdb)




def plotClusters(fields1, fields2, crit1, crit2, output, png=False):
    labels = ["cluster_{}".format(i) for i in np.arange(len(fields1))]
    fig, ax = plt.subplots()
    ax.scatter(fields1, fields2, label=labels)
    ax.set_title('RMSD Cluster {} vs {}'.format(crit1, crit2))
    ax.set_xlabel(crit1)
    ax.set_ylabel(crit2)
    print("Plotting")
    if png:
        fig.savefig(os.path.join(output, "ClusterMap.png"))
    else:
        fig.savefig(os.path.join(output, "ClusterMap.pdf"), format='pdf', dpi=1200)

def writePDB(pmf_xyzg, title="clusters.pdb"):
    templateLine = "HETATM%s  H%sCLT L 502    %s%s%s  0.75%s           H\n"

    content = ""
    for j, line in enumerate(pmf_xyzg):
        number = str(j).rjust(5)
        number3 = str(j).ljust(3)
        x = ("%.3f" % line[0]).rjust(8)
        y = ("%.3f" % line[1]).rjust(8)
        z = ("%.3f" % line[2]).rjust(8)
        g = 0
        content += templateLine % (number, number3, x, y, z, g)

    with open(title, 'w') as f:
        f.write(content)


def writeInitialStructures(field1, field2, crit1, crit2, centers_info, filename_template, traj, topology=None, use_pdb=False):
    for cluster_num, field1, field2 in zip(centers_info, field1, field2):
        epoch_num, traj_num, snap_num = map(int, centers_info[cluster_num]['structure'])
        trajectory = "{}/{}{}.xtc".format(epoch_num, traj, traj_num) if topology else "{}/{}{}.pdb".format(epoch_num, traj, traj_num)
        snapshots = utilities.getSnapshots(trajectory, topology=topology, use_pdb=use_pdb)
        filename = filename_template.format(cluster_num, crit1, field1, crit2, field2)
        if not topology:
            with open(filename, "w") as fw:
                fw.write(snapshots[snap_num])
        else:
            splitTrajectory.main("", [trajectory, ], topology, [snap_num+1,],template=filename, use_pdb=use_pdb)


def get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters):
    centersInfo = {x: {"structure": None, "minDist": 1e6, "center": None} for x in range(num_clusters)}

    trajFiles = glob.glob(os.path.join(trajectoryFolder, trajectoryBasename))
    for traj in trajFiles:
        _, epoch, iTraj = os.path.splitext(traj)[0].split("_", 3)
        trajCoords = np.loadtxt(traj)
        if len(trajCoords.shape) < 2:
            trajCoords = [trajCoords]
        for snapshot in trajCoords:
            nSnap = snapshot[0]
            snapshotCoords = snapshot[1:]
            dist = np.sqrt(np.sum((clusterCenters-snapshotCoords)**2, axis=1))
            for clusterInd in range(num_clusters):
                if dist[clusterInd] < centersInfo[clusterInd]['minDist']:
                    centersInfo[clusterInd]['minDist'] = dist[clusterInd]
                    centersInfo[clusterInd]['structure'] = (epoch, int(iTraj), nSnap)
                    centersInfo[clusterInd]['center'] = snapshotCoords
    return centersInfo


def get_metric(criteria, epoch_num, traj_num, snap_num, report):
    report = os.path.join(str(epoch_num), "{}{}".format(report, traj_num))
    report_data = pd.read_csv(report, sep='    ', engine='python')
    print(snap_num, type(snap_num))
    print(type(report_data["Step"].tolist()[0]))
    print(report_data)
    print(report_data["numberOfAcceptedPeleSteps"] == (snap_num-1))
    value = report_data[(report_data["numberOfAcceptedPeleSteps"] == (snap_num-1))].values[0][criteria-1]
    header = list(report_data)[criteria-1]
    return value, header


def assesClusterConvergence(df, num_clusters, traj_name = "trajectory_", topology=None):
    for i in range(num_clusters):
        path = "ClustersSummary/Cluster{}".format(i)
        if not os.path.exists(path):
            os.makedirs(path)
        selected_clust = df[df["Cluster"] == float(i)].astype(float)
        maximum = selected_clust.nlargest(10, "Binding Energy")
        epochs = maximum[simulationToCsv.EPOCH].values
        trajs = maximum[simulationToCsv.TRAJ].values
        steps = maximum[simulationToCsv.STEPS].values
        for j, (epoch, traj, step) in enumerate(zip(epochs, trajs, steps)):
            trajectory = "{}/{}{}.xtc".format(int(epoch), traj_name, int(traj)) if topology else "{}/{}{}.pdb".format(int(epoch), traj_name, int(traj))
            write_snapshot(int(step), trajectory, os.path.join(path, "Cluster_{}_{}.pdb".format(i, j)), topology=topology)
        

def save_to_df(input):
    df, traject, dtraj = input
    dfs_tmp = []
    epoch, traj = traject.strip(".dat").split("_")[-2:]
    for i, d in enumerate(dtraj):
        df_tmp = df[(df[simulationToCsv.EPOCH] == float(epoch)) & (df[simulationToCsv.TRAJ] == float(traj)) & (df[simulationToCsv.STEPS] == float(i))]
        df_tmp["Cluster"] = d
        dfs_tmp.append(df_tmp)
        #df.update(df_tmp)
    return dfs_tmp

def main(num_clusters, criteria1, criteria2, ligand_resname, output_folder = "ClusterCentroids", atom_ids="", cpus=2, topology=None, report="report_", traj="trajectory_", use_pdb=False, png=False, CA=0, sidechains=0, restart="all"):
    #Create multiprocess pool
    if cpus>1:
        pool = mp.Pool(cpus)
    else:
        pool=mp.Pool(1)
    #Extract COM ligand for each snapshot
    if not glob.glob("allTrajs/traj*"):
    	extractCoords.main(lig_resname=ligand_resname, non_Repeat=True, atom_Ids=atom_ids, nProcessors=cpus, parallelize=True, topology=topology, use_pdb=use_pdb, protein_CA=CA, sidechains=sidechains)

    print("Clusterize trajectories by RMSD of COM")
    trajectoryFolder = "allTrajs"
    trajectoryBasename = "*traj*"
    stride = 1
    clusterCountsThreshold = 0
    folders = utilities.get_epoch_folders(".")
    folders.sort(key=int)
    if not restart:

        clusteringObject = cluster.Cluster(num_clusters, trajectoryFolder,
                                           trajectoryBasename, alwaysCluster=True,
                                           stride=stride)
        clusteringObject.clusterTrajectories()
        clusteringObject.eliminateLowPopulatedClusters(clusterCountsThreshold)
        clusterCenters = clusteringObject.clusterCenters
        np.savetxt("clustercenters.dat", clusterCenters)
        dtrajs = clusteringObject.dtrajs

        print("Extract metrics for each snapshot")
        min_metric_trajs = {}
        epochs = [folder for folder in glob.glob("./*/") if folder.isdigit()]
        reports = simulationToCsv.gather_reports()
        fields = simulationToCsv.retrieve_fields(reports[0])
        df = simulationToCsv.init_df(fields)
        df = simulationToCsv.fill_data(reports, df, pool)

        print("Update data with metrics and clusters")
        df.index = range(df.shape[0])
        df["Cluster"] = [None]*df.shape[0]
        input_list = [ [df, Traj, d] for d, Traj in zip(dtrajs, clusteringObject.trajFilenames) ]
        results = pool.map(save_to_df, input_list)
        for data in results:
            for df_tmp in data:
                df.update(df_tmp)
        df.to_csv("Simulation.csv", index=False) 
    if restart:
        df = pd.read_csv("Simulation.csv")    
        clusterCenters = np.loadtxt("clustercenters.dat")
        print(clusterCenters)
    centersInfo = get_centers_info(trajectoryFolder, trajectoryBasename, num_clusters, clusterCenters)
    COMArray = [centersInfo[i]['center'] for i in range(num_clusters)]

    print("Retrieve clusters and metric")
    fields1 = []
    fields2 = []
    print(centersInfo)
    for cluster_num in centersInfo:
        epoch_num, traj_num, snap_num = map(int, centersInfo[cluster_num]['structure'])
        field1, crit1_name = get_metric(criteria1, epoch_num, traj_num, snap_num, report)
        field2, crit2_name = get_metric(criteria2, epoch_num, traj_num, snap_num, report)
        fields1.append(field1)	
        fields2.append(field2)	

    if output_folder is not None:
        outputFolder = os.path.join(output_folder, "")
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)
    else:
        outputFolder = ""
    print("Output structures")
    writePDB(COMArray, outputFolder+"clusters_%d_KMeans_allSnapshots.pdb" % num_clusters)
    writeInitialStructures(fields1, fields2, crit1_name, crit2_name, centersInfo, outputFolder+"cluster_{}_{}_{}_{}_{}.pdb", traj, topology=topology, use_pdb=use_pdb) 
    plotClusters(fields1, fields2, crit1_name, crit2_name, outputFolder, png=png)
    assesClusterConvergence(df, num_clusters, traj, topology)
    return 

if __name__ == "__main__":
    n_clusters, criteria1, criteria2, lig_name, atom_id, output, top, cpus, report, traj, use_pdb, png, CA, sidechains, restart = parseArgs()
    main(n_clusters, criteria1, criteria2, lig_name, output, atom_id, cpus, top, report, traj, use_pdb, png, CA, sidechains, restart)
