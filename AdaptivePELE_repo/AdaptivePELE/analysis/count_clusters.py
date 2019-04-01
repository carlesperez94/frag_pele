import numpy as np
import os
import mdtraj as md
import glob
import matplotlib.pyplot as plt
import argparse


def obtainLigandIndexes(trajectory, ligand):
    """
    Extract the indexes for the ligand in the trajectory
    :param trajectory: mdtraj trajectory
    :param ligand: name of the ligand
    :return: list of the atom indexes of the heavy atoms of the ligand
    """
    residueIndexes = []
    for residue in trajectory.topology.residues:
        if residue.name == ligand:
            for atom in residue.atoms:
                if "H" not in atom.name:
                    residueIndexes.append(atom.index)
    return residueIndexes


def calculate_norm(member1, member2, clust):
    """
    Calculates the euclidian distance of two mass centers and adds the cluster information
    :param member1: coords of the first mass center
    :param member2: coords of the second mass center
    :param clust: cluster that is beign processed
    :return:
    """
    return (np.linalg.norm(member1 - member2), clust)


def extrac_traj_num(traj):
    """
    Extracts the number of the trajectory
    :param traj: name of the trajectory (str)
    :return: number of the trajectory (int)
    """
    num = traj.split(".")[0].split("_")[-1]
    return int(num)


def load_cluster_data(clusters_file):
    """
    Loads into  a dictionary the coords of each cluster center
    :param clusters_file: name of the MSM clusters file
    :return: dict with the clusters as keys and the coords as values
    """
    clusters_data = {}
    with open(clusters_file, "r") as inputfile:
        for line in inputfile:
            line = line.split()
            clusters_data[int(line[1])] = [float(line[-6]), float(line[-5]), float(line[-4])]
    return clusters_data


def parseArguments():
    """
        Parse command line arguments

        :returns: str, str, str, str, str, str, bool, bool -- path with the MSM clusters.pdb file,
        path to the reference file, name of the ligand, path to the trajectories,
        string that matches the trajectories names, path to the topology if needed,
        whether to plot the clustering distances, wether to plot the trajectory data
    """
    desc = """ Script that plots the MSM clusters that each trajectory has and the distance 
    of the clusters to the reference position """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-clust", required=True, help="Path with the MSM clusters.pdb file")
    parser.add_argument("-ref", required=True,  help="Path to the reference file")
    parser.add_argument("-lig", required=True,  help="Name of the Ligand")
    parser.add_argument("-dat", required=True,  help="Path to the trajectory data")
    parser.add_argument("-name", default="traj*", help="Template that matches the trajectory names")
    parser.add_argument("-top", default=None, help="Path to the topology file if needed")
    parser.add_argument("-dists", action="store_true", default=False, help="whether to plot the clustering distances or not")
    parser.add_argument("-clustdat", action="store_true", default=False, help="whether to plot the Trajectory data or not")

    args = parser.parse_args()
    return args.clust, args.ref, args.lig, args.dat, args.name, args.top, args.dists, args.clustdat


def main(cluster_file, reference, ligand, data_folder, traj_template, topology, plot_clust_dist, plot_traj):

    # Variable declaration
    clusters_map = {}      # Dict with the cluster equivalences
    clusters_counts = {}   # Dict with the number of times that each cluster appears
    newClusters = {}       # Dict with the data of the new clusters
    dataFramedistances = {"Clusters": [], "Distances": []}
    trajdataFrame = {"Frame": [], "Trajectory": [], "Cluster": [], "Trajectory_names": []}
    bar_dataframe = {"Cluster": [], "Counts": []}
    # Loading data from cluster
    clusters_data = load_cluster_data(cluster_file)
    # Sort the trajectories to process
    trajectories = glob.glob(os.path.join(data_folder, traj_template))
    trajectories.sort(key=lambda x: extrac_traj_num(x))
    # Load the data from the reference
    reference_traj = md.load(reference)
    reference_ligand = reference_traj.atom_slice(obtainLigandIndexes(reference_traj, ligand))
    ligand_center = md.compute_center_of_mass(reference_ligand)*10  # multiplies x10 to change from nm to A
    # Sort Clusters according to the distance to the reference
    clusters_distance =[a for a in map(lambda x: calculate_norm(np.asarray(clusters_data[x]), ligand_center, x), clusters_data)]
    clusters_distance.sort()
    # Initialize the dictionaries
    for i, element in enumerate(clusters_distance):
        newClusters[i] = clusters_data[element[1]]
        clusters_map[i] = element[1]
        clusters_counts[i] = 0
        dataFramedistances["Clusters"].append(i)
        dataFramedistances["Distances"].append(element[0])
    # Load trajectory data if necessary
    if plot_traj:
        for i, traj in enumerate(trajectories):
            traj_num = extrac_traj_num(traj)
            trajdataFrame["Frame"].append([])
            trajdataFrame["Trajectory"].append([])
            trajdataFrame["Cluster"].append([])
            print("Starting with trajectory: %s" % traj_num)
            if topology:
                traj_obj = md.load(traj, top=topology)
            else:
                traj_obj = md.load(traj)
            ligandtraj = traj_obj.atom_slice(obtainLigandIndexes(traj_obj, ligand))
            for j, frame in enumerate(ligandtraj):
                center = md.compute_center_of_mass(frame)*10  # multiplies x10 to change from nm to A
                closestClust = min(map(lambda x: calculate_norm(center[0], newClusters[x], x), newClusters))[1]
                clusters_counts[closestClust] += 1
                trajdataFrame["Frame"][i].append(j)
                trajdataFrame["Trajectory"][i].append(int(traj_num))
                trajdataFrame["Cluster"][i].append(int(closestClust))
            trajdataFrame["Trajectory_names"].append("Traj_%s" % traj_num)
        for key in clusters_counts:
            bar_dataframe["Cluster"].append(key)
            bar_dataframe["Counts"].append(clusters_counts[key])
    plt.style.use('ggplot')
    with open("Cluster_equivalences.txt", "w") as clustequi:
        clustequi.write("File with the equivalences with the original clusters and the new clusters\nNew Cluster    Old Cluster\n")
        for key in clusters_map:
            clustequi.write("%s %s\n" % (key, clusters_map[key]))
    if plot_clust_dist:
        plt.plot(dataFramedistances["Clusters"], dataFramedistances["Distances"])
        plt.xlabel("Clusters")
        plt.ylabel("Distance to Reference")
        plt.savefig("clusters_distances.png")

    if plot_traj:
        plt.figure(1)
        plt.subplot(1, 2, 1)
        cmap = plt.get_cmap("viridis")
        data = zip(trajdataFrame["Trajectory"], cmap.colors[::int(len(cmap.colors)/len(trajdataFrame["Trajectory"]))])
        traj_count = 0
        for trajectory, color in data:
            plt.scatter(x=trajdataFrame["Cluster"][traj_count], y=trajectory, c=color)
            traj_count += 1
        plt.xlabel('Clusters')
        plt.ylabel('Trajectory')
        # plt.yticks([i for i in range(len(trajdataFrame["Trajectory_names"]))], trajdataFrame["Trajectory_names"])
        plt.subplot(1, 2, 2)
        plt.bar(x=bar_dataframe["Cluster"], height=bar_dataframe["Counts"], color=cmap.colors)
        plt.xlabel('Clusters')
        plt.ylabel('Counts')
        plt.savefig("clusters_trajectory.png")


if __name__ == "__main__":
    cluster_file, reference, ligand, data_folder, traj_template, topology, plot_clust_dist, plot_traj = parseArguments()
    main(cluster_file, reference, ligand, data_folder, traj_template, topology, plot_clust_dist, plot_traj)
