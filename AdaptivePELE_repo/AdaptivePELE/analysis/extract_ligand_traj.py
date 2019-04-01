from __future__ import print_function
import argparse
import mdtraj as md
import numpy as np
import multiprocessing as mp
import AdaptivePELE.analysis.trajectory_processing as tp


def parseArguments():
    desc = "Program that extracts the center of mass of the ligand and puts them in a pdb."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--processors", type=int, default=4, help="Number of cpus to use")
    parser.add_argument("resname", help="Ligand resname")
    parser.add_argument("output_file", help="Ligand resname")
    parser.add_argument("topology", help="Glob string for the topology")
    parser.add_argument("trajectories", help="Glob string for the trajectories")
    args = parser.parse_args()
    return args.trajectories, args.resname, args.topology, args.output_file, args.processors


def extract_trajectory(traj, ligand_name, topology, output_file):
    matrix = np.ndarray(shape=(0, 3))
    md_traj = md.load(traj, top=topology)
    ligan_traj = md_traj.atom_slice(tp.extract_ligand_indexes(md_traj, ligand_name))
    mass_centers = md.compute_center_of_mass(ligan_traj) * 10
    matrix = np.append(matrix, mass_centers, axis=0)
    return matrix


def main(trajectory_template, ligand_name, topology, output_file, processors):
    pool = mp.Pool(processors)
    workers = []
    matrix_list = []
    # matrix = np.ndarray(shape=(0,3))
    num = 0
    for traj, top in tp.load_trajs(trajectory_template, topology, PELE_order=False):
        print("Procesing %s num %s" % (traj, num))
        num = num + 1
        workers.append(pool.apply_async(extract_trajectory, args=(traj, ligand_name, top, output_file)))
    for worker in workers:
        matrix_list.append(worker.get())
    print("SAVING")
    # for m in matrix_list:
    #     matrix = np.concatenate((matrix, m))
    tp.save_clusters_to_pdb(matrix_list, output_file)
    print("FINISHED")

if __name__ == "__main__":
    trajectory_template, ligand_name, topology, output_file, processors = parseArguments()
    main(trajectory_template, ligand_name, topology, output_file, processors)
