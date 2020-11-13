import glob
import os
import argparse
import joblib
import mdtraj as md
import pandas as pd

def parseArguments():
    """
        Parse the command-line options
        :returns: str, int, int --  path to file to results folder,
            index of the first atom,
            index of the second atom
    """
    desc = "It includes the atom-atom distance of the specified ones to report files\n"
    parser = argparse.ArgumentParser(description=desc)
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("sim_folder", type=str, help="Path to the simulation results.")

    required_named.add_argument("-a", "--atoms", type=int, nargs='+', action='append', 
                                help="List of pairs of indexes with atoms to compute the distance between them.")
    parser.add_argument("-t", "--traj", default="trajectory_",
                        help="Trajectory file prefix.")
    parser.add_argument("-r", "--rep", default="report_",
                        help="Report file prefix.")
    parser.add_argument("-p", "--proc", type=int, default=4,
                        help="Number of processors to paralellize the computation.")
    args = parser.parse_args()
    return args.sim_folder, args.atoms, args.traj, args.rep, args.proc


def compute_atom_atom_dist(infile, atoms_list):
    distances = []
    names = []
    traj = md.load_pdb(infile)
    print(atoms_list)
    for at_pair in atoms_list:
        name ="{}-{}".format(traj.topology.atom(at_pair[0]), traj.topology.atom(at_pair[1]))
        distance = md.compute_distances(traj, [at_pair])
        distances.append(distance)
        names.append(name)
    return distances, names


def compute_distances_from_report(atomlist, report, trajectory):
    distances, colnames = compute_atom_atom_dist(trajectory, atomlist)
    new_lines = []
    with open(report) as rep:
        rep_lines = rep.readlines()
        rep_lines = [x.strip("\n") for x in rep_lines]
        for ind, line in enumerate(rep_lines):
            new_content = list(line.split("    "))
            if new_content[-1] == '':
                new_content = new_content[:-1]
            if ind == 0:
                for colname in colnames:
                    new_content.append(colname)
            else:
                for dist in distances:
                    value = "{:.3f}".format(dist[ind-1][0]*10)
                    new_content.append(value)
            new_line = "    ".join(new_content)
            new_lines.append(new_line)
    new_report = "\n".join(new_lines)
    new_rep_name = report.split("/")
    new_rep_name[-1] = "dist" + new_rep_name[-1]
    new_rep_name = "/".join(new_rep_name)
    with open(new_rep_name, "w") as out:
        out.write(new_report)
    print("{} completed".format(new_rep_name))

def compute_simulation_distance(sim_folder, atomlist, traj_pref="trajectory_", report_pref="report_", processors=4):
    trajectories = sorted(glob.glob("{}*".format(os.path.join(sim_folder, traj_pref))))
    reports = sorted(glob.glob("{}*".format(os.path.join(sim_folder, report_pref))))
    joblib.Parallel(n_jobs=processors)(joblib.delayed(compute_distances_from_report)(atomlist, report, traj) for report, traj in zip(reports, trajectories))

if __name__ == '__main__':
    sim_fold, atom_list, traj, report, processors = parseArguments()
    compute_simulation_distance(sim_fold, atom_list, traj, report, processors)
