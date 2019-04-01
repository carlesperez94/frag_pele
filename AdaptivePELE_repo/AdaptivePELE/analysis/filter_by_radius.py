from __future__ import print_function
import argparse
import numpy as np
import mdtraj as md
import multiprocessing as mp
import AdaptivePELE.analysis.trajectory_processing as tp
try:
    basestring
except NameError:
    basestring = str


def parseArguments():
    desc = "Program that filters the trajectories outside or inside an sphere defined by a radi and a point"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--append_results", action="store_true", help="Flag to set whether all the filtered points should be merged together")
    parser.add_argument("--filter_larger", action="store_false", help="Flag to set whether we should filter the points outside or inside the radi (default, ie not specified, is outside)")
    parser.add_argument("--processors", type=int, default=4, help="Number of cpus to use")
    parser.add_argument("--minimum_length", type=int, default=300, help="Minimum steps for the filtered traj to be considered valid")
    parser.add_argument("radi", type=int, default=2, help="Number of cpus to use")
    parser.add_argument("--center_point", type=int, nargs='+', help="Point to filter", required=True)
    parser.add_argument("resname", help="Ligand resname")
    parser.add_argument("topology", help="Glob string for the topology")
    parser.add_argument("trajectories", help="Glob string for the trajectories")
    args = parser.parse_args()
    return args.trajectories, args.topology, args.resname, args.radi, args.center_point, args.processors, args.filter_larger, args.append_results


def radifilter(num, traj, top, ligand, radi, center_point, minimum_length, Larger, Append_Results):
    valid_frames = []
    NOT_CORD = False
    interval_start = None
    interval_final = None
    final_traj = None
    off_count = 0
    md_traj = md.load(traj, top=top)
    if isinstance(center_point, basestring):
        NOT_CORD = True
        DUM_traj = md_traj.atom_slice(tp.extract_ligand_indexes(md_traj, center_point))
        DUM_center = md.compute_center_of_mass(DUM_traj) * 10
    else:
        center_point = np.array(center_point)
    ligand_traj = md_traj.atom_slice(tp.extract_ligand_indexes(md_traj, ligand))
    center_traj = md.compute_center_of_mass(ligand_traj) * 10
    for i, frame in enumerate(center_traj):
        if NOT_CORD:
            center_point = DUM_center[i]
        distance = np.linalg.norm(frame-center_point)
        if (distance <= radi and not Larger) or (distance >= radi and Larger):
            if interval_start is None:
                interval_start = i
            interval_final = i
        elif interval_start is not None:
            valid_frames.append((interval_start, interval_final))
            interval_start = None
            interval_final = None
    if interval_start is not None:
        valid_frames.append((interval_start, interval_final))
    if len(valid_frames) == 0:
        print("traj %s does not have any frame inside the defined sphere" % num)
    for interval in valid_frames:
        if not Append_Results:
            if (interval[1] - interval[0]) < minimum_length:
                print("Interval from %s %s of traj %s not long enough" % (interval[0], interval[1], num))
            else:
                if off_count:
                    name = "%s.%s" % (num, off_count)
                else:
                    name = num
                md_traj[interval[0]:interval[1]+1].save_xtc("trajectory_radi_%s_filtered_%s.xtc" % (radi, name))
                print("trajectory_radi_%s_filtered_%s.xtc with length %s" % (radi, name, (interval[1] - interval[0])))
                off_count += 1
        else:
            if off_count:
                print("%s.%s" % (num, off_count))
            if final_traj:
                final_traj = final_traj + md_traj[interval[0]:interval[1]+1]
            else:
                final_traj = md_traj[interval[0]:interval[1]+1]
            off_count += 1
    if Append_Results:
        final_traj.save_xtc("trajectory_radi_%s_filtered_%s.xtc" % (radi, num))


def main(trajectory_template, topology, ligand, radi, center_point, minimum_length, processors, Larger, Append_Results):
    pool = mp.Pool(processors)
    workers = []
    num = 0
    for traj, top in tp.load_trajs(trajectory_template, topology, PELE_order=False):
        print("Procesing %s num %s" % (traj, num))
        workers.append(pool.apply_async(radifilter, args=(num, traj, top, ligand, radi, center_point, minimum_length, Larger, Append_Results)))
        num = num + 1
    for worker in workers:
        worker.get()
    print("FINISHED")

if __name__ == "__main__":
    trajectory_template, topology, ligand, radi, center_point, minimum_length, processors, Larger, Append_Results = parseArguments()
    main(trajectory_template, topology, ligand, radi, center_point, minimum_length, processors, Larger, Append_Results)
