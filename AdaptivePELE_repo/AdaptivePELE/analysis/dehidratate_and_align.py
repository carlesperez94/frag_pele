from __future__ import print_function
import argparse
import mdtraj as md
import multiprocessing as mp
import AdaptivePELE.analysis.trajectory_processing as tp


def parseArguments():
    desc = "Program that extracts residue coordinates for a posterior MSM analysis."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--dont_image", action="store_false", help="Flag to set whether trajectories should be imaged before the alignment (if not specfied performs the imaging)")
    parser.add_argument("--offset", type=int, default=0, help="Offset to add to trajectory number")
    parser.add_argument("--processors", type=int, default=4, help="Number of cpus to use")
    parser.add_argument("resname", help="Ligand resname")
    parser.add_argument("reference", help="Ligand resname")
    parser.add_argument("topology", help="Glob string for the topology")
    parser.add_argument("trajectories", help="Glob string for the trajectories")
    args = parser.parse_args()
    return args.trajectories, args.resname, args.topology, args.reference, args.processors, args.offset, args.dont_image


def process_traj(traj, top, ligand_name, reference, num, image=True):
    reference = md.load(reference)
    reference = tp.dehidratate(reference)
    md_traj = md.load(traj, top=top)
    if image:
        md_traj = md_traj.image_molecules()
    nowat_traj = tp.dehidratate(md_traj)
    aligned_traj = nowat_traj.superpose(reference, frame=0, atom_indices=tp.extract_heavyatom_indexes(nowat_traj), ref_atom_indices=tp.extract_heavyatom_indexes(reference))
    aligned_traj.save_xtc("trajectory_aligned_%s.xtc" % num)
    if num == 0:
        aligned_traj[0].save_pdb("top%s.pdb" % ligand_name)


def main(trajectory_template, ligand_name, topology, reference, processors, off_set, image):
    pool = mp.Pool(processors)
    workers = []
    num = off_set
    for traj, top in tp.load_trajs(trajectory_template, topology, PELE_order=True):
        print("Procesing %s num %s with top %s" % (traj, num, top))
        workers.append(pool.apply_async(process_traj, args=(traj, top, ligand_name, reference, num, image)))
        num = num + 1
    for worker in workers:
        worker.get()


if __name__ == "__main__":
    trajectory_template, ligand_name, topology, reference, processors, off_set, image = parseArguments()
    main(trajectory_template, ligand_name, topology, reference, processors, off_set, image)
