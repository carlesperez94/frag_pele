from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.utilities import utilities
from AdaptivePELE.atomset import atomset
import argparse
import os
try:
    basestring
except NameError:
    basestring = str


def parseArguments():
    desc = "Program that writes a trajectory into separate pdbs."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("files", nargs='+', help="Trajectory files to split")
    parser.add_argument("--structs", nargs='*', type=int, default=None, help="Snapshots to write, start with 1")
    parser.add_argument("-o", type=str, default=".", help="Output dir")
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories")
    args = parser.parse_args()
    return args.files, args.o, args.top, args.structs


def main(outputDir, files, topology, structs, template=None, use_pdb=False):
    found=False
    if outputDir:
        utilities.makeFolder(outputDir)
    if topology is not None:
        topology_contents = utilities.getTopologyFile(topology)
    else:
        topology_contents = None
    if structs is not None:
        structs = set(structs)
    for f in files:
        name = os.path.split(f)[-1]
        templateName = os.path.join(outputDir,template) if template else os.path.join(outputDir, name[:-4] + "_%d.pdb")
        snapshots = utilities.getSnapshots(f, topology=topology, use_pdb=use_pdb)
        for i, snapshot in enumerate(snapshots):
            if structs is not None and i+1 not in structs:
                continue
            if not isinstance(snapshot, basestring):
                PDB = atomset.PDB()
                PDB.initialise(snapshot, topology=topology_contents)
                snapshot = PDB.get_pdb_string(model_num=i+1)
            if template:
                with open(templateName, 'w') as of:
                    of.write(snapshot)
                found=True
            else:
                with open(templateName % i, 'w') as of:
                    of.write(snapshot)
                found=True
    return found

if __name__ == "__main__":
    traj_files, output_dir, top, conf = parseArguments()
    main(output_dir, traj_files, top, conf)
