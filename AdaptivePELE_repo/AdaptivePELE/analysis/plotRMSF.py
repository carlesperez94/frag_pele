from __future__ import print_function, unicode_literals
import argparse
import numpy as np
from AdaptivePELE.atomset import atomset
from AdaptivePELE.utilities import utilities
import matplotlib.pyplot as plt


def parseArguments():
    desc = "Calculate the aggregate RMSF for every residue of the protein\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("trajectory", type=str, nargs='+', help="Path to the trajectory or pdbs to analyse")
    parser.add_argument("-ref", default=None, help="Reference structure, if not specified use average of the trajectories")
    parser.add_argument("-nRes", type=int, default=10, help="Number of top RMSF residues to display")
    parser.add_argument("-top", type=str, default=None, help="Topology file for non-pdb trajectories")
    args = parser.parse_args()
    return args.trajectory, args.ref, args.nRes, args.top


def extractAvgPDB(trajs, topology, topology_content):
    nSnapshots = 0
    avgStruct = {}
    snapshotsTot = []
    for traj in trajs:
        snapshots = utilities.getSnapshots(traj, topology=topology)
        for snapshot in snapshots:
            nSnapshots += 1
            PDB = atomset.PDB()
            PDB.initialise(snapshot, type="PROTEIN", topology=topology_content)
            snapshotsTot.append(PDB)
            for atomID, atom in PDB.atoms.items():
                if atomID in avgStruct:
                    avgStruct[atomID] += (atom.getAtomCoords()-avgStruct[atomID])/nSnapshots
                else:
                    avgStruct[atomID] = atom.getAtomCoords()
    return avgStruct, snapshotsTot


def mapReference(ref, trajs, topology, topology_content):
    refPDB = atomset.PDB()
    refPDB.initialise(ref, type="PROTEIN", topology=topology_content)
    avgStruct = {atomID: refPDB.atoms[atomID].getAtomCoords() for atomID in refPDB.atoms}
    snapshotsTot = []
    for traj in trajs:
        snapshots = utilities.getSnapshots(traj, topology=topology)
        for snapshot in snapshots:
            PDB = atomset.PDB()
            PDB.initialise(snapshot, type="PROTEIN", topology=topology_content)
            snapshotsTot.append(PDB)

    return avgStruct, snapshotsTot

if __name__ == "__main__":
    trajs, ref, nResidues, top = parseArguments()
    if top is None:
        top_content = None
    else:
        top_content = utilities.getTopologyFile(top)
    if ref is None:
        avgPDB, totPDBs = extractAvgPDB(trajs, top, top_content)
    else:
        avgPDB, totPDBs = mapReference(ref, trajs, top, top_content)
    RMSF = {atom: 0.0 for atom in avgPDB}
    residueMapping = {}
    # TODO: Handle multiple chains and insertion residues in PDB
    for PDBobj in totPDBs:
        for atomID, atom in PDBobj.atoms.items():
            RMSF[atomID] += np.sum((atom.getAtomCoords()-avgPDB[atomID])**2)
    for atomID, atom in PDBobj.atoms.items():
        if atom.resnum not in residueMapping:
            residueMapping[atom.resnum] = {atomID}
        else:
            residueMapping[atom.resnum].add(atomID)
    RMSFresidue = {}
    for residue, atoms in residueMapping.items():
        RMSFresidue[residue] = sum([RMSF[atom] for atom in atoms])
        RMSFresidue[residue] /= len(atoms)
        RMSFresidue[residue] = np.sqrt(RMSFresidue[residue])

    print("Residue\tRMSF")
    for res in sorted(RMSFresidue, key=lambda x: RMSFresidue[x], reverse=True)[:nResidues]:
        print("%s\t%.4f" % (res, RMSFresidue[res]))

    plt.plot(list(RMSFresidue.keys()), RMSFresidue.values(), 'x')
    plt.xlabel("Residue number")
    plt.ylabel("RMSF")
    plt.savefig("RMSF-residue.png")
    plt.show()
