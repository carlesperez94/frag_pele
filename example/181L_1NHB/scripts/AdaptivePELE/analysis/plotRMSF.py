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
    args = parser.parse_args()
    return args.trajectory, args.ref, args.nRes


def extractAvgPDB(trajs):
    nSnapshots = 0
    avgStruct = {}
    snapshotsTot = []
    for traj in trajs:
        snapshots = utilities.getSnapshots(traj)
        for snapshot in snapshots:
            nSnapshots += 1
            PDB = atomset.PDB()
            PDB.initialise(snapshot, type="PROTEIN")
            snapshotsTot.append(PDB)
            for atomID, atom in PDB.atoms.iteritems():
                if atomID in avgStruct:
                    avgStruct[atomID] += (atom.getAtomCoords()-avgStruct[atomID])/nSnapshots
                else:
                    avgStruct[atomID] = atom.getAtomCoords()
    return avgStruct, snapshotsTot


def mapReference(ref, trajs):
    refPDB = atomset.PDB()
    refPDB.initialise(ref, type="PROTEIN")
    avgStruct = {atomID: refPDB.atoms[atomID].getAtomCoords() for atomID in refPDB.atoms}
    snapshotsTot = []
    for traj in trajs:
        snapshots = utilities.getSnapshots(traj)
        for snapshot in snapshots:
            PDB = atomset.PDB()
            PDB.initialise(snapshot, type="PROTEIN")
            snapshotsTot.append(PDB)

    return avgStruct, snapshotsTot

if __name__ == "__main__":
    trajs, ref, nResidues = parseArguments()
    if ref is None:
        avgPDB, totPDBs = extractAvgPDB(trajs)
    else:
        avgPDB, totPDBs = mapReference(ref, trajs)
    RMSF = {atom: 0.0 for atom in avgPDB}
    residueMapping = {}
    for PDBobj in totPDBs:
        for atomID, atom in PDBobj.atoms.iteritems():
            RMSF[atomID] += np.sum((atom.getAtomCoords()-avgPDB[atomID])**2)
    for atomID, atom in PDBobj.atoms.iteritems():
        if atom.resnum not in residueMapping:
            residueMapping[atom.resnum] = set([atomID])
        else:
            residueMapping[atom.resnum].add(atomID)
    RMSFresidue = {}
    for residue, atoms in residueMapping.iteritems():
        RMSFresidue[residue] = sum([RMSF[atom] for atom in atoms])
        RMSFresidue[residue] /= len(atoms)
        RMSFresidue[residue] = np.sqrt(RMSFresidue[residue])

    print "Residue\tRMSF"
    for res in sorted(RMSFresidue, key=lambda x: RMSFresidue[x], reverse=True)[:nResidues]:
        print "%s\t%.4f" % (res, RMSFresidue[res])

    plt.plot(RMSFresidue.keys(), RMSFresidue.values(), 'x')
    plt.xlabel("Residue number")
    plt.ylabel("RMSF")
    plt.savefig("RMSF-residue.png")
    plt.show()
