from __future__ import absolute_import, division, print_function, unicode_literals
import argparse
from AdaptivePELE.atomset import atomset, SymmetryContactMapEvaluator
from AdaptivePELE.utilities import utilities
import matplotlib.pyplot as plt


def parseArguments():
    desc = "Calculate the histogram of ContactMap contacts over a trajectory, either as pdb files or from a clustering\n"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("resname", type=str, help="Ligand resname in pdb")
    parser.add_argument("contactThreshold", type=int, help="Contact threshold to calculate contactMap")
    parser.add_argument("-trajectory", type=str, nargs='+', help="Path to the trajectory or pdbs to analyse")
    parser.add_argument("-clustering", help="Path to the clustering object to analyse")
    parser.add_argument("-nRes", type=int, default=10, help="Number of top residues to display")
    parser.add_argument("--top", type=str, default=None, help="Topology file needed for non-pdb trajectories")
    args = parser.parse_args()
    return args.trajectory, args.clustering, args.nRes, args.resname, args.contactThreshold, args.top


def generateConformations(resname, clAcc, trajectory, topology):
    if topology is None:
        topology_contents = None
    else:
        topology_contents = utilities.getTopologyFile(topology)
    if clAcc is None:
        for traj in trajectory:
            snapshots = utilities.getSnapshots(traj, topology=topology)
            for snapshot in snapshots:
                PDBobj = atomset.PDB()
                PDBobj.initialise(snapshot, resname=resname, topology=topology_contents)
                yield PDBobj
    else:
        for cluster in clAcc.clusters.clusters:
            yield cluster.pdb

if __name__ == "__main__":
    traj_name, clustering, nRes, lig_resname, contactThreshold, top = parseArguments()

    if clustering is None:
        clusterAcc = None
    else:
        clusetrAcc = utilities.readClusteringObject(clustering)

    totalAcc = []
    symEval = SymmetryContactMapEvaluator.SymmetryContactMapEvaluator()
    refPDB = None

    for pdb in generateConformations(lig_resname, clusetrAcc, traj_name, top):
        if refPDB is None:
            refPDB = pdb
        contactMap, foo = symEval.createContactMap(pdb, lig_resname, contactThreshold)
        if len(totalAcc):
            totalAcc += contactMap.sum(axis=0, dtype=bool).astype(int)
        else:
            totalAcc = contactMap.sum(axis=0, dtype=bool).astype(int)

    proteinList = symEval.proteinList
    residueCounts = {}
    totCounts = 0
    for atomID, counts in zip(proteinList, totalAcc):
        res = refPDB.atoms[atomID].resnum
        totCounts += counts
        if res in residueCounts:
            residueCounts[res] += counts
        else:
            residueCounts[res] = counts

    for res in residueCounts:
        residueCounts[res] /= float(totCounts)

    print("Residue\tResidue frequency")
    for res in sorted(residueCounts, key=lambda x: residueCounts[x], reverse=True)[:nRes]:
        print("%s\t%.4f" % (res, residueCounts[res]))

    plt.figure()
    plt.ylabel("Contact frequency")
    plt.xlabel("Residue number")
    plt.bar(list(residueCounts.keys()), residueCounts.values())
    plt.savefig("hist_CM.png")
    plt.show()
