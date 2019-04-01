import glob
import itertools


class CustomError(Exception):
    __module__ = Exception.__module__


def load_trajs(trajects, topology, PELE_order=True):
    """
    Generator that yields tuples with the trajectory name in  and its topology in numerical order.
    This method requieres that the trajectories and topologies follow this pattern X_num.extension
    Where X can be anything, num must be an integer and extension can be anyithing.
    The method acepts string templates with wildcards such as trajectory_*xtc and system_*prmtop.
    If the topology string only matches one file the method will asume that all the trajectories have the same topology.
    If multiple files match the topology, the method will use the PELE naming convention if PELE_order =  True
    or the topologies in numerical order otherwise.
    PARAMS:
    trajects: string that matches the trajectories to be procesed
    topology: string that matches the topologies to use
    PELE_order: Boolean to follow PELE naming convention or not.
    OUTPUT:
    Tuples with the trajectory name in the first position and its topology in the second one.
    """
    trajectories = [x for x in glob.glob(trajects)]
    topologies = [x for x in glob.glob(topology)]
    if len(topologies) > 1:
        topology_extension = len(topologies[0].split(".")[-1])+1
        topologies = sorted(topologies, key=lambda x: int(x.split("_")[-1][:-topology_extension]))
        if PELE_order:
            topologies = topologies[1:] + topologies[:1]
    extension_len = len(trajectories[0].split(".")[-1])+1
    trajectories = sorted(trajectories, key=lambda x: int(x.split("_")[-1][:-extension_len]))
    for file_pair in zip(trajectories, itertools.cycle(topologies)):
        yield file_pair


def extract_ligand_indexes(traj, ligand):
    """
    Function that extracts the indexes corresponding to the ligand
    PARAMS:
    traj: mdtraj trajectory object, with topology
    ligand: string with the name that the ligand has.
    return a list with the indexes of the ligand or an exceptionn if the ligand is not found.
    """
    ligand_indexes = []
    for residue in traj.topology.residues:
        if residue.name == ligand:
            for atom in residue.atoms:
                if atom.element != "H":
                    ligand_indexes.append(atom.index)
    if len(ligand_indexes) == 0:
        raise CustomError("Choosed Ligan %s does not apper in the trajectory" % ligand)
    return ligand_indexes


def save_clusters_to_pdb(matrix_list, output):
    # Method that given a matrix of points, shape (x,3) returns a PDB of dummy atoms in that positions
    try:
        matrix_list[0][0]
    except:
        tmp_list = matrix_list
        matrix_list = []
        matrix_list.append(tmp_list)
    chains_list = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    off_chain = 0
    off_res = 0
    resnames = ["DU%s" % x for x in range(10)]
    with open(output, "w") as out:
        for matrix in matrix_list:
            i = 0
            if off_chain >= len(chains_list):
                off_chain = 0
                off_res += 1
                if off_res > 9:
                    off_res = 0
            chain = chains_list[off_chain]
            resname = resnames[off_res]
            off_chain += 1
            for frame in matrix:
                i += 1
                out.write("%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n" % ("HETATM", i, "H", resname, chain, i, round(frame[0], 3), round(frame[1], 3), round(frame[2], 3)))


def dehidratate(traj, other=("Na+", "Cl-")):
    """
    Method that removes the waters and ions from the trajectorie.
    PARAMS:
    traj: mdtraj trajectory with topology
    other: tuple with additional atom names to remove
    system_indexes: list with the indexes of the system, this is to save time, if provided,
    the function will not calculate them again.
    RETURNS:
    mdtraj trajectory without waters and the indexes of the system.
    """
    system_indexes = []
    if not system_indexes:
        for residue in traj.topology.residues:
            if residue.name not in ("HOH", "WAT") and residue.name not in other:
                for atom in residue.atoms:
                    system_indexes.append(atom.index)
    traj = traj.atom_slice(system_indexes)
    return traj


def extract_heavyatom_indexes(traj):
    # Method that extracts the heavy atoms indexes for aligment porpouses, Hydrogens have different atom names depending on the forcefield used
    heavy_indexes = []
    for atom in traj.topology.atoms:
        if atom.element.name not in ("H", "hydrogen") and atom.name not in ("OXT", "DUM") and atom.residue.name != "DUM":
            heavy_indexes.append(atom.index)
    return heavy_indexes
