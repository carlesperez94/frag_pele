# coding: utf-8
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import range
import os
import re
import sys
import glob
import shutil
import argparse
import itertools
import numpy as np
import mdtraj as md
from AdaptivePELE.atomset import atomset
from AdaptivePELE.freeEnergies import utils
from AdaptivePELE.utilities import utilities
PARALELLIZATION = True
try:
    import multiprocessing as mp
except ImportError:
    PARALELLIZATION = False
PRODY = True
try:
    import prody as pd
except ImportError:
    PRODY = False


MDTRAJ_FORMATS = set(['.xtc', '.dcd', '.dtr', '.trr', 'mdcrd', 'nc'])
VALID_CM_MODES = ["p-lig", "p-p"]
# consider most extreme atoms
EXTRA_ATOMS = {u"ALA": u"empty", u"VAL": u"CG1", u"LEU": u"CD1", u"ILE": u"CD1",
               u"MET": u"CE", u"PRO": u"empty", u"PHE": u"CE1", u"TYR": u"CE1",
               u"TRP": u"CZ3", u"SER": u"OG", u"THR": u"OG1", u"CYS": u"SG",
               u"ASN": u"OD1", u"GLN": u"OE1", u"LYS": u"NZ", u"HIS": u"CE1",
               u"HIE": u"CE1", u"HID": u"CE1", u"HIP": u"CE1", u"ARG": u"NH1",
               u"ASP": u"OD1", u"GLU": u"OE1", u"GLY": u"empty"}
# consider more central atoms in the side-chain
EXTRA_ATOMS_CENTRAL = {u"ALA": u"empty", u"VAL": u"empty", u"LEU": u"CG", u"ILE": u"CG2",
                       u"MET": u"CG", u"PRO": u"empty", u"PHE": u"CZ", u"TYR": u"CZ",
                       u"TRP": u"CE2", u"SER": u"OG", u"THR": u"OG1", u"CYS": u"SG",
                       u"ASN": u"CG", u"GLN": u"CG", u"LYS": u"CG", u"HIS": u"CG",
                       u"HIE": u"CG", u"HID": u"CG", u"HIP": u"CG", u"ARG": u"CD",
                       u"ASP": u"CG", u"GLU": u"CG", u"GLY": u"empty"}


class Constants(object):
    def __init__(self):
        self.extractedTrajectoryFolder = "%s/extractedCoordinates"
        self.baseExtractedTrajectoryName = "coord_"
        self.reportName = '*report_'
        self.baseGatheredFilename = "traj_*.dat"
        self.outputTrajectoryFolder = "%s/repeatedExtractedCoordinates"
        self.ligandTrajectoryFolder = "ligand_trajs"
        self.ligandTrajectoryBasename = os.path.join(self.ligandTrajectoryFolder, "traj_ligand_%s.pdb")
        self.gatherTrajsFolder = "allTrajs"
        self.gatherTrajsFilename = os.path.join(self.gatherTrajsFolder, "traj_%s_%s.dat")
        self.gatherNonRepeatedFolder = os.path.join(self.gatherTrajsFolder, "extractedCoordinates")
        self.gatherNonRepeatedTrajsFilename = os.path.join(self.gatherNonRepeatedFolder, "traj_%s_%s.dat")


class TopologyCompat(object):
    def __init__(self, pdb_file):
        self.topologyFiles = os.path.abspath(pdb_file)
        self.path = os.path.split(self.topologyFiles)[0]
        self.topologies = [utilities.getTopologyFile(self.topologyFiles)]

    def getTopologyFile(self, epoch, trajectory_number):
        return self.topologyFiles

    def topologyFilesIterator(self):
        yield self.topologyFiles

    def getTopologyIndex(self, epoch, trajectory_number):
        """
            Get the topology index for a particular epoch and trajectory number

            :param epoch: Epoch of the trajectory of interest
            :type epoch: int
            :param trajectory_number: Number of the trajectory to select
            :type trajectory_number: int

            :returns: int -- Index of the corresponding topology
        """
        return 0

    def getTopology(self, epoch, trajectory_number):
        """
            Get the topology for a particular epoch and trajectory number

            :param epoch: Epoch of the trajectory of interest
            :type epoch: int
            :param trajectory_number: Number of the trajectory to select
            :type trajectory_number: int

            :returns: list -- List with topology information
        """
        return self.topologies[0]


class ParamsHandler(object):
    def __init__(self, folderWithTrajs, atom_id, lig_name, total_steps, sequential, writeLigandTrajectory, set_number, protein_CA, noRepeat, numProcessors, parallelize, topol, sidechains, sidechains_folder, CM, use_extra_atoms, CM_mode, dihedrals, dihedrals_projection):
        self.folder_name = folderWithTrajs
        self.atomIds = atom_id
        self.lig_resname = lig_name
        self.numtotalSteps = total_steps
        self.enforceSequential_run = sequential
        self.writeLigandTrajectory = writeLigandTrajectory
        self.setNumber = set_number
        self.protein_CA = protein_CA
        self.non_Repeat = noRepeat
        self.nProcessors = numProcessors
        self.parallelize = parallelize
        self.topology = topol
        self.sidechains = sidechains
        self.sidechain_folder = sidechains_folder
        self.contact_map = CM
        self.extra_atoms = use_extra_atoms
        self.cm_mode = CM_mode
        self.dihedrals = dihedrals
        self.dihedrals_projection = dihedrals_projection
        if self.contact_map and self.cm_mode == "p-lig" and self.lig_resname == "":
            raise ValueError("Ligand resname needed for protein-ligand contact map")
        if self.contact_map and self.cm_mode not in VALID_CM_MODES:
            raise ValueError("Unrecognized type of contact map, valids are: %s" " ".join(VALID_CM_MODES))
        self.com = not self.protein_CA and (self.atomIds is None or len(self.atomIds) == 0) and not self.sidechains and not self.contact_map and not self.dihedrals


def parseArguments():
    desc = "Program that extracts residue coordinates for a posterior MSM analysis.\
            It either extracts the resname COM coordinates or those of an atomId, depending on the input.\
            It then fills the rejected steps, which is not done by PELE.\
            Finally, trajectories are gathered together in the same allTrajs folder.\
            It automatically detects whether it is an adaptive or a sequential PELE run by looking for folders\
            with numeric names."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-f", "--folderWithTrajs", default=".",
                        help="Folder with trajectories (or epochs)")
    # parser.add_argument("-atomId", action=AtomIdAction, help="serial:atomName:resname, e.g. 2048:C1:AIN")
    parser.add_argument("-atomIds", nargs='*', help="serial:atomName:resname, e.g. 2048:C1:AIN. May contain more than one atomId")
    parser.add_argument("-resname", default="", help="Ligand resname")
    parser.add_argument("-CA", "--proteinCA", action="store_true", help="Extract protein alpha carbons coordinates")
    parser.add_argument("-s", "--enforceSequential", action="store_true", help="Force the consideration as sequential run (non-adaptive)")
    parser.add_argument("--setNum", type=int, default=0, help="Sets the number to appear in gathered trajectory in order to avoid clashes between different sequential runs. Ignored in adaptive runs.")
    parser.add_argument("-w", "--writeLigandTrajectory", action="store_true", help="It writes a traj_ligand_XXX.pdb file with the ligand coordinates. The user must delete the original trajectory (if wanted)")
    parser.add_argument("-t", "--totalSteps", type=int, default=0, help="Total number of steps in traj. Equivalent to epoch length in adaptive runs")
    parser.add_argument("-nR", "--noRepeat", action="store_true", help="Flag to avoid repeating the rejected steps")
    parser.add_argument("-n", "--numProcessors", type=int, default=None, help="Number of cpus to use")
    parser.add_argument("--top", type=str, default=None, help="Topology file for non-pdb trajectories or path to Adaptive topology object")
    parser.add_argument("--sidechains", action="store_true", help="Flag to extract sidechain coordinates")
    parser.add_argument("-sf", "--sidechains_folder", default=".", type=str, help="Folder with the structures to obtain the sidechains to extract")
    parser.add_argument("--serial", action="store_true", help="Flag to deactivate parallelization")
    parser.add_argument("--contact_map", action="store_true", help="Flag to activate contact map creation")
    parser.add_argument("--extra_atoms", action="store_true", help="Flag to use extra atoms in contact map creation (in addition to alpha carbons)")
    parser.add_argument("--dihedrals", action="store_true", help="Flag to activate dihedral angles calculations")
    parser.add_argument("--dihedrals_projection", action="store_true", help="Flag to project dihedral angles calculations into their cos and sin")
    parser.add_argument("--cm_mode", default="p-lig", help="Type of contact map to create (p-lig for protein-ligand or p-p protein-protein)")
    args = parser.parse_args()

    return args.folderWithTrajs, args.atomIds, args.resname, args.proteinCA, args.enforceSequential, args.writeLigandTrajectory, args.totalSteps, args.setNum, args.noRepeat, args.numProcessors, args.top, args.sidechains, args.sidechains_folder, args.serial, args.contact_map, args.extra_atoms, args.cm_mode, args.dihedrals, args.dihedrals_projection


def loadAllResnameAtomsInPdb(filename, params):
    prunedFileContent = []
    sidechains_bool = bool(params.sidechains)
    with open(filename) as f:
        prunedSnapshot = []
        for line in f:
            if utils.is_model(line):
                prunedFileContent.append("".join(prunedSnapshot))
                prunedSnapshot = []
            elif utils.is_end(line) or utils.is_remark(line) or utils.is_cryst(line):
                continue
            elif line[17:20] == params.lig_resname or utils.isAlphaCarbon(line, params.protein_CA or params.contact_map) or utils.isSidechain(line, sidechains_bool, params.sidechains) or (params.contact_map and params.extra_atoms and utils.extraAtomCheck(line, EXTRA_ATOMS)):
                prunedSnapshot.append(line)
        if prunedSnapshot:
            prunedFileContent.append("".join(prunedSnapshot))
    return prunedFileContent


def extractFilenumber(filename):
    last = filename.rfind('.')
    first = filename.rfind('_')
    number = re.sub("[^0-9]", "", filename[first+1:last])
    return number


def getOutputFilename(directory, filename, baseOutputFilename):
    filenumber = extractFilenumber(filename)
    return os.path.join(directory, baseOutputFilename+filenumber+".dat")


def extractContactMapCoordinatesPDB(allCoordinates, params):
    trajCoords = []
    for coordinates in allCoordinates:
        if params.cm_mode == "p-lig":
            PDB = atomset.PDB()
            PDB.initialise(coordinates, resname=params.lig_resname)
            snapshotCoords = [coord for at in PDB.atomList for coord in PDB.atoms[at].getAtomCoords()]
        else:
            snapshotCoords = []
        PDBCA = atomset.PDB()
        if params.extra_atoms:
            PDBCA.initialise(coordinates, type=u"PROTEIN", extra_atoms=EXTRA_ATOMS)
        else:
            PDBCA.initialise(coordinates, type=u"PROTEIN")
        snapshotCoords.extend([coord for at in PDBCA.atomList for coord in PDBCA.atoms[at].getAtomCoords()])
        trajCoords.append(snapshotCoords)
    return trajCoords


def getLigandAlphaCarbonsCoords(allCoordinates, lig_resname, sidechains=False):
    trajCoords = []
    for coordinates in allCoordinates:
        PDB = atomset.PDB()
        PDB.initialise(coordinates, resname=lig_resname)
        snapshotCoords = [coord for at in PDB.atomList for coord in PDB.atoms[at].getAtomCoords()]
        PDBCA = atomset.PDB()
        if not sidechains:
            PDBCA.initialise(coordinates, type="PROTEIN")
        else:
            PDBCA.initialise(coordinates, type="PROTEIN", heavyAtoms=True)
        snapshotCoords.extend([coord for at in PDBCA.atomList for coord in PDBCA.atoms[at].getAtomCoords()])
        trajCoords.append(snapshotCoords)
    return trajCoords


def getPDBCOM(allCoordinates, lig_resname):
    COMs = []
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=lig_resname, heavyAtoms=True)
        COMs.append(pdb.extractCOM())
    return COMs


def getAtomCoord(allCoordinates, lig_resname, atom_Ids):
    coords = []
    # If ever need to speed this up, build a Trajectory class that inherits from PDB
    # and loads the atom according to the position in the snapshot, rather than looking
    # for the atom
    for coordinates in allCoordinates:
        pdb = atomset.PDB()
        pdb.initialise(coordinates, resname=lig_resname, heavyAtoms=True)
        snapshotcoords = []
        for atomId in atom_Ids:
            snapshotcoords.extend(pdb.getAtom(atomId).getAtomCoords())
        coords.append(snapshotcoords)
    return coords


def writeToFile(COMs, outputFilename):
    with open(outputFilename, 'w') as f:
        for i, line in enumerate(COMs):
            f.write(str(i) + ' ')
            for i in range(len(line) - 1):
                f.write(str(line[i]) + ' ')
            f.write(str(line[-1]) + '\n')


def extractIndexesTopology_CM(topology, lig_resname, CM_mode, use_extra_atoms):
    selection = []
    iline = 0
    if CM_mode == "p-p":
        check_ligand = False
    else:
        check_ligand = True
    with open(topology) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if line[76:80].strip().upper() != "H" and (check_ligand and line[17:20] == lig_resname or utils.isAlphaCarbon(line, True) or use_extra_atoms and utils.extraAtomCheck(line, EXTRA_ATOMS)):
                selection.append(iline)
            iline += 1
    return selection


def extractIndexesTopology(topology, lig_resname, atoms, writeCA, sidechains):
    selection = []
    if atoms:
        atoms_set = set(atoms)
    template = "%s:%s:%s"
    iline = 0
    bool_sidechains = bool(sidechains)
    with open(topology) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if atoms:
                serial_num = line[6:11].strip()
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                if template % (serial_num, atom_name, residue_name) in atoms_set:
                    selection.append(iline)

            elif (line[17:20] == lig_resname or utils.isAlphaCarbon(line, writeCA) or utils.isSidechain(line, bool_sidechains, sidechains)) and line[76:80].strip().upper() != "H":
                selection.append(iline)
            iline += 1
    return selection


def contactMapNonPDB(file_name, params, topology, selected_indices):
    trajectory = md.load(file_name, top=topology)
    atom_pairs = list(itertools.combinations(selected_indices, 2))
    return 10*md.compute_distances(trajectory, atom_pairs, periodic=True)


def calculateDihedrals(file_name, params, topology):
    trajectory = md.load(file_name, top=topology)
    _, psi_angles = md.compute_psi(trajectory, periodic=True)
    _, phi_angles = md.compute_phi(trajectory, periodic=True)
    return np.hstack((psi_angles, phi_angles))


def projectDihedrals(dihedrals):
    cos_proj = np.cos(dihedrals)
    sin_proj = np.sin(dihedrals)
    return np.hstack((cos_proj, sin_proj))


def extractCoordinatesXTCFile(file_name, params, topology, selected_indices):
    trajectory = md.load(file_name, top=topology)
    if params.com:
        # getCOM case
        # convert nm to A
        coordinates = 10*md.compute_center_of_mass(trajectory.atom_slice(selected_indices))
    else:
        coordinates = 10*trajectory.xyz[:, selected_indices, :].reshape((trajectory.n_frames, -1))
    return coordinates


def writeFilenameExtractedCoordinates(filename, params, pathFolder, constants, topology, indexes=None):
    """
        Process the coordinates of a trajectory

    """
    if params.dihedrals:
        coords = calculateDihedrals(filename, params, topology)
        if params.dihedrals_projection:
            coords = projectDihedrals(coords)
        outputFilename = getOutputFilename(constants.extractedTrajectoryFolder, filename,
                                           constants.baseExtractedTrajectoryName)
        writeToFile(coords, outputFilename % pathFolder)
        return
    ext = utilities.getFileExtension(filename)
    if ext == ".pdb":
        allCoordinates = loadAllResnameAtomsInPdb(filename, params)
        if params.writeLigandTrajectory:
            outputFilename = os.path.join(pathFolder, constants.ligandTrajectoryBasename % extractFilenumber(filename))
            with open(outputFilename, 'w') as f:
                f.write("\nENDMDL\n".join(allCoordinates))
        if params.protein_CA:
            coords = getLigandAlphaCarbonsCoords(allCoordinates, params.lig_resname)
        elif params.sidechains:
            coords = getLigandAlphaCarbonsCoords(allCoordinates, params.lig_resname, sidechains=params.sidechains)
        else:
            if params.com:
                coords = getPDBCOM(allCoordinates, params.lig_resname)
            elif params.contact_map:
                coords = np.array(extractContactMapCoordinatesPDB(allCoordinates, params))
                coords = utils.contactMap(coords, int(coords.shape[1]/3))
            else:
                coords = getAtomCoord(allCoordinates, params.lig_resname, params.atomIds)
    elif ext in MDTRAJ_FORMATS:
        if params.contact_map:
            coords = contactMapNonPDB(filename, params, topology, indexes)
        else:
            coords = extractCoordinatesXTCFile(filename, params, topology, indexes)
    else:
        raise ValueError("Unrecongnized file extension for %s" % filename)

    outputFilename = getOutputFilename(constants.extractedTrajectoryFolder, filename,
                                       constants.baseExtractedTrajectoryName)
    writeToFile(coords, outputFilename % pathFolder)


def writeFilenamesExtractedCoordinates(pathFolder, params, constants, pool=None):
    if not os.path.exists(constants.extractedTrajectoryFolder % pathFolder):
        os.makedirs(constants.extractedTrajectoryFolder % pathFolder)

    originalPDBfiles = glob.glob(os.path.join(pathFolder, '*traj*.*'))
    ext = os.path.splitext(originalPDBfiles[0])[1]
    if ext in MDTRAJ_FORMATS:
        if params.topology is None:
            raise ValueError("Necessary topology not provided!")
        # get topology for the first trajectory
        top_file = params.topology.getTopologyFile(0, 1)
        if params.contact_map:
            indexes = extractIndexesTopology_CM(top_file, params.lig_resname, params.cm_mode, params.extra_atoms)
        else:
            indexes = extractIndexesTopology(top_file, params.lig_resname, params.atomIds, params.protein_CA, params.sidechains)
    else:
        indexes = None
    workers = []
    for filename in originalPDBfiles:
        if params.topology is not None:
            epoch, traj_num = get_epoch_traj_num(filename)
            topology_file = params.topology.getTopologyFile(epoch, traj_num)
        else:
            topology_file = None
        if pool is None:
            # serial version
            writeFilenameExtractedCoordinates(filename, params, pathFolder, constants, topology_file, indexes=indexes)
        else:
            # multiprocessor version
            workers.append(pool.apply_async(writeFilenameExtractedCoordinates, args=(filename, params, pathFolder, constants, topology_file, indexes)))
    for w in workers:
        w.get()


def parseResname(atom_Ids, lig_resname, CM, CM_mode, dihedrals):
    if dihedrals:
        return ""
    if atom_Ids is not None and len(atom_Ids) > 0:
        differentResnames = {atomId.split(":")[-1] for atomId in atom_Ids}
        if len(differentResnames) > 1:
            sys.exit("Error! Different resnames provided in atomIds!")
        elif len(differentResnames) == 1:
            extractedResname = differentResnames.pop()

    if CM:
        if CM_mode == "p-lig":
            if lig_resname == "":
                sys.exit("Ligand resname should be provided for the protein-ligand contact map")
            else:
                return lig_resname
        else:
            return ""
    if (atom_Ids is None or len(atom_Ids) == 0) and lig_resname == "":
        sys.exit("Either resname or atomId should be provided")
    elif lig_resname == "":
        lig_resname = extractedResname  # the atom Id last element is the resname
    elif atom_Ids is not None and len(atom_Ids) > 0:
        if extractedResname != lig_resname:
            sys.exit("Residue name in resname and atomId do not match!")
    return lig_resname


def buildFullTrajectory(steps, trajectory, numtotalSteps, inputTrajectory):
    completeTrajectory = []
    counter = 0
    if len(trajectory) > 0:
        sthWrongInTraj = False
        for i in range(len(steps) - 1):
            repeated = steps[i+1, 0] - steps[i, 0]
            for _ in range(repeated):
                try:
                    snapshot = trajectory[steps[i, 1]].split()
                except IndexError:
                    print("sth wrong in trajectory %s. This is likely to disagreement between report and trajectory files. Please, fix it manually" % inputTrajectory)
                    sthWrongInTraj = True
                    break
                snapshot[0] = str(counter)
                snapshot = ' '.join(snapshot)
                completeTrajectory.append(snapshot)
                counter += 1

        if sthWrongInTraj:
            return completeTrajectory

        if numtotalSteps == 0:
            iterations = list(range(1))
        else:
            iterations = list(range(numtotalSteps + 1 - counter))

        for i in iterations:
            snapshot = trajectory[-1].split()
            snapshot[0] = str(counter)
            snapshot = ' '.join(snapshot)
            completeTrajectory.append(snapshot)
            counter += 1

    return completeTrajectory


def repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants, numtotalSteps):
    extractedTrajFolder, trajFilename = os.path.split(inputTrajectory)
    trajectoryNumber = re.sub(r'\.dat$', '', trajFilename)
    trajectoryNumber = re.sub(constants.baseExtractedTrajectoryName, '', trajectoryNumber)

    origDataFolder = re.sub(constants.extractedTrajectoryFolder % "", "", extractedTrajFolder)
    try:
        reportFile = glob.glob(os.path.join(origDataFolder, constants.reportName + trajectoryNumber))[0]
    except IndexError:
        sys.exit("Couldn't find file that matches: %s" % os.path.join(origDataFolder, constants.reportName + trajectoryNumber))

    with open(inputTrajectory) as f:
        trajectory = f.read().splitlines()

    acceptedSteps = np.loadtxt(reportFile, dtype='int', comments='#', usecols=(1, 2))
    if len(acceptedSteps.shape) < 2:
        acceptedSteps = acceptedSteps[np.newaxis, :]

    fullTrajectory = buildFullTrajectory(acceptedSteps, trajectory, numtotalSteps, inputTrajectory)

    if len(fullTrajectory) > 0:
        outputFilename = os.path.join(constants.outputTrajectoryFolder % origDataFolder, constants.baseExtractedTrajectoryName + trajectoryNumber + '.dat')
        with open(outputFilename, "w") as outputFile:
            for snapshot in fullTrajectory:
                outputFile.write("%s\n" % snapshot)


def repeatExtractedSnapshotsInFolder(folder_name, constants, numtotalSteps, pool=None):
    inputTrajectoryFolder = constants.extractedTrajectoryFolder % folder_name
    outputTrajectoryFolder = constants.outputTrajectoryFolder % folder_name

    if not os.path.exists(outputTrajectoryFolder):
        os.makedirs(outputTrajectoryFolder)

    inputTrajectories = glob.glob(os.path.join(inputTrajectoryFolder, constants.baseExtractedTrajectoryName + '*'))
    workers = []
    for inputTrajectory in inputTrajectories:
        if pool is None:
            # serial version
            repeatExtractedSnapshotsInTrajectory(inputTrajectory, constants, numtotalSteps)
        else:
            # multiprocessor version
            workers.append(pool.apply_async(repeatExtractedSnapshotsInTrajectory, args=(inputTrajectory, constants, numtotalSteps)))
    for w in workers:
        w.get()


def makeGatheredTrajsFolder(constants):
    if not os.path.exists(constants.gatherTrajsFolder):
        os.makedirs(constants.gatherTrajsFolder)
    if not os.path.exists(constants.gatherNonRepeatedFolder):
        os.makedirs(constants.gatherNonRepeatedFolder)


def copyTrajectories(traj_names, destFolderTempletized, folderName, setNumber=0, epochNum=None):
    for inputTrajectory in traj_names:
        trajectoryNumber = extractFilenumber(os.path.split(inputTrajectory)[1])
        if folderName != ".":  # if not sequential
            setNumber = folderName
        if epochNum is not None:
            setNumber = epochNum
        shutil.copyfile(inputTrajectory, destFolderTempletized % (setNumber, trajectoryNumber))


def gatherTrajs(constants, folder_name, setNumber, non_Repeat, epochNum=None):
    nonRepeatedTrajs = glob.glob(os.path.join(constants.extractedTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*"))
    copyTrajectories(nonRepeatedTrajs, constants.gatherNonRepeatedTrajsFilename, folder_name, setNumber, epochNum=epochNum)
    if not non_Repeat:
        # copy the repeated coordinates to the allTrajs folder
        trajectoriesFilenames = os.path.join(constants.outputTrajectoryFolder % folder_name, constants.baseExtractedTrajectoryName + "*")
        trajectories = glob.glob(trajectoriesFilenames)
        copyTrajectories(trajectories, constants.gatherTrajsFilename, folder_name, setNumber)
    else:
        # if we ask to not repeat trajectories, copy the non-repeated to the
        # allTrajs folder
        copyTrajectories(nonRepeatedTrajs, constants.gatherTrajsFilename, folder_name, setNumber, epochNum=epochNum)


def extractSidechainIndexes_prody(traj, ligand_resname, topology=None):
    if not PRODY:
        raise utilities.UnsatisfiedDependencyException("Prody module not found, will not be able to extract sidechain coordinates")
    atoms = pd.parsePDB(traj)
    sidechains = atoms.select("protein within 5 of resname {}".format(ligand_resname))
    return [atom.getIndex() for atom in sidechains]


def extractSidechainIndexes_mdtraj(traj, lig_resname, topology=None):
    atoms = md.load(traj, top=topology)
    ligand_indices = atoms.top.select("resname '{lig}'".format(lig=lig_resname))
    water_indices = set(atoms.top.select("not protein or not resname '{lig}'".format(lig=lig_resname)))
    # the distance is specified in nm
    sidechains = md.compute_neighbors(atoms, 0.5, ligand_indices)
    sidechains_trajs = []
    for _, sidechain in enumerate(sidechains):
        sidechains_trajs.extend(list(set(sidechain.tolist())-water_indices))
    return sidechains_trajs


def extractSidechainIndexes(params, pool=None):
    trajs = glob.glob(params.sidechain_folder)
    sidechains_trajs = []
    workers = []
    for traj in trajs:
        ext = utilities.getFileExtension(traj)

        if ext == ".pdb":
            if PRODY:
                if pool is None:
                    sidechains_trajs.extend(extractSidechainIndexes_prody(traj, params.lig_resname))
                else:
                    workers.append(pool.apply_async(extractSidechainIndexes_prody, args=(traj, params.lig_resname)))
            else:
                if pool is None:
                    sidechains_trajs.extend(extractSidechainIndexes_mdtraj(traj, params.lig_resname))
                else:
                    workers.append(pool.apply_async(extractSidechainIndexes_mdtraj, args=(traj, params.lig_resname)))
        elif ext in MDTRAJ_FORMATS:
            epoch, traj_num = get_epoch_traj_num(traj)
            if pool is None:
                sidechains_trajs.extend(extractSidechainIndexes_mdtraj(traj, params.lig_resname, topology=params.topology.getTopologyFile(epoch, traj_num)))
            else:
                workers.append(pool.apply_async(extractSidechainIndexes_mdtraj(traj, params.lig_resname, params.topology)))
        else:
            raise ValueError("Unrecongnized file extension for %s" % traj)
    for w in workers:
        sidechains_trajs.extend(w.get())
    return list(set(sidechains_trajs))


def get_epoch_traj_num(filename):
    # assumes trajectories come from an Adaptive simulation
    path, traj_name = os.path.split(filename)
    try:
        epoch = int(os.path.split(path)[-1])
    except ValueError:
        # if for some reason epoch number can't be inferred, assume first
        # epoch
        epoch = 0
    try:
        traj_num = utilities.getTrajNum(traj_name)
    except ValueError:
        # if for some reason trajectory number can't be inferred, assume
        # first trajectory
        traj_num = 1
    return epoch, traj_num


def getTopologyObject(topology_file):
    ext = utilities.getFileExtension(topology_file)
    if ext == ".pdb":
        return TopologyCompat(topology_file)
    elif ext == ".pkl":
        return utilities.readClusteringObject(topology_file)
    else:
        raise ValueError("The topology parameter needs to be the path to a pickled Topology object or a pdb!")


def main(folder_name=".", atom_Ids="", lig_resname="", numtotalSteps=0, enforceSequential_run=0, writeLigandTrajectory=True, setNumber=0, protein_CA=0, non_Repeat=False, nProcessors=None, parallelize=True, topology=None, sidechains=False, sidechain_folder=".", cm=False, use_extra_atoms=False, CM_mode="p-lig", calc_dihedrals=False, dihedrals_projection=False):
    params = ParamsHandler(folder_name, atom_Ids, lig_resname, numtotalSteps, enforceSequential_run, writeLigandTrajectory, setNumber, protein_CA, non_Repeat, nProcessors, parallelize, topology, sidechains, sidechain_folder, cm, use_extra_atoms, CM_mode, calc_dihedrals, dihedrals_projection)
    constants = Constants()

    if params.topology is not None:
        params.topology = getTopologyObject(params.topology)

    params.lig_resname = parseResname(params.atomIds, params.lig_resname, params.contact_map, params.cm_mode, params.dihedrals)

    folderWithTrajs = params.folder_name

    makeGatheredTrajsFolder(constants)

    if params.enforceSequential_run:
        folders = ["."]
    else:
        allFolders = os.listdir(folderWithTrajs)
        folders = [epoch for epoch in allFolders if epoch.isdigit()]
        if len(folders) == 0:
            folders = ["."]

    # if multiprocess is not available, turn off parallelization
    params.parallelize &= PARALELLIZATION

    if params.parallelize:
        if params.nProcessors is None:
            params.nProcessors = utilities.getCpuCount()
        params.nProcessors = max(1, params.nProcessors)

        print("Running extractCoords with %d cores" % (params.nProcessors))
        pool = mp.Pool(params.nProcessors)
    else:
        pool = None

    params.sidechains = extractSidechainIndexes(params, pool=pool) if params.sidechains else []

    for folder_it in folders:
        pathFolder = os.path.join(folderWithTrajs, folder_it)
        print("Extracting coords from folder %s" % folder_it)
        ligand_trajs_folder = os.path.join(pathFolder, constants.ligandTrajectoryFolder)
        if params.writeLigandTrajectory and not os.path.exists(ligand_trajs_folder):
            os.makedirs(ligand_trajs_folder)
        writeFilenamesExtractedCoordinates(pathFolder, params, constants, pool=pool)
        if not params.non_Repeat:
            print("Repeating snapshots from folder %s" % folder_it)
            repeatExtractedSnapshotsInFolder(pathFolder, constants, params.numtotalSteps, pool=None)
        print("Gathering trajs in %s" % constants.gatherTrajsFolder)
        gatherTrajs(constants, folder_it, params.setNumber, params.non_Repeat)


if __name__ == "__main__":
    folder, atomIds, resname, proteinCA, enforceSequential, writeLigandTraj, totalSteps, setNum, nonRepeat, n_processors, top, side_chains, sideChain_folder, serial, contact_map, extra_atoms, cm_mode, dihedral_angles, dihedrals_proj = parseArguments()
    main(folder, atomIds, resname, totalSteps, enforceSequential, writeLigandTraj, setNum, proteinCA, nonRepeat, n_processors, topology=top, sidechains=side_chains, sidechain_folder=sideChain_folder, parallelize=(not serial), cm=contact_map, use_extra_atoms=extra_atoms, CM_mode=cm_mode, calc_dihedrals=dihedral_angles, dihedrals_projection=dihedrals_proj)
