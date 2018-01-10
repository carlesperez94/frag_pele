import os
import shutil
import numpy as np
# import pickle
import cPickle as pickle
from AdaptivePELE.atomset import RMSDCalculator
import AdaptivePELE.atomset.atomset as atomset


def cleanup(tmpFolder):
    """
        Remove folder if exists

        :param tmpFolder: Folder to remove
        :type tmpFolder: str
    """
    if os.path.exists(tmpFolder):
        shutil.rmtree(tmpFolder)


def makeFolder(outputDir):
    """
        Makes folder

        :param outputDir: Folder filename
        :type outputDir: str
    """
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)


def getSnapshots(trajectoryFile, verbose=False):
    """
        Gets the snapshots

        :param trajectoryFile: Trajectory filename
        :type trajectoryFile: str
        :param verbose: Add verbose to snapshots
        :type verbose: bool

        :returns: Snapshots with information
        :rtype: str
    """
    inputFile = open(trajectoryFile, "r")
    inputFileContent = inputFile.read()

    snapshots = inputFileContent.split("ENDMDL")[:-1]
    if not verbose:
        return snapshots

    remarkInfo = "REMARK 000 File created using PELE++\nREMARK source            : %s\nREMARK original model nr : %d\nREMARK First snapshot is 1, not 0 (as opposed to report)\n"
    snapshotsWithInfo = [remarkInfo % (trajectoryFile, i+1)+snapshot for i, snapshot in enumerate(snapshots)]
    return snapshotsWithInfo


def getTrajNum(trajFilename):
    """
        Gets the trajectory number

        :param trajFilename: Trajectory filename
        :type trajFilename: str

        :returns: Trajectory number
        :rtype: int
    """
    return int(trajFilename.split("_")[-1][:-4])


def calculateContactMapEigen(contactMap):
    """
        Calculates eigenvectors and values of an extended contact map

        :param contactMap: Contact map
        :type contactMap: np.array

        :returns: (eiv, eic) 
        **eiv** (numpy.ndarray) eigenvalues
        **eic** (numpy.ndarray) eigenvectors
    """
    nLig, nCA = contactMap.shape
    extendedCM = np.zeros((nLig+nCA, nLig+nCA))
    extendedCM[nCA:, :nCA] = contactMap
    extendedCM[:nCA, nCA:] = contactMap.T
    assert (extendedCM == extendedCM.T).all(), "Extended ContactMap not symmetric"
    eiv, eic = np.linalg.eigh(extendedCM)
    return eiv, eic


def assertSymmetriesDict(symmetries, PDB):
    """
        Asserts the symmetry list in a PDB

        :param symmetries: List of symmetry groups
        :type symmetries: list of str
        :param PDB: PDB object to check symmetry list against
        :type PDB: :py:class:`.PDB`

        :raise AssertionError: If an atom is not found in the structure
    """
    for group in symmetries:
        for key in group:
            assert key in PDB.atoms, "Symmetry atom %s not found in initial structure" % key
    if symmetries:
        print "Symmetry dictionary correctly defined!"


def getRMSD(traj, nativePDB, resname, symmetries):
    """
        Computes the RMSD of a trajectory, given a native and symmetries

        :param traj: Trajecotry filename
        :type traj: str
        :param nativePDB:  Native PDB object
        :type native PDB: :py:class:`.PDB`
        :param resname: Resname to compute its RMSD
        :type resname: str
        :param symmetries: Symmetries dictionary list with independent symmetry groups
        :type symmetries: list of dict

        :return: rmsds
        :rtype: np.array
    """

    snapshots = getSnapshots(traj)
    rmsds = np.zeros(len(snapshots))
    RMSDCalc = RMSDCalculator.RMSDCalculator(symmetries)
    for i, snapshot in enumerate(snapshots):
        snapshotPDB = atomset.PDB()
        snapshotPDB.initialise(snapshot, resname=resname)

        rmsds[i] = RMSDCalc.computeRMSD(nativePDB, snapshotPDB)

    return rmsds


def readClusteringObject(clusteringObjectPath):
    """
        Reads and returns a clustering object

        :param clusteringObjectPath: Clustering object path
        :type clusteringObjectPath: str

        :raise EOFError: If the object is empty

        :return: clusteringObject
        :rtype: :py:class:`.Clustering`
    """
    with open(clusteringObjectPath, 'rb') as f:
        try:
            return pickle.load(f)
        except EOFError:
            raise EOFError("Empty clustering object!")
