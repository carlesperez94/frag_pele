from __future__ import absolute_import, division, print_function, unicode_literals
import os
import socket
machine = socket.getfqdn()

print("MACHINE", machine)
if "bsccv" in machine:
    PELE_EXECUTABLE = "/data/EAPM/PELE/PELE++/bin/rev12360/Pele_rev12360_mpi"
    DATA_FOLDER = "/data/EAPM/PELE/PELE++/data/rev12360/Data"
    DOCUMENTS_FOLDER = "/data/EAPM/PELE/PELE++/Documents/rev12360"

    PYTHON = "/data2/apps/PYTHON/2.7.5/bin/python2.7"

elif "mn.bsc" in machine:
    PELE_EXECUTABLE = "/gpfs/projects/bsc72/PELE++/nord/rev090518/bin/PELE-1.5_mpi"
    DATA_FOLDER = "/gpfs/projects/bsc72/PELE++/nord/rev090518/Data"
    DOCUMENTS_FOLDER = "/gpfs/projects/bsc72/PELE++/nord/rev090518/Documents"
    PYTHON = "python"


elif "bsc.mn" in machine:
    PELE_EXECUTABLE = "/gpfs/projects/bsc72/PELE++/mniv/rev090518/bin/PELE-1.5_mpi"
    DATA_FOLDER = "/gpfs/projects/bsc72/PELE++/mniv/rev090518/Data"
    DOCUMENTS_FOLDER = "/gpfs/projects/bsc72/PELE++/mniv/rev090518/Documents"

elif "bullx" in machine:
    # this values are not correct for the minoTauro hardware, just leaving it
    # here as a placeholder
    PELE_EXECUTABLE = "/gpfs/projects/bsc72/PELE++/nord/rev090518/bin/PELE-1.5_mpi"
    DATA_FOLDER = "/gpfs/projects/bsc72/PELE++/nord/rev090518/Data"
    DOCUMENTS_FOLDER = "/gpfs/projects/bsc72/PELE++/nord/rev090518/Documents"

elif machine == "bscls309":
    PELE_EXECUTABLE = "/home/jgilaber/PELE/PELE-1.5/bin/PELE-1.5_mpi"
    DATA_FOLDER = "/home/jgilaber/PELE/PELE-1.5/Data"
    DOCUMENTS_FOLDER = "/home/jgilaber/PELE/PELE-1.5/Documents"

else:
    PELE_EXECUTABLE = None
    DATA_FOLDER = None
    DOCUMENTS_FOLDER = None


inputFileTemplate = "{ \"files\" : [ { \"path\" : \"%s\" } ] }"
trajectoryBasename = "*traj*"


class AmberTemplates:
    forcefields = {"ff99SB": "oldff/leaprc.ff99SB", "ff14SB": "leaprc.protein.ff14SB"}
    antechamberTemplate = "antechamber -i $LIGAND -fi pdb -o $OUTPUT -fo mol2 -c bcc -pf y -nc $CHARGE"
    parmchk2Template = "parmchk2 -i $MOL2 -f mol2 -o $OUTPUT"
    tleapTemplate = "source $FORCEFIELD\n" \
                    "source leaprc.gaff\n" \
                    "source leaprc.water.tip3p\n" \
                    "$MODIFIED_RES " \
                    "$RESNAME = loadmol2 $MOL2\n" \
                    "loadamberparams $FRCMOD\n" \
                    "$DUM " \
                    "COMPLX = loadpdb $COMPLEX\n" \
                    "$BONDS " \
                    "addions COMPLX Cl- 0\n" \
                    "addions COMPLX Na+ 0\n" \
                    "solvatebox COMPLX TIP3PBOX $BOXSIZE\n" \
                    "saveamberparm COMPLX $PRMTOP $INPCRD\n" \
                    "savepdb COMPLX $SOLVATED_PDB\n" \
                    "quit"
    tleapTemplatenoLigand = "source $FORCEFIELD\n" \
                    "source leaprc.gaff\n" \
                    "source leaprc.water.tip3p\n" \
                    "$MODIFIED_RES " \
                    "$DUM " \
                    "COMPLX = loadpdb $COMPLEX\n" \
                    "$BONDS " \
                    "addions COMPLX Cl- 0\n" \
                    "addions COMPLX Na+ 0\n" \
                    "solvatebox COMPLX TIP3PBOX $BOXSIZE\n" \
                    "saveamberparm COMPLX $PRMTOP $INPCRD\n" \
                    "savepdb COMPLX $SOLVATED_PDB\n" \
                    "quit"
    DUM_atom = "DUM"
    DUM_res = "DUM"
    DUM_prep = "  0  0  0\n" \
               "\n" \
               "------%s--------------\n" \
               "%s\n" \
               "%s   INT  0\n" \
               "CHANGE OMIT DU  BEG\n" \
               " 0.0\n" \
               "   1 DUMM  DU  M  0  -1  -2  0.000     0.000     0.000    0.000\n" \
               "   2 DUMM  DU  M  1   0  -1  1.0000    0.0000    0.0000   0.000\n" \
               "   3 DUMM  DU  M  2   1   0  1.0000   90.0000    0.0000   0.000\n" \
               "   4 %s    C  E  0.00 0.00 0.00  0.00\n" \
               "   5 %s    C  E  0.00 0.00 0.00  0.00\n" \
               "   6 %s    C  E  0.00 0.00 0.00  0.00\n" \
               "\n" \
               "\n" \
               "DONE\n" \
               "STOP\n" \
               "\n" % (DUM_res, DUM_res, DUM_res, DUM_atom, DUM_atom+"B", DUM_atom+"T")
    DUM_frcmod = "invented MM atom\n" \
                 "MASS\n" \
                 "%s     0.00     0.00\n" \
                 "%s     0.00     0.00\n" \
                 "%s     0.00     0.00\n" \
                 "\n" \
                 "NONB\n" \
                 "  %s        0.00     0.00\n" \
                 "  %s        0.00     0.00\n" \
                 "  %s        0.00     0.00\n" % (DUM_atom, DUM_atom+"B", DUM_atom+"T", DUM_atom, DUM_atom+"B", DUM_atom+"T")

    trajectoryTemplate = "trajectory_%d.%s"
    CheckPointReporterTemplate = "checkpoint_%d.chk"


class OutputPathConstants():
    """
        Class with constants that depend on the outputPath
    """
    def __init__(self, outputPath):
        self.originalControlFile = ""
        self.epochOutputPathTempletized = ""
        self.clusteringOutputDir = ""
        self.clusteringOutputObject = ""
        self.equilibrationDir = ""
        self.tmpInitialStructuresTemplate = ""
        self.tmpControlFilename = ""
        self.tmpInitialStructuresEquilibrationTemplate = ""
        self.tmpControlFilenameEqulibration = ""
        self.topologies = ""
        self.allTrajsPath = ""
        self.MSMObjectEpoch = ""
        self.buildConstants(outputPath)

    def buildConstants(self, outputPath):
        self.buildOutputPathConstants(outputPath)

        self.tmpFolder = "tmp_" + outputPath.replace("/", "_")

        self.buildTmpFolderConstants(self.tmpFolder)

    def buildOutputPathConstants(self, outputPath):
        self.originalControlFile = os.path.join(outputPath, "originalControlFile.conf")
        self.epochOutputPathTempletized = os.path.join(outputPath, "%d")
        self.clusteringOutputDir = os.path.join(self.epochOutputPathTempletized, "clustering")
        self.clusteringOutputObject = os.path.join(self.clusteringOutputDir, "object.pkl")
        self.MSMObjectEpoch = os.path.join(self.epochOutputPathTempletized, "MSM_object.pkl")
        self.topologies = os.path.join(outputPath, "topologies")
        self.equilibrationDir = os.path.join(outputPath, "equilibration")
        self.allTrajsPath = os.path.join(outputPath, "allTrajs")

    def buildTmpFolderConstants(self, tmpFolder):
        self.tmpInitialStructuresTemplate = tmpFolder+"/initial_%d_%d.pdb"
        self.tmpInitialStructuresEquilibrationTemplate = tmpFolder+"/initial_equilibration_%d.pdb"
        self.tmpControlFilename = tmpFolder+"/controlFile%d.conf"
        self.tmpControlFilenameEqulibration = tmpFolder+"/controlFile_equilibration_%d.conf"

md_supported_formats = set(["xtc", "dcd"])
formats_md_string = ", ".join(md_supported_formats)
