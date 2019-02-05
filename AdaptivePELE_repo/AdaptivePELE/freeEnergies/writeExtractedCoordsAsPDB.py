"""
    Program that writes extracted coordinates as pdb files.

    Execution: python generateFiles.py outputDir rawData/traj_*

"""
from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.freeEnergies import computeDeltaG as dg
import sys
import numpy as np
import os


def getNum(trajFilename):
    left = trajFilename.rfind("_")
    return int(trajFilename[left+1:-4])

outputFolder = sys.argv[1]
files = sys.argv[2:]

print("files", files)

for filename in files:
    content = np.loadtxt(filename)

    length = content.shape[0]

    dummyContent = np.ones((length, 1))

    paddedContent = np.delete(np.hstack((content, dummyContent)), 0, axis=1)

    outputFilename = os.path.join(outputFolder, os.path.split(filename)[1])

    dg.writePDB(paddedContent, outputFilename)

    with open(outputFilename, 'r') as f:
        content = "MODEL\n" + "ENDMDL\nMODEL\n".join(f.readlines()) + "ENDMDL\n"
    with open(outputFilename, 'w') as f:
        f.write(content)
