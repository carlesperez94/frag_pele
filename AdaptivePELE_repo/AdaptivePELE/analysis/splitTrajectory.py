from __future__ import absolute_import, division, print_function, unicode_literals
from AdaptivePELE.utilities import utilities
import argparse
import os


def parseArguments():
    desc = "Program that writes a trajectory into separate pdbs."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("files", nargs='+', help="Trajectory files to split")
    parser.add_argument("-o", type=str, default=".", help="Output dir")
    args = parser.parse_args()
    return args.files, args.o


files, outputDir = parseArguments()
utilities.makeFolder(outputDir)

for f in files:
    name = os.path.split(f)[-1]
    templateName = os.path.join(outputDir, name[:-4] + "_%d.pdb")
    snapshots = utilities.getSnapshots(f)
    print(len(snapshots))
    for i, snapshot in enumerate(snapshots):
        print(templateName % i)
        with open(templateName % i, 'w') as of:
            of.write(snapshot)
