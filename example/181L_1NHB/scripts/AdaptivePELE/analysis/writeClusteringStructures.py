from AdaptivePELE.utilities import clusteringUtilities
import argparse
import math


def parseArgs():
    parser = argparse.ArgumentParser(description="Write the requested cluster "
                                     "structures from a clustering object")
    parser.add_argument('clObject', type=str)
    parser.add_argument('outputPath', type=str,
                        help="Path where to write the structures, including "
                        "name of the files, i.e output/path/cluster.pdb")
    parser.add_argument("structures", nargs='*', type=list, default=None,
                        help="Structures to write")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Only print those structures with mathcing threshold")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parseArgs()
    if not args.threshold is None:
        condition = lambda x: abs(x.threshold-args.threshold) < 0.01
    else:
        condition = None
    clusteringUtilities.writeStructures(args.clObject, args.structures, condition,
                                        args.outputPath)
