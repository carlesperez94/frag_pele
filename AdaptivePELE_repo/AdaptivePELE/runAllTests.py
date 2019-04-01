from __future__ import absolute_import, division, print_function, unicode_literals
import unittest
import argparse
from AdaptivePELE.tests import testSpawning as tSpawning
from AdaptivePELE.tests import testAtomset as tAtomset
from AdaptivePELE.tests import testClustering as tClustering
from AdaptivePELE.tests import testAdaptiveSampling as tAdaptive
from AdaptivePELE.tests import testThresholdcalculator as tThreshold
from AdaptivePELE.tests import testDensityCalculator as tDensity
from AdaptivePELE.tests import testMD as tMD
from AdaptivePELE.tests import testMD_CUDA as tMD_CUDA
from AdaptivePELE.tests import testReporter as tR


def parse_args():
    desc = ("Run testing suite. Possible options are:\na  -- Run all tests\n"
            "at -- Run atomset tests\ns  -- Run spawning tests\nth -- Run threshold "
            "calculator tests\nd  -- Run density tests\nc  -- Run clustering tests\n"
            "Ad -- Run adaptive integration tests\nMD -- Run adaptive MD tests\nMD_CUDA"
            " -- Run adaptive MD tests with CUDA\nR -- Run reporter tests\n")
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--run", default=None, nargs="*", help="Tests to run")
    parser.add_argument("--exclude", default=[], nargs="*", help="Tests to exclude")
    args = parser.parse_args()

    return args.run, args.exclude


def main(run, exclude):
    testSuite = unittest.TestSuite()
    if run is None:
        run = ["at", "s", "th", "d", "c", "Ad", "MD", "MD_CUDA", "R"]
    to_run = set(run)-set(exclude)

    if "at" in to_run or "a" in to_run:
        print("Will run atomset tests")
        testSuite.addTest(unittest.makeSuite(tAtomset.atomsetTest))
    if "s" in to_run or "a" in to_run:
        print("Will run spawning tests")
        testSuite.addTest(unittest.makeSuite(tSpawning.TestSpawningCalculator))
    if "th" in to_run or "a" in to_run:
        print("Will run threshold tests")
        testSuite.addTest(unittest.makeSuite(tThreshold.thresholdCalculatorTest))
    if "d" in to_run or "a" in to_run:
        print("Will run denstity tests")
        testSuite.addTest(unittest.makeSuite(tDensity.densityCalculatorTest))
    if "c" in to_run or "a" in to_run:
        print("Will run clustering tests")
        testSuite.addTest(unittest.makeSuite(tClustering.clusteringTest))
    if "Ad" in to_run or "a" in to_run:
        print("Will run integration tests")
        testSuite.addTest(unittest.makeSuite(tAdaptive.TestadaptiveSampling))
    if "MD" in to_run or "a" in to_run:
        print("Will run integration tests with md")
        testSuite.addTest(unittest.makeSuite(tMD.TestMD))
    if "MD_CUDA" in to_run or "a" in to_run:
        print("Will run integration tests with md in CUDA")
        testSuite.addTest(unittest.makeSuite(tMD_CUDA.TestMD_CUDA))
    if "R" in to_run or "a" in to_run:
        print("Will run repoter tests for OpenMM")
        testSuite.addTest(unittest.makeSuite(tR.TestReporter))

    runner = unittest.TextTestRunner()
    runner.run(testSuite)

if __name__ == "__main__":
    run_list, exclude_list = parse_args()
    main(run_list, exclude_list)
