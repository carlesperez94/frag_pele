from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import time
import shutil
import unittest
import AdaptivePELE.adaptiveSampling as adaptiveSampling
from AdaptivePELE.utilities import utilities


class TestMD_CUDA(unittest.TestCase):

    def check_succesful_simulation(self, output, epochs, nTrajs):
        for epoch in range(epochs):
            self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "summary.txt")))
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "trajectory*"))), nTrajs)
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "report*"))), nTrajs)
        self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "object.pkl")))

    def testOpenMM3ptb_2replicas_CUDA(self):
        output_path = "tests/data/openmm_3ptb_2replica_CUDA"
        controlFile = "tests/data/templetized_controlFile_3ptb_CUDA_2replicas_md.conf"

        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # ensure all replicas have enough time to check the output before
        # deleting it
        time.sleep(10)
        # cleanup
        utilities.cleanup(output_path)

    def testOpenMM1ab1_CUDA(self):
        output_path = "tests/data/openmm_1ab1_CUDA"
        controlFile = "tests/data/templetized_controlFile_1ab1_md_CUDA.conf"
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # ensure all replicas have enough time to check the output before
        # deleting it
        time.sleep(10)
        # cleanup
        utilities.cleanup(output_path)
