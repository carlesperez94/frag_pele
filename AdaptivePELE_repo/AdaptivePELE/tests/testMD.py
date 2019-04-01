from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import shutil
import unittest
import numpy as np
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj as md
import AdaptivePELE.adaptiveSampling as adaptiveSampling
from AdaptivePELE.simulation.openmm_simulations import XTCReporter


class TestMD(unittest.TestCase):

    def check_succesful_simulation(self, output, epochs, nTrajs):
        for epoch in range(epochs):
            self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "summary.txt")))
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "trajectory*"))), nTrajs)
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "report*"))), nTrajs)
        self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "object.pkl")))

    def testOpenMM3ptb(self):
        output_path = "tests/data/openmm_3ptb"
        controlFile = "tests/data/templetized_controlFile_3ptb_md.conf"

        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testOpenMM1ab1(self):
        output_path = "tests/data/openmm_1ab1"
        controlFile = "tests/data/templetized_controlFile_1ab1_md.conf"
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testRestartAt0(self):
        output_path = "tests/data/openmm_restart_0"
        controlFile = "tests/data/templetized_controlFile_restart_0_md.conf"
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        shutil.copytree("tests/data/restart_0", output_path)
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testRestartAt1(self):
        output_path = "tests/data/openmm_restart_1"
        controlFile = "tests/data/templetized_controlFile_restart_1_md.conf"
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        shutil.copytree("tests/data/restart_1", output_path)
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)
