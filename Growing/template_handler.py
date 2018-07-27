import sys
import re
import os
import logging

# Getting the name of the module for the log system
logger = logging.getLogger(__name__)

# Definition of reggex patterns
ATOM_PATTERN = "\s+(\d+)\s+(\d+)\s+(\w)\s+(\w*)\s+(\w{4,})\s+(\d*)\s+(-?[0-9]*\.[-]?[0-9]*)\s+(-?[0-9]*\.[0-9]*)\s+(-?[0-9]*\.[0-9]*)"
NBON_PATTERN = "\s+(\d+)\s+(\d+\.\d{4})\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
BOND_PATTERN = "\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)"
WRITE_NBON_PATTERN = " {:5d}   {:3.4f}   {:3.4f}  {: 3.6f}   {:3.4f}   {:3.4f}   {:3.9f}  {: 3.9f}\n"
WRITE_BOND_PATTERN = " {:5d} {:5d}   {:5.3f} {: 2.3f}\n"

# Class definitions
def