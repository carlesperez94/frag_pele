import sys
sys.settrace
import os
import socket
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import argparse


DATAFILE = 'pmf_xyz.dat'

data = np.loadtxt(DATAFILE)

#x =  data[:,0]

templateLine = "HETATM%s  C1  CLT L 502    %s%s%s  0.75%s           C"

#x0 = data[0][0]
#a = "%.3f"%x0
#print a.rjust(8)
#print templateLine%(x0,x0,x0,x0,x0)
for i,line in enumerate(data):
    number = str(i).rjust(5)
    x = ("%.3f"%line[0]).rjust(8)
    y = ("%.3f"%line[1]).rjust(8)
    z = ("%.3f"%line[2]).rjust(8)
    g = ("%.3f"%line[3]).rjust(8)

    print templateLine%(number, x, y, z, g)
