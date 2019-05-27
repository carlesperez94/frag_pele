import pyemma
import pickle
import numpy as np


def inout(i):
    print "in", C[:,i].sum() - C[i,i]
    print "out", C[i,:].sum() - C[i,i]
    print "self", C[i,i]

def updown(i, down):
    print "from states", C[down,i].sum()
    print "to states", C[i,down].sum()
    print "***"
    print "from conjugate", C[:,i].sum() - C[down,i].sum() - C[i,i]
    print "to conjugate", C[i,:].sum() - C[i,down].sum() - C[i,i]

with open("MSM_object_0.pkl", 'rb') as f:
    msmObj = pickle.load(f)

C = msmObj.count_matrix_full

notFinish = True
while notFinish:
    i = int(raw_input('Enter state: '))
    inout(i)

    print "\nLeaving guys..."
    print C[i,:].argsort()[-10:]
    print C[i,C[i,:].argsort()[-10:]]

    print "\nEntering guys..."
    print C[:,i].argsort()[-10:]
    print C[C[:,i].argsort()[-10:],i]

    #compute transitions up & transitions down 
    computeTransitionsUpDown = (str(raw_input('\nEnter "y" to compute transitions to a set of states and the conjugate: ')).lower() == 'y')
    if computeTransitionsUpDown == True:
            down = np.array([int(j) for j in str(raw_input('\nEnter set of states:')).split()])
            print down
            updown(i,down)
    
    notFinish = not bool(raw_input("\nPress a key other than ENTER to finish "))
