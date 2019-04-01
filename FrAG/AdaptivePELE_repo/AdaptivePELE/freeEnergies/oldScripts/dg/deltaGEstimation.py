import os
import numpy as np
import firstSnapshots
from simulation import simulationrunner
import shutil
import computeDG
import glob
import checkDetailedBalance

def rm(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

def rmTrajFiles(trajWildcard):
    allfiles = glob.glob(trajWildcard)
    for f in allfiles:
        rm(f)

def copyAllTrajectories(trajWildcard, folderWithTraj):
    allfiles = glob.glob(os.path.join(folderWithTraj, trajWildcard))
    for f in allfiles:
        shutil.copy(f, ".")

ntrajs = 50
length = 2500
lagtime = 200
nclusters = 100
nruns = 100
useAllTrajInFirstRun = True

computeDetailedBalance = True

controlFile = "templetized_control_MSM.conf" #unused, but remember it needs to be defined
trajWildcard = "traj*"
folderWithTraj = "rawData"
folderTraj = "."
deltaGs = {}
detailedBalance = {}
allFolders = os.listdir(folderTraj)
Epochs = [int(epoch) for epoch in allFolders if epoch.isdigit()]
[shutil.rmtree(fold) for fold in glob.glob("coord_*")]
i = -1
for folder in sorted(Epochs):
    print "Epoch", folder
    # rmTrajFiles(trajWildcard)
    # if useAllTrajInFirstRun and i == 0: #this uses full length trajectories
    #     rm("clustering_object.pkl") 
    #     copyAllTrajectories(trajWildcard, folderWithTraj)
    # else:
    #     firstSnapshots.main(length, ntrajs)
    i += 1
    folder = str(folder)
    shutil.copytree(os.path.join(folderTraj, folder, "repeatedExtractedCoordinates"), "coord_"+folder)
    rm("clustering_object.pkl") 
    rm("MSM_object.pkl") 
    rmTrajFiles("discretized/traj_*")
    if i < 1:
        nclusters = 10
    elif i < 3:
        nclusters = 50
    elif i < 5:
        nclusters = 100  
    else:
        nclusters = 200
    deltaG = computeDG.computeDG(lagtime, nclusters, os.path.join("coord_%s"%folder,"coord_*"))
    try:
        deltaGs[folder].append(deltaG)
    except KeyError:
        deltaGs[folder] = [deltaG]

    if computeDetailedBalance:
        frobeniusAvg, relativeEntropyT, unused = checkDetailedBalance.main("discretized", 0, lagtime)
        try:
            detailedBalance['frobenius'].append(frobeniusAvg)
            detailedBalance['relativeEntropy'].append(relativeEntropyT)
        except KeyError:
            detailedBalance['frobenius'] = [frobeniusAvg]
            detailedBalance['relativeEntropy'] = [relativeEntropyT]

    shutil.copyfile("its.png", "its_%d.png"%i)
    shutil.copyfile("volumeOfClusters.dat", "volumeOfClusters_%d.dat"%i)
    shutil.copyfile("clusters.pdb", "clusters_%d.pdb"%i)
    shutil.copyfile("pmf_xyzg.dat", "pmf_xyzg_%d.dat"%i)
    if i == 0:
        try:
            shutil.copyfile("db_frobenius.eps", "db_frobenius_%d.eps"%i)
            shutil.copyfile("db_abs_diff.eps", "db_abs_diff_%d.eps"%i)
            shutil.copyfile("db_flux.eps", "db_flux_%d.eps"%i)
        except:
            pass

#PLOT RESULTS
print "clusters: %d, ntrajs: %d, trajLength: %d"%(nclusters, ntrajs, length)
for key, val in deltaGs.iteritems():
    print key
    print "====="
    dGs = []
    for element in val:
        print element
        dG = element.split()[1]
        dGs.append(float(dG))
    try:
        print "dG = %f +- %f"%(np.mean(dGs), np.std(dGs))
    except:
        print "dG = %f"%(np.mean(dGs))

if computeDetailedBalance:
    for key, val in detailedBalance.iteritems():
        print key
        print "====="
        values = []
        for element in val:
            print element
            values.append(float(element))
        print "%s = %f +- %f"%(key, np.mean(values), np.std(values))
