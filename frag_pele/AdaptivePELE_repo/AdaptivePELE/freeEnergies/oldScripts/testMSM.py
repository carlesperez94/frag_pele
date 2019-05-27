import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyemma.plots as plots
import trajectories
import msm
import tpt
import helper
"""
############# This script is a copy of the buildMSM.py script, this is intended
to executed inside a python terminal to dump all variables in memory and debug
or test different pyemma methods
#############
"""

    #TODO: Define blocks with tasks to make the program more modular

    ### parameters
trajectoryFolder = "test/MSM2"
trajectoryBasename = "*traj_*"

numClusters = 100

lagtimes = [1, 2, 5, 10, 20, 50, 100, 200, 500]
itsOutput = "its.png"
numberOfITS = -1
itsErrors=None #'bayes'
nsetsCK = 2
error_estimationCK=False
membershipsCK=None
numPCCA = 4
state_labels = 'auto' # default value of the labels in the flux diagram of
# the TPT
outfile_fluxTPT = None # file to store the flux diagram of the TPT, default
# not saving it(None)

### constants
discretizedFolder = "discretized"
clusterCentersFile = os.path.join(discretizedFolder, "clusterCenters.dat")
discTraj = os.path.join(discretizedFolder, "%s.disctraj")
clusteringFile="clustering_object.pkl"
MSMFile="MSM_object.pkl"

#program

if os.path.exists(MSMFile) and os.path.exists(clusteringFile):
    MSM_object = helper.loadMSM(MSMFile)
    cl = helper.loadClustering(clusteringFile)
else:
    print "Loading trajectories..."
    x = trajectories.loadCOMFiles(trajectoryFolder, trajectoryBasename)

    #cluster & assign
    print "Clustering data..."
    cl = trajectories.clusterTrajectories(x, numClusters)

    #write output
    print "Writing clustering data..."
    helper.makeFolder(discretizedFolder)
    helper.writeClusterCenters(cl, clusterCentersFile)

    """
    fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    ax.scatter(clX, clY, clZ)
    plt.show()
    """

    plot_its = None #just a placeholder so it doesn't go out of scope
    is_converged = False
    #its
    print ("Calculating implied time-scales, when it's done will prompt for "
            "confirmation on the validity of the lagtimes...")
    while not is_converged:

        its_object = msm.calculateITS(cl.dtrajs, lagtimes, itsErrors)
        plot_its = msm.plotITS(its_object, itsOutput, numberOfITS)
        plt.show()
        while True:
            convergence_answer = raw_input("Has the ITS plot converged?[y/n] ")
            convergence_answer.rstrip()
            convergence_answer = convergence_answer or "y" #Making yes the default
            #answer
            if convergence_answer.lower() == "y" or convergence_answer.lower() == "yes":
                is_converged = True
                lagtime_str = raw_input("Please input the lagtime to construct the MSM: ")
                lagtime = int(lagtime_str.rstrip())
                break
            elif convergence_answer.lower() == "n" or convergence_answer.lower() == "no":
                break
            else:
                print "Answer not valid. Please answer yes or no"
            if not is_converged:
                new_lagtimes = raw_input("Do you want to define new lagtimes or add to the previous?[add(a)/new(n)] ")
                new_lagtimes.rstrip()
                if new_lagtimes.lower() == "add" or new_lagtimes.lower() == "a":
                    lag_list = raw_input("Please input the lagtimes you want to add separated by a space: ")
                    lag_list.rstrip()
                    lagtimes.extend(map(int,lag_list.split(" ")))
                elif new_lagtimes.lower() == "new" or new_lagtimes.lower() == "n":
                    lag_list = raw_input("Please input the new lagtimes separated by a space: ")
                    lag_list.rstrip()
                    lagtimes = map(int,lag_list.split(" "))
                    lagtimes.sort()

    #estimation
    print "Estimating MSM with lagtime %d..."%lagtime
    MSM_object = msm.estimateMSM(cl.dtrajs, lagtime)


    #PCCA
    print "Calculating PCCA cluster with %d sets..."%numPCCA
    MSM_object = msm.calculatePCCA(MSM_object, numPCCA)

    print "Saving MSM and clustering objects..."
    helper.saveMSM(MSM_object, cl)

#Chapman-Kolgomorov validation
print ("Performing Chapman-Kolmogorov validation with the %d sets from the "
        "PCCA, when it's done will prompt for the validity of the model...")%numPCCA
nsetsCK = len(MSM_object.metastable_sets)
membershipsCK = MSM_object.metastable_memberships
CKObject = msm.ChapmanKolmogorovTest(MSM_object,
                                        nsetsCK,memberships=membershipsCK,
                                        error_estimation=error_estimationCK)
msm.plotChapmanKolmogorovTest(CKObject)
plt.show()

#Identify relevant sets for TPT
print ("Plotting PCCA sets, identify the sets that will serve as source and sink "
       "for TPT...")
msm.plot_PCCA_clusters(cl, MSM_object)
plt.show()
SetA_index = int(raw_input("Please input index of Set A(source):"))
SetB_index = int(raw_input("Please input index of Set B(sink):"))
SetA, SetB = tpt.selectTPTSets(MSM_object, SetA_index, SetB_index)
print "Creating TPT object..."
TPT_object = tpt.createTPT(MSM_object, SetA, SetB)
print "Coarsing TPT for visulization..."
coarseTPT_object = tpt.coarseTPT(TPT_object, MSM_object)
print "Plotting TPT flux diagram..."
flux_figure = tpt.plotTPT(coarseTPT_object, state_labels=state_labels,
                          outfile=outfile_fluxTPT)
plt.show()
print "Writing the main properties of the TPT in the tpt/ folder..."
tpt.writeTPTOutput(coarseTPT_object)

#Free energy estimation
print "Calculating free energies..."
kbt = 0.0019872041*300
pi = MSM_object.stationary_distribution


"""
for each cluster
for each disctraj
these_frames = np.argwhere(disc_traj==i)
theta_i = cont_traj[these_frames]
rho_i, _ = np.histogram(theta_i, range=myrange, normed=True)

pi_total = rho_i * pi_i, donde pi son los counts relativos de cada centro.

Si el mismo "myrange" para para las distintas bases (regspace y kmeans),
la Delta_G deberia ser igual!

_total = \Sigma_i rho_i * pi_i,  donde pi son los counts relativos de cada
Ojo, que para que esto funcione, cada rho_i tiene que ser una densidad
verdadera (rho_i.sum() * dx = 1) y pi tambien \Sigma_i pi_i = 0).
"""

Gpmf = -kbt*np.log(pi/np.max(pi))

Gpmf = np.expand_dims(Gpmf, axis=1)
output = np.hstack((cl.clustercenters, Gpmf))
np.savetxt("output.txt", output)


