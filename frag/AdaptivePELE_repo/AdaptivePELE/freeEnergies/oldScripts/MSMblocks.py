import os
import trajectories
import helper
import msm
import tpt
import json
import matplotlib.pyplot as plt


class PrepareMSM:
    def __init__(self, numClusters, trajectoryFolder, trajectoryBasename, stride=1, cluster=True):

        self.discretizedFolder = "discretized"
        self.clusterCentersFile = os.path.join(self.discretizedFolder,
                                               "clusterCenters.dat")
        self.dTrajTemplateName = os.path.join(self.discretizedFolder, "%s.disctraj")
        self.clusteringFile = "clustering_object.pkl"
        self.stride = stride
        self.trajFilenames = []
        self.dtrajs = []
        self.cl = None

        if os.path.exists(self.clusteringFile):
            self.cl = helper.loadClustering(self.clusteringFile)

        self.clusterTrajectories(numClusters, trajectoryFolder, trajectoryBasename)


    def clusterTrajectories(self, numClusters, trajectoryFolder, trajectoryBasename):
        print "Loading trajectories..."
        self.x, self.trajFilenames = trajectories.loadCOMFiles(trajectoryFolder, trajectoryBasename)

        # cluster & assign
        if self.cl is None:
            print "Clustering data..."
            self.cl = trajectories.clusterTrajectories(self.x, numClusters, stride=self.stride)
            print "Assigning data..."
            self.dtrajs = self.cl.dtrajs[:]
        else:
            print "Assigning data..."
            #THIS DOES NOT WORK AS SAID IN PYEMMA DOCUMENTATION
            #IT KEEPS PREVIOUS DTRAJS, AND I COULDNT MANAGE TO REMOVE THEM
            self.dtrajs = self.cl.assign(self.x)
        # self.cl = trajectories.clusterRegularSpace(self.x, 3, stride=self.stride)
        # write output
        print "Writing clustering data..."
        helper.makeFolder(self.discretizedFolder)
        helper.writeClusterCenters(self.cl, self.clusterCentersFile)
        helper.writeDtrajs(self.trajFilenames, self.dtrajs, self.dTrajTemplateName)
        print "Saving clustering object..."
        helper.saveClustering(self.cl, self.clusteringFile)

    def getClusteringObject(self):
        return self.cl


class MSM:
    def __init__(self, cl, lagtimes, numPCCA, itsOutput=None, numberOfITS=-1,
                 itsErrors=None, error_estimationCK=None, mlags=2, lagtime=None, dtrajs=[]):
        self.MSMFile = "MSM_object.pkl"
        self.MSM_object = None
        if os.path.exists(self.MSMFile):
            self.MSM_object = helper.loadMSM(self.MSMFile)

        self.cl = cl #only used for pcca
        self.lagtimes = lagtimes
        self.numPCCA = numPCCA
        self.itsOutput = itsOutput
        self.numberOfITS = numberOfITS
        self.itsErrors = itsErrors
        self.error_estimationCK = error_estimationCK
        self.mlags = mlags
        self.lagtime = lagtime
        self.dtrajs = dtrajs

    def estimate(self):
        if self.lagtime is None:
            self.lagtime = self.calculateITS() #keep calculating until convergence is reached
        else:
            its_object = msm.calculateITS(self.dtrajs, self.lagtimes, self.itsErrors)
            plot_its = msm.plotITS(its_object, self.itsOutput, self.numberOfITS)
        if self.MSM_object is None:
            self.createMSM(self.lagtime)
        #self.PCCA(self.numPCCA)
        print "Saving MSM object..."
        helper.saveMSM(self.MSM_object)
        #self.performCKTest(self.error_estimationCK)

    def getMSM_object(self):
        return self.MSM_object

    def performCKTest(self, error_estimationCK=None):
        # Chapman-Kolgomorov validation
        nsetsCK = len(self.MSM_object.metastable_sets)
        print ("Performing Chapman-Kolmogorov validation with the %d sets from"
               " the PCCA, when it's done will prompt for the validity of "
               "the model...") % nsetsCK
        membershipsCK = self.MSM_object.metastable_memberships
        CKObject = msm.ChapmanKolmogorovTest(self.MSM_object,
                                             nsetsCK, memberships=membershipsCK,
                                             error_estimation=error_estimationCK,
                                             mlags=self.mlags)
        msm.plotChapmanKolmogorovTest(CKObject)
        plt.show()

    def PCCA(self, numPCCA):
        # PCCA
        print "Calculating PCCA cluster with %d sets..." % numPCCA
        self.MSM_object = msm.calculatePCCA(self.MSM_object, numPCCA)

    def createMSM(self, lagtime):
        # estimation
        print "Estimating MSM with lagtime %d..." % lagtime
        error_est = self.itsErrors or self.error_estimationCK
        self.MSM_object = msm.estimateMSM(self.dtrajs, lagtime, error_est)

    def calculateITS(self):
        is_converged = False
        # its
        print ("Calculating implied time-scales, when it's done will prompt "
               "for confirmation on the validity of the lagtimes...")
        while not is_converged:
            its_object = msm.calculateITS(self.dtrajs, self.lagtimes,
                                          self.itsErrors)
            plot_its = msm.plotITS(its_object, self.itsOutput, self.numberOfITS)
            plt.show()
            while True:
                convergence_answer = raw_input("Has the ITS plot converged?[y/n] ")
                convergence_answer.rstrip()
                convergence_answer = convergence_answer or "y"  # Making yes the default answer
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
                    self.lagtimes.extend(map(int, lag_list.split(" ")))
                elif new_lagtimes.lower() == "new" or new_lagtimes.lower() == "n":
                    lag_list = raw_input("Please input the new lagtimes separated by a space: ")
                    lag_list.rstrip()
                    self.lagtimes = map(int, lag_list.split(" "))
                self.lagtimes.sort()
        return lagtime

    def writeClustersForWMD(self, outfile="clusters.pdb"):
        tempStr = "HETATM  {:>5d} H1 CLT L {:>4d} {:>8.3f}{:>8.3f}{:>8.3f} 0.75  {:<6f}           H\n"
        with open(outfile, "w") as f:
            for i, coords in enumerate(self.cl.clustercenters):
                f.write(tempStr.format(i+1, i+1, coords[0], coords[1], coords[2], self.MSM_object.stationary_distribution[i]))


class TPT:
    def __init__(self, MSM_object, cl, outfile_fluxTPT, state_labels):
        # Identify relevant sets for TPT
        print ("Plotting PCCA sets, identify the sets that will serve as "
               "source and sink for TPT...")
        msm.plot_PCCA_clusters(cl, MSM_object)
        plt.show()
        SetA_index = int(raw_input("Please input index of Set A(source):"))
        SetB_index = int(raw_input("Please input index of Set B(sink):"))
        SetA, SetB = tpt.selectTPTSets(MSM_object, SetA_index, SetB_index)
        print "Creating TPT object..."
        self.TPT_object = tpt.createTPT(MSM_object, SetA, SetB)
        print "Coarsing TPT for visualization..."
        self.coarseTPT_object = tpt.coarseTPT(self.TPT_object, MSM_object)
        print "Plotting TPT flux diagram..."
        flux_figure = tpt.plotTPT(self.coarseTPT_object,
                                  state_labels=state_labels,
                                  outfile=outfile_fluxTPT)
        plt.show()
        print "Writing the main properties of the TPT in the tpt/ folder..."
        tpt.writeTPTOutput(self.coarseTPT_object)

    def getTPTObject(self):
        return self.TPT_object

    def getCoarseTPTObject(self):
        return self.coarseTPT_object


def readParams(control_file):
    """Reads a JSON control file for the paramaters.
    List of parameters:
        trajectoryFolder : string, folder where the trajectory data is
        trajectoryBasename : string, names of the trajectory files is
        numClusters : int, number of clusters to partition the data
        lagtimes: list of ints, lagtimes to calculate the implied time-scales
        numPCCA : int, number of sets to perform PCCA
        itsOutput : string (optional), name of the file where to store the implied time-scales plot
        numberOfITS : int (optional), number of eigenvalues to include in the ITS plot
        itsErrors : string (optional), 'bayes' to include error estimation in ITS plot
        error_estimationCK : bool (optional), wether to include error estimation in CK test
        state_labels : list of strings (optional), list of labels for the sets in the coarse TPT flux plot
        outfile_fluxTPT : string, name of the file to store the coarse TPT flux plots
        """
    with open(control_file, "r") as f:
        paramsJSON = json.load(f)
    return paramsJSON
