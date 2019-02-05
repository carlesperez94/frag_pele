from __future__ import absolute_import, division, print_function, unicode_literals
import os
try:
    import cPickle
except ImportError:
    import pickle as cPickle
import matplotlib.pyplot as plt
import pyemma.msm as msm
import pyemma.plots as mplt
import numpy as np


class MSM:
    """
    def __init__(self, lagtimes, numPCCA, itsOutput=None, numberOfITS=-1,
                 itsErrors=None, error_estimationCK=None, mlags=2, lagtime=None, dtrajs=[]):
    """
    def __init__(self, error=False, dtrajs=None):
        """
            If MSM_object.pkl exists, and we call estimate, does it override whatever was before?
        """

        self.MSMFile = "MSM_object.pkl"
        self.MSM_object = None
        if os.path.exists(self.MSMFile):
            self.MSM_object = self.loadMSM(self.MSMFile)

        self.lagtimes = []
        self.itsOutput = None
        self.numberOfITS = -1
        self.error = error
        self.lagtime = 0
        if dtrajs is None:
            dtrajs = []
        self.dtrajs = dtrajs
        self.stationaryDistributionFilename = "stationaryDistribution.dat"
        self.numPCCA = None

    def estimate(self, lagtime=None, lagtimes=None, numberOfITS=-1):
        self.lagtime = lagtime
        if lagtimes is None:
            lagtimes = []
        self.lagtimes = lagtimes
        self.numberOfITS = numberOfITS
        print("LAGTIME", self.lagtime)
        try:
            self.lagtime = self._calculateITS()  # keep calculating until convergence is reached
        except Exception as err:
            # This error should happen only when running in the MN computing queues
            # with the tkinter matplotlib backend, to avoid it set the default
            # backend to pdf
            if "connect to display" in err.message:
                print("ITS plots not saved because of MN error")
            else:
                raise err
        print("Using lagtime = ", self.lagtime)
        self.buildMSM()
        np.savetxt(self.stationaryDistributionFilename, self.MSM_object.pi)
        self.saveMSM(self.MSM_object)

    def computePCCA(self, numberOfPCCAclusters):
        self.numPCCA = numberOfPCCAclusters
        self.MSM_object.pcca(self.numPCCA)

    def runCKTest(self, mlags):
        self._performCKTest(self.error)

    def getMSM_object(self):
        return self.MSM_object

    def _performCKTest(self, mlags):
        # Chapman-Kolgomorov validation
        nsetsCK = len(self.MSM_object.metastable_sets)
        print(("Performing Chapman-Kolmogorov validation with the %d sets from"
               " the PCCA, when it's done will prompt for the validity of "
               "the model...") % nsetsCK)
        membershipsCK = self.MSM_object.metastable_memberships
        CKObject = self.ChapmanKolmogorovTest(self.MSM_object,
                                              nsetsCK, memberships=membershipsCK,
                                              error_estimation=self.error,
                                              mlags=mlags)
        self.plotChapmanKolmogorovTest(CKObject)
        plt.show()

    def _calculateITS(self):
        is_converged = False
        # its
        print(("Calculating implied time-scales, when it's done will prompt "
               "for confirmation on the validity of the lagtimes..."))
        while not is_converged:
            if not self.error:
                itsErrors = None
            elif self.error:
                itsErrors = "bayes"
            if self.lagtimes and self.lagtimes is not None:
                # workaround to get new its plot at each iteration, the
                # plot_implied_timescales function is calling plt.gca() and
                # recovers the previous plot's axes, by creating a new figure
                # gca gets a set of empty axes and plots are fine
                plt.figure()
                its_object = msm.its(self.dtrajs, lags=self.lagtimes, errors=itsErrors)
                mplt.plot_implied_timescales(its_object, outfile=self.itsOutput, nits=self.numberOfITS)
                plt.savefig("its.png")
            if self.lagtime is not None:
                return self.lagtime
            while True:
                plt.show()
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
                    print("Answer not valid. Please answer yes or no")
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

    def plotITS(self, its_object, its_plot_file=None, nits=-1):
        its_plot = mplt.plot_implied_timescales(its_object, outfile=its_plot_file, nits=nits)
        plt.savefig("its.eps")
        return its_plot

    def buildMSM(self):
        """ Estimate a MSM from the trajectories using a provided lagtime that
        should be big enough so that the relevant processes have converged.

        self.error: whether to estimate errors or not
        """
        if self.error:
            self.MSM_object = msm.bayesian_markov_model(self.dtrajs, self.lagtime)
        else:
            self.MSM_object = msm.estimate_markov_model(self.dtrajs, self.lagtime)

    def ChapmanKolmogorovTest(self, MSM_object, nsets, memberships=None, error_estimation=False, mlags=2):
        """ Perform the ChapmanKolmogorov test to validate the MSM"""
        return MSM_object.cktest(nsets, memberships=memberships, err_est=error_estimation, mlags=mlags)

    def plotChapmanKolmogorovTest(self, CKObject, layout=None, padding_between=0.1,
                                  padding_top=0.075):
        """ Plot the results of the Chapman-Kolgomorov tests"""
        mplt.plot_cktest(CKObject, layout=layout, padding_between=padding_between,
                         padding_top=padding_top)
        plt.savefig("CK.eps")

    def loadMSM(self, MSMFile):
        with open(MSMFile) as MSMfile:
            MSM_object = cPickle.load(MSMfile)
        return MSM_object

    def saveMSM(self, MSM_object):
        """Save the MSM object to avoid having to run again
        the more computationally expensive part"""
        with open("MSM_object.pkl", "wb") as MSMfile:
            cPickle.dump(MSM_object, MSMfile, 2)
