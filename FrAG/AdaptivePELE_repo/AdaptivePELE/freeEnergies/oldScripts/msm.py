import pyemma.msm as MSM
import pyemma.plots as mplt
import matplotlib.pyplot as plt

def calculateITS(trajectories, lagtimes, errors = None):
    """ Calulate the implied time-scales at the given lagtimes"""
    its_object = MSM.its(trajectories, lags=lagtimes, errors=errors)
    return its_object

def plotITS(its_object, its_plot_file=None, nits=-1):
    its_plot = mplt.plot_implied_timescales(its_object, outfile=its_plot_file, nits=nits)
    plt.savefig("its.eps")
    return its_plot

def estimateMSM(trajectories,lagtime, error_est=False):
    """ Estimate a MSM from the trajectories using a provided lagtime that
    should be big enough so that the relevant processes have converged.
    Return a MaximumLikelihoodMSM object"""
    if error_est:
        print "Computing msm with bayes error calc"
        MSM_object = MSM.bayesian_markov_model(trajectories, lagtime)
    else:
        print "Computing msm with no error calc"
        MSM_object = MSM.estimate_markov_model(trajectories, lagtime, count_mode='sliding')
    return MSM_object

def calculatePCCA(MSM_object, numPCCA):
    """ Coarse-cluster the MSM usin numPCCA clusters.
    Return a PCCA object"""
    MSM_object.pcca(numPCCA)
    return MSM_object


def ChapmanKolmogorovTest(MSM_object, nsets,memberships=None, error_estimation=False, mlags=2):
    """ Perform the ChapmanKolmogorov test to validate the MSM"""
    return MSM_object.cktest(nsets,memberships=memberships,err_est=error_estimation, mlags=mlags)

def plotChapmanKolmogorovTest(CKObject, layout=None, padding_between=0.1,
                              padding_top=0.075):
    """ Plot the results of the Chapman-Kolgomorov tests"""
    mplt.plot_cktest(CKObject,layout=layout, padding_between=padding_between,
                     padding_top=padding_top)
    plt.savefig("CK.eps")

def plot_PCCA_clusters(cluster_object, MSM_object):
    cols = ['orange', 'magenta', 'red', 'black', 'blue', 'green',]
    ccx = cluster_object.clustercenters[:,0]
    ccy = cluster_object.clustercenters[:,1]
    ccz = cluster_object.clustercenters[:,2]
    fig,axes = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row')
    for index,set in enumerate(MSM_object.metastable_sets):
        axes[0][0].scatter(ccx[set],ccy[set],color=cols[index],label="Set %d"%index)
        axes[0][1].scatter(ccz[set],ccy[set],color=cols[index],label="Set %d"%index)
        axes[1][0].scatter(ccx[set],ccz[set],color=cols[index],label="Set %d"%index)
    axes[1][0].legend(loc='center right', bbox_to_anchor=[1.8,0.5])
    axes[0][0].set_ylabel('y')
    axes[1][0].set_ylabel('z')
    axes[1][0].set_xlabel('x')
    axes[1][1].axis('off')
    axes[0][1].set_xticks(axes[1][1].get_xticks())
    axes[0][1].set_xticklabels(axes[1][1].get_xticklabels())
    axes[0][1].set_xlabel('z')
    return fig


