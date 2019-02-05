import pyemma.msm as msm
import pyemma.plots as mplt
import helper
import numpy as np

def selectTPTSets(MSM_object, indexA, indexB):
    """ Extract from the sets of the PCCA clustering the sets that will serve
    as the extrems of the TPT
    """
    PCCAsets = MSM_object.metastable_sets
    SetA = PCCAsets[indexA]
    SetB = PCCAsets[indexB]
    return SetA, SetB

def createTPT(MSM_object, A, B):
    """ Calculate the reactive flux between sets A and B.
    Return a ReactiveFlux object"""
    return msm.tpt(MSM_object, A, B)

def coarseTPT(TPT_object, MSM_object):
    """Coarse the TPT object into one with as many states as the PCCA to ease
    visualization """
    (foo, coarseTPT) = TPT_object.coarse_grain(MSM_object.metastable_sets)
    return coarseTPT

def plotTPT(TPT_object, state_labels='auto', outfile=None):
    """ Plot the flux network from a Reactive Flux object"""
    flux_figure = mplt.plot_flux(TPT_object,state_labels=state_labels)
    if not outfile is None:
        flux_figure.savefig(outfile)
    return flux_figure

def writeTPTOutput(TPT_object):
    """ Write the main properties of the TPT calculated"""
    helper.makeFolder("tpt")
    with open("tpt/sets.dat","w") as f:
        f.write("TPT created with PCCA sets: %d and %d"%(TPT_object.A[0],TPT_object.B[0]))
    np.savetxt("tpt/backward_committor.dat", TPT_object.backward_committor)
    np.savetxt("tpt/forward_committor.dat", TPT_object.forward_committor)
    np.savetxt("tpt/gross_flux.dat", TPT_object.gross_flux)
    np.savetxt("tpt/net_flux.dat", TPT_object.net_flux)

