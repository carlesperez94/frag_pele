AdaptivePELE
============

AdaptivePELE aims at enhancing sampling of molecular simulations.

----

AdaptivePELE uses cython in order to speedup internal processes
Please, compile AdaptivePELE with:

python setup.py build_ext --inplace

----

AdaptivePELE is called with a control file as a
parameters. The control file is a json document that contains 4 sections:
general parameters, simulation parameters, clustering parameters and spawning
parameters. The first block refers to general parameters of the adaptive run,
while the other three blocks configure the three steps of an adaptive sampling
run, first run a propagation algorithm (simulation), then cluster the
trajectories obtained (clustering) and finally select the best point to start
the next iteration (spawning).

An example of usage:

python -m AdaptivePELE.adaptiveSampling controlFile.conf
