============
AdaptivePELE
============


|MIT license| |GitHub release| |DOI|

AdaptivePELE is a Python module to perform enhancing sampling of molecular
simulation built around the Protein Energy Landscape Exploration method (`PELE <https://pele.bsc.es/pele.wt>`_) developed in the Electronic and Atomic Protein Modelling grop (`EAPM <https://www.bsc.es/discover-bsc/organisation/scientific-structure/electronic-and-atomic-protein-modeling-eapm>`_) at the Barcelona Supercomputing Center (`BSC <https://www.bsc.es>`_).

Usage
-----

AdaptivePELE is called with a control file as a
parameters. The control file is a json document that contains 4 sections:
general parameters, simulation parameters, clustering parameters and spawning
parameters. The first block refers to general parameters of the adaptive run,
while the other three blocks configure the three steps of an adaptive sampling
run, first run a propagation algorithm (simulation), then cluster the
trajectories obtained (clustering) and finally select the best point to start
the next iteration (spawning).

An example of usage::

    python -m AdaptivePELE.adaptiveSampling controlFile.conf

Installation
------------

In order to have a running copy of AdaptivePELE, you need to install and compile cython files in the base folder with::

    git clone https://github.com/AdaptivePELE/AdaptivePELE.git
    cd AdaptivePELE
    python setup.py build_ext --inplace

Also, if AdaptivePELE was not installed in a typical library directory, a common option is to add it to your local PYTHONPATH::

    export PYTHONPATH="/location/of/AdaptivePELE:$PYTHONPATH"

Documentation
-------------

The documentation for AdaptivePELE can be found `here <https://adaptivepele.github.io/AdaptivePELE/>`_


Citation 
--------

AdaptivePELE is research software. If you make use of AdaptivePELE in scientific publications, please cite it. The BibTeX reference is::

    @article{Lecina2017,
    author = {Lecina, Daniel and Gilabert, Joan Francesc and Guallar, Victor},
    doi = {10.1038/s41598-017-08445-5},
    issn = {2045-2322},
    journal = {Scientific Reports},
    number = {1},
    pages = {8466},
    pmid = {28814780},
    title = {{Adaptive simulations, towards interactive protein-ligand modeling}},
    url = {http://www.nature.com/articles/s41598-017-08445-5},
    volume = {7},
    year = {2017}
    }


.. |MIT license| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://lbesson.mit-license.org/


.. |GitHub release| image:: https://img.shields.io/github/release/AdaptivePELE/AdaptivePELE.svg
    :target: https://github.com/AdaptivePELE/AdaptivePELE/releases/

.. |DOI| image:: https://zenodo.org/badge/DOI/10.1038/s41598-017-08445-5.svg
  :target: https://doi.org/10.1038/s41598-017-08445-5
