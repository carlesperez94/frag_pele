.. FrAG-PELE documentation master file, created by
   sphinx-quickstart on Mon Jan 21 17:55:18 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FrAG-PELE. Fragment-based Automatic Growing.
============================================

FrAG a new tool for in silico hit-to-lead drug design, capable of growing a fragment into a core while exploring the protein-ligand conformational space.

The general workflow starts with the system preparation, where core and fragment are put together, followed by several growing stages, where various parallel Monte Carlo simulations are performed while increasing linearly the fragment characteristics. In each stage a clustering and spawning is done to trigger the next growing stage. Finally, a simulation analysis package scores best poses and  retrieved its PDBs.

FrAG is build on top of Protein Energy Landscape Exploration (PELE) Software.

.. image:: img/growing.gif
   :scale: 80 %
   :align: center

Requirements (you can use other versions of all modules, but they have not been tested with the newest version):

- Python >= 3.7.6

- AdaptivePELE >= 1.6.1

- Schrödinger's Python = 2.7.13

- Biopython = 1.70

- Prody >= 1.10.10

- Pandas >= 1.0.1

- rdkit = 2019.09.3

Installation
============
.. toctree::
   installation
   first_steps

Protocols
==========
.. toctree::
   protocols

Analysis
========
.. toctree::
   analysis_tools

Documentation
=============
.. toctree::
   documentation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`