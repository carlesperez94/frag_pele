********
Overview
********

A PELE simulation consists of several steps of the following 3-part algorithm:

1. Perform perturbation

  a. Ligand (if present): Translation and rotation. Basic side chain relocation to avoid clashes.
  b. Protein: Get direction vectors with ANM. Perform movement through minimization.

2. Relaxation.

  a. Side Chain Prediction

    i. Around the ligand (if present).
    ii. For the top side chains among those whose energy has changed the most in the ANM step.

  b. Global minimization.

3. Accept or reject the step based on the Metropolis criterium.


The following image shows the previous steps:

.. image:: ../images/pele_scheme.png
