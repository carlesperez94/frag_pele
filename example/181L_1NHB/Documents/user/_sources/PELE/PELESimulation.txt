.. _sec-PELESimulation:

PeleSimulation block in PELE++
==============================

commandType
-----------

For a whole PELE simulation, it must be set to:

.. code-block:: json

   "commandType": "peleSimulation"

Algorithms
----------

Pele must have the following algorithms defined:
:ref:`Perturbation <sec-perturbation>`,
:ref:`ANM <sec-anm>`, :ref:`Side Chain
Prediction <sec-sideChainPrediction>` and
:ref:`Minimization <sec-minimization>`. If an algorithm is
defined at peleSimulation level it is going to be used throughout the
whole simulation, unless it is overriden in a task by defining a new
one.

RandomGenerator
---------------

For more information, check :ref:`random generator <sec-randomGenerator>`. 

Example
^^^^^^^

.. code-block:: json

       "RandomGenerator": {
           "seed": 100586
       }

selectionToPerturb
------------------

Selects the ligands to be perturbed in the ligand perturbation. If no perturbation is performed (that is, the ``Perturbation`` block is empty), then the selection to perturb is empty for all purposes, and this parameter will be ignored.

If the ``Perturbation`` block (:ref:`sec-perturbation`) is non-empty (that is, perturbation is active), then you must provide this parameter, since it has no default value.

The selection is done for whole chains.

Example
^^^^^^^

.. code-block:: json

           "selectionToPerturb": {
               "chains": {
                   "names": [ "L" ]
               }
           }

Perturbation of multiple chains.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This feature does not have complete support in PELE (for example, you
cannot use different options for different perturb chains; besides, the
steering vector is shared between all perturb chains). By knowing which
chain is perturbed at each PELE step (see below), you can change the
perturbation parameters using the "parametersChanges" block (see
:ref:`PELEDynamicChangesInSimulationParameters <sec-PELEDynamicChangesInSimulationParameters>`)
and a counter-based condition (see :ref:`example 9 of
Conditions <sec-conditionsByExample-example9>`).

If you specify more than one chain, then PELE will cycle through them
perturbing only one of them at each PELE step. The order in which
multiple chains are considered for perturbation is:

-  first, protein chains in the order they appear in the input PDB.
-  second, DNA chains in the order they appear in the input PDB.
-  third, RNA chains in the order they appear in the input PDB.
-  fourth, all other chains (ligands) in the order they appear in the
   input PDB.

Therefore, if you select chains "L", "M", both of ligand type, and "M"
appears first in the PDB file, then the first PELE step will perturb
chain "M", next step will perturb chain "L", the step after that will
perturb "M", and so on.

NonBonding
----------

Configuration of the algorithm for non-bonding energy-related terms
calculation.

For more details see
:ref:`NonBonding <sec-nonBondingAlgorithms>`.

Pele_parameters
----------------

Parameters used by the Pele algorithm. They can be defined at
PeleSimulation or PeleTask level.

For more details see :ref:`Pele Parameters <sec-peleParameters>`.

Pele_Output
------------

Parameters related to the output files generated in the simulation. They
must be defined at the PeleSimulation level.

For more details see :ref:`Pele Output Parameters <sec-peleOutput>`.

Constraints
-----------

Constraints applied in all minimizations during a Pele simulation. Check
the :ref:`Constraints section <sec-permanentConstraints>`.

PeleTasks
---------

Serie of tasks that Pele has to carry out before a command (and usually
a simulation) is considered to be finished.

For more information check the :ref:`Pele Tasks section <sec-peleTasks>`.

