.. _sec-minimization:

Minimization in Pele++ control file
===================================

The ``"Minimizer"`` block in a ``peleSimulation`` command configures the final minimization of the relaxation phase. This block also configures the minimization in a ``minimization`` command. In a ``peleSimulation``, this phase will be active if either the side chain prediction phase is active, or if the ``minimizationFrequency`` PELE parameter says so (see :ref:`sec-peleParameters-minimizationFrequency`). In particular, every time a side chain prediction is performed, a final minimization also happens; besides, every ``minimizationFrequency`` steps, a minimization is performed (if it is also a step where the side chain prediction is performed, only one minimization, after the prediction, is executed).

There also exist minimization blocks for the minimizers used in the ANM phase (:ref:`sec-anm`) and in the Side Chain Prediction phase (:ref:`sec-sideChainPrediction`), where they are respectively called ``"anmMinimizer"`` and ``"sideChainMinimizer"``. Unless noted, all parameters shown below can also be used when configuring those minimizers. The default values shown below correspond to the ``"Minimizer"`` block; see the ANM and Side Chain Prediction manual pages for the defaults used in their corresponding minimizers.

In the minimization phase of a PELE simulation, the following will be minimized:

- if ``minimizationRegionRadius`` (:ref:`sec-peleParameters-minimizationRegionRadius`) is not zero, then the perturbed links and the links that are at most ``minimizationRegionRadius`` :math:`\AA{}` away from them.
- the *top side links*, which are those links that worsen their energy the most between the beginning of a PELE step and the end of the perturbation and ANM phase. The number of these links is selected by the PELE parameter ``numberOfTopSideSelectedLinks`` (:ref:`sec-peleParameters-numberOfTopSideSelectedLinks`). Also the links that are at most ``topSideRadius`` :math:`\AA{}` away from the *top side links* (:ref:`sec-peleParameters-topSideRadius`) will be included.
- those atoms selected in ``includeInMinimization`` (:ref:`sec-minimization-includeInMinimization`).
- removing, from all the atoms previously selected, those atoms in ``doNotIncludeInMinimization`` (:ref:`sec-minimization-doNotIncludeInMinimization`).

Besides the permanent constraints (:ref:`sec-permanentConstraints`), the following constraints are considered during the minimization phase of a PELE simulation:

- Constraints on ANM nodes (see :ref:`sec-anm-relaxationSpringConstant`).
- If the perturbation is active, constraints on the perturbed atoms to their center of mass (see :ref:`sec-peleParameters-perturbationCOMConstraintConstant`).

Algorithm selection
-------------------

algorithm
^^^^^^^^^

Use: It sets the minimization algorithm that will be used.

Possible options:

-  "TruncatedNewton": Newton's method minimizes :math:`Q(x_1,\ldots,x_n) \simeq F(x_1,\ldots,x_n)`, where Q is a quadratic function. Truncated Newton
   consists of two sub-algorithms, the first one controlling the entire
   minimization and the second one iteratively minimizing Q.
-  "ConjugatedGradient": Minimizes the energy using the
   conjugated-gradient algorithm.
-  "SteepestDescent": Minimizes the energy using the steepest-descent
   method.

Default value: "TruncatedNewton"

Selection of links to minimize
------------------------------

These selection options appear at the peleSimulation level. Example:

.. code-block:: json

      "commands" : [
         {
            "commandType" : "peleSimulation",
            "selectionToPerturb" : {
               "chains" : {
                  "names" : [ "L" ]
               }
            },
            "PELE_Output" : {
            ...
            },
            "PELE_Parameters" : {
            ...
            },
            "includeInMinimization": { "chains" : "all" },
            "doNotIncludeInMinimization": { "chains": { "names" : ["A"] }},
            "Minimizer" : {
               ...
            }
            ...
         }
       ]

.. _sec-minimization-includeInMinimization:

includeInMinimization
^^^^^^^^^^^^^^^^^^^^^

Use: Selection string used to select the atoms to which the minimization
is going to be performed. Make sure this parameter does not interfere
with the minimizationRegionRadius PELE parameter
(:ref:`PeleParameters <sec-peleParameters>`).

You can find examples of selections in (:ref:`Selection
Examples <sec-selectionExamples>`).

Parameter: Selector \* doMinimizeSelector

Default value: Select the whole system (equivalent to '{"chains":
"all"}').

.. _sec-minimization-doNotIncludeInMinimization:

doNotIncludeInMinimization
^^^^^^^^^^^^^^^^^^^^^^^^^^

Use: Selection string used to select the atoms that will be frozen
during the minimization.

Notice that doNotIncludeInMinimization prevails over
includeInMinimization, i.e. if an atom or link is selected by both
selection, it will be omitted. Therefore, check it does not deselects
the whole system.

You can find examples of selections in (:ref:`Selection Examples <sec-selectionExamples>`).

Parameter: Selector \* doNotMinimizeSelector

Default value: The empty selection (equivalent to '"empty"').

.. _sec-minimization-parameters:

Parameters
----------

EnergyDifference
^^^^^^^^^^^^^^^^

Parameter: double energyTol

Use: Used as a convergence criterion in the different minimization
algorithms.

Note for developers: Used in ConjugatedGradient::minimize,
SteepestDescent::minimize and
TruncatedNewton::runMinimizationIterations.

Units: kCal/mol

Range: [0, inf)

Default value: 1.0

Plop info: Not accessible from control file (hardcoded).

MinimumRMS
^^^^^^^^^^

Parameter: double rmsMin

Use: Used as a convergence criterion in the different minimization
algorithms.

Note for developers: Used in ConjugatedGradient::minimize,
SteepestDescent::minimize and TruncatedNewton::hasOuterLoopConverged.

Units: :math:`(\text{kCal}/\text{mol}/\AA{})^2`

Range: (0, inf)

Default value: 0.1

Plop info:

-  Plop control file name: 'rmsg'
-  Plop parameter name: min\_params%rmsmin
-  Plop default value: 0.001

MaximumMinimizationIterations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameter: unsigned int iterMax

Use: Used as a convergence criterion in the different minimization
algorithms.

Note for developers: Used in ConjugatedGradient::minimize,
SteepestDescent::minimize and
TruncatedNewton::runMinimizationIterations.

Units: Dimensionless

Range: [0, inf)

Default value: 3

Plop info:

-  Plop control file name: 'iter'
-  Plop parameter name: min\_params%iter\_max
-  Plop default value: 3

Truncated Newton Parameters
---------------------------

MaximumNewtonIterations
^^^^^^^^^^^^^^^^^^^^^^^

Parameter: unsigned int maximumNewtonIterations

Use: Number of times that the outer loop is going to be run.

Note for developers: Used in TruncatedNewton::runOuterLoop

Units: Dimensionless

Range: [0, inf)

Default value: 65

Plop info:

-  Plop control file name: 'mxitn'
-  Plop parameter name: min\_params%mxitn
-  Plop default value: 65

nonBondingListUpdatedEachMinStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameter: bool updateNblist

Use: Controls whether the NB list is updated after each minimization
step or not.

Note for developers: Used in
TruncatedNewton::revertToPreviousValidStateInCaseOfMinimizationFailure,
TruncatedNewton::prepareOuterLoopNextIteration and
TruncatedNewton::runOneIteration.

Units: Dimensionless

Range: True or false

Default value: true

Plop info:

-  Plop control file name: 'nbup'
-  Plop parameter name: min\_params%update\_nblist
-  Plop default value: True

alphaUpdated
^^^^^^^^^^^^

Parameter: bool updateAlpha

Use: Controls whether the Born alpha radius is updated after each
minimization step or not.

Note for developers: Used in TruncatedNewton::runOneIteration and
TruncatedNewton::revertToPreviousValidStateInCaseOfMinimizationFailure.

Units: Dimensionless

Range: True or false

Default value: false

Plop info:

-  Plop control file name: 'alphaup'
-  Plop parameter name: min\_params%update\_alpha
-  Plop default value: False

sgbUpdated
^^^^^^^^^^

Parameter: bool updateSgb

Use: Controls whether SGB cell list is updated after each minimization
step or not.

Note for developers: Used in TruncatedNewton::runOneIteration and
TruncatedNewton::revertToPreviousValidStateInCaseOfMinimizationFailure.

Units: Dimensionless

Range: True or false

Default value: true

Plop info:

-  Plop control file name: 'gbup'
-  Plop parameter name: min\_params%update\_sgb
-  Plop default value: True

iterationsBetweenNBlistLongUpdate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This parameter is currently not assignable through the control file.

Parameter: unsigned int iterationsBetweenNonBondingLongUpdate

Use: The list of non bonded long is updated every
iterationsBetweenNBlistLongUpdate.

Note for developers: Used, indirectly, in
OuterLoop::prepareForNextIteration().

Units: Dimensionless

Range: [1,inf)

Default value: 2

Plop info:

-  Plop control file name: 'niterup'
-  Plop parameter name: min\_params%nfull
-  Plop default value: 2

