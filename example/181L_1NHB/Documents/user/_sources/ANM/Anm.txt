.. _sec-anm:

ANM Block in Pele++ control file
================================

This block configures the ANM phase in a ``peleSimulation`` command, or the ANM execution of an ``anmComputation`` command. In a ``peleSimulation``, this phase will be active depending on the ``anmFrequency`` PELE parameter (:ref:`sec-peleParameters-anmFrequency`).

.. warning::

  Even if ANM is deactivated in your simulation, its configuration may affect the relaxation phase minimization in a PELE simulation, since the ANM nodes are constrained according to the `relaxationSpringConstant` parameter (:ref:`sec-anm-relaxationSpringConstant`).

Location
--------

The ANM block can be defined, as all algorithms, at two different
levels:

Inside the PeleSimulation command block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parameters defined at PeleSimulation level will be common for all
the tasks.

.. code-block:: json

   ...
   "commands":
   [
       {
           "commandType": "peleSimulation",

           ...
           ...

          
           "ANM": {
               "algorithm": "ALPHACARBONS",
               "options": {
                   "moveMagnitudeGeneration": "noRandom",
                   "modesMixingOption": "doNotMixModes",
                   "pickingCase": "LOWEST_MODE",
                   "directionGeneration": "useAlwaysTheSame"
               },
               "parameters": {
                   "eigenUpdateFrequency": 1,
                   "modesChangeFrequency": 1,
                   "numberOfModes": 8,
                   "thermalScaling": true,
                   "displacementFactor": 1.5,
                   "cutoff": 12,
                   "normalMode": 0,
                   "bondedConstantForHessian": 1.667,
                   "relaxationSpringConstant": 40,
                   "equilibriumParameter": 0,
                   "steeringForce": 30
               },
               "Native": {
                   "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
               },
               "anmMinimizer": {
                   "algorithm": "TruncatedNewton",
                   "parameters": {
                       "EnergyDifference": 1,
                       "MaximumMinimizationIterations": 1,
                       "MaximumNewtonIterations": 100,
                       "MinimumRMS": 0.01,
                       "nonBondingListUpdatedEachMinStep": false,
                       "alphaUpdated": false,
                       "sgbUpdated": false
                   }
               }
           }
           ...
           ]
       }
   ]

Inside a PeleTask block
^^^^^^^^^^^^^^^^^^^^^^^

The parameters defined at PeleTask level will be used only by the
PeleTask in which they are defined.

.. code-block:: json

   ...
   "commands":
   [
       {
           "commandType": "peleSimulation",
           ....
           ....
           "PeleTasks": [
               {
             ....
                   ....
                   "ANM": {
                       "algorithm": "ALPHACARBONS",
                       "options": {
                           "moveMagnitudeGeneration": "noRandom",
                           "modesMixingOption": "doNotMixModes",
                           "pickingCase": "LOWEST_MODE",
                           "directionGeneration": "useAlwaysTheSame"
                       },
                       "parameters": {
                           "eigenUpdateFrequency": 1,
                           "modesChangeFrequency": 1,
                           "numberOfModes": 8,
                           "thermalScaling": true,
                           "displacementFactor": 1.5,
                           "cutoff": 12,
                           "normalMode": 0,
                           "bondedConstantForHessian": 1.667,
                           "relaxationSpringConstant": 40,
                           "equilibriumParameter": 0,
                           "steeringForce": 30
                       },
                       "Native": {
                           "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
                       },
                       "anmMinimizer": {
                           "algorithm": "TruncatedNewton",
                           "parameters": {
                               "EnergyDifference": 1,
                               "MaximumMinimizationIterations": 1,
                               "MaximumNewtonIterations": 100,
                               "MinimumRMS": 0.01,
                               "nonBondingListUpdatedEachMinStep": false,
                               "alphaUpdated": false,
                               "sgbUpdated": false
                           }
                       }
                   }
             ...
             ...
               }
           ]
       }
   ]

In both PeleSimulation command block and PeleTask block
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, the parameters defined at PeleSimulation level will be
common for all the tasks, while the parameters defined at PeleTask level
will be used only by the PeleTask in which they are defined.

If a given parameter is defined in both blocks, the value defiened in
the PeleTask level has preference.

Algorithm selection
-------------------

algorithm
^^^^^^^^^

Use: It sets the ANM algorithm that will be used.

Possible options:

-  "CARTESIANS" An ANM using cartesian coordinates on the selected nodes
   (see "nodes").
-  "ALPHACARBONS" Special case of "CARTESIANS" algorithm where the nodes
   are only Cα.

When using either ``CARTESIANS`` or ``ALPHACARBONS`` options, you can specify to
load the eigenvalues and eigenvectors from a file. See "Modes
computation" below.

Default value: "ALPHACARBONS"

Nodes selection
---------------

ANM works on a set of nodes. First, a set of links is created by removing the ``linksToOmit`` (:ref:`sec-anm-linksToOmit`) from the ``linksToInclude`` (:ref:`sec-anm-linksToInclude`) set. Then:

- for the ``ALPHACARBONS`` algorithm, all Cα atoms of that set of links will be used as nodes for ANM.
- for the ``CARTESIANS`` algorithm, the ``nodes`` selection is applied to that set of links, and will be used as nodes for ANM.

Notice that the ``nodes`` are disposed linearly, as they appear in the input PDB, and that consecutive ``nodes`` are considered bonded if they are interacting (that is, under the given :ref:`sec-anm-cutoff`). This makes sense in the general case; also notice that this is only significant when the elastic network model used is "powerContactModel" (see :ref:`sec-anm-elasticNetworkModel`) and the :ref:`sec-anm-bondedConstantForHessian` is greater than zero. For multichain systems, this implies that the last node from a chain and the first node from the next chain may be considered bonded. Also, if the :ref:`sec-initialization-general-allowMissingTerminals` option of the :ref:`sec-initialization-general` section is active, it considers as bonded those nodes opening and closing a gap, depending on :ref:`sec-anm-cutoff`.

nodes
^^^^^

Use: String used to select the nodes that the ANM will use (see :ref:`sec-selectionExamples`). This parameter is not available for "ALPHACARBONS", since in this case it is equivalent to selecting all :math:`C\alpha`, as shown in example 2 below.

Default value: All atoms for ``CARTESIANS``.

Example 1
"""""""""

This is an example of an ANM simulation with DNA:

.. code-block:: json

   "ANM":
   {
       "algorithm": "CARTESIANS",
       "nodes": {
           "atoms": {
               "names": [
                   "_P__",
                   "_C2_",
                   "_C4'"
               ]
           }
       },
       ...
   }

Example 2
"""""""""

.. code-block:: json

    "algorithm": "ALPHACARBONS"

is equivalent to using:

.. code-block:: json

   "ANM":
   {
       "algorithm": "CARTESIANS",
       "nodes": {
           "atoms": {
               "names": [
                   "_CA_"
               ]
           }
       },
       ...
   }

Modes computation
-----------------

.. _sec-anm-preloadedModesIn:

preloadedModesIn
^^^^^^^^^^^^^^^^

Parameter: string preloadedModesIn

Use: By default, PELE computes the eigenvalues and eigenvectors.
However, when the algorithm specified is CARTESIANS or ALPHACARBONS,
then they can be loaded from an external file. In order to do so, set a
value for this parameter, indicating the path to that file. The path may
be absolute, or relative to the execution directory.

.. warning::

   If a value is set for this parameter, it will introduce a redundancy in
   the input of PELE. For one side, the nodes to use when computing ANM
   will be indicated in the control file, and they also will be indicated
   in the nmd input file. PELE checks this, and complains if the list of
   nodes is not the same in both files. You must ensure that your nmd file
   has the same nodes as the ones that you indicate in your control file.

Units: Path to a nmd/prody file

Default value: No value.

More information about the input format at :ref:`sec-fileFormats-NMD`.


Example
"""""""

As an example of the warning above, if you specify the following block
in your control file, then your nmd file must contain the information
for all the alphacarbons of the system, except for those from link 71 to
76 of chain A (both included).

.. code-block:: json

   "ANM":
   {
       "algorithm": "ALPHACARBONS",
       "linksToOmit":
       {
           "links":
           {
               "ranges": [ "A:71 A:76" ]
           }
       },
       "prelodadedModesIn" : "/path/to/some/nmd/file",
       "parameters": {
           ...
       },
       ...
   }

Options
-------

Main mode selection
^^^^^^^^^^^^^^^^^^^

pickingCase
"""""""""""

Parameter: PickingType pickingType

Use: It is used to choose the picking method to select a mode.

Note for developers: Used in AnmParametersBuilder::setParameters
function.

Possible options:

-  "LOWEST\_MODE" Selects the first mode.
-  "RANDOM\_MODE" Selects a random mode.
-  "DEFINED\_MODE" Selects a specific mode, that is the one indicated by
   "definedMode" value.
-  "MODE\_USING\_NATIVE" Currently being implemented

Units: No units

Default value: "RANDOM\_MODE"

Plop info:

-  Plop control file name: 'anm\_altm\_type'
-  Plop parameter name: liga\_params%anm\_altm\_type
-  Plop default value: 3 ("RANDOM\_MODE")

There is no mapping between Plop and Pele++ because direction and
picking are mixed in Plop producing a "case explosion".

definedMode
"""""""""""

Parameter: unsigned int definedMode

Use: Index of the chosen mode. It is only used when the "DEFINED\_MODE"
picker option is selected.

Note for developers: Used in PickerBuilder::createPicker function.

Range: An integer [1, inf)

Note for developers: To avoid confusing the users, in the control file
the definedMode index starts counting from 1 instead of 0. However, in
the code the definedMode parameter is set to "definedMode" - 1, so that,
for instance:

        "definedMode": 1 sets the definedMode parameter to 0, or

        "definedMode": 5 sets the definedMode parameter to 4.

Units: No units

Default value: 1 in control file

Plop info:

-  Plop control file name: 'mode'
-  Plop parameter name: anm\_params%mode
-  Plop default value: 1

Direction generation
^^^^^^^^^^^^^^^^^^^^

directionGeneration
"""""""""""""""""""

Parameter: DirectionCase directionCase

Use: It selects the method to generate the direction of the move vector.

Note for developers: Used in DirectionSelector::getDirection function.

Possible options:

-  "useAlwaysTheSame" It always uses the same direction. This direction
   is the value indicated by "initialDirection"
-  "oscillate" The new direction is equal to minus the last direction.
   The initial direction is indicated by "initialDirection"
-  "random" The new direction is a random one.

Units: No units

Default value: "random"

Plop info: In Plop the direction generation is mixed with other features
producing a "case explosion".

initialDirection
""""""""""""""""

Parameter: double initialDirection

Use: Initial direction of the move vector. It is required when you use
"useAlwaysTheSame" or "oscillate" direction generation options.

Note for developers: Used in DirectionSelector class constructor.

Range: It can be 1.0 or -1.0.

Units: No units

Default value: 1.0

Plop info: There is no equivalence in Plop (see "directionGeneration"
Plop info)

Modes mixing
^^^^^^^^^^^^


modesMixingOption
"""""""""""""""""

Parameter: MixCase mixCase

Use: Selects the method to mix the modes.

Note for developers: Used in AnmMoveVectorCalculator::computeMoveVector
function.

Possible options:

-  "mixAllModesEquallyRandom": All modes are mixed randomly with the
   same weight.
-  "mixMainModeWithOthersModes": Some contribution of the non selected
   modes is added to the chosen mode to create the move vector. The
   value of the parameter "mainModeWeightForMixModes" is the weight of
   the chosen mode and (1.0 - "mainModeWeightForMixModes") is the weight
   of the contribution of the rest of the modes.
-  "doNotMixModes" : The modes are not mixed, i. e. only the chosen mode
   is used to compute the movement vector.

Units: No units

Default value: "mixMainModeWithOthersModes"

Plop info: There is not an equivalent parameter in Plop. Plop selects
the mixing behaviour according to the value of several parameters which
is not a good idea. Mapping between Plop and Pele++: > if
anm\_params%xmoddir is true -> "mixAllModesEquallyRandom"

    if anm\_params%xmoddir is false and anm\_params%mix\_modes > 0.0 ->
    "mixMainModeWithOthersModes"

    default behaviour -> "doNotMixModes"

Movement magnitude computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

moveMagnitudeGeneration
"""""""""""""""""""""""

Parameter: MoveMagnitudeGenerationOption moveMagnitudeGenerationOption

Use: Selects the method to generate the movement magnitude that will be
used in ANM.

Note for developers: Used in
AnmMagnitudeCalculator::computeDisplacementMagnitude function.

Possible options:

-  "noRandom": The movement magnitude is equal to the value indicated in
   the control file by "displacementFactor".
-  "biasedRandom": The movement magnitude is equal to
   "displacementFactor" \* random number in [0.625, 1.0)
-  "scaledBiasedRandom": The movement magnitude is equal to
   "displacementFactor"/(0.6 \* chosen eigenvalue/ biggest eigenvalue +
   0.4) \* random number in [0.625, 1.0)
-  "random": The movement magnitude is equal to "displacementFactor" \*
   random number in [0.0, 1.0)

Units: No units

Default value: "noRandom"

Plop info:

-  Plop control file name: 'rand\_move'
-  Plop parameter name: anm\_params%rand\_move
-  Plop default value: 0 -> equivalent to "noRandom"

Mapping between Plop and Pele++:

-  0 -> "noRandom"
-  1 -> "biasedRandom"
-  2 -> "scaledBiasedRandom"
-  3 -> "random"

.. _sec-anm-elasticNetworkModel:

Elastic network model
^^^^^^^^^^^^^^^^^^^^^

Parameter: elasticNetworkModel

Use: Sets the elastic network model that will be applied for the ANM
computation. Affects the calculation of the constant of force applied to
the hessian matrix within the ANM algorithm.

Note for developers: Check
ModesCalculatorBuilder::getConstForceCalculator

Possible options:

-  "powerContactModel": A power-based function will be applied to
   calculate the constant of force
-  "exponentialContactModel": An exponential function will be applied to
   calculate the constant of force. Check the parameter
   "exponentialDistanceDecay", that can be applied in this case to
   adjust the exponent's value.

Units: It's a string

Default value: "powerContactModel"

Plop info:

Mapping between Plop and Pele++:

.. _sec-anm-parameters:

Parameters
----------

Movement magnitude computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

displacementFactor
""""""""""""""""""

Parameter: double displacement

Use: Nodes displacement factor.

Note for developers: Used in
AnmMagnitudeCalculator::computeDisplacementMagnitude.

Units: :math:`\AA{}`

Default value: 1.0

Range: (-inf, inf)

Plop info:

-  Plop control file name: 'move\_ca'
-  Plop parameter name: anm\_params%move\_ca
-  Plop default value: 0.5

Hessian matrix computation and diagonalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

numberOfModes
"""""""""""""

Parameter: unsigned int numberOfModes

Use: Number of desired eigen values.

Units: No units

Default value: 6

Range: 1,..., 3N-6 (N = number of atoms)

Plop info:

-  Plop control file name: 'neig'
-  Plop parameter name:anm\_params%neig
-  Plop default value: 6

thermalScaling
""""""""""""""

Parameter: bool useThermalScalingInPicking

Use: If true, the eigenvectors are thermally scaled, otherwise, they are
not.

Note for developers: Used in LapackCartesianModesCalculator to choose
the kind of normalization to apply to the eigenvectors.

Range: True or false

Default value: True

Plop info:

-  Plop control file name: 'thermal'
-  Plop parameter name: anm\_params%thermal
-  Plop default value: true

.. _sec-anm-cutoff:

cutoff
""""""

Parameter: double cutoff

Use: Used to determine the pairs of nodes that have a harmonic
interaction between them in the ANM Hessian matrix generation.

Note for developers: Used in LapackCartesianModesCalculator

Range: [0, inf)

Units: :math:`\AA{}`

Default value: 12.0

Plop info:

-  Plop control file name: 'cutoff'
-  Plop parameter name: anm\_params%cutoff
-  Plop default value: 12.0

.. _sec-anm-bondedConstantForHessian:

bondedConstantForHessian
""""""""""""""""""""""""

Parameter: double bondedConstantForHessian

Use: Spring constant used to compute the forces between nodes to
generate the ANM Hessian matrix. It is only used in the "powerContactModel" (:ref:`sec-anm-elasticNetworkModel`), and only if it is greater than zero, otherwise a different parameter (:ref:`sec-anm-constantForHessian`) is used.

Note for developers: Used in
LapackCartesianModesCalculator::calculateEigValuesAndEigVectors
function.

Units: :math:`kT/\AA{}^2`

Range: [0, inf)

Default value: 0.0

Plop info:

-  Plop control file name: 'bkcons'
-  Plop parameter name: anm\_params%bkcons
-  Plop default value: 0.0

decayConstant
"""""""""""""

Parameter: double springConstantDecay

Use: Exponent of the decay (kx2./springConstantDecay) with the distance
of the spring energy. If it is equal to 1, the energy follows Hook's law
(grows with the square of the distance).

Note for developers: Used in
LapackCartesianModesCalculator::calculateEigValuesAndEigVectors
function.

Units: No units

Default value: 1.0

Range: (0, inf)

Plop info:

-  Plop control file name: 'kdecay'
-  Plop parameter name: anm\_params%kdecay
-  Plop default value: 1.0

.. _sec-anm-constantForHessian:

constantForHessian
""""""""""""""""""

Parameter: double constantForHessian

Use: Spring constant used to compute the forces between nodes to
generate the ANM Hessian matrix. In the powerContactModel (see :ref:`sec-anm-elasticNetworkModel`), it is used for non consecutive nodes when :ref:`sec-anm-bondedConstantForHessian` is greater than zero, and for all nodes when :ref:`sec-anm-bondedConstantForHessian` is zero. In the "exponentialContactModel", this constant is used for all nodes.

Note for developers: Used in
LapackCartesianModesCalculator::calculateEigValuesAndEigVectors
function.

Units: :math:`kT/\AA{}^2`

Range: [0, inf)

Default value: 1.667

Plop info:

-  Plop control file name: 'kcons'
-  Plop parameter name: anm\_params%kcons
-  Plop default value: 1.667

.. _sec-anm-constraints:

Constraints
^^^^^^^^^^^

steeringForce
"""""""""""""

Parameter: double steeringForce

Use: Force for the constraints used in the ANM inner minimization. It is
used to compute the elastic constant for the constraints depending on
the magnitude of the displacement: > springConstant =
steeringForce/displacementMagnitude

Note for developers: Used in AnmCalculator::minimize.

Units: :math:`kT/\AA{}` (force units)

Range: [0, inf)

Default value: 20.0

Plop info:

-  Plop control file name: 'fsteer'
-  Plop parameter name: anm\_params%fsteer
-  Plop default value: 20.0

.. _sec-anm-relaxationSpringConstant:

relaxationSpringConstant
""""""""""""""""""""""""

Parameter: double relaxationSpringConstant

Use: Anm constraints spring constant for relaxation step minimization.
If it is 0, no constraints are applied in the relaxation step
minimization.

Note for developers: Used in AnmCalculator::constrainCurrentAnmNodes
function.

Units: :math:`kT/\AA{}^2`

Range: (0, inf)

Default value: 0. By default, no constraints are applied in the
relaxation step minimization.

Plop info:

-  Plop control file name: 'caconst'
-  Plop parameter name: anm\_params%caconst
-  Plop default value: 0.0

equilibriumParameter
""""""""""""""""""""

Parameter: double equilibriumParameter

Use: Elastic constraints equilibrium position. It might be a distance or
an angle depending on the algorithm you choose.

Note for developers: Used in AnmCalculator::constrainCurrentAnmNodes and
AnmCalculator::minimize functions.

Units: :math:`\AA{}`

Range: (-inf, inf)

Default value: 0.0

Plop info: No equivalent in Plop

Frequencies
^^^^^^^^^^^

eigenUpdateFrequency
""""""""""""""""""""

Parameter: unsigned int eigenUpdateFrequency

Use: It is the frequency with which the eigen values and eigenvectors of
the ANM Hessian matrix are computed.

Note for developers: Used in AnmCalculator::isTimeToUpdateEigenVectors
function.

Range: An integer in [1, inf).

Units: No units

Default value: 1000000 (it means: calculate it the first time, never calculate it again).

Plop info:

-  Plop control file name: 'anm\_eig\_freq'
-  Plop parameter name: liga\_params%anm\_eig\_freq
-  Plop default value: 1

modesChangeFrequency
""""""""""""""""""""

Parameter: unsigned int modesChangeFrequency

Use: It is the frequency with which a new mode is selected.

Note for developers: Used in AnmCalculator::isTimeToPickNewMode
function.

Range: An integer in [1, inf).

Units: No units

Default value: 6

Plop info:

-  Plop control file name: 'anm\_altm\_freq'
-  Plop parameter name: liga\_params%anm\_altm\_freq
-  Plop default value: 1

Note: The eigenUpdateFrequency must be greater than or equal to the
modesChangeFrequency.

Modes mixing
^^^^^^^^^^^^

mainModeWeightForMixModes
"""""""""""""""""""""""""

Parameter: double mainModeWeightForMixModes

Use: It is used to add some contribution of the non selected modes to
the move vector. It represents the weight of the chosen mode.

Note for developers: Used in AnmMoveVectorCalculator::computeMoveVector.

Range: (0, 1].

Units: No units

Default value: 0.6

Plop info:

-  Plop control file name: 'mix\_modes'
-  Plop parameter name: anm\_params%mix\_modes
-  Plop default value: -1.0d0

exponentialDistanceDecay
^^^^^^^^^^^^^^^^^^^^^^^^

Parameter: double exponentialDistanceDecay

Use: It applies only if the option "elasticNetworkModel" has been set to
"exponentialNetworkModel". Then, given a value "d" for this parameter,
then the computation of the constant of force of each pair of atoms for
the hessian matrix will be done as follows:

.. code-block:: c++

   constant = parameter_constantForHessian/exp(pow(sqrt(distance_between_atoms)/d, 2.0))

   // (The above line comes directly from the code)

Note for developers: Used in
ExponentialConstForceCalculator::calculateConstForce.

Range: (-inf, inf).

Units: No units

Default value: 5.0

Plop info:

.. _sec-anm-linksToInclude:

linksToInclude
--------------

Selection of links to include in the ANM (see :ref:`sec-selectionExamples`).

Default: Everything but ligands (``"chains": "everythingButLigands"``).

.. _sec-anm-linksToOmit:

linksToOmit
-----------

Selection of links not to include in the ANM (see :ref:`sec-selectionExamples`).

Default: No selection (``"empty"``).

Example:

.. code-block:: json

   "ANM":
   {
       "algorithm": "ALPHACARBONS",
       "linksToOmit":
       {
           "links":
           {
               "ranges": [ "A:71 A:76" ]
           }
       },
       "parameters": {
           ...
       },
       ...
   }

.. _sec-anm-anmMinimizer:

anmMinimizer
------------

Minimizer block.

Check the minimization block documentation (:ref:`sec-minimization`) for more information.

The default minimization type is "TruncatedNewton", with the following default parameters (other parameters have the same default as shown in the minimization block documentation):

- "EnergyDifference": 1.0
- "MinimumRMS": 0.1
- "MaximumMinimizationIterations": 1
- "MaximumNewtonIterations": 30,
- "nonBondingListUpdatedEachMinStep": false
- "alphaUpdated": false
- "sgbUpdated": false
- "iterationsBetweenNBlistLongUpdate": 2

The same atoms as those shown in the minimization block documentation are considered for minimization, except that the top side links and nearby links are not added to the list of atoms/links to minimize.

This minimization makes the actual ANM movement by constraining the ANM nodes to target coordinates (see :ref:`sec-anm-constraints`). If the perturbation is active, also the perturbation constraints will be on (see :ref:`sec-peleParameters-perturbationCOMConstraintConstant`). The permanent constraints (:ref:`sec-permanentConstraints`) are also on.

