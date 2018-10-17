.. _sec-PELEDynamicChangesInSimulationParameters:

PELE: Dynamic Changes in simulation parameters
==============================================

Some PELE simulation parameters can be changed on the fly based on the
evaluation of logical conditions that use values measured by PELE
sensors (metrics and state trackers, see
:ref:`Sensors <sec-sensors>`).

This document explains how to write those logical conditions in a
configuration file. To see which simulation parameters can be changed,
see :ref:`ParametersThatCanChange <sec-parametersThatCanChange>`.

The key of the block to dynamically change simulation parameters is:
"parametersChanges".

This block is optional, but if you use it, it must be included inside a
PeleTask block.

Let's see an example:

.. code-block:: json

   "parametersChanges": 
   [
       {
           "ifAnyIsTrue": [
               "rmsd > 1.0 and rmsd < 1.2",
               "rmsd_CA > 1"
           ],
           "doThesechanges": {
               "Minimizer::parameters": {
                   "alphaUpdated": false,
                   "sgbUpdated": false
               }
           }
       },
       {
           "ifAnyIsTrue": [
               "rmsd < 1.0 or rmsd <= 0.5",
               "not (rmsd_CA < 1 or currentEnergy > 3.0)"
           ],
           "doThesechanges": {
               "LigandPerturbation::parameters": {
                   "rotationScalingFactor": 0.4,
                   "translationRange": 3.5
               },
               "ANM::parameters": {
                   "thermalScaling": false,
                   "displacementFactor": 1.5
               }
           },
           "otherwise": {
               "Minimizer::parameters": {
                   "alphaUpdated": true,
                   "sgbUpdated": true
               }
           }
       }
   ]

In this case "parametersChanges" corresponds to a list (square brackets
define lists in JSON) of two parameter changes (each defined between
curly brackets).

The first parameters change block of the list in this example is:

.. code-block:: json

   {
       "ifAnyIsTrue": [
           "rmsd > 1.0 and rmsd < 1.2",
           "rmsd_CA > 1"
       ],
       "doThesechanges": {
           "Minimizer::parameters": {
               "alphaUpdated": false,
               "sgbUpdated": false
           }
       }
   }

and the second one is:

.. code-block:: json

   {
       "ifAnyIsTrue": [
           "rmsd < 1.0 or rmsd <= 0.5",
           "not (rmsd_CA < 1 or currentEnergy > 3.0)"
       ],
       "doThesechanges": {
           "LigandPerturbation::parameters": {
               "rotationScalingFactor": 0.4,
               "translationRange": 3.5
           },
           "ANM::parameters": {
               "thermalScaling": false,
               "displacementFactor": 1.5
           }
       },
       "otherwise": {
           "Minimizer::parameters": {
               "alphaUpdated": true,
               "sgbUpdated": true
           }
       }
   }

Each parameters change block is comprised of three blocks:

-  A list of conditions for the parameters changes to happen.
-  The changes that must be done if any of the conditions evaluates to
   true.
-  The changes that must be done otherwise, i.e., if all the conditions
   evaluate to false.

Let's have a closer look at the second parameters change block.

In this case, the list of conditions for the changes to happen is given
by:

.. code-block:: json

   "ifAnyIsTrue": 
   [
       "rmsd<1.0 or rmsd<=0.5",
       "not(rmsd_CA<1 or currentEnergy>3.0)"
   ]

As it says, if any of the conditions in the list (remember in JSON lists
are defined between square brackets) evaluates to true, the changes will
take place. To know how to write conditions check :ref:`Conditions By Example <sec-conditionsByExample>`.

The block that defines the changes that take place when any of the
previous conditions is true is:

.. code-block:: json

   "doThesechanges": {
       "LigandPerturbation::parameters": {
           "rotationScalingFactor": 0.4,
           "translationRange": 3.5
       },
       "ANM::parameters": {
           "thermalScaling": false,
           "displacementFactor": 1.5
       }
   }

As its name implies the "doThesechanges" block contains the sets of
parameter changes that will take place when the conditions for the
change are true. Inside the "doThesechanges" block, you should define
the different set of changes by indicating which parameters to change
and which are the new values of the parameters.

The different sets of parameters changes are separated by commas.

For instance, the first set of changes is:

.. code-block:: json

   "LigandPerturbation::parameters": {
       "rotationScalingFactor": 0.4,
       "translationRange": 3.5,
   }

where we are indicating that the parameters that will change are the
"rotationScalingFactor" and the "translationRange" of the
LigandPerturbation, being their new values 0.4 and 3.5, respectively.

Similarly, the second set of changes:

.. code-block:: json

   "ANM::parameters": {
       "thermalScaling": false,
       "displacementFactor": 1.5,
   }

indicates that the parameters that will change are the "thermalScaling"
and the "displacementFactor" of the ANM with new values of false and
1.5, respectively.

There is a list of which parameters can be changed available in:
:ref:`Parameters That Can Change <sec-parametersThatCanChange>`.

Finally, there is a block, "otherwise", that defines the changes that
take place when none of the previous conditions is true:

.. code-block:: json

   "otherwise": {
       "Minimizer::parameters": {
           "alphaUpdated": true,
           "sgbUpdated": true
       }
   }

The "otherwise" block contains the sets of parameter changes that take
place when all the conditions for the change are false. Inside it, you
should define the different set of changes by indicating which
parameters to change and which are the new values of the parameters.

.. code-block:: json

   "Minimizer::parameters": {
       "alphaUpdated": true,
       "sgbUpdated": true
   }

In this case, we are changing two parameters of the minimizer.

