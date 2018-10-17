.. _sec-nonBondingAlgorithms:

NonBonding in Pele++ control file
=================================

Configuration of the algorithm for non-bonding energy-related terms
calculation.

It appears in the "peleSimulation", "minimization", "anmComputation",
"sideChainSimulation" commands.

NonBonding v2 algorithm
-----------------------

In this algorithm you only define one cutoff

Example
^^^^^^^

.. code-block:: json

   "NonBonding": 
       {
           "nonBondingVersion": "NB_GEN_2",
           "cutoffList": 
           [
               {
                   "cutoff_v2": 50 
               } 
           ] 
       }

NonBonding multiscale algorithm
-------------------------------

This is the default algorithm chosen if you do not include this section in the control file.

In this algorithm we have six different cutoffs. It is the default
algorithm for the NonBonding computations. 

Example
^^^^^^^

Example with the default cutoffs.

.. code-block:: json

   "NonBonding": 
       {
                  "nonBondingVersion": "NB_GEN_MSCALE",
                  "cutoffList": 
                  [
                      {
                          "type": "shortNeutrumNeutrum",
                          "value": 10 
                      },
                      {
                          "type": "shortChargedNeutrum",
                          "value": 10 
                      },
                      {
                          "type": "shortChargedCharged",
                          "value": 15 
                      },
                      {
                          "type": "longNeutrumNeutrum",
                          "value": 15 
                      },
                      {
                          "type": "longChargedNeutrum",
                          "value": 20 
                      },
                      {
                          "type": "longChargedCharged",
                          "value": 30 
                      } 
                  ] 
               }

