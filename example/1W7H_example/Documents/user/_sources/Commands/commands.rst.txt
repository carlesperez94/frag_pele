.. _sec-commands:

************************
Commands block in Pele++
************************

A simulation is made up of a set of different commands. Once they all
finish, the simulation also finishes.

commands
========

Array containing the list of all the commands.

commandType
-----------

Every command must contain this field. We can launch different kinds of
command:

-  "peleSimulation": Most important block. It performs a :ref:`Pele simulation <sec-PELESimulation>`.
-  "energyComputation": Performs an energy computation.
-  "minimization": Performs a minimization.
-  "anmComputation": Performs an ANM computation.
-  "sideChainSimulation": Performs a Side Chain Prediction simulation.

Example 1
---------

Energy computation:

.. code-block:: json

   {
      "controlFileSavingPath" : "../../simulations/ain/originalControlFile.conf",
      "Initialization" : {
         "Complex" : {
            "files" : [
               {
                  "path": "../proteins/ain_fixed.pdb"
               }
            ]
         },
         "ForceField" : "OPLS2005",
         "Solvent" : {
            "ionicStrength" : 0.250,
            "solventType" : "VDGBNP",
            "useDebyeLength" : true
         }
      },
      "commands" : [
         {
          "commandType":"energyComputation"
         }
      ]
   }

Example 2
---------

We can link together different commands. For example, we can minimize
the structure before running a Pele simulation.

.. code-block:: json

   {
       ...
       "commands": [
           {
               "commandType": "minimization",
               "Minimizer": 
               {
                   "algorithm": "TruncatedNewton",
                   "parameters": 
                   { 
                       "MinimumRMS": 0.05
                   }
               }
           },
           {
               "commandType": "peleSimulation",
               ...
           }
       ]
   }

