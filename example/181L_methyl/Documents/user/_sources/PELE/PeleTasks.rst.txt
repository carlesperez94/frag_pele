.. _sec-peleTasks:

PeleTasks
=========

A peleSimulation has a "PeleTasks" block. It is an array comprised of all
the tasks that the program has to carry out before a command (usually a
simulation) is considered to be finished [#f1]_. Each task is a block, so it
should be enclosed in curly braces.

Contents:
---------

Pele_parameters and algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It can redefine both Pele\_parameters and algorithms (e.g.
ligandPerturbation, ANM, sideChainPrediction, Minimization...).

If an algorithm is defined at PeleTask level, it uses this one instead
of the one at peleSimulation level.

Example
^^^^^^^

Side chain prediction is defined at peleSimulation level, so, unless
overriden, this algorithm is going to be used thorughout all the
simulation. However, the first PeleTask redefines it, so the latter one
is going to be used for that specific task.

.. code-block:: json

   "commands": [
           {
               "commandType": "peleSimulation",
               ...
                "SideChainPrediction": {
                   "algorithm": "zhexin"
                   },
                   "parameters": {
                       "randomize": false,
                       "numberOfIterations": 1
                   }
               },
               ...
               "PeleTasks": [
                   {
                       ...
                        "SideChainPrediction": {
                           "algorithm": "zhexin"
                           },
                           "parameters": {
                               "randomize": true,
                               "numberOfIterations": 3
                           }
                       },
                       ...
                   },
                   ...
               ]
   ]

Constraints
-----------

A PELE Task can include additional constraints to those already defined
for the complete PELE simulation. Check the :ref:`Constraints section <sec-permanentConstraints>`.

metrics
-------

It is an array that contains all the metrics to be measured in that
task.

Check the :ref:`metrics documentation <sec-metricTypes>` for
more information.

jumpif/jumpController block
---------------------------

Makes explorer jump if certain conditions are met.

Check the :ref:`jumpController documentation <sec-jumpController>` for more
information.

parametersChanges
-----------------

Array of conditions that make some defined parameters to be changed.

Check the :ref:`dynamic changes section <sec-PELEDynamicChangesInSimulationParameters>` for more information.

exitCondition
-------------

When the condition is met, it makes the task finish.

Check the :ref:`exit conditions section <sec-exitConditions>` for more information.

.. rubric:: Footnotes

.. [#f1] It can also finish because numberOfPeleSteps is reached.

