.. _sec-exitConditions:

Exit Conditions
===============

If any of the conditions in "exitConditions" list evaluates to true, the
current PeleTask ends.

Let's see an example of a PeleTask that uses an "exitConditions" block:

.. code-block:: json

   "PeleTasks":
   [
       {
           "metrics": [
               {
                   "type": "rmsd",
                   "tag": "rmsd_2",
                   "Native": {
                       "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
                   },
                   "selection": {
                       "atoms": {
                           "names": [ "_CA_", "_C4'", "_C2_", "_P__"]
                       }
                   }
               }
           ],
           "exitConditions": [
               "rmsd_2 < 1.1"
           ]
       }
   ]

In this example we have only one PeleTask block in which we have defined
a Metrics block and an "exitConditions" block. In the Metrics block a
RMSD metric is defined and tagged as "rmsd_2".

Now let's take a closer look to the "exitConditions" block:

.. code-block:: html

   "exitConditions": [
       "rmsd_2 < 1.1"
   ]

The "exitConditions" block is a list of conditions (Json lists are
defined between square brackets), therefore it only contains one
condition in this case.

Notice that the condition uses a Metric that has been previously tagged
when it was defined in the :ref:`metrics block <sec-metricTypes>`. Conditions can only work with tagged Metrics or StateTrackers names. To learn more about conditions, check :ref:`Conditions By Example <sec-conditionsByExample>`.

Another example using StateTrackers is:

.. code-block:: json

   "exitConditions": [
           "numberOfAcceptedPeleSteps > 1"
   ]

The "exitConditions" block can contain as many conditions as you like,
as long as they are separated by commas and each of them bounded between
double quotation marks ("). Of course the
:ref:`Sensors <sec-sensors>` (Metrics and/or StateTrackers)
must have been defined and tagged before using them.

