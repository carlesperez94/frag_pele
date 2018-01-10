.. _sec-stateTrackerTypes:

StateTracker types in Pele++ control file
=========================================

The StateTrackers are defined inside the "track" block as is explained
in :ref:`Sensors <sec-sensors>`.

A general state tracker is define as follows:

.. code-block:: json

   {
       "variable": "variableName",
       "logIt": false
   }


-  The "variable" parameter indicates the Pele++ variable that is
   followed. It's different an unique for each type of StateTracker and
   it's automatically used as a tag in logical conditions. This fact is
   very important because the names used to define StateTracker types
   can not be used to tag Metrics.

-  The "logIt" parameter is optional and it indicates if we'd like the
   variable values to be logged or not. By default is set to true, so
   they are always logged. If you want them not to be logged, set
   "logIt" to false.

Currently available trackers are:

-  currentEnergy (available and logged by default; cannot be
   configured).
-  numberOfAcceptedPeleSteps (available and logged by default; cannot be
   configured).
-  currentStep.
-  stepsSinceLastAccepted.

Following are the types of state trackers:

Current energy
--------------

It tracks the current energy variable in PeleTask. Its variable name is
currentEnergy.

This tracker is active and logged by default, and its logging property
cannot be changed.

Steps since last accepted
-------------------------

It tracks the steps of the simulation that has passed since the last
accepted step. ###Example.

.. code-block:: json

   {
       "variable":"stepsSinceLastAccepted",
       "logIt": false
   }

In this case, the values of the current energy are not logged.

Current step
------------

It tracks the steps of the simulation. When evaluated, it considers the
number of the step just about to start.

Example
^^^^^^^

.. code-block:: json

   {
       "variable":"currentStep",
       "logIt": false
   }

In this case, the values of the current step are not logged (actually,
the step number is logged under column "Step", but it refers to the step
just finished; if logIt is set to true, an additional "currentStep"
column would be added, with the value of the step about to start).

Number of accepted PELE steps
-----------------------------

It tracks the number of accepted PELE steps during the simulation. Its
variable name is numberOfAcceptedPeleSteps.

This tracker is active and logged by default, and its logging property
cannot be changed.

StateTrackers block example
---------------------------

.. code-block:: json

   "track":
   [
       {
           "variable": "currentEnergy"
       },
       {
           "variable": "stepsSinceLastAccepted",
           "logIt": false
       },
   ]

In this second example, we show a "track" configuration block which is a
list that contains the definition of two StateTrackers.

The first StateTracker block of the list is:

.. code-block:: json

   {
       "variable": "currentEnergy"
   }

The StateTrackers track variables of the current state of a
PeleSimulation.

You only need to indicate the variable you'd like to track. For a list
of the state variables that can be tracked check this section, above.

In this case, we are tracking the "currentEnergy" state variable (which
is the energy of the current Pele simulation step).

Each StateTracker has a "variable name" which identifies it univocally
and which often coincides with the variable name in our code, (in
previous sections you'll
also find the corresponding variable name for each StateTracker).

In the first example we saw how Metrics need to be tagged in order to be
used in logical conditions. StateTrackers don't need to be tagged,
though, because its variable name acts as a unique identifier that is
then used in logical conditions.

The second StateTracker block is:

.. code-block:: json

   {
       "variable": "stepsSinceLastAccepted",
       "logIt": false
   }

Here we are tracking the "stepsSinceLastAccepted" state variable (which
is the number of steps since the last Metropolis acceptance event).

Notice how we have disabled it being logged in the PeleSimulation
reports by setting "logIt" to false. As in the case of Metrics the
"logIt" default value is set to true, so all StateTrackers values are
logged unless you explicitly indicate otherwise.
