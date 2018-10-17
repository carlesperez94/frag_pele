.. _sec-sensors:

Sensors Block in Pele++ control file
====================================

Pele++ sensors measure the value of different observable variables in
the landscape that is being explored.

Depending on what they observe, Pele++ sensors can be classified in two
groups:

-  Metrics: They measure the value of external observable variables,
   such as, binding energy, rmsd, distance, etc.

-  StateTrackers: They track the value of internal simulation variables,
   such as, steps, accepted steps, current energy, etc.

Each type of sensors has its own configuration block in the
configuration file which are defined at Pele task level, so you can
observe, if you wish, different variables in each different Pele task.

The values measured by any Sensor can be used in any logical condition
as long as the Sensors have been previously tagged (check the following
examples). To refer to the sensors, the tag name is used (in StateTrackers, this corresponds to the variable name). Notice that the tag name is also used in the metrics file as the header for that metric, unless it is not specified, in which case a standard description for the metric is used.

To learn how to write logical conditions based on Sensor measurements
see :ref:`Conditions By Example <sec-conditionsByExample>`,

Contents:

.. toctree::
    :maxdepth: 2

    MetricTypes
    StateTrackerTypes

