.. _sec-metricTypes:

Metrics block in Pele++ control file
====================================

Pele++ has different types of Metrics. Metrics are a type of Sensor,
(see :ref:`Sensors <sec-sensors>`).

In general a Metric declaration in the control file is comprised of
three parts:

-  The type of Metric.
-  The specific parameters for that type of Metric (each metric has its
   own parameters).
-  A tag which is a name to identify the Metric in logical conditions
   (see :ref:`Conditions By Example <sec-conditionsByExample>`). The tag 
   cannot contain spaces. This part is optional; if present, it is used
   to identify the metric in the report file; otherwise, a generic 
   description is used.
-  A boolean value indicating if you want the Metric to be logged or
   not. This is also optional. By default all Metrics are logged (it's
   true by default).

This would be a general example of a Metric definition block:

.. code-block:: json

   {
       "type": "MetricType",
       "parameter1": "xxx",
       ...
       "parameterN": "xxx",
       "tag": "whateverYouWant",
       "logIt": false
   }

Metrics block example
---------------------

.. code-block:: json

   "metrics": [
       {
           "type": "rmsd",
           "Native": {
               "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
           }
       },
       {
           "type": "rmsd",
           "logIt": false,
           "Native": {
               "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
           },
           "selection": {
               "atoms": {
                   "names": [
                       "_CA_",
                       "_C4'",
                       "_C2_",
                       "_P__"
                   ]
               }
           },
           "tag": "selectionRmsd"
       }
   ]

In this first example, we show a "metrics" configuration block which is
a list (in `JSON <http://www.json.org/>`__ lists elements are surrounded
by []) that contains two Metrics.

The first metrics block of the list is:

.. code-block:: json

   {
       "type": "rmsd",
       "Native": {
           "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
       }
   }

It measures the RMSD between the current configuration and a given
"Native" and its values are logged each step.

The second metrics block of the list is:

.. code-block:: json

   {
       "type": "rmsd",
       "logIt": false,
       "Native": {
           "path": "src/Molecules/Test/Data/plop_isr_out.pdb"
       },
       "selection": {
           "atoms": {
               "names": [
                   "_CA_",
                   "_C4'",
                   "_C2_",
                   "_P__"
               ]
           }
       },
       "tag": "selectionRmsd"
   }

It measures the RMSD between a group of selected atoms in the current
configuration and their counterparts in a given "Native".

In this case its values are not logged because "logIt" is set to false.
By default "logIt" is set to true.

It's important to notice that this second metric has been tagged ("tag":
"selectionRmsd").

Doing so allows it to be used in logical conditions, being selectionRmsd
the variable that you'll have to use inside any logical conditions that
uses this metrics values.

There are many other Metric types that you can use. For a comprehensive
list, check below.


Metric types
------------

Now we will show you an example of each type of metrics so you can see
their names and specific parameters.

The currently available metrics are:

-  :ref:`Binding energy <sec-metricTypes-bindingEnergy>`
-  :ref:`Distance bewteen centers of mass <sec-metricTypes-comDistance>`
-  :ref:`SASA <sec-metricTypes-SASA>`
-  :ref:`Distance to point <sec-metricTypes-distanceToPoint>`
-  :ref:`Random number <sec-metricTypes-randomNumber>`
-  :ref:`RMSD <sec-metricTypes-RMSD>`

.. _sec-metricTypes-bindingEnergy:

Binding energy
^^^^^^^^^^^^^^


It computes the binding energy of a ligand.

Example
"""""""

.. code-block:: json

   {
       "type":"bindingEnergy",
       "boundPartSelection": {
           "chains": { "names": [ "L" ] }
           }
   }

In this example, we are indicating that we want to compute the binding
energy of the ligand "L". The "boundPartSelection" parameter is a
selection block, see :ref:`Selections By Example <sec-selectionExamples>` to see different
examples of selection blocks. The calculation is done by dividing the
system into two parts: that selected by the "boundPartSelection"
selection, and the rest. The "boundPartSelection" must be a single
chain, and it should correspond to the ligand, though it does not have
to.

.. _sec-metricTypes-comDistance:

Distance between centers of mass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


It computes the distance between the center of masses of two groups of
atoms. Each group of atoms is defined by a selection block, see
:ref:`Selections By Example <sec-selectionExamples>` to see
different examples of selection blocks. 

Example
"""""""

.. code-block:: json

   {
       "type":"com_distance",
       "selection_group_1":{
           "links": { "ids":["A:5"]}
       },
       "selection_group_2":{
           "links": { "ids":["A:3"]}
       }
   }

In this example, we are indicating that we want to compute the distance
between the center of masses of two Links from chain A, the ones with
number 5 and 3. Computing the distance between two Atoms is just a
special case of this Metric type. You just need to select just one Atom
in each selection block.

.. _sec-metricTypes-SASA:

SASA
^^^^

Computes the SASA for the indicated set of atoms. The SASA metric is a
normalized value between 0 and 1, and is calculated as follows:
SASAOfSelectedAtoms/SASAOfSelectedAtomsAlone. SASAOfSelectedAtomsAlone is the value calculated as if the selected
atoms were alone in the system, so it is the maximum value for these
atoms.

Example
"""""""

.. code-block:: json

   {
       "type":"sasa",
       "selection" : {
           "chains": {
                       "names": [ "L" ]
                   }
       }
   }

Computes the SASA for the chain named "L".

.. _sec-metricTypes-distanceToPoint:

Distance to point
^^^^^^^^^^^^^^^^^

It computes the distance from a fixed point to the center of mass of a
group of atoms. The group of atoms is defined used a selection block,
see :ref:`Selections By Example <sec-selectionExamples>` to
see different examples of selection blocks, whereas the point is given
as a list of coordinates. 

Example
"""""""

.. code-block:: json

   {
       "type":"distanceToPoint",
       "point": [0.5, 0.1, 0.0],
       "atoms": {
           "links" : { "ids":["A:2"]} 
       }
   }

In this example, we are indicating that we want to compute the distance
between the center of masses of the second Link of chain A and the point
with coordinates (0.5, 0.1, 0.0). Again to compute the distance between
an Atom and a given point, you need to select just one Atom in the
selection block.

.. _sec-metricTypes-randomNumber:

Random number
^^^^^^^^^^^^^


This is a special Metric. It doesn't actually compute anything. It just
returns a random number. We've added it, because random numbers are
often used in Pele++ simulations as part of logical conditions in
dynamic parameter changes.

Example
"""""""

.. code-block:: json

   {
       "type" : "random"
   }

.. _sec-metricTypes-RMSD:

RMSD
^^^^

It computes the RMSD between the Complex and a Native.

Notice that due to the numerical nature of the RMSD calculation, you may
sometimes obtain an almost-zero result (about 1e-15) for Complex and
Native having the same conformation when the real value should be zero.
To all extent, both results are the same and the metric calculation is
right.

You have to provide the following parameters:

Native
""""""

A section with a "path" element with the path to the PDB file containing
the native structure.

See :ref:`example <sec-metricTypes-RMSD-example1>`.

selection
"""""""""

A selection of the atoms to use for doing the RMSD calculation (see
:ref:`Selections By Example <sec-selectionExamples>`).

If the "doSuperposition" parameter is true, then an RMSD calculation
with superposition is performed.

If "includeHydrogens" is false, then only non-hydrogen atoms of the
selection are used for the RMSD calculation and, if "doSuperposition" is
true, for the superposition.

Notice that an atom, to be actually used, must have the same chain ID,
residue name, residue sequence number and atom name in both the native
and the non-native structure.

It is an error to have an actual selection of no atoms (since at least
one atom is needed for the metric calculation).

See :ref:`example <sec-metricTypes-RMSD-example2>`.

Default value: all atoms comprise the selection.

superpositionSelection
""""""""""""""""""""""

A selection of the atoms to use for doing the superposition before doing
the actual RMSD calculation (see :ref:`Selections By Example <sec-selectionExamples>`).

This parameter will only be used if the "doSuperposition" parameter is
true.

If "includeHydrogens" is false, then only non-hydrogen atoms of the
selection are used for the superposition.

Notice that an atom, to be actually used, must have the same chain ID,
residue name, residue sequence number and atom name in both the native
and the non-native structure.

It is an error to have an actual superposition of less than 3 atoms
(since otherwise there will be rotation ambiguity).

See :ref:`example <sec-metricTypes-RMSD-example2>`.

Default value: the same atoms as in selection.

doSuperposition
"""""""""""""""

This states whether a superposition must be done before calculation of
the RMSD. Notice that this superposition is only done for the purpose of
the RMSD calculation and does not affect the simulation coordinates.

The atoms of "superpositionSelection" (removing hydrogens if
"includeHydrogens" is true) are used to perform the superposition;
notice that, by default, "superpositionSelection" is the same as
"selection".

See :ref:`example <sec-metricTypes-RMSD-example3>`.

Range: true or false

Default value: true

includeHydrogens
""""""""""""""""

This states whether hydrogen atoms should be used for the RMSD
calculation and the optional superposition.

Range: true or false

Default value: false

.. _sec-metricTypes-RMSD-example1:

Example 1
"""""""""

.. code-block:: json

   {
       "type": "rmsd",
       "Native": {
           "path": "src/Molecules/Test/Data/AIN/native_fixed.pdb"
       }
   }

In this case we'd be computing the RMSD between the Complex and the
Native in the PDB stored in the indicated path. If you want to compute
the RMSD to the initial position you just have to make the "path"
parameter equal to the path of the initial PDB.

.. _sec-metricTypes-RMSD-example2:

Example 2
"""""""""

The RMSD computation can also be limited to a given selection of the
Complex, as we show in the following example:

.. code-block:: json

   {
       "type": "rmsd",
       "Native": {
           "path": "src/Molecules/Test/Data/AIN/native_fixed.pdb"
       },
       "selection": { 
           "chains": {
               "names": [ "L" ] 
           }
       }
       "superpositionSelection": { 
           "atoms": {
               "names": [ "_CA_" ] 
           }
       }
   }

In this case we'd be just computing the RMSD between the selected atoms
of the Complex and Native, (the atoms in the L chain), first doing a
superposition of the CA atoms of the system (which are in the protein
and not in the ligand).

.. _sec-metricTypes-RMSD-example3:

Example 3
"""""""""

By default the RMSD is computed after superposing the selected atoms of
the Complex and the Native. In the next example this superposition is
disabled by setting the "doSuperposition" parameter to false.

.. code-block:: json

   {
       "type": "rmsd",
       "savingFrequency": 1,
       "Native": {
           "path": "src/Molecules/Test/Data/AIN/native_fixed.pdb"
       },
       "selection": {
           "chains": { 
               "names": ["L"]
           }
       },
       "doSuperposition": false
   }

