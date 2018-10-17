.. _sec-permanentConstraints:

Permanent constraints in Pele++ control file
============================================

Constraints are applied in all minimizations during a Pele simulation. A
constraints list is specified within an array (a "[ ]" delimited block)
called "constraints", where each element of the array is a constraint's
specific configuration. e.g:

.. code-block:: json

   "constraints" : [constraint_configuration_1, constraint_configuration_2, ..., constraint_configuration_n] 

Constraints appear in the top level section of a
:ref:`peleSimulation <sec-PELESimulation>` command, but they can
also appear when dealing with a specific
:ref:`peleTask <sec-peleTasks>`, as an addition to the top level
constraints. The following snippet shows how to include them:

.. code-block:: json

   ...
       "commands": [
           {
               "commandType": "peleSimulation",
               "selectionToPerturb": {
                   "chains": {
                       "names": [
                           "L"
                       ]
                   }
               },
               "constraints":[
                   {
                       "type": "constrainAtomsDistance",
                       "springConstant": 200,
                       "equilibriumDistance": 2.55,
                       "constrainThisAtom":  "B:1:CA__",
                       "toThisOtherAtom": "A:48:_OD1"
                   },
                   {
                       "type": "constrainAtomsDistance",
                       "springConstant": 200,
                       "equilibriumDistance": 2.55,
                       "constrainThisAtom":  "B:1:CA__",
                       "toThisOtherAtom": "A:48:_OD2"
                   }
               ],
               ...
           }
       ]
   ...    

Below are examples of specific contraint's configuration

Constraint types
----------------

Harmonic constraint on atom
^^^^^^^^^^^^^^^^^^^^^^^^^^^

It's a harmonic constraint on an Atom using a target position

Example:

.. code-block:: json

   {
       "type": "constrainAtomToPosition",
       "springConstant": 2.2,
       "equilibriumDistance": 0.0,
       "constrainThisAtom":  "A:1:_CB_",
       "toThisPoint": [3.3, 2.2, 1.1]
   }   

Notice that if the "toThisPoint" key is not specified, by default PELE
will constrain the atom to its original position

Harmonic constraint on distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It's a harmonic constraint on the distance between two Atoms

Example:

.. code-block:: json

   {
       "type": "constrainAtomsDistance",
       "springConstant": 2.2,
       "equilibriumDistance": 3.3,
       "constrainThisAtom":  "A:1:_CB_",
       "toThisOtherAtom": "A:2:_CA_"
   }

Harmonic constraint on distance to center of masses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It's a harmonic constraint on the distance between an Atom and the
center of masses of a group of Atoms

Example:

.. code-block:: json

   {
       "type": "constrainDistanceToCOM",
       "springConstant": 2.2,
       "equilibriumDistance": 3.3,
       "constrainThisAtom":  "A:1:_CB_",
       "toCenterOfMassesOf": {
           "chains": {
               "names": [ "B" ]
           }
       }
   }

Notes: "toCenterOfMassesOf" can be any valid selection block, see
:ref:`Selections By Example <sec-selectionExamples>`.

