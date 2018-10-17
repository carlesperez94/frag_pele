.. _sec-boxes:

Boxes
-----

In the control file we can define Boxes in which the ligand will move.

These Boxes are added in the LigandPerturbation block, see
:ref:`Perturbation <sec-perturbation>`.

This would be an example of a spherical box definition inside a
Perturbation block:

.. code-block:: json

   ...
   ...
   "Perturbation": {
      "perturbationType": "naive",
      "parameters": {
          "rotationScalingFactor": 0.05,
          "translationRange": 2.5,
          "overlapFactor": 0.7
      },
      "Box": {
          "type": "sphericalBox",
          "radius": 45,
          "selectionOfAtomsForCOM": {
              "chains": "everythingButLigands"
          }
      }
   }, 
   ...
   ...   

If you are using boxes and steering, and if a steric trail gets out of
the box, then the steering vector (if it was active) will be changed,
and the internal counter used to decide when to change of steering
vector (based on steeringUpdateFrequency), is reset.

Also, and independently of steering, if for whatever reason (for
example, during the ANM or relaxation minimization) the ligand gets
moved outside of the box, the PELE step will be rejected.

There are three type of boxes: 

Spherical Boxes
^^^^^^^^^^^^^^^

In this case the
exploration will take place inside a sphere of a given radius which is
centered either in a fixed point or in the center of masses of a given
selection of atoms.

If no selection or fixed center were indicated, the center of the Box
would be the center of masses of all the atoms.

You can find more examples of selections in :ref:`Selections By Example <sec-selectionExamples>`.

Example centered in the center of masses of all atoms:

.. code-block:: json

   "Box": {
      "type": "sphericalBox",
      "radius": 45
   }

Example centered in the center of masses of a selection of atoms:

.. code-block:: json

   "Box": {
      "type": "sphericalBox",
      "radius": 45,
      "selectionOfAtomsForCOM": {
          "chains": "everythingButLigands"
      }
   }

Example centered in a fixed point:

.. code-block:: json

   "Box": {
      "type": "sphericalBox",
      "radius": 45,
      "fixedCenter": [0.0, -5.0, 3.0]
   }

Prismatic Boxes
^^^^^^^^^^^^^^^

In this case the exploration will take place inside a prismatic box of
given dimensions which is centered either in a fixed point or in the
center of masses of a given selection of atoms.

Example centered in the center of masses of all atoms:

.. code-block:: json

   "Box": { 
      "type": "prismaticBox", 
      "x_y_z": [2, 4, 6] 
   }

Example centered in the center of masses of a selection of atoms:

.. code-block:: json

   "Box": {
      "type": "prismaticBox", 
      "x_y_z": [2, 4, 6],
      "selectionOfAtomsForCOM": {
          "chains": "everythingButLigands"
      }
   }

Example centered in a fixed point:

.. code-block:: json

   "Box": {
      "type": "prismaticBox", 
      "x_y_z": [2, 4, 6],
      "fixedCenter": [0.0, -5.0, 3.0]
   }

Spherical Shell Boxes
^^^^^^^^^^^^^^^^^^^^^

In this case the exploration will take place inside a shell between to
concentric spheres of given radii which is centered either in a fixed
point or in the center of masses of a given selection of atoms.

Example centered in the center of masses of all atoms:

.. code-block:: json

   "Box": { 
      "type": "sphericalShellBox", 
      "innerRadius": 2.0,
      "outerRadius": 5.0
   }

Example centered in the center of masses of a selection of atoms:

.. code-block:: json

   "Box": {
      "type": "sphericalShellBox", 
      "innerRadius": 2.0,
      "outerRadius": 5.0,
      "selectionOfAtomsForCOM": {
          "chains": "everythingButLigands"
      }
   }

Example centered in a fixed point:

.. code-block:: json

   "Box": {
      "type": "sphericalShellBox", 
      "innerRadius": 2.0,
      "outerRadius": 5.0,
      "fixedCenter": [0.0, -5.0, 3.0]
   }

