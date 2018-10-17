.. _sec-initialization-general:

General description
===================

Location
--------

Outermost layer, typically at the beginning of the document.

Example:

.. code-block:: json

   {
      "Initialization" : {
	     "fixTrees": true,
   	     "fixChainNames": false,
         "allowMissingTerminals": false,
         "Complex" : {
            "files" : [
               {
                  "path": "pdbs/ain_fixed.pdb"
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
      "commands" : [ ... ]
   }

Options and parameters
----------------------

fixTrees
^^^^^^^^

Use: Bolean. When it is set to true, for any complex being built,
it builds the trees equal to Plop.

Note for developers: Used in generateAtomsTrees

Default value: true

fixChainNames
^^^^^^^^^^^^^

Use: Boolean. When true, for any complex being built,
it renames all chains starting at 'A' through 'Z' and then '1' through '0'.

Notice that if you have many molecules (such as water molecules) you may
run out of chain identifiers. In such a case, chain ids will get reused
(it will restart from 'A' again).

Note for developers: Used in
ComplexBuilder::createComplexBuilder.

Default value: false

.. _sec-initialization-general-allowMissingTerminals:

allowMissingTerminals
^^^^^^^^^^^^^^^^^^^^^

Use: Boolean. When true, when reading a complex, if a terminal residue in a chain does not have the atoms of the terminal type for that residue, then it is treated as a normal internal/middle residue. If this does not work either, then PELE fails as usual, since it is not able to assign the right template. Notice that this is valid for any chain, though only proteins and nucleic acids actually have terminal types of residues.

This option is useful if you have gaps in your structure. In such a case, you mark the gaps with TER records, and instruct PELE to interpret the terminal residues at the gaps (which must not have terminal atoms such as OXT) as normal internal residues.

Notice that, being internal residues, they will be treated as such in all parts of PELE. For example, in the VDGBNP case, those residues will be treated as internal as for the dielectric assignment.

Default value: false


Complex
^^^^^^^

See :ref:`Complex <sec-complex>`

MultipleComplex
^^^^^^^^^^^^^^^

See
:ref:`AdditionalFeaturesUsingMPI <sec-additionalFeaturesUsingMPI>`


.. _sec-initialization-general-solvent:

Solvent
^^^^^^^

See :ref:`Solvent <sec-solvent>`

.. _sec-initialization-general-forcefield:

ForceField
^^^^^^^^^^

Parameter: enum PotentialParameterization

Use: Sets the force field to be used throughout the simulation

Note for developers: Used in SgbAlphaSasaUpdater,
NonBondingTermCalculator, EnergyAlgorithm and
BindingEnergyMetricsBuilder.

Note: The force field parameter files are found under :file:`Data/Templates/`

Possible options:

-  "OPLS2005"

-  "AMBER99sb" (not thouroghly tested in production, but it should work)

-  "AMBER99" (not thoroughly tested in production, but it should work)

-  "AMBER99sbBSC0"

Defalut value: None, it must be specified

Plop info: No correspondence, it just uses OPLS2005.

useGPU
^^^^^^

Use: Boolean that enables the use of GPUs. We remind that its in an
experimental stage.

