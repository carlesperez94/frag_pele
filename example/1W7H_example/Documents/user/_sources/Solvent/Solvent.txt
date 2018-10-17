.. _sec-solvent:

Solvent Block in Pele++ control file
====================================

Location
--------

The Solvent block must be defined inside the initialization block as
shown in the following example:

.. code-block:: json

   {
       ...
       
       "Initialization":
       {
           ...
          
           "Solvent": 
           {
               "solventType": "SGBNP",
               "useDebyeLength": true,
               "ionicStrength": 0.15 
           }
       },
       ...
   }

Options and parameters
----------------------

.. _sec-solvent-solventType:

solventType
^^^^^^^^^^^

Uses: It sets the solvent type to be used throughout the whole
simulation.

Possible options: 

- "VACUUM": Vacuum.
- "VDGBNP": Variable-Dielectric Generalized Born + Non-Polar model. 
- "SGBNP": Surface Generalized Born + Non-Polar.
- "OBC": OBC + ACE Non-Polar model.

Default value: "SGBNP", unless the "Solvent" block is not defined, in which case the "VACUUM" solvent is used.

useDebyeLength
^^^^^^^^^^^^^^

Uses: Permits the use of ionic strength. This parameter only makes sense for non-vacuum solvent types. For vacuum, the value of this parameter is ignored.

Default value: true

Plop info: No direct correspondence, if 'energy\_params%kappa' > 0.001
it is set to True.

ionicStrength
^^^^^^^^^^^^^

Use: If useDebyeLength is true, the ionic strenght to be used

Default value: 0.15

Plop info:

-  Plop control file name: 'ionic'
-  Plop parameter name: 'energy\_params%ionic\_strength'
-  Plop default value: 0

