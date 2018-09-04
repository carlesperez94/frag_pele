.. _sec-complex:

Complex Block in Pele++ control file
====================================

It sets the complex to be used throughout the simulation.

Plop info: "load pdb" command.

Location
--------

The Complex block must be defined inside the initialization block as
shown in the following example:

.. raw:: html

   <pre>
   {
       ...
       
       "Initialization":
       {
           ...
          
           "Complex": 
           {
             "files": [{"path": "complex.pdb"}],
           }
       },
       ...
   }
   </pre>

Options and parameters
----------------------

files
^^^^^

An array of paths, since the PDB can be split in different files. It
will define the complex to be used throughout the simulation.

