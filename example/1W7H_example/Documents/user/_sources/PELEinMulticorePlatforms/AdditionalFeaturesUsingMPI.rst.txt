.. _sec-additionalFeaturesUsingMPI:

Additional features in control file using MPI
=============================================

Loading multiple PDBs
---------------------

In a serial execution we only have the option to load one complex.

.. code-block:: json

   "Complex": {
               "files": [
                   {
                       "path": "complex.pdb"
                   }
               ]
           }

However when we use the MPI implementation, we can load a different PDB
to each explorer using the "MultipleComplex" option instead of
"Complex".

.. code-block:: json

   "MultipleComplex":
           [
            {
               "files":
               [
                   {
                       "path":"complex1.pdb"
                   }
               ]
           },
           {
               "files":
               [
                   {
                       "path":"complex2.pdb"
                   }
               ]
           }       
          ]

If we set less PDBs than explorers, there will be explorers using the
same PDBs. The PDBs are equitably distributed.

