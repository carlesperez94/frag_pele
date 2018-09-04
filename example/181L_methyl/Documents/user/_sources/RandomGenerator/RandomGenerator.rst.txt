.. _sec-randomGenerator:

Random Generator in Pele++ control file
=======================================

Random numbers are generated using a Mersenne Twister pseudo random number generator. The generator can be initialized with a user-provided seed, or it may start at the state stored in a file, or, if none of the above options are given, it uses a seed created from the current time.

If working with several trajectories, then you can only provide a seed or use the default seed based on the current time. Notice that each MPI process will use a different seed: rank 1 process gets the original seed (either provided in the configuration file, or obtained by using the current time), all other processes increment the seed according to their rank; for example, rank 3 process gets the original seed + 2. The MPI controller process does not require any random number generator.

Example using a seed
--------------------

Seed range: [0, inf)

.. code-block:: json

   "RandomGenerator": {
       "seed": 30786
       }

Example using a previous random generator stored in a file
----------------------------------------------------------

.. code-block:: json

   "RandomGenerator": {
       "randomGeneratorFile": "../example.txt"
       }

