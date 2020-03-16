Protocols
===========


High Throughput
---------------

To run in HT mode use the flag **-HT** that will perform: 3 growing steps of 3 pelesteps
and a final exploration of 20 pele steps

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf -HT




Standard Precision
---------------------

To run in SP mode use the default values that will perform: 6 growing steps of 6 pelesteps
and a final exploration of 20 pele steps

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf


Explorative
-----------

To run in EX mode use the flag **-EX** that will perform an standar growing simulation with a sampling simulation of: 2 steering, 0.5/0.3 translation, 0.4/0.15 rotation and 25 of box radius. Select the number of steps in this sampling simulation with the flag **-es**.

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf -EX -es 200


Personalized Sampling Simulation
--------------------------------

If you want to use a different templatized control file in the sampling simulation, use the flag **-sc** to set its path

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf -sc /path/to/control_personalized.conf

Only prepare
------------

If for any reason you would need to prepare all PDB files, templates and configuration files but without growing you can
use this mode setting the flag -op. FragPELE will generate one folder for each fragment-core connection specified in the
serie_file.conf.

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf -op

Only grow
---------
After using the previous mode, to only execute the growing part of the code set the flag -og. Take into account that you
will use the configuration specified in the preparation step.

::

    python frag_pele/main.py -cp core.pdb -sef serie_file.conf -og
