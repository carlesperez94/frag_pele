Protocols
===========


HT
------

To run in HT mode use the flag -HT that will perform: 3 growing steps of 3 pelesteps
and a final exploration of 20 pele steps

::
 
    python frag/grow_for_pele.py -cp core.pdb -sef serie_file.conf -HT




Standard Precision
---------------------

To run in SP mode use the dafault values that will perform: 6 growing steps of 6 pelesteps
and a final exploration of 20 pele steps

::

    python frag/grow_for_pele.py -cp core.pdb -sef serie_file.conf

