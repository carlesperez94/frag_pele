===============
Getting Started
===============

.. toctree::
   :maxdepth: 2

Here you will learn how to run your first FrAG-PELE simulation.

How to do it?
-------------

**Previous Requisites**
++++++++++++++++++++++++

- **Complex PDB**: PDB processed, with the protein and the ligand prepared to run with PELE.
  Remember to rename the ligand chain to "L".
- **Fragment PDB**: PDB with the desired fragment. Remember to rename the chain to "L".
- **Serie file**: The instructions that explain how the growing has to be done are stored in this file. It must have the
  format described bellow.

**Growing of a fragment to a core** 
+++++++++++++++++++++++++++++++++++++

Grow 1 fragment to 1 cores

+---------------+-------------------+-----------------------+
| fragment PDB 1| heavy atom core 1 | heavy atom fragment 1 |
+---------------+-------------------+-----------------------+

ie::

    amino.pdb C1 N1

**Growing of two fragment to the same core**  
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Add 1 fragment to an initial core. Then, add another fragment to the
same core on the output of the 1st growing. Example: Grow two radicals onto a phenyl

+---------------+------------------------+-----------------------+---------------+------------------------+-----------------------+
| fragment PDB 1| heavy atom core 1 (C1) | heavy atom fragment 1 | fragment PDB 2| heavy atom core 1 (C2) | heavy atom fragment 2 |
+---------------+------------------------+-----------------------+---------------+------------------------+-----------------------+

i.e::

    amino.pdb C1 N1 phenyl C1 C1


**Growing of a fragment to a previously grown fragment**  
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Add 1 fragment to an initial core. Then, add another fragment to the
previously grown fragment. Example: Grow one radical onto a phenyl and
then another fragment into the previously grown radical.

+-----------+-------------------+-----------------------+---------------+--------------------------------------------+-----------------------+
| fragment 1| heavy atom core 1 | heavy atom fragment 1 | fragment PDB 2| fragment atom core 1 (C2)*fragment_number* | heavy atom fragment 2 |
+-----------+-------------------+-----------------------+---------------+--------------------------------------------+-----------------------+

i.e::

    amino.pdb   C1  N1  phenyl.pdb  N1*1*    C1


**All together** 
++++++++++++++++++

Add 2 fragments onto different position of the same core. Then, perform another growing
with a 3rd fragment onto a different core and a 4th into the previosuly grown fragment.


+---------------+------------------------+-----------------------+---------------+------------------------+-----------------------+
| fragment PDB 1| heavy atom core 1 (C1) | heavy atom fragment 1 | fragment PDB 2| heavy atom core 1 (C2) | heavy atom fragment 2 |
+---------------+------------------------+-----------------------+---------------+------------------------+-----------------------+
| fragment PDB 3| heavy atom core 2      | heavy atom fragment 3 | fragment PDB 4| heavy atom fragment 3  | heavy atom fragmnet 4 |
+---------------+------------------------+-----------------------+---------------+------------------------+-----------------------+

i.e::

    amino.pdb C1 N1 phenyl C3 C1
    ethyl C1  C2 methyl C1*1* C1


**Launch simple FrAG-PELE simulation**
+++++++++++++++++++++++++++++++++++++++

Default: 10 growing steps of 6 pele steps, and a final equilibration simulation of 20 pele steps

::

    python3.X grow_for_pele.py -cp path_to_complex_pdb -sef serie_file



