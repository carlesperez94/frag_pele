.. _sec-selectionExamples:

*********************
Selections By Example
*********************

A selection allows selecting specific atoms, complete residues (links) and complete chains. It is possible to create a selection that is a union of atomic, link and chain selections (however, it is not currently possible to create a selection by intersection). It is also possible to remove some parts of a selection through the "omit" clause.

A selection that includes all possibilities is of the form:

.. code-block:: json

    "selection" : {
      "atoms" : {
	 "ids": [...],
	 "names": [...],
	 "elements": [...]
      },
      "links" : {
	 "ids": [...],
	 "names": [...],
	 "ranges": [...]
      },
      "chains" : {
	 "names": [...]
      },
      "omit": {
	  "atoms" : {
	     "ids": [...],
	     "names": [...],
	     "elements": [...]
	  },
	  "links" : {
	     "ids": [...],
	     "names": [...],
	     "ranges": [...]
	  },
	  "chains" : {
	     "names": [...]
	  }
      }
    }

There are also special selections like the `"empty"` string (see :ref:`sec-selectionExamples-example15`), or the special keywords used in chain and link selection (see :ref:`sec-selectionExamples-example7`).

The following examples show how the specific sections work. Notice that the atoms can be selected, individually, by id, name or chemical element. Links can be selected by id, name or through a range. Finally, chains can be selected by name.
    
The omit section works only at the atom, link or chain level. Therefore, all atoms selected under the omit "atoms" section will be removed from those atoms in the selection "atoms" section, but it won't affect the atoms selected through the selection "links" or "chains" sections. Similarly happens for the omit "links" and "chains" sections. This behavior can be understood by studying the examples :ref:`sec-selectionExamples-example13`, :ref:`sec-selectionExamples-example14` and :ref:`sec-selectionExamples-example16`).

Example 1
=========

.. code-block:: json

   "selection": {
       "atoms": {
           "names": [ "_N__" ]
       }
   }

Explanation: This is selecting all the atoms in PDB the names of which
are " N ". The underscores must be there because atom names in PDBs have
4 characters.

Example 2
=========

.. code-block:: json

   "selection": {
       "atoms": {
           "names": [ "_N__", "_C__" ]
       }
   }

Explanation: This is selecting all the atoms in PDB the names of which
are " N " or " C ".

Example 3
=========

.. code-block:: json

   "selection": {
       "atoms": {
           "names": [ "_N__", "_CA_", "_C__" ]
       }
   }

Explanation: This is selecting all the atoms in PDB the names of which
are " N ", " CA " or " C ".

Example 4
=========

.. code-block:: json

   "selection": {
       "atoms": {
           "ids": [ "A:8:_CB_", "B:32:_H2'" ]
       }
   }

Explanation: This is selecting two atoms in the PDB indicated by its
unique identification (id) -> "chain\_name:residue\_number:atom\_name".

Note that the different components of the identification are separated
by colons, ":".

In this case the atoms are:

1. An Atom named " CB " in the 8th residue of the "A" chain.
2. An Atom named " H2'" in the 32nd residue of the "B" chain.

Example 5
=========

.. code-block:: json

   "selection": {
       "chains": {
           "names": [ "A" ]
       }
   }

Explanation: This is selecting a chain indicating it in the PDB
indicated its unique identification (id) -> "chain\_name".

In this case its the "A" chain.

Example 6
=========

.. code-block:: json

   "selection": {
       "chains": {
           "names": [ "B" ]
       }
   }

Explanation: This is also selecting a chain. In this case it is the "B"
chain.

.. _sec-selectionExamples-example7:

Example 7
=========

.. code-block:: json

   "selection": {
       "chains": "all"
   }

Explanation: This is selecting all the chains of the PDB.

The complete set of special strings for chain selection is:

- "all": Selects all chains.
- "everythingButLigands": Selects all chains that aren't of ligand type.
- "ligands": Selects all chains of ligand type.

Notice that a chain is considered a ligand if it has an hetero atom.

For links, you can also select all links by using:

.. code-block:: json

   "selection": {
       "links": "all"
   }

Currently, only "all" is allowed for links.

Example 8
=========

.. code-block:: json

   "selection": {
       "links": {
           "names": [ "ALA_" ]
       }
   }

Explanation: This is selecting all the links (residues in this case) in
the PDB which have "ALA " for a name.

Again the underscores must be there because link names in PDBs have 4
characters.

.. _sec-selectionExamples-example9:

Example 9
=========

.. code-block:: json

   "selection": {
       "links": {
           "names": [ "ASN_", "MET_" ]
       }
   }

Explanation: This is selecting all the links (residues in this case) in
the PDB which have "ALA " or "MET " for a name.

Example 10
==========

.. code-block:: json

   "selection": {
       "links": {
           "ids": [ "_:5", "_:26", "_:31" ]
       }
   }

Explanation: This is selecting three links in the PDB indicated by its
unique identification (id) -> "chain\_name:residue\_number".

Note that the different components of the identification are separated
by colons, ":".

Note that whenever the chain does not have a name, you must indicate it
by using an underscore.

.. _sec-selectionExamples-example11:

Example 11
==========

.. code-block:: json

   "selection": {
       "links": {
           "ranges": [ "A:1 A:4", "B:3 B:5" ]
       }
   }

Explanation: This is selecting two ranges of links in the PDB.

Link selection ranges are indicated in the following way:
initial_link_id last_link_id.

Note that the initial and last link ids are separated by a space, " ".

The only limitation is that the initial and last links must be in the
same chain.

This selection is selecting the links in the "A" chain from the 1st to
the 4th and the links in the "B" chain from the 3rd to the 5th.

Example 12
==========

.. code-block:: json

   "selection": {
       "links": {
           "ranges": [ "A:1 A:4", "B:19 B:23" ]
       }
   }

Explanation: This is selecting the links in the "A" chain from the 1st
to the 4th and the links in the "B" chain from the 19th to the 23th.

.. _sec-selectionExamples-example13:

Example 13
==========

.. code-block:: json

   "selection": {
       "atoms" : { 
           "ids":["A:8:_CB_", "B:32:_H2'"],
           "names":["_CA_", "_H__"]
       },
       "links" : {
           "names":["ALA_","GLY_"], 
           "ids":["_:5", "_:26", "_:31"],
           "ranges":["A:1 A:10", "B:3 B:8"]
       },
       "chains" : { 
           "names":["A","B","C"]
       },
       "omit":{
           "atoms" : { 
               "ids":["B:32:_H2'"],
               "names":["_H__"]
           },
           "links" : {
               "names":["ALA_","GLY_"], 
               "ids":["_:26", "_:31"],
               "ranges":["A:3 A:7", "B:3 B:8"]
           },
           "chains" : { 
               "names":["B","C"]
           }"
       } 
   }

Explanation: Complete example of use of omit.

The 'omit' clause is a selection of the things that will be eliminated
from the current selection.

This means that something HAS TO BE EXPLICITLY DECLARED in order to omit
it or a part of it ( see example below ).

It automatically handles residue ranges and some light errors
(repetitions, overlappings, etc...), but it does not work (YET) when
using the special words which describe "allChains", "allLigands",
"everythingButLigand" or "allLinks".

To eliminate generalizations from selection like atom name or link name,
another generalization has to be written.

.. _sec-selectionExamples-example14:

Example 14
==========

Our purpose is to select CA atoms of the first 100 residues of a protein
with 120 residues using omit (a version without is also possible, of
course). This cannot be currently represented with the selection language.

Imagine, however, that it could be done. In such a case, one could write:

.. code-block:: json

   "selection": {
       "atoms" : { 
           "names":["_CA_"]
       },
       "omit":{
           "links" : {
               "ranges":["A:101 A:120"]
           }
       }
   }

But it won't work because the range containing all the links wasn't
defined. We should write instead:

.. code-block:: json

   "selection": {
       "atoms" : { 
           "names":["_CA_"]
       },
       "links" : {
           "ranges":["A:1 A:120"]
       },
       "omit":{
           "links" : {
               "ranges":["A:101 A:120"]
           }
       }
   }

Since the selection in PELE is additive, the previous example does not obtain the alpha carbons of residues 1-100 in chain A but, instead, it obtains all alpha carbons, as well as all atoms in residues 1-100 of chain A.

.. _sec-selectionExamples-example15:

Example 15
==========

The empty selector can be expressed with the "empty" keyword:

.. code-block:: json

   "selection": "empty"

.. _sec-selectionExamples-example16:

Example 16
==========

The omit clause works at the level of "atoms", "links" and "chains", and it doesn't matter how those atoms, links or chains were selected (either by id, name, element, or range). However, remember that an omit for links only affects selection of "links", and similarly for "atoms" and "chains".

The following example, where residue 2 is of type GLY, would end up selecting residues 1,3 and 4:

.. code-block:: json

  "selection": {
    "links": "ranges":["A:1 A:4"],
    "omit":{
       "links" : {"names":["GLY"]}
    }
  }


