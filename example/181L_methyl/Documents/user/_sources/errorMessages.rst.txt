.. _sec-errorMessages:

*******************
PELE Error Messages
*******************

This section explains the error messages that may appear during the execution of PELE, and provides hints as to how to fix the errors.

General errors in the control file
==================================

Error in the validation
-----------------------

Error message:

.. code-block:: console

  validation: FAIL

Cause:

There is a problem with the syntax of the control file. You may probably have received furthur information about the specific error in the validation, so please check more detailed error messages in this same section.

Suggested fix:

Check the more specific error. Review the control file syntax.

Missing configuration file
--------------------------

Error message:

.. code-block:: console

  /!\ Configuration file doesn't exist!! /!\

Cause:

  The name of the configuration file (also known as control file) that you provided does not match to an existing file in the file system.

Suggested fix:

Provide the right path to your configuration file for the corresponding command line argument of PELE.

Missing PeleTasks section, or missing comma before the PeleTasks keyword in the control file
--------------------------------------------------------------------------------------------

Error message:

.. code-block:: console

    PeleTasks failed validation
	    required item 'PeleTasks' not found.
    validation: FAIL
    Configuration file not valid. Macro won't start

Cause:

Either there is no PeleTasks section in the control file, or there is a missing comma separating the PeleTasks section from the previous section (or even in some, unrelated, part of the control file), so it has not been recognized.

Suggested fix:

Add the PeleTasks section if it is missing, or make sure there is a comma separating this section from the previous one (or that all commas are rightly placed in the complete control file). If the problem is not with the PeleTasks section, try to place this section at the beginning of the command, since in those cases the error message is usually more meaningful.

Missing comma (``,``)  or curly brace (``}``) in the control file
-----------------------------------------------------------------

Error message:

.. code-block:: console

    terminate called after throwing an instance of 'JsonParserException'
    PRINTS THE CONF FILE
    Errors: * Line 141, Column 6
      Missing ',' or '}' in object declaration

Cause:

There is a syntax error in the control file. It generally detects the error and the place where the missing curly brace or, less frequently, the comma, should be placed.
    
Suggested fix:
    
First, review if there is a missing comma or curly brace at the shown position. If not, review your control file, since you have problably missed a comma, curly brace, or other separator, before that position.

Missing and opening square bracket (``[``) in the commands section
------------------------------------------------------------------
    
Error message:
    
.. code-block:: console

    commands failed validation
	    Complex failed validation
		    files failed validation
			    files is of incorrect type
	    required item 'commands' not found.
    validation: FAIL
    Configuration file not valid. Macro won't start    

Cause:
    
The Complex section is missing the opening square bracket in the files argument.
    
Suggested fix:
    
Make the Complex files argument have a value enclosed by square brackets.
    
Property found in the wrong place, or mistyped property name
------------------------------------------------------------
    
Error message:
    
.. code-block:: console
    
    (root): additional property 'trajectoryPath' found.
    validation: FAIL
    Configuration file not valid. Macro won't start


.. code-block:: console

    ANM failed validation
	    anmMinimizer failed validation
		    anmMinimizer: additional property 'MaximumMinimizationIterations' found.
    validation: FAIL
    Configuration file not valid. Macro won't start
    

	
Cause:

The given configuration keyword (trajectoryPath in the example) could not be recognized because it is not at the right level (it has to be included in the PELE_Output section, but it was placed in a different level). The ``(root)`` part mentions that the property was found at the top level of the command checking. In the second example, the property was found on the anmMinimizer section of the ANM section.

Suggested fix:

Make sure the configuration keyword is at the right level in the control file.

Missing value after keyword
---------------------------

Error message:

.. code-block:: console

    PeleTasks failed validation
	    anmFrequency failed validation
		    anmFrequency is of incorrect type
	    required item 'PeleTasks' not found.
    validation: FAIL
    Configuration file not valid. Macro won't start
    
Cause:

There is no value for a given configuration property, so it is understood as of being of the wrong type. It may also happen that the value is of incorrect type (for example, a string where it expects a number).
    
Suggested fix:

Add the missing value with the correct type for the given configuration property.

Missing command type value
--------------------------

Error message:
    
.. code-block:: console
    
  terminate called after throwing an instance of 'controlFileValidation::ValidationException'
    what():  Error when validating control file: basic_string::_S_construct NULL not valid

Cause:

Each command must have a "commandType" property (for example: "commandType": "peleSimulation"), if not, this error appears.

Suggested fix:

Add the commandType property to all your commands.

Template and parameter related errors
=====================================

Missing atom in the PDB file, or missing TER record
---------------------------------------------------

Error message:

.. code-block:: console

    terminate called after throwing an instance of 'PeleBuildException'
      what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in LinkBuilder::createNewLink: Trying template LEUE:
    Template atom 'HD23' in template 'LEUE' not found in residue 'LEU:296:A'
    !!!!
    
Cause:
    
Your PDB file is missing an atom that is specified in the template file for that residue; or you missed a TER record separating molecules, and your N- or C-terminal residue is understood as a normal internal residue.
    
Suggested fix:
    
Make sure all the atoms in your PDB file for that residue match the names in the template: you may have to rename some of them. Take care with spaces, since they are significant. Also, for protein terminal residues, you have to check special templates: those ending in B (for example, "LEUB") for N-terminal types, and those ending in E (for example, "LEUE") for C-terminal types. You may also want to check you placed a TER record separating your molecules: if this residue is a C-terminal residue, but you forgot to separate it from the next molecule with a TER (or your residue is an N-terminal, and is not separated from the preceding molecule with a TER), PELE will understand the residue as a non-terminal residue, and it will probably complain about non-matching atoms.
    
Mismatching atoms between template and residue, or extra atoms, or missing TER
------------------------------------------------------------------------------

Error message:

.. code-block:: console

  terminate called after throwing an instance of 'PeleBuildException'
    what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in LinkBuilder::createNewLink: Trying template LIGZ: 
  Atom name not found in template 'LIGZ'; atom is ' C16' in residue 'LIG:1:Z

Cause:

Your PDB file has atoms that are not found in the template file for the atom's residue; or you missed a TER record separating molecules, and your N- or C-terminal residue is understood as a normal internal residue.

Suggested fix:

Make sure all the atoms in your PDB file have the right names (as included in the residue template). Also, make sure you are using as residue name the name corresponding to the right template. Take care with spaces, since they are significant. Also, for protein terminal residues, you have to check special templates: those ending in B (for example, "LEUB") for N-terminal types, and those ending in E (for example, "LEUE") for C-terminal types. You may also want to check you placed a TER record separating your molecules: if this residue is a C-terminal residue, but you forgot to separate it from the next molecule with a TER (or your residue is an N-terminal, and is not separated from the preceding molecule with a TER), PELE will understand the residue as a non-terminal residue, and it will probably complain about non-matching atoms.

Missing template
----------------

Error message:

.. code-block:: console

    terminate called after throwing an instance of 'PeleBuildException'
      what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in LinkBuilder::createNewLink: Trying template FE Z:
    Error getting structural template: /!\ Error in ForceField::getStructuralTemplate: Warning BAD template name!! -> *FE Z*!!
    !!!!

Cause:

One of your residues (it may be the ligand or a biomolecular residue) does not have a corresponding template.
    
Suggested fix:
    
Either make sure your residue is named exactly as your template (in the case of ligands, the template has an extra `z` character; for biomolecules, if the residue is an N-terminal residue, the template ends in `b`; if it is a C-terminal residue, the template ends in 'e'). Make sure your template is in the template directory (:file:`Templates/XXX/YYY/` where `XXX` is the solvent model name, and YYY is the type of biomolecule (DNA, Protein or RNA) or, for ligands, ions and cofactors, HeteroAtoms).
    
Missing OBC parameters
----------------------
    
Error message:
    
.. code-block:: console
    
  terminate called after throwing an instance of 'PeleBuildException'
    what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in ImplementationObcAlphaSasaUpdater::checkAtomSolventId: Error atom parameters NOT found in template for AINZ_O1!!!!

Cause:

The OBC parameters for a residue (usually a ligand, ion or cofactor) are missing. 

Suggested fix:

For a ligand, ion or cofactor, you should generate the parameters with :program:`solventOBCParamsGenerator.py` (see :ref:`sec-molecularParameters-obc`).

For a residue in a biomolecule, PELE should provide this data, so please contact PELE support for getting this issue fixed.
    
    
Selection errors
================

Wrong range specification
-------------------------

Error message:

.. code-block:: console

     what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in SelectionBuilder::getRange: Expecting 2 elements in range specification, but found 1. Offending link range string: 'A:35:A:37'.!!!!

Cause:

When specifying a range, you used the wrong format.

Suggested fix:

Use the range specification format: two selection elements, separated by a space. For example, ``"A:1 A:4"``. See :ref:`sec-selectionExamples-example11`.


Miscellanea
===========

Couldn't create the given file
------------------------------
    
Error message:
    
.. code-block:: console
    
    terminate called after throwing an instance of 'PeleBuildException'
      what():  /!\ Error in MacroBuilder::createMacro: /!\ Error in PeleOutPutBuilder::createPeleReport: Error: Could't create PeleReport file! (/gpfs/scratch/bsc72/bsc72156/WORK/VS2015/Dataset_OK/TESTset/PL_opls/22_3F80/output/22_3F80_test_report.txt)!!!!
    Aborted
    
Cause:
    
PELE couldn't create the given file, either because the parent directory does not exist, or because it has no write permissions in the directory, or because the disk is full.
    
Suggested fix:
    
Make sure the parent directory exists (or create it if not). Also, make sure you have the write permission for that directory, and the execute permission for that directory and all parent directories in the path. Finally, make sure there is enogh disk space.
    
