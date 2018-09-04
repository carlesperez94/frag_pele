.. _sec-generalStructure:

*********************************************
General Structure of a control file in Pele++
*********************************************

It contains 3 blocks.

The top-level configuration block
=================================

This contains a few global configuration variables.

licenseDirectoryPath
--------------------

Path of the directory where the license file (and its digest file) are
located. Default value is "/home/user". It expects to find a
"peleLicense.txt" licese file and its corresponding digest
"peleLicense.txt.digest".

.. code-block:: json

   "licenseDirectoryPath" : "/home/user/mylicenses"

Default value: Either the command-line provided directory path or, if not provided, "/home/user".


controlFileSavingPath
---------------------

Saves the current control file in the given path.

.. code-block:: json

  "controlFileSavingPath" : "../../simulations/ain/originalControlFile.conf"

Default value: "" (no control file is saved)


simulationLogPath
-----------------

Use: Path of the log file where all the commands' important information
is logged.

Parameter: std::string simulationLogPath

Note for developers: It's defined in SystemVars class.

Default value: "" (no logging).

.. _sec-generalStructure-commandFilePath:

commandFilePath
---------------

Use: Path of the file where commands are programmatically sent to the application for Real Time Control. This should not be configured by the user. It is an error if the file exists before running the application.

Default value: No default value. If commandFilePath and commandOutputFilePath are not present, Real Time Control is not activated. It is an error if only one of them is configured.


.. _sec-generalStructure-commandOutputFilePath:

commandOutputFilePath
---------------------

Use: Path of the file where the program will store the results of Real Time Control commands. This variables should be programmatically used by the GUI, and not configured by the user. It is an error if the file exists before running the application.

Default value: No default value. If commandFilePath and commandOutputFilePath are not present, Real Time Control is not activated. It is an error if only one of them is configured.

Initialization
==============

Contains the setup of the system (typically pdb file(s), force field and
solvent). To read more, check the :ref:`initialization
document <sec-initialization>`.

commands
========

Array containing the set of commands to perform throughout the
simulation. To read more, check the :ref:`commands
document <sec-commands>`.
