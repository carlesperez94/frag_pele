.. _sec-peleOutput:

PELE output parameters in Pele++ control file
=============================================

The "PELE_Output" block defines the output options and the paths for
all the output files.

The "PELE_Output" block must be defined inside the PeleSimulation
command as this parameters must remain the same during the entire Pele
simulation.

Output options and parameters.
------------------------------

savingMode
^^^^^^^^^^

Use: It tells Pele++ how the accepted models must be saved to disk.

Parameter: ConfigurationSavingMode savingMode

Note for developers: It's defined in PeleOutputParameters class and its
options in the ConfigurationSavingMode enum.

Options:

-  "savingTrajectory": It saves the accepted models in a given
   trajectory file. Model 1 always corresponds to the initial state.
-  "savingSteps": It saves the accepted models in different files in a
   given folder. The file corresponding to accepted step 0 is the
   initial state.
-  "savingBoth": It does both previous things.

Default value: "savingTrajectory"

trajectoryPath
^^^^^^^^^^^^^^

Use: Path of the file where all the models are saved to form a
trajectory. Model 1 always corresponds to the initial state. Model 2
corresponds to the first recorded accepted step (see
saveFrequencyForAcceptedSteps), Model 3 to the second recorded accepted
step, and so on. It's only used in "savingTrajectory" and "savingBoth"
options. If the route is invalid, the program cannot write the
trajectory, but the execution of the program proceeds.

Parameter: std::string trajectoryFilePath

Note for developers: It's defined in PeleOutputParameters class.

Default value: "../../trajectory.pdb"

stepsFolder
^^^^^^^^^^^

Use: Path of the folder where all the accepted models are saved in
separated files. The file numbered as 0 is the initial state. For all
other files, the number is the accepted step (thus, n corresponds to the
nth accepted step; if savingFrequencyForAcceptedSteps is 2, the values
for n will be 2, 4, 6, etc.). It's only used in "savingSteps" and
"savingBoth" options.

Parameter: std::string stepsFolder

Note for developers: It's defined in PeleOutputParameters class.

Default value: ""

savingFrequencyForAcceptedSteps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use: Accepted steps are saved only if they are multiples of this
parameter. Thus, besides the initial state, which is always saved, if
savingFrequenceyForAcceptedSteps is n, only accepted steps n, 2\*n, etc.
are saved.

Parameter: unsigned int savingFrequency

Note for developers: It's defined in PeleOutputParameters class.

Range: (1, inf)

Default value: 1, which means that every accepted step is saved.

reportPath
^^^^^^^^^^

Use: Path of the file where the :ref:`Sensors <sec-sensors>`
(metrics and tracked variables) values are stored. The first data line
corresponds to the initial state (marked with a value of 0 for both Step
and AcceptedSteps). Then, there is a data line per accepted step (both
Step and AcceptedSteps start counting at 1). If the route is invalid,
the execution of the program is aborted.

Parameter: std::string fileReportPath

Note for developers: It's defined in PeleOutputParameters class.

Default value: "../../test\_report.txt".

initialPdbPath
^^^^^^^^^^^^^^

Use: Path of the file where the initial PDB is stored. If the route is
invalid, the program cannot write the initial structure, but the
execution of the program proceeds. If it is an empty string, no file is
created nor PDB is stored.

Parameter: std::string initialPDBPath

Note for developers: It's defined in PeleOutputParameters class.

Default value: "".

finalPdbPath
^^^^^^^^^^^^

Use: Path of the file where the final PDB is stored. If the route is
invalid, the program cannot write the final structure, but the program
finishes the execution normally. If it is an empty string, no file is
created nor PDB is stored.

Parameter: std::string finalPDBPath

Note for developers: It's defined in PeleOutputParameters class.

Default value: "".

flushingFrequency
^^^^^^^^^^^^^^^^^

Use: When the number of report data lines is a multiple of
"flushingFrequency" the report file buffer is flushed to the file. There
are as many report data lines as number of accepted steps plus one (the
initial state report line). The lower this parameter, the faster you
will see your results in the belonging log.

Parameter: unsigned int flushingFrequency

Note for developers: It's defined in PeleOutputParameters class.

Range: (1, inf)

Default value: 50.

.. _sec-peleOutput-controllerEventsLogPath:

controllerEventsLogPath
^^^^^^^^^^^^^^^^^^^^^^^

Use: Path of the log file where all ExplorationController events are
logged. Examples of this are: parameter changes, new best explorer
election, explorer jump (both forced and requested), explorer
finalization... It's only needed when a Controller for jump events is
defined in the control file.

Parameter: std::string controllerEventsLogPath

Note for developers: It's defined in PeleOutputParameters class.

Default value: "".

Example
-------

.. code-block:: json

   "PELE_Output":
   {
       "reportPath": "../../test_report",
       "trajectoryPath": "../../trajectory.pdb",
       "stepsFolder": "../../step_",
       "savingMode": "savingBoth",
       "savingFrequencyForAcceptedSteps": 1,
       "initialPdbPath": "../../theInitial.pdb",
       "finalPdbPath": "../../theFinal.pdb",
       "simulationLogPath": "../../simulationLog.txt",
       "controllerEventsLogPath": "../../controllerEventsLog.txt",
       "flushingFrequency": 5
   }           

Simulations using MPI
---------------------

In simulations using MPI the same parameters are used. Pele++ adds a
suffix to all file and folder names to separate the output of each
Explorer. The controller process (rank 0) does not output any of these
files. Explorers start being numbered at 1 (following their MPI rank).

