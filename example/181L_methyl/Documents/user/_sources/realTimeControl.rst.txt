.. _sec-realTimeControl:

****************************************************
Real Time Control of PELE simulations (experimental)
****************************************************

This is an experimental feature to allow real time control of an already started PELE simulation, in a similar way to :ref:`sec-parametersThatCanChange`, but this time performing the change of parameters at the moment the user decides to issue the command. This feature is intended to be used indirectly by the user, with a middleware that will understand the user control commands and transform them into the low-level way of interacting with the PELE simulation.

This feature is currently only available for MPI simulations; in those cases, only PELE simulation commands can be executed.

The following commands are available:

- PAUSE: Pause the simulation.
- RESUME: Resume a paused simulation.
- TERMINATE: Terminate a simulation.
- CHANGE_PARAMETERS: Change simulation parameters (the same that can be changed in :ref:`sec-parametersThatCanChange`).
- NOOP: No-operation (command without effect).

To use real time control of a PELE simulation, you must enable this through the control file. Then, you must write commands to the command input file, while you can see how the commands are processed at the command output file.

Configure real time control of a PELE simulation
================================================

To enable real time control of a PELE simulation, you must add the ``commandFilePath`` and ``commandOutputFilePath`` configuration variables to the top level of the control file (see :ref:`sec-generalStructure-commandFilePath` and :ref:`sec-generalStructure-commandOutputFilePath`).

For example:

.. code-block:: json

    {
     "commandFilePath": "commands.txt", 
     "commandOutputFilePath": "commandsOutput.txt", 

     // ...
       "commands" : [
	  {
	     "commandType" : "peleSimulation",
	  // ...
	  }
      ]
    }

You must make sure that those files do not exist when you start the simulation, since PELE++ creates those files itself and, to make sure you write your commands in the right file, it will terminate execution if those files already existed.

Write commands to the command input file
========================================

Once the application has started and created the command files, you can write to the command stream (``commandFilePath``), where each command is a JSON object with three fields, and to let the program know that the command is finished, a final line ``*END`` is added to terminate the command. The fields of the JSON command object are:

- commandId: A numeric identifier to match command output in the ``commandOutputFilePath`` to a given input command.
- commandType: The type of command. One of the allowed ones (see above, in the introduction section).
- commandData: The actual contents of the command, as a string. It will contain a JSON object (codified as a string, so quote symbols must be escaped, and backslashes also must be escaped). It will contain a set of parameter change instructions, as in :ref:`sec-parametersThatCanChange`.

A PELE simulation, regarding real time control, can be in one of two states:

- Running the simulation, in which case it accepts both PAUSE and TERMINATE commands. Commands are only processed just before starting a new simulation step.
- In a paused state, where it can accept a CHANGE_PARAMETERS command, or it can receive a RESUME command to return to the running state, or a TERMINATE command to end with the simulation.

It always can receive a NOOP command, which has no effect.

For example, the following commands will pause the simulation, make a change of parameters, and then resume the simulation:

.. code-block:: text

    {"commandId": 1,
     "commandType": "PAUSE",
     "commandData": ""
    }
    *END
    {"commandId": 5,
     "commandType": "CHANGE_PARAMETERS",
     "commandData": "{
	    \"Perturbation::parameters\": {
		    \"numberOfSteps\" : 40,
		    \"numberOfTrials\" : 1,
		    \"rotationScalingFactor\": 0.4,
		    \"translationRange\": 3.5,
		    \"overlapFactor\": 1.3, 
		    \"temperature\": 81,
		    \"numberOfStericTrials\": 1
	    }
	   }"
    }
    *END
    {"commandId": 10,
     "commandType": "RESUME",
     "commandData": ""
    }
    *END

The command output file will register all results of a command execution while in the running state. Therefore, for the previous commands, only one result will be seen, corresponding to the pause command. Notice that the output is only written once all active explorers have acknowledged the command; so if some explorers are doing long calculations at a given step, you will have to wait for those calculations (and the step) to finish, before the command is actually effective.

The contents of the result are a JSON object with two fields:

- commandId: Matches the commandId of the command this result corresponds to.
- output: A string with the output for the given command.

The output is a string codifying a list of steps where the different explorers are after executing the command, which corresponds to the next step they plan to run. The format is: ``X@Y:Z`` where X is the explorer id, Y is the task id, and Z is the step number.

An example output, corresponding to the previous input commands:

.. code-block:: text

    {"commandId": 1, "output": "1@1:2 2@1:2 3@1:2"}
    *END

Notice that each reply is also terminated by a ``*END`` line. In the example, explorers 1, 2 and 3 were active when they received the command, and all of them were ready to start step 2 of task 1 at that time.
    

