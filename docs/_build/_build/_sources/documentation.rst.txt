==============
Documentation
==============

.. toctree::
   :maxdepth: 2

Before starting, take into account that only the Basic parameters are mandatory. Otherwise, Optional parameters have
predefined default values. You can find and modify them in ``constants.py``.

Required parameters
-------------------

- **-cp**, **---complex_pdb**: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be used
  as the core structure. Remember to rename the ligand chain with a different character in order to detect it.
- **-sef**, **---serie_file**:  Name of the tabular file which must contain the instructions required to perform several
  successive growings, using different fragments or different growing positions.

Optional parameters
-------------------

Protocol
----------

- **-HT**, **---highthroughput**: High throughput protocol.
- **-EX**, **---explorative**: Explorative protocol.
- **-sc**, **---sampling_control**: Templatized control file that will be used in the sampling simulation instead of the default one.
- **default option**: Sp protocol.

For more information refer to the protocol tab.

Frag configuration
``````````````````
This are the options directly related to Frag itself.

- **-x**, **---growing_steps**: Number of Growing Steps (GS). By default: 10. Set 6 for small ligands and 2 to perform HT mode.
- **-cr**, **---criteria**: Name of the column of the report file to select the structures that will spawn in the
  next GS. Additionally, this parameter will be the selection criteria to extract the best structure after completing
  the growing. By default: 'Binding Energy'.
- **-rst**, **---restart**: If set it will continue from the last GS detected. If all GS are finished it will restart
  in the equilibration phase.
- **-cc**, **---c_chain**: Chain name of the core. By default: 'L'.
- **-fc**, **---f_chain**: Chain name of the fragment. By default: 'L'.
- **-tc**, **---clash_thr**: Threshold distance to detect intramolecular clashes between the fragment and the core.

PlopRotTemp configuration
`````````````````````````
- **-pl**, **--plop_path**: Absolute path to PlopRotTemp.py. Its included in FragPELE's library. It's default is defined in the constants.py.
- **-sp**, **---sch_python**: Absolute path to Schrodinger's python. It's default is defined in the constants.py
- **-rot**, **---rotamers**: Degrees of rotation to select the rotamers in the library to do the side-chain prediction. By default: "30.0".

PELE configuration
``````````````````
Most of the following arguments are only focused in localize PELE's directories in the user machine. If you want to
modify parameters directly related to PELE configuration (control_file) it is highly recommendable to read PELE's
documentation.

- **-d**, **---pele_dir**: Complete path to Pele_serial.
- **-c**, **---contrl**: Path to PELE's control file templatized.
- **-l**, **---license**: Absolute path to PELE's licenses folder.
- **-r**, **---resfold**: Name for PELE's results folder. By default: 'growing_output'.
- **-rp**, **---report**: Suffix name of the report file from PELE.  By default: 'report'
- **-tj**, **---traject**: Suffix name of the trajectory file from PELE. By default: 'trajectory'
- **-cs**, **---cpus**: Number of cores (computational power) to paralellize PELE's simulations. By default: 48.
- **-stp**, **---steps**: Number of simulation steps inside each GS. By default: 6.
- **-es**, **---pele_eq_steps**: Number of PELE steps in sampling simulation. By default: 20.
- **-miov**, **---min_overlap**: Minimum value of overlapping factor used in the control_file of PELE. By default: 0.5.
- **-maov**, **---max_overlap**: Maximum value of overlapping factor used in the control_file of PELE. By default: 0.7.
- **-tmp**, **---temperature**: Temperature value to add in the control file. If the temperature is high more steps of
  PELE will be accepted when applying the Metropolis Criteria. By default: 1000.
- **-sd**, **---seed**: Seed to get the random numbers of the PELE's simulation.
- **-st**, **---steering**: Steering to use in the PELE's simulation. By default 0.
- **-trh**, **---translation_high**: Translation to use in the PELE's simulation to displace the ligand. High value. By default: 0.05.
- **-trl**, **---translation_low**: Translation to use in the PELE's simulation to displace the ligand. Low value. By default: 0.02.
- **-roth**, **---rotation_high**: Rotation (radians) to use in the PELE's simulation to perturb the ligand. High value. By default: 0.10.
- **-rotl**, **---rotation_low**: Rotation (radians) to use in the PELE's simulation to perturb the ligand. Low value. By default: 0.05.
- **-rad**, **---radius_box**: Size of the radius to define the box in the PELE's simulation where the ligand will be perturbed. By default: 4.

Clustering and spawning configuration
`````````````````````````````````````
The following parameters are extracted from `AdaptivePELE <https://github.com/AdaptivePELE/AdaptivePELE/>`_. It is not
recommendable to modify default parameters.

- **-dis**, **---distcont**: Distance used to determine which amino acids are in contact with the ligand to generate
  different clusters of structures to initialize the next GS. By default: 4.
- **-ct**, **---threshold***: Threshold distance used in the clustering. By default: 0'3.
- **-e**, **---epsilon**: An epsilon fraction of processors are distributed proportionally to the value of a metric,
  and the rest are inverselyProportional distributed. A param n can be specified to only consider the n clusters with
  best metric. By default: 0'5.
- **-cn**, **---condition**: Selects whether to take into account maximum or minimum values in epsilon related spawning,
  values are min or max. By default: min.
- **-mw**, **---metricweights**: Selects how to distribute the weights of the cluster according to its metric, two
  options: linear (proportional to metric) or Boltzmann weigths (proportional to exp(-metric/T). Needs to define the
  temperature T. By default: 'linear'.
- **-ncl**, **--.nclusters**: Number of initial structures that we want to use in each new GS. By default: 5.

The next ones have been created specially for FragPELE's simulations:

- **-pdbf**, **---pdbout**: Folder where PDBs selected to spawn in the next GS will be stored. By default: 'PDBs_growing'.
- **-ban**, **---banned**: List of lists of quartets of PDB atom names that form the banned dihedrals.
  These dihedrals will be examined after each growing step and if they cross the limit they will not be selected to
  next initialize the new growing step.
- **-lim**, **---limit**: Maximum degrees that can accept the banned dihedral.

Output configuration
``````````````````````
- **--mae**: Output mae files with binding energy trajectory as property inside the file.
- **---rename**: To avoid core renaming.











