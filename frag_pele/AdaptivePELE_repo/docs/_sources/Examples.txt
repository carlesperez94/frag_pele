User Manual
===========


Installation
------------

In order to have a running copy of AdaptivePELE, you need to install and compile cython files in the base folder with::

    cd AdaptivePELE
    python setup.py build_ext --inplace

Also, if AdaptivePELE was not installed in a typical library directory, a common option is to add it to your local PYTHONPATH::

    export PYTHONPATH=<location_of_AdaptivePELE>


Overview
--------

This page tries to give an overview of AdaptivePELE, enough to get it
up and running. The program aims to enhance the exploration of standard
molecular simulations by iteratively running short simulations, 
assessing the exploration with a clustering, and spawning new trajectories
in *interesting* regions using the multi-armed bandit as a framework.
The algorithm can be summarized in the flow diagram shown below:

.. image:: adaptiveDiagram.png
    :align: center
    :alt: image trouble


Control file outline
--------------------

AdaptivePELE is a Python program called with a parameter: the control 
file.

::

    python -m AdaptivePELE.adaptiveSampling controlFile.conf

The control file is a JSON document that permits tuning its behavior. It contains 4 blocks: 
the general parameters, simulation parameters, clustering
parameters and spawning parameters. The first block refers to general
parameters such as the output path, initial structures, or whether to restart (resume) a simulation. 
The other three blocks configure the three iterative steps of an adaptive run, and have a 
``{"type":str, "params":{}, "optionalClasses":{}}`` structure. Namely, the **simulation** block 
chooses and tunes the propagation algorithm, typically PELE, the **clustering** block tunes the clustering method,
cluster sizes... and the **spawning** block selects the strategy to choose the most interesting
structures. 

::

        {
            "generalParams" : {},
            "simulation": {"type":str, "params":{}},
            "clustering": {"type":str, "params":{}},
            "spawning": {"type":str, "params":{}}
        }


generalParams block
-------------------

These parameters control the general aspects of the simulation, such as the initial structures or the output path.
For example, if you consider that a simulation did not explore the energy landscape sufficiently, 
it can be automatically resumed. In this case, if you change the clustering parameters from the previous run,
it will need the trajectory files to recluster all snapshots. Otherwise, it is enough with the automatically 
handled clustering objects.
If you want to test the effect of using a particular cluster threshold or spawning parameter (see below),
using debug=true may be a great option.



Parameters
..........

The general parameters block has five possible fields:

* **restart** (*boolean*, default=True): This parameter specifies whether you want to
  continue a previous simulation or not. It requires the last epoch's clustering object, 
  which is stored as a binary file to optimize disk usage.

* **debug** (*boolean*, default=False): Run adaptive in debug mode, without calling
  any propagation algorithm. Typically used to tinker with the clustering/spawning.

* **outputPath** (*string*, mandatory): The path where the results of the simulation will be
  written

* **initialStructures** (*list*, mandatory): The path(s) to the intial structure(s)

* **writeAllClusteringStructures** (*boolean*, default=False): Whether to write all the cluster 
  center structures as pdbs. Setting it to True is inneficient, and cluster center structures 
  can still be recovered from the binary clustering object.

Additionaly, it can also have a nativeStructure parameter, a string containing
the path to the native structure. This structure will only be used to correct
the RMSD in case of symmetries. The symmetries will also have to be specified
(see `clustering block`_).

Example::

    "generalParams" : { 
        "restart": false,
        "debug" : false,
        "outputPath":"simulationOutput",
        "nativeStructure" : "nativeStructure.pdb",
        "initialStructures" : ["initial1.pdb", "initial2.pdb"]
    }

simulation block
-----------------

Currently, there are two implemented simulation types: 

* **pele**. `PELE <https://pele.bsc.es/pele.wt>`_ is a great tool to efficiently explore the energy landscape. Parameters have been optimized for its use.

* **test**. The test type has no real use outside of testing. 

We plan to implement an MD type in future versions.


Templetized PELE control file
.............................

In order to run adaptivePELE, PELE control file needs to be templetized. In particular:

* **MultipleComplexes"**: AdaptivePELE requieres the use of multiple complexes: ``"MultipleComplex": [ $COMPLEXES ]``

* **seed**: The seed needs to be templetized: ``"seed": $SEED``

* **outputPath**: The output path needs to be changed for: ``"reportPath": "$OUTPUT_PATH/report"``

* **numberOfPeleSteps**: The number of pele steps needs to be templetized as ``"numberOfPeleSteps": $PELE_STEPS``

parameters
..........

When using PELE as a propagator, the following parameters are mandatory:

* **iterations** (*integer*, mandatory): Number of adaptive sampling iterations to run
* **processors** (*integer*, mandatory): Number of processors to use with PELE
* **peleSteps** (*integer*, mandatory):  Number of PELE steps in a epoch (iteration)
* **seed** (*integer*, mandatory): Seed for the random number generator
* **controlFile** (*string*, mandatory): Path to the templetized PELE control file (see below)

Optionally, you can also use the following parameters:

* **data** (*string*, default=MareNostrum or Life cluster path): Path to the Data folder needed for PELE
* **documents** (*string*, default=MareNostrum or Life cluster path): Path to the Documents folder needed for PELE
* **executable** (*string*, default=MareNostrum or Life cluster path): Path to the Pele executable folder, it is already

Additionally, the block may have an exit condition that stops the execution:

* **exitCondition** (*dict*, default=None): Block that specifies an exit condition for the simulation.
  Currently only the metric type is implemented, this type accepts a
  *metricCol* which represents a column in the report file and a *exitValue*
  which represents a value. The simulation will terminate after the metric
  written in the *metricCol* reaches a value smaller than *exitValue*. An example of the exit condition block 
  that would terminate the program after a trajectory reaches a value of less than 2 for the fifth column 
  (4th starting to count from 0) of the report file would look like::

    "exitCondition" : {
        "type" : "metric",
        "params" : {
            "metricCol" : 4,
            "exitValue" : 2.0
        }
    }

Example of simulation block::

    "simulation": {
        "type" : "pele",
        "params" : { 
            "iterations" : 25,
            "processors" : 128,
            "peleSteps" : 4,
            "seed": 30689,
            "controlFile" : "templetizedPELEControlFile.conf"
        }   
    }


clustering block
----------------

Currently there are two functional types of clustering:

* **rmsd**, which solely uses the ligand rmsd

* **contactMap**, which uses a protein-ligand contact map matrix

These clusterings are based on the leader algorithm, an extremely fast clustering method that in the 
worst case makes *kN* comparisons, where *N* is the number of snapshots to cluster and *k* the number of existing clusters.
The procedure is as follows. Given some clusters, a conformation is said to belong to a cluster 
when it differs in less than a certain metric threshold (e.g. ligand RMSD)
to the corresponding cluster center. Cluster centers are always compared in the same order, and,
if there is no similar cluster, it generates a new one. 

Aside form the speed, a big advantage of using this method 
is that it permits the user to define different criteria in different regions.
This way, we can optimize the number of clusters, giving more importance to regions with more interactions,
potentially being more metastable.

In order to measure the potential metastability, 
we use the ratio of the number of protein-ligand heavy atom contacts over the number of ligand heavy atoms, *r*. 
Two atoms are considered to be in contact if the distance between 
them is less than a certain **contactThreshold** (8Å by default). Although these values depend on the particular
protein-ligand geometry and ligand size, this measure is more ligand-independent compared to the overall
number of contacts and a value of 1 typically indicates that the ligand is on the surface entering a protein pocket.
We encourage the use of default parameters with very few exceptions such as in the study 
of the diffusion of ions or tiny molecules (e.g. a oxygen molecule).


thresholdCalculator
...................

* **constant**, where all clusters have the same threshold. A sound value may be 3Å.

* **heaviside** (default), where thesholds (values) are assigned according to a set of step functions that 
  vary according to a ratio of protein-ligand contacts and ligand size , *r*, (conditions, see below). The values and conditions 
  of change are defined with two lists. The condition list is iterated until *r* > condition[i], and the used
  threshold is values[i]. If r <= conditions[i] for all i, it returns the last element in values. 
  Thresholds typically vary from 5Å in the bulk to 2Å in protein pockets. This method is preferred, as it
  optimizes the number of clusters, giving more importance to regions with more contacts and interactions, 
  where metastability occurs. Default values: [2,3,4,5], default conditions: [1, 0.75, 0.5].

rmsd clustering
...............

In the **rmsd** clustering, if the RMSD between two ligand conformations is less than 
a certain threshold, the conformation is added to the cluster, and otherwise, a new cluster 
is generated.


contactMap clustering
.....................

The **contactMap** uses the similarity between protein-ligand contact maps.
The contact map is a boolean matrix with the protein
atoms (or a subset of them, typically one or two per residue) as columns and 
ligand atoms (typically only heavy atoms) as rows, and a value of True indicates a contact.
There are currently three implemented methods to evaluate the similarity of contactMaps:

* **Jaccard**, which calulates the Jaccard Index (`Wikipedia page <https://en.wikipedia.org/wiki/Jaccard_index>`_). The recommended values using the heaviside threshold calculator are [0.375, 0.5, 0.55, 0.7] for the conditions [1, 0.75 , 0.5].

* **correlation**, which calculates the correlation between the two matrices and

* **distance**, which evaluates the similarity of two contactMaps by calculating the ratio of the number of differences over the average of elements in the contacts maps.

parameters
..........

* **ligandResname** (*string*, default=""): Ligand residue name in the PDB (if necessary)
* **ligandChain** (*string*, default=""): Ligand chain (if necessary)
* **ligandResnum** (*int*, default=0): Ligand residue number (if necessary). If 0 or not specified, it is ignored. The ligand ought to be univoquely identified with any combination of this and the two former parameters
* **contactThresholdDistance** (*float*, default=8): Maximum distance at which two atoms have to
  be separated to be considered in contact
* **symmetries** (*list*, default=[]) List of symmetry groups of key:value maps with the names of atoms
  that are symmetrical in the ligand
* **similarityEvaluator** (*string*, mandatory)  Name of the method to evaluate the similarity
  of contactMaps, only available and mandatory in the contactMap clustering
* **alternativeStructure** (*bool*, default=False): It stores alternative spawning structures within each cluster to be used in the spawning (see below). Any two pairs of alternative structures within a cluster are separated a minimum distance of cluster_threshold_distance/2.

Example
.......

A typical setting of the rmsd clustering is::

    "clustering" : { 
        "type" : "rmsd",
        "params" : { 
            "ligandResname" : "AEN",
            "contactThresholdDistance" : 8,  
            "symmetries": [{"3225:C3:AEN":"3227:C5:AEN","3224:C2:AEN":"3228:C6:AEN"}, {"3230:N1:AEN": "3231:N2:AEN"}]
        },  
        "thresholdCalculator" : { 
            "type" : "heaviside",
            "params" : { 
                "values" : [2, 3, 4, 5], 
                "conditions": [1.0, 0.75, 0.5]
            }   
        }   
    }

which, given the default options, is equivalent to::

    "clustering" : { 
        "type" : "rmsd",
        "params" : { 
            "ligandResname" : "AEN",
            "symmetries": [{"3225:C3:AEN":"3227:C5:AEN","3224:C2:AEN":"3228:C6:AEN"}, {"3230:N1:AEN": "3231:N2:AEN"}]
        }  
    }


In this exemple, clusters having a contacts ration greater than 1 have a
treshold of 2, those with contacts ratio between 1 and 0.75 have a treshold of
3, between 0.75 and 0.5 a threshold of 4 and the rest have a threshold size of
5. This means that for greater contacts ratio, typically closer to the binding site,
the cluster size will be smaller and therefore those regions will be more
finely discretized.


spawning block
---------------

Finally, trajectories are spawned in different *interesting* clusters, according to a reward function.
There are several implemented strategies:

* **sameWeight**: Uniformly distributes the processors over all clusters

* **inverselyProportional**: Distributes the processors with a weight that is inversely proportional to the cluster population.

* **epsilon**: An *epsilon* fraction of processors are distributed proportionally to the value of a metric, and the rest are inverselyProportional distributed.  A param **n** can be specified to only consider the *n* clusters with best metric.

* **variableEpsilon**: Equivalent to epsilon, with an epsilon value changing over time

* **independent**: Trajectories are run independently, as in the original PELE. It may be useful to restart simulations or to use the analysis scripts built for AdaptivePELE.

* **UCB**: Upper confidence bound.

* **FAST**: FAST strategy (see J. Chem. Theory Comput., 2015, 11 (12), pp 5747–5757).

According to our experience, the best strategies are **inverselyProportional** and **epsilon**, guided with either PELE binding energy or the RMSD to the bound pose if available.


density calculator
..................

Each cluster is assigned a relative density of points compared to other clusters.
Again, and in analogy to the threshold calculator, the aim is to give more emphasis to interesting regions.
There are two types of density calculators:

* **constant** (or **null**, default), which assigns the same density to all the clusters regardless of the number of contacts

* **heaviside**, which assigns different densities using a heaviside function, much like the thresholdCalculator (values and conditions are mandatory)

* **continuous**, which assings increasing densities for an increasing number of contacts. Default values, if **r** > 1, density = 8, otherwise, density = 64.0/(-4 **r** + 6)^3


parameters
..........

* **reportFilename** (*string*, mandatory): Basename to match the report file with metrics. E.g. "report". 

* **metricColumnInReport** (*integer*, mandatory): Column of the report file that contains the metric of interest (zero indexed)

* **epsilon** (*float*, mandatory in **epsilon** spawning): The fraction of the processors that will be assigned according to the selected metric

* **metricWeights** (*string*, default=linear): Selects how to distribute the weights of the cluster according to its metric, two options: linear (proportional to metric) or Boltzmann weigths (proportional to exp(-metric/T). Needs to define the temperature **T**.

* **T** (*float*, default=1000): Temperature, only used for Boltzmann weights

The following parameters are mandatory for **variableEpsilon**:

* **varEpsilonType** (*string*,default=linear): Selects the type of variation for the epsilon value. At the moment only a linear variation is implemented
* **maxEpsilon** (*float*): Maximum value for epsilon
* **minEpsilon** (*float*): Minimum value for epsilon
* **variationWindow** (*integer*): Last iteration over which to change the epsilon value
* **maxEpsilonWindow** (*integer*): Number of iteration with epsilon=maxEpsilon
* **period** (*integer*) Variation period (in number of iterations)


Examples
..........

::

    "spawning" : {
        "type" : "inverselyProportional",
        "params" : {
            "reportFilename" : "report"
        }
    }


::

    "spawning" : {
        "type" : "epsilon",
        "params" : {
            "reportFilename" : "report",
            "metricColumnInReport" : 5,
            "epsilon" : 0.25
        },
        "density" : {
            "type" : "continuous"
        }
    }


Control File Examples
---------------------

Example 1
.........

The first example makes use of default parameters (used in the AdaptivePELE paper).

::

    {
        "generalParams" : {
            "restart": false,
            "outputPath":"example1",
            "nativeStructure" : "native.pdb",
            "initialStructures" : ["initial1.pdb", "initial2.pdb"]
        },

        "simulation": {
            "type" : "pele",
            "params" : {
                "iterations" : 25,
                "processors" : 128,
                "peleSteps" : 4,
                "seed": 30689,
                "controlFile" : "templetizedPELEControlFile.conf"
                
            }
        },

        "clustering" : {
            "type" : "rmsd",
            "params" : {
                "ligandResname" : "AEN"
            }
        },

        "spawning" : {
            "type" : "inverselyProportional",
            "params" : {
                "reportFilename" : "report"
            }
        }
    }


To summarize, below there is a screenshot of a simple functional control file:

Example 2
.........

A more complete (although not so comprehensible) example::

    {
        "generalParams" : {
            "restart": true,
            "debug" : false,
            "outputPath":"example2",
            "writeAllClusteringStructures": false,
            "nativeStructure" : "native.pdb",
            "initialStructures" : ["initial1.pdb", "initial2.pdb"]
        },

        "spawning" : {
            "type" : "epsilon",
            "params" : {
                "reportFilename" : "report",
                "metricColumnInReport" : 5,
                "epsilon":0.1
            },
            "density" : {
                "type" : "null"
            }
        },

        "simulation": {
            "type" : "pele",
            "params" : {
                "executable" : "PELE++/bin/rev12025/Pele_rev12025_mpi",
                "data" : "PELE++/data/rev12025/Data",
                "documents" : "PELE++/Documents/rev12025",
                "iterations" : 25,
                "processors" : 51,
                "peleSteps" : 4,
                "seed": 30689,
                "controlFile" : "/gpfs/scratch/bsc72/bsc72755/adaptiveSampling/data/3ptb/3ptb_a_1000.conf"
                
            }
        },

        "clustering" : {
            "type" : "rmsd",
            "params" : {
                "ligandResname" : "AEN",
                "contactThresholdDistance" : 8, 
                "symmetries": [{"3225:C3:AEN":"3227:C5:AEN","3224:C2:AEN":"3228:C6:AEN"}, {"3230:N1:AEN": "3231:N2:AEN"}]
            },
            "thresholdCalculator" : {
                "type" : "heaviside",
                "params" : {
                    "values" : [2, 3, 4, 5],
                    "conditions": [1.0, 0.75, 0.5]
                }
            }
        }

    }


Output
------

The output for each epoch is redirected to a different folder, with a name that matches the epoch number. For example, if we run three epochs, we will have three folders named
0, 1, and 2.
Aside from the regular simulation program output each directory contains a clustering subdirectory with the clustering summary information, and 
eventually, the cluster center pdb files and the clustering object. This clustering object is used to restart simulations, and only that of the last
finished epoch is kept for disk usage optimization. 
If we change a clustering parameter in a restart run, AdaptivePELE will recluster all the snapshots, which will fail if previous trajectories are not present.


Analysis
........

In order to analyse simulation results, a bunch of scripts are provided in ``AdaptivePELE/analysis``. Get help to run them with: ``python <script> -h``

Example to print column 5 evolution with gnuplot::

    python -m AdaptivePELE.analysis.plotAdaptive 4 2 5 report_ -rmsd | gnuplot -persist

It prints the evolution of column 5 (e.g. RMSD) in report_* files with lines (-rmsd) in epochs of 4 steps.

Example to print BE against RMSD with gnuplot::

    python -m AdaptivePELE.analysis.plotAdaptive 4 5 6 report_ -be | gnuplot -persist

It prints the column 6 against column 5 with points (-be). Epoch length is ignored in this case

To plot the evolution of the number of clusters along the simulation::

    python -m AdaptivePELE.analysis.numberOfClusters -filename "plot"

It shows the evolution of the total number of clusters, and the number of clusters divided in different densities and cluster thresholds.
It also prints a histogram with the ratio of counts *r* (see above). When ``-filename`` is provided, it saves the plots as png files.
