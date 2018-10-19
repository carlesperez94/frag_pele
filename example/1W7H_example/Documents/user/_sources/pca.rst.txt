.. _sec-pca:

**************************************
PCA for the perturbation phase of PELE
**************************************

You can do a PCA (Principal Component Analysis) of a set of structures, so that the eigenvalues and eigenvectors representing the principal components of the variation in the conformation of the macromolecular structure in the set can be used as the eigenvalues and eigenvectors directing the perturbation moves for the macromolecule in PELE (see :ref:`sec-anm`).

.. warning::

  When running the script, you input PDB files may be modified, since they are rewritten in a PDB format consistent to that used by ProDy, and the HIE and HID residues are renamed as HIS residues.

To do so, you need an NMD file (see :ref:`sec-fileFormats-NMD`) codifying the eigenvalues and eigenvectors. To ease the generation of this file, you can use the script :file:`scripts/pca/calculatePCA4PELE.py`, which has been provided by Christoph Grebner, from AstraZeneca. You will need Python 2.7 installed, with the ProDy Python package also installed (http://prody.csb.pitt.edu/index.html).

To see how to use the script, you will perform some calculations with the data in :file:`samples/aspirin/`. This example contains a subdirectory :file:`normalModes/`, with a generated trajectory from which to do the PCA study.

Go to that directory, and run the :program:`calculatePCA4PELE.py` script:

.. code-block:: console

    $ cd samples/aspirin/normalModes
    $ python /path/to/scripts/pca/calculatePCA4PELE.py --pdb ain_trajectory_1.pdb
    @> ProDy is configured: verbosity='info'
    Create Ensemble
    Found 1 Chain(s) in ain_trajectory_1_ca
    None
    @> Starting iterative superposition:
    @> Step #1: RMSD difference = 7.2761e-01
    @> Step #2: RMSD difference = 1.5657e-04
    @> Step #3: RMSD difference = 1.2805e-07
    @> Final superposition to calculate transformations.
    Calculate PCA
    @> Covariance is calculated using 82 coordinate sets.
    PCA
    <PCA: ain_trajectory_1_ca (20 modes; 119 atoms)>
    PCA is saved in: ain_trajectory_1_ca_pca_modes.nmd
    Calculate ANM

    ANM is saved in: ain_trajectory_1_ca_m1_anm_modes.nmd

You will then obtain two NMD files, one with the PCA components (calculated with respect to the first structure in the trajectory), and another with the ANM modes for the first structure in the trajectory. For this example, you can copy the PCA NMD file to the main directory of the example, and execute the simulation configured in the :file:`pca_control_file`.
    
.. code-block:: console
    
    $ cp ain_trajectory_1_ca_pca_modes.nmd ../pca_modes.nmd
    $ cd ..
    $ /path/to/PELE pca_control_file

You will observe the report and trajectory generated under the :file:`results/` subdirectory (you must previously create it if it doesn't exist).

For more information on how to configure the control file of PELE to include precalculated normal modes, see :ref:`sec-anm-preloadedModesIn`.

Usage
=====

.. code-block:: console

    $ python /path/to/scripts/pca/calculatePCA4PELE.py --help
    usage: calculatePCA4PELE.py [-h] [--pdb PDB] [--selection SELECTION]
				[--compare] [--vmd] [--debug] [--ref REF]
				[--print_mode] [--print_diff]
				[--print_file PRINT_FILE] [--anm_only]

    Calculate PCA and ANM for given structures

    optional arguments:
      -h, --help            show this help message and exit
      --pdb PDB             specify the pdbs which should be used.
					If not specified, all PDBs in the current folder are taken by default.
					Possibilities:
					use all PDBs in directory: --pdb all
					give a list of PDBs, either single or trajectory file: --pdb "name1.pdb name2.pdb name3.pdb name4.pdb"

      --selection SELECTION
			    Specify the selection for normal mode analysis
					valid options:
					--selection "calpha": 
					--selection "backbone"
					--selection "all"

      --compare             compare PCA and ANM
      --vmd                 Visualize result in VMD
      --debug               Print debug information
      --ref REF             Select a reference structure for PCA and ANM
      --print_mode          Print the square fluctuations of the modes
      --print_diff          Print difference between ANM and PCA modes
      --print_file PRINT_FILE
			    Print modes from a NMD file. This will exit the
			    program after printing
      --anm_only            Only perform ANM analysis

As in the example, you can provide a trajectory file with the ``--pdb`` option. Alternatively, you can give a list of PDB files (``--pdb "name1.pdb name2.pdb name3.pdb name4.pdb"``) and mark one of them as a reference (for example, ``--ref "name2.pdb") if you do want to use a reference different from the first PDB (in the case of a trajectory, it is always the first model the one used as a reference). By default, both the PCA and ANM studies are done regarding the alpha carbons (as if ``--selection "calpha"``), but you can specify a different selection with the ``--selection`` command line option.
    
