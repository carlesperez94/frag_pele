Frequently Asked Questions
==========================

This section compiles some of the most often encountered issues running an
AdaptivePELE simulation.

    - Not having AdaptivePELE in your PYTHONPATH variable. If you get an error like::

        ImportError: No module named AdaptivePELE

      python is not able to find the adaptive module. You need to add the path as::

        export PYTHONPATH="/path/to/repo/AdaptivePELE:PYTHONPATH"

    - Similar to the previous message, if the module affected is one of atomset,
      freeEnergies.utils, RMSDCalculator or SymmetryContactMapEvaluator, the issue
      might be that the cython extension of AdaptivePELE are not compiled, to do
      so find the setup.py script in the AdaptivePELE repository and run::

        python setup.py build_ext --inplace

    - If you get an error like::

        ValueError: The input pdb file/string was empty, no atoms loaded!

      tipically one of two things happened, either the pdb file passed is not
      correctly formed or empty or the resname introduced in the control file does
      not match the one in the control file.

    - Simulation dies with no error. This happen sometimes in a PELE simulation,
      almost always is related to a PELE error that was not handled properly. One
      frequent source of this issue in simulation running in the life or nord
      clusters is a bug of mpi that requieres to declare the following environment
      variables::

        export OMPI_MCA_coll_hcoll_enable=0
        export OMPI_MCA_mtl=^mxm

    - If you get an error including the message::

        No trajectories to cluster! Matching path:

      adaptivePELE has not been able to find trajectories to cluster. This is
      tipically because the simulation died before it could produce any
      trajectories, in such case check for PELE errors.

    - If you get an error like::

        UnsatisfiedDependencyException: Missing package....

      it means that an optional dependence is missing. That is a package that
      is not absolutely required to run AdaptivePELE but it is needed for some
      particular feature that you are using cannot be found. To solve it ensure
      that the package in question is installed and accessible to the python
      interpreter (e.g check the PYTHONPATH value)

    - There is a weird atom named *DUM* in my trajectories! When using MD with
      an spherical box for the ligand, a dummy atom is introduced to act as the
      center of the box and is listed in the trajectories as *DUM* atom in
      a residue also named *DUM*. This atom is massless so it will not be moved
      during the simulation.
