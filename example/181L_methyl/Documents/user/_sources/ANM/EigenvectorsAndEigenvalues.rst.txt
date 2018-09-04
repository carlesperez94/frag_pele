Treatment of eigenvector and eigenvalues by the PELE ANM module
===============================================================

The eigenvectors of the hessian matrix are used by the ANM algorithm in
order to calculate the displacement vector for each node. The sequence
of steps that PELE follows when dealing with eigenvectors in order to
calculate the final displacement vector is the following:

-  The eigenvectors, and they corresponding eigenvalues, are computed
   both at the same time.
-  The vectors are normalized.
-  A mode is chosen among the set of eigenvectors
-  For the selected mode, a new direction is set.
-  The chosen mode is scaled to a certain modulus.
-  A final direction vector is obtained mixing the selected mode with
   the other ones.
-  For each atom, a target position is calculated using the previously
   calculated move vector.
-  A set of harmonic constraints is added to energy computation
   function. The equilibrium position for each constraint (i.e. for each
   atom) is the target position obtained at the previous step, while the
   "k" constant is the same for all the constraints.

The previous sequence of steps describes the algorithm, but the actual
procedure that is followed at each step is parameterizable by the user.
So the user can choose how the vector will be normalized or how the
modulus will be calculated for the selected mode. All the options are
intercheangeable, leading to a wide range of combinations. Information
about each specific option can be obtained in PELE's documentation.

The purpose of the eigenvectors in the previous algorithm is clear: a
prefered normal mode will be chosed among them. The eigenvalues,
however, may or may not be used, depending on the algorithms chosen by
the user for each step.

The following section explains which are the places where eigenvalues
are used, and for each one, indicates in an algorithmic way the math
that is performed, and also which parameteres have to be enabled in
order to perform the calculation:

From now on, ``values`` refers to the array of eigenvalues, and
``vectors`` refers to the array of eigenvectors (both with dimension
``numberOfModes``). Functions like ``inverseEuclideanNorm``,
``vectorByScalar``, etc. represent the mathematical functions they are
called after. Also ``sqrt`` means squareroot, and ``rand(arg, arg)``
returns a random number between the given interval.

-  At the normalization step:

   -  If the parameter "thermalScaling" has been set to true, the
      eigenvalue associated with each vector will be used in order to
      normalize the vector.
   -  Calculation:

      .. code-block:: text

        For i from 1 to numberOfModes
	  val = values[i]
          vec = vectors[i]
          factorNormalization = inverseEuclideanNorm(vec) * (1 / sqrt(hessianConstant * val)) //"hessianConstant" is a constant introduced by the user
          vec = vectorByScalar(vec, factorNormalization)

-  When calculating the modulus for the chosen mode:

   -  If the option "moveMagnitudeGeneration" has been set to
      "scaledBiasedRandom" the eigenvalues of the chosen mode and the
      first mode will be used in order to set an scaling factor for the
      modulus.
   -  Calculation:

      .. code-block:: text

        val = values[chosenMode]
        base = values[0]
        scaledFactor = (0.6 * val) / base + 0.4
        scaledMove = displacement / scaledFactor //"displacement" is a constant introduced by the user
        modulus = scaledMove * rand(0.625, 1.0)

-  When mixing the normal modes:

   -  Eigenvalues are also used for mixing the chosen mode with the
      another ones. However, PELE by default does not mix the modes, so
      in order to do so, the option "modesMixingOption" must be set to
      either "mixAllModesEquallyRandom" or "mixMainModeWithOthersModes",
      although the only one which will use the eigenvalues is the latter
      one.
   -  Calculation (for "mixMainModeWithOthersModes" option):

      .. code-block:: text

        For i from 1 to numberOfModes
          if i != from chosenModeIndex
            vec = vectors[i]
            magnitude = rand(-1.0, 1.0) / sqrt(values[i]) * inverseEuclideanNorm(vec)
	    moveVector = linearCombination(magnitude, vec, 1.0, moveVector)

	chosenEigenVectorWeight = mainModeWeight * inverseEuclideanNorm(chosenEigenVector)
	otherModesWeight = (1.0 - mainModeWeight) * inverseEuclideanNorm(moveVector)
	moveVector = linearCombination(chosenEigenVectorWeight, chosenEigenVector, otherModesWeight, moveVector);

	moveVector = normalizeVector(moveVector);

