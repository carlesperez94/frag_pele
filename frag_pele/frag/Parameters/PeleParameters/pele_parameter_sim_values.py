# Python Imports

# Third-Party Imports

# Project Imports


class PeleParameterSimulationValues:

    def __init__(self, cpus: int, steps: int, pele_eq_steps: int, min_overlap: float, max_overlap: float,
                 temperature: int, seed: int, steering: int, translation_high: float, rotation_high: float,
                 translation_low: float, rotation_low: float, radius_box: float):
        """
        :param cpus: Number of CPUs that will be used to perform the simulation.
        :param steps: PELE steps to do in each GS.
        :param pele_eq_steps: Number of PELE steps that we want to perform in the equilibration.
        :param min_overlap: Minimun value of overlapping factor that will be replaced in the control file.
        :param max_overlap: Maximum value of overlapping factor that will be replaced in the control file.
        :param temperature: Temperature to add in the control file.
        :param seed: Seed in the PELE's control file.
        :param steering: Steering parameter of the PELE's control file.
        :param translation_high: Translation (high value) of the PELE's control file. In anstrongs.
        :param translation_low: Translation (low value) of the PELE's control file. In anstrongs.
        :param rotation_high: Rotation (high value) of PELE's control file. In rad.
        :param rotation_low: Rotation (low value) of PELE's control file. In rad.
        :param radius_box: Radius box size of PELE's control file. In anstrongs.
        """
        self._cpus = cpus
        self._steps = steps
        self._pele_eq_steps = pele_eq_steps
        self._min_overlap = min_overlap
        self._max_overlap = max_overlap
        self._temperature = temperature
        self._seed = seed
        self._steering = steering
        self._translation_high = translation_high
        self._rotation_high = rotation_high
        self._translation_low = translation_low
        self._rotation_low = rotation_low
        self._radius_box = radius_box

    # Properties
    @property
    def cpus(self):
        return self._cpus

    @property
    def steps(self):
        return self._steps

    @property
    def pele_eq_steps(self):
        return self._pele_eq_steps

    @property
    def min_overlap(self):
        return self._min_overlap

    @property
    def max_overlap(self):
        return self._max_overlap

    @property
    def temperature(self):
        return self._temperature

    @property
    def seed(self):
        return self._seed

    @property
    def steering(self):
        return self._steering

    @property
    def translation_high(self):
        return self._translation_high

    @property
    def rotation_high(self):
        return self._rotation_high

    @property
    def translation_low(self):
        return self._translation_low

    @property
    def rotation_low(self):
        return self._rotation_low

    @property
    def radius_box(self):
        return self._radius_box
