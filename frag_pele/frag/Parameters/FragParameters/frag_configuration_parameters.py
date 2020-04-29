# Python Imports

# Third-Party Imports

# Project Imports


class FragConfigurationParameters:

    def __init__(self, growing_steps: int, criteria: str, sampling_control: str):
        """
        :param growing_steps: Number of Growing Steps (GS).
        :param criteria: Name of the column of the report file to select the structures that will spawn in the
        next GS. Additionally, this parameter will be the selection criteria to extract the best
        structure after completing the growing.
        :param sampling_control: templatized control file to be used in the sampling simulation.
        """
        self._growing_steps = growing_steps
        self._criteria = criteria
        self._sampling_control = sampling_control

    # Properties
    @property
    def growing_steps(self):
        return self._growing_steps

    @property
    def criteria(self):
        return self._criteria

    @property
    def sampling_control(self):
        return self._sampling_control
