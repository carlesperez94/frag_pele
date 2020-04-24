# Python Imports

# Third-Party Imports

# Project Imports


class FragRunningModes:

    def __init__(self, restart: bool, no_check: bool):
        """
        :param restart: If set FrAG will find your last iteration completed and will restart from this point.
        :param no_check = todo description
        """
        self._restart = restart
        self._no_check = no_check

    # Properties
    @property
    def restart(self):
        return self._restart

    @property
    def no_check(self):
        return self._no_check
