# Python Imports

# Third-Party Imports

# Project Imports


class PlopParameters:

    def __init__(self, plop_path: str, sch_python: str, rotamers: str):
        self._plop_path = plop_path
        self._sch_python = sch_python
        self._rotamers = rotamers

    # Properties
    @property
    def plop_path(self):
        return self._plop_path

    @property
    def sch_python(self):
        return self._sch_python

    @property
    def rotamers(self):
        return self._rotamers
