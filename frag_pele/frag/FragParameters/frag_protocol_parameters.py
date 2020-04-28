# Python Imports

# Third-Party Imports

# Project Imports


class FragProtocolParameters:

    def __init__(self, only_prepare: bool, only_grow: bool, explorative: bool):
        """
        :param only_prepare = todo description
        :param only_grow = todo description
        :param explorative = todo description
        """
        # Protocol
        self._only_prepare = only_prepare
        self._only_grow = only_grow
        self._explorative = explorative

    # Properties
    @property
    def only_prepare(self):
        return self._only_prepare

    @property
    def only_grow(self):
        return self._only_grow

    @property
    def explorative(self):
        return self._explorative
