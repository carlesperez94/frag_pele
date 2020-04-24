# Python Imports

# Third-Party Imports

# Project Imports


class FragProtocolParameters:

    def __init__(self, only_prepare: bool, only_grow: bool):
        """
        :param only_prepare = todo description
        :param only_grow = todo description
        """
        # Protocol
        self._only_prepare = only_prepare
        self._only_grow = only_grow

    # Properties
    @property
    def only_prepare(self):
        return self._only_prepare

    @property
    def only_grow(self):
        return self._only_grow
