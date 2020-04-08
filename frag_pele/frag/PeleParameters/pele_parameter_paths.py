# Python Imports

# Third-Party Imports

# Project Imports


class PeleParameterPaths:

    def __init__(self, pele_dir: str, pele_license: str, data: str, documents: str):
        """
        :param pele_dir: Absolute path to PELE executables.
        :param pele_license: Absolute path to PELE's license file.
        :param data: Path to PELE Data folder.
        :param documents: Path to PELE Documents folder.
        """
        self._pele_dir = pele_dir
        self._pele_license = pele_license
        self._data = data
        self._documents = documents

    # Properties
    @property
    def pele_dir(self):
        return self._pele_dir

    @property
    def pele_license(self):
        return self._pele_license

    @property
    def data(self):
        return self._data

    @property
    def documents(self):
        return self._documents
