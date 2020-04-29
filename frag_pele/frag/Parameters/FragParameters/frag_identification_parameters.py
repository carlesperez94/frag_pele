# Python Imports

# Third-Party Imports

# Project Imports


class FragIdentificationParameters:

    def __init__(self, ligand_id: str):
        """
        :param ligand_id: Name to identify the new ligand in certain folder.
        """
        # Identification Params
        self._ligand_id = ligand_id

    @property
    def ligand_id(self):
        return self._ligand_id
