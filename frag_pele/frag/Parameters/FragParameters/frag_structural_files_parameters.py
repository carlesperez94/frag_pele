# Python Imports

# Third-Party Imports

# Project Imports


class FragStructuralFilesParameters:

    def __init__(self, complex_pdb: str, fragment_pdb: str):
        """
        :param complex_pdb: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
        used as the core structure. Remember to rename the ligand chain with a different character in
        order to detect it.
        :param fragment_pdb: Path to the PDB file containing the fragment.
        """
        self._complex_pdb = complex_pdb
        self._fragment_pdb = fragment_pdb

    # Properties
    @property
    def complex_pdb(self):
        return self._complex_pdb

    @property
    def fragment_pdb(self):
        return self._fragment_pdb
