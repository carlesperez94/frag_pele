

class PathHandler:

    def __init__(self, complex_pdb: str, fragment_pdb: str, plop_path: str, sch_python: str, pele_dir: str,
                 pele_license: str):
        """
        :param complex_pdb: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
        used as the core structure. Remember to rename the ligand chain with a different character in order to detect it
        :param fragment_pdb: Path to the PDB file containing the fragment.
        :param plop_path: Absolute path to PlopRotTemp.py.
        :param sch_python: Absolute path to Schrodinger's python.
        :param pele_dir: Absolute path to PELE executables.
        :param pele_license: Absolute path to PELE's license file.
        """
        self._complex_pdb = complex_pdb
        self._fragment_pdb = fragment_pdb
        self._plop_path = plop_path
        self._sch_python = sch_python
        self._pele_dir = pele_dir
        self._pele_license = pele_license

    # Properties
    @propertyc
    def complex_pdb(self):
        return self._complex_pdb

    @property
    def fragment_pdb(self):
        return self._fragment_pdb

    @property
    def plop_path(self):
        return self._plop_path

    @property
    def sch_python(self):
        return self._sch_python

    @property
    def pele_dir(self):
        return self._pele_dir

    @property
    def pele_license(self):
        return self._pele_license
