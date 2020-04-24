# Python Imports

# Third-Party Imports

# Project Imports


class FragStructuralConfiguration:

    def __init__(self, core_atom: str, fragment_atom: str, h_core: str, h_frag: str, c_chain: str, f_chain: str,
                 threshold_clash: float):
        """
        :param core_atom: PDB-atom-name of the atom of the core that will be used as starting point to grow the fragment
        :param fragment_atom: PDB-atom-name of the atom of the fragment that will bond the core.
        :param h_core: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the core that will be replaced
        for the fragment, being the initial point of the growing.
        :param h_frag: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the fragment that will removed
        to bond the fragment with the core.
        :param c_chain: Chain name for the ligand in the complex_pdb.
        :param f_chain: Chain name for the ligand in the fragment_pdb.
        :param threshold_clash: Distance that will be used to detect contacts between atoms of the fragment and atoms
        of the core.
        """
        # Structural configuration
        self._core_atom = core_atom
        self._fragment_atom = fragment_atom
        self._h_core = h_core
        self._h_frag = h_frag
        self._c_chain = c_chain
        self._f_chain = f_chain
        self._threshold_clash = threshold_clash

    # Properties
    @property
    def core_atom(self):
        return self._core_atom

    @property
    def fragment_atom(self):
        return self._fragment_atom

    @property
    def h_core(self):
        return self._h_core

    @property
    def h_frag(self):
        return self._h_frag

    @property
    def c_chain(self):
        return self._c_chain

    @property
    def f_chain(self):
        return self._f_chain

    @property
    def threshold_clash(self):
        return self._threshold_clash
