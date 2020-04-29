# Python Imports

# Third-Party Imports

# Project Imports
from frag_pele.frag.Parameters.PeleParameters import PeleParameterArchives
from frag_pele.frag.Parameters.PeleParameters import PeleParameterPaths
from frag_pele.frag.Parameters.PeleParameters import PeleParameterSimulationValues


class PeleParameters:

    def __init__(self, pele_params_path: PeleParameterPaths, pele_params_archives: PeleParameterArchives,
                 pele_params_sim_values: PeleParameterSimulationValues):
        """
        :param pele_params_path: Path related arguments of PELE.
        :param pele_params_archives: Archives related arguments of PELE.
        :param pele_params_sim_values: Simulation values (control_file) related arguments of PELE.
        """
        self._pele_params_path = pele_params_path
        self._pele_params_archives = pele_params_archives
        self._pele_params_sim_values = pele_params_sim_values

    # Properties
    @property
    def pele_params_path(self):
        return self._pele_params_path

    @property
    def pele_params_archives(self):
        return self._pele_params_archives

    @property
    def pele_params_sim_values(self):
        return self._pele_params_sim_values

    # Methods
    def extract_parameters(self):
        return self.pele_params_path, self.pele_params_archives, self.pele_params_sim_values
