# Python Imports

# Third-Party Imports

# Project Imports
from frag_pele.frag.FragParameters.frag_configuration_parameters import FragConfigurationParameters
from frag_pele.frag.FragParameters.frag_identification_parameters import FragIdentificationParameters
from frag_pele.frag.FragParameters.frag_protocol_parameters import FragProtocolParameters
from frag_pele.frag.FragParameters.frag_running_modes_parameters import FragRunningModesParameters
from frag_pele.frag.FragParameters.frag_structural_configuration_parameters import FragStructuralConfigurationParameters
from frag_pele.frag.FragParameters.frag_structural_files_parameters import FragStructuralFilesParameters


class FragParameters:

    def __init__(self, structural_files: FragStructuralFilesParameters,
                 structural_config: FragStructuralConfigurationParameters,
                 identification_params: FragIdentificationParameters,
                 configuration_params: FragConfigurationParameters,
                 protocol_params: FragProtocolParameters,
                 running_modes_params: FragRunningModesParameters
                 ):
        """
        :param structural_files: todo definition
        :param structural_config: todo definition
        :param identification_params: todo definition
        :param configuration_params: todo definition
        :param protocol_params: todo definition
        :param running_modes_params: todo definition
        """
        self._structural_files = structural_files
        self._structural_configuration = structural_config
        self._identification_parameters = identification_params
        self._configuration_parameters = configuration_params
        self._protocol_parameters = protocol_params
        self._running_modes = running_modes_params

    # Properties
    @property
    def structural_files(self):
        return self._structural_files

    @property
    def structural_configuration(self):
        return self._structural_configuration

    @property
    def identification_parameters(self):
        return self._identification_parameters

    @property
    def configuration_parameters(self):
        return self._configuration_parameters

    @property
    def protocol_parameters(self):
        return self._protocol_parameters

    @property
    def running_modes(self):
        return self._running_modes
