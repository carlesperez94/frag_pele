# Python Imports

# Third-Party Imports

# Project Imports


class FragParameters:

    def __init__(self):
        pass
    """
    pele_parameters: PeleParameters, complex_pdb, fragment_pdb, core_atom, fragment_atom, growing_steps, criteria,
    plop_path, sch_python, pdbout, distance_contact, clusterThreshold,
    epsilon, condition, metricweights, nclusters, restart, ID,
    h_core = None, h_frag = None, c_chain = "L", f_chain = "L", rotamers = "30.0",
    banned = None, limit = None, mae = False, rename = False, threshold_clash = 1.7, explorative = False,
    sampling_control = None, 
    
    only_prepare = False, 
    only_grow = False,
    no_check = False
    """

    """
    :param complex_pdb: Path to the PDB file which must contain a protein-ligand complex. Its ligand will be
    used as the core structure. Remember to rename the ligand chain with a different character in
    order to detect it.
    :type complex_pdb: str
    :param fragment_pdb: Path to the PDB file containing the fragment.
    :type fragment_pdb: str
    :param core_atom: PDB-atom-name of the atom of the core that will be used as starting point to grow the fragment.
    :type core_atom: str (max length: 4)
    :param fragment_atom: PDB-atom-name of the atom of the fragment that will bond the core.
    :type fragment_atom: str (max length: 4)
    :param growing_steps: Number of Growing Steps (GS).
    :type growing_steps: int
    :param criteria: Name of the column of the report file to select the structures that will spawn in the
    next GS. Additionally, this parameter will be the selection criteria to extract the best
    structure after completing the growing.
    :type criteria: str
    :param plop_path: Absolute path to PlopRotTemp.py.
    :type plop_path: str
    :param sch_python: Absolute path to Schrodinger's python.
    :type sch_python: str
    :param pdbout: Prefix name of the output PDB extracted from FrAG in each iteration.
    :type pdbout: str
    :param distance_contact: This distance will be used to determine which amino acids are in contact with the ligand.
    Then, this information will be useful to generate different clusters of structures to initialize the next iteration.
    :type distance_contact: float
    :param clusterThreshold: Threshold that will be used in the clustering step.
    :type clusterThreshold: float (from 0 to 1)
    :param epsilon: Epsilon parameter used in the clustering.
    :type epsilon: float (from 0 to 1)
    :param condition: Condition to select best values to cluster: minimum or maximum.
    :type condition: str
    :param metricweights: Selects how to distribute the weights of the cluster according to its metric,
    two options: linear (proportional to metric) or Boltzmann weigths (proportional to exp(-metric/T).
    Needs to define the temperature T.
    :type metricweights: str
    :param nclusters: Number of cluster that we would like to generate.
    :type nclusters: int
    :param restart: If set FrAG will find your last iteration completed and will restart from this point.
    :type restart: bool
    :param ID: Name to identify the new ligand in certain folder.
    :type ID: str
    
    :param h_core: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the core that will be replaced for
     the fragment, being the initial point of the growing.
    :type h_core: str (max length: 4)
    :param h_frag: PDB-atom-name of the hydrogen bonded to the selected heavy atom of the fragment that will removed to
    bond the fragment with the core.
    :type h_frag: str (max length: 4)
    :param c_chain: Chain name for the ligand in the complex_pdb.
    :type c_chain: str (max length: 1)
    :param f_chain: Chain name for the ligand in the fragment_pdb.
    :type f_chain: str (max length: 1)

    :param rotamers: Degrees of rotation used to select the rotamer's library. Higher value, faster but less accuracy.
    :type rotamers: str
    :param banned: If set, list of tuples with the atom names of the angles that will not be modified.
    :type banned: list
    :param limit: Limit angle in degrees of the banned angles. Any structure with a higher value will be discarted.
    :type limit: int
    :param mae: If set, the output structures will be saved in "mae." format.
    :type mae: bool
    :param rename: If set, the pdb atom names of the fragment will be renamed with "G+atom_number".
    :type rename: bool
    :param threshold_clash: Distance that will be used to detect contacts between atoms of the fragment and atoms
    of the core.
    :type threshold_clash: float

    :param explorative: Set all parameters to perform a more explorative simulation after the growing, in the sampling
    simulation.
    :type explorative: bool

    :param sampling_control: templatized control file to be used in the sampling simulation.
    :type sampling_control: str
    
    only_prepare = False, 
    only_grow = False,
    no_check = False
    :return:
    """
