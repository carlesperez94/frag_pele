# Python Imports
from typing import List

# Third-Party Imports

# Project Imports


class ClusterParameters:

    def __init__(self, distance_contact: float, cluster_threshold: float, epsilon: float, condition: str,
                 metricweights: str, nclusters: int, pdbout: str, banned: List, limit: int):
        """
        :param distance_contact: This distance will be used to determine which amino acids are in contact with the ligand.
        Then, this information will be useful to generate different clusters of structures to initialize the next iteration.
        :param cluster_threshold: Threshold that will be used in the clustering step.
        :param epsilon: Epsilon parameter used in the clustering.
        :param condition: Condition to select best values to cluster: minimum or maximum.
        :param metricweights: Selects how to distribute the weights of the cluster according to its metric,
        two options: linear (proportional to metric) or Boltzmann weigths (proportional to exp(-metric/T).
        Needs to define the temperature T.
        :param nclusters: Number of cluster that we would like to generate.
        :param pdbout: Prefix name of the output PDB extracted from FrAG in each iteration.
        :param banned: If set, list of tuples with the atom names of the angles that will not be modified.
        :param limit: Limit angle in degrees of the banned angles. Any structure with a higher value will be discarted.
        """
        self._distance_contact = distance_contact
        self._cluster_threshold = cluster_threshold
        self._epsilon = epsilon
        self._condition = condition
        self._metricweights = metricweights
        self._nclusters = nclusters
        self._pdbout = pdbout
        self._banned = banned
        self._limit = limit

    # Properties
    @property
    def distance_contact(self):
        return self._distance_contact

    @property
    def cluster_threshold(self):
        return self._cluster_threshold

    @property
    def epsilon(self):
        return self._epsilon

    @property
    def condition(self):
        return self._condition

    @property
    def metric_weights(self):
        return self._metricweights

    @property
    def number_clusters(self):
        return self._nclusters

    @property
    def pdnout(self):
        return self._pdbout

    @property
    def banned_list(self):
        return self._banned

    @property
    def limit(self):
        return self._limit
