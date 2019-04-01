from AdaptivePELE.constants import blockNames


class CLUSTERING_TYPES:
    rmsd, contactMap, lastSnapshot, null, MSMClustering = list(range(5))

CLUSTERING_TYPE_TO_STRING_DICTIONARY = {
    CLUSTERING_TYPES.rmsd: blockNames.ClusteringTypes.rmsd,
    CLUSTERING_TYPES.contactMap: blockNames.ClusteringTypes.contactMap,
    CLUSTERING_TYPES.lastSnapshot: blockNames.ClusteringTypes.lastSnapshot,
    CLUSTERING_TYPES.null: blockNames.ClusteringTypes.null,
    CLUSTERING_TYPES.MSMClustering: blockNames.ClusteringTypes.MSMClustering
}


class SIMILARITY_TYPES:
    differenceDistance, Jaccard, correlation = list(range(3))

SIMILARITY_TYPES_TO_STRING_DICTIONARY = {
    SIMILARITY_TYPES.differenceDistance: blockNames.ClusteringTypes.differenceDistance,
    SIMILARITY_TYPES.Jaccard: blockNames.ClusteringTypes.Jaccard,
    SIMILARITY_TYPES.correlation: blockNames.ClusteringTypes.correlation
}

SIMILARITY_TYPES_NAMES = {blockNames.ClusteringTypes.differenceDistance,
                          blockNames.ClusteringTypes.Jaccard,
                          blockNames.ClusteringTypes.correlation}
