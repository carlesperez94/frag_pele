from AdaptivePELE.constants import blockNames


class SPAWNING_TYPES:
    sameWeight, inverselyProportional, epsilon, simulatedAnnealing, FAST, variableEpsilon, UCB, independent = range(8)

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    SPAWNING_TYPES.independent: blockNames.StringSpawningTypes.independent,
    SPAWNING_TYPES.sameWeight: blockNames.StringSpawningTypes.sameWeight,
    SPAWNING_TYPES.inverselyProportional: blockNames.StringSpawningTypes.inverselyProportional,
    SPAWNING_TYPES.epsilon: blockNames.StringSpawningTypes.epsilon,
    SPAWNING_TYPES.FAST: blockNames.StringSpawningTypes.fast,
    SPAWNING_TYPES.variableEpsilon: blockNames.StringSpawningTypes.variableEpsilon,
    SPAWNING_TYPES.simulatedAnnealing: blockNames.StringSpawningTypes.simulatedAnnealing,
    SPAWNING_TYPES.UCB: blockNames.StringSpawningTypes.UCB
}


class EPSILON_VARIATION_TYPES:
    linearVariation, contactsVariation = range(2)

EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY = {
    EPSILON_VARIATION_TYPES.linearVariation: blockNames.VariableEpsilonTypes.linearVariation,
    EPSILON_VARIATION_TYPES.contactsVariation: blockNames.VariableEpsilonTypes.contactsVariation,
}
