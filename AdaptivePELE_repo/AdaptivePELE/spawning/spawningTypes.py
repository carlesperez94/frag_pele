from AdaptivePELE.constants import blockNames


class SPAWNING_TYPES:
    sameWeight, inverselyProportional, epsilon, simulatedAnnealing, FAST, variableEpsilon, UCB, independent, REAP, null, independentMetric, ProbabilityMSMCalculator, MetastabilityMSMCalculator, UncertaintyMSMCalculator, IndependentMSMCalculator = list(range(15))

SPAWNING_TYPE_TO_STRING_DICTIONARY = {
    SPAWNING_TYPES.independent: blockNames.StringSpawningTypes.independent,
    SPAWNING_TYPES.sameWeight: blockNames.StringSpawningTypes.sameWeight,
    SPAWNING_TYPES.inverselyProportional: blockNames.StringSpawningTypes.inverselyProportional,
    SPAWNING_TYPES.epsilon: blockNames.StringSpawningTypes.epsilon,
    SPAWNING_TYPES.FAST: blockNames.StringSpawningTypes.fast,
    SPAWNING_TYPES.variableEpsilon: blockNames.StringSpawningTypes.variableEpsilon,
    SPAWNING_TYPES.simulatedAnnealing: blockNames.StringSpawningTypes.simulatedAnnealing,
    SPAWNING_TYPES.UCB: blockNames.StringSpawningTypes.UCB,
    SPAWNING_TYPES.REAP: blockNames.StringSpawningTypes.REAP,
    SPAWNING_TYPES.null: blockNames.StringSpawningTypes.null,
    SPAWNING_TYPES.independentMetric: blockNames.StringSpawningTypes.independentMetric,
    SPAWNING_TYPES.ProbabilityMSMCalculator: blockNames.StringSpawningTypes.ProbabilityMSMCalculator,
    SPAWNING_TYPES.MetastabilityMSMCalculator: blockNames.StringSpawningTypes.MetastabilityMSMCalculator,
    SPAWNING_TYPES.UncertaintyMSMCalculator: blockNames.StringSpawningTypes.UncertaintyMSMCalculator,
    SPAWNING_TYPES.IndependentMSMCalculator: blockNames.StringSpawningTypes.IndependentMSMCalculator

}

MSMSpawning = set([blockNames.StringSpawningTypes.ProbabilityMSMCalculator, blockNames.StringSpawningTypes.MetastabilityMSMCalculator, blockNames.StringSpawningTypes.UncertaintyMSMCalculator, blockNames.StringSpawningTypes.IndependentMSMCalculator])


class EPSILON_VARIATION_TYPES:
    linearVariation, contactsVariation = list(range(2))

EPSILON_VARIATION_TYPE_TO_STRING_DICTIONARY = {
    EPSILON_VARIATION_TYPES.linearVariation: blockNames.VariableEpsilonTypes.linearVariation,
    EPSILON_VARIATION_TYPES.contactsVariation: blockNames.VariableEpsilonTypes.contactsVariation,
}

SPAWNING_NO_DEGENERACY_TYPES = set([SPAWNING_TYPES.null, SPAWNING_TYPES.independent, SPAWNING_TYPES.independentMetric, SPAWNING_TYPES.IndependentMSMCalculator])
