using Catalyst

"""
    DiscreteStochasticSystem

A structure representing a discrete stochastic system defined by its stoichiometric
vectors and propensity functions.

# Fields
- `stoichvecs::Vector{ElementType}`: A vector of stoichiometric change vectors.
- `propensities::Vector{Function}`: A vector of functions computing reaction propensities.
"""
struct DiscreteStochasticSystem{ElementType}
    stoichvecs::Vector{ElementType}
    propensities::Vector{Function}
end

"""
    DiscreteStochasticSystem(reaction_system::ReactionSystem)

Constructs a `DiscreteStochasticSystem` from a given `reaction_system` (a Catalyst `ReactionSystem`).

This constructor performs the following steps:
1. Extracts species `X` and parameters `P` from the reaction system.
2. Computes the stoichiometric change vectors by extracting the net stoichiometry matrix 
   and converting each column into a `CartesianIndex`.
3. Builds the propensity functions for each reaction using Catalyst's `jumpratelaw` and 
   `build_function`.

# Arguments
- `reaction_system::ReactionSystem`: The reaction system defining the discrete stochastic model.

# Returns
A `DiscreteStochasticSystem{CartesianIndex}` instance.
"""
function DiscreteStochasticSystem(reaction_system::ReactionSystem)
    # Extract species and parameters.
    X = Catalyst.get_species(reaction_system)
    P = Catalyst.get_ps(reaction_system)
    
    # Compute stoichiometric vectors: convert each column of the net stoichiometry matrix to a CartesianIndex.
    stoichvecs = map(eachcol(netstoichmat(reaction_system))) do stoichvec
        CartesianIndex(stoichvec...)
    end
    
    # Build propensity functions for each reaction.
    propensities = map(Catalyst.get_rxs(reaction_system)) do reaction
        eqn = jumpratelaw(reaction; combinatoric_ratelaw=true)
        build_function(eqn, X, P, expression=Val{false})
    end
    
    return DiscreteStochasticSystem{CartesianIndex}(stoichvecs, propensities)
end

# Getter functions
stoichvecs(system::DiscreteStochasticSystem) = system.stoichvecs
propensities(system::DiscreteStochasticSystem) = system.propensities

export DiscreteStochasticSystem

