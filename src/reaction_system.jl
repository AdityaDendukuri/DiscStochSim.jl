"""
    DiscreteStochasticSystem{E}

A discrete stochastic system with stoichiometric vectors and propensity functions,
constructed from a Catalyst `ReactionSystem`.

# Fields
- `stoichvecs::Vector{E}`: Net stoichiometric change vectors (one per reaction).
- `propensities::Vector{Function}`: Propensity functions `(state, rates, t) -> α`.
"""
struct DiscreteStochasticSystem{E}
    stoichvecs::Vector{E}
    propensities::Vector{Function}
end

"""
    DiscreteStochasticSystem(rn::ReactionSystem)

Construct a `DiscreteStochasticSystem` from a Catalyst `ReactionSystem`.
Stoichiometric vectors are stored as `CartesianIndex`, and propensity functions
are built via `Catalyst.jumpratelaw` + `build_function`.
"""
function DiscreteStochasticSystem(rn::ReactionSystem)
    X = Catalyst.get_species(rn)
    P = Catalyst.get_ps(rn)

    stoichvecs = map(eachcol(netstoichmat(rn))) do col
        CartesianIndex(col...)
    end

    propensities = map(Catalyst.get_rxs(rn)) do rx
        eqn = jumpratelaw(rx; combinatoric_ratelaw=true)
        build_function(eqn, X, P, expression=Val{false})
    end

    DiscreteStochasticSystem{CartesianIndex}(stoichvecs, propensities)
end
