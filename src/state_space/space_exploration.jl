"""
    expand_forward(x, model, boundary_condition)

Computes the set of connected states from `x` by adding each stoichiometric vector in `model.stoichvecs`.
Only states satisfying `boundary_condition` are included.
"""
function expand_forward(x::Element, model::Model, boundary_condition::Function) where {Element, Model}
    # For each stoichiometric vector S, compute the candidate state and return it if it satisfies the boundary condition.
    states = map(model.stoichvecs) do S
        x_new = x + S
        boundary_condition(x_new) ? x_new : nothing
    end
    # Filter out any `nothing` values and convert the result to a Set.
    return Set(filter(!isnothing, states))
end

"""
    expand_backward(x, model, boundary_condition)

Computes the set of source states from `x` by subtracting each stoichiometric vector in `model.stoichvecs`.
Only states satisfying `boundary_condition` are included.
"""
function expand_backward(x::Element, model::Model, boundary_condition::Function) where {Element, Model}
    states = map(model.stoichvecs) do S
        x_new = x - S
        boundary_condition(x_new) ? x_new : nothing
    end
    return Set(filter(!isnothing, states))
end

"""
    ssa_step(x, model, rates, t)

Performs one step of the stochastic simulation algorithm (SSA) starting from state `x`.
- Computes propensities using each function in `model.propensities`.
- Chooses a reaction index based on the cumulative sum of propensities.
- Returns the new state obtained by adding the corresponding stoichiometric vector.
"""
function ssa_step(x, model, rates, t)
    # Calculate propensities for the current state.
    propensities = map(model.propensities) do α
        α(x, rates, t)
    end
    # Compute the cumulative sum of propensities.
    cumsum_propensities = cumsum(propensities)
    total = cumsum_propensities[end]
    # Choose a random threshold.
    threshold = rand() * total
    # Find the first index where the cumulative sum exceeds the threshold.
    idx = findfirst(p -> threshold < p, cumsum_propensities)
    # Return the new state.
    return x + model.stoichvecs[idx]
end

