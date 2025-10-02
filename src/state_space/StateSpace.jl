using LinearAlgebra


"""
    StateSpace{E,T}()

Hybrid active set with stable global indexing for reconstruction.
- `active_states::Set{E}`: O(1) membership
- `state_vec::Vector{E}`: stable iteration order for current active set
- `state_to_global::Dict{E,Int}` / `global_to_state::Vector{E}`: stable global ids
- `next_global_idx::Int`: next available global id (1-based)
- `probabilities::Vector{T}`: p aligned with `state_vec`
- `global_active::BitVector`: which global ids are active
"""
mutable struct StateSpace{E,T<:Real}
    active_states::Set{E}
    state_vec::Vector{E}

    state_to_global::Dict{E,Int}
    global_to_state::Vector{E}
    next_global_idx::Int

    probabilities::Vector{T}
    global_active::BitVector
end

function StateSpace{E,T}() where {E,T}
    StateSpace{E,T}(
        Set{E}(),
        E[],
        Dict{E,Int}(),
        E[],
        1,
        T[],
        BitVector()
    )
end

"""
    add_state!(space, state, prob)

Add (or update) a state with given probability. Assigns a permanent global id on first insertion.
"""
function add_state!(space::StateSpace{E,T}, state::E, prob::T) where {E,T}
    if state ∉ space.active_states
        push!(space.active_states, state)
        push!(space.state_vec, state)
        push!(space.probabilities, prob)

        if !haskey(space.state_to_global, state)
            g = space.next_global_idx
            space.state_to_global[state] = g
            push!(space.global_to_state, state)
            push!(space.global_active, true)
            space.next_global_idx += 1
        else
            g = space.state_to_global[state]
            space.global_active[g] = true
        end
    else
        # Update existing probability
        idx = findfirst(==(state), space.state_vec)
        @assert idx !== nothing "Internal error: state present but not in state_vec."
        space.probabilities[idx] = prob
    end
    return space
end

"""
    remove_states!(space, states_to_remove)

Remove a batch of states from the active set (keeps their global ids for later reactivation).
"""
function remove_states!(space::StateSpace, states_to_remove)
    keep_mask = [s ∉ states_to_remove for s in space.state_vec]

    # Flip global active flags
    for state in states_to_remove
        if haskey(space.state_to_global, state)
            g = space.state_to_global[state]
            space.global_active[g] = false
        end
    end

    # Remove from active set and compact arrays
    setdiff!(space.active_states, states_to_remove)
    space.state_vec       = space.state_vec[keep_mask]
    space.probabilities   = space.probabilities[keep_mask]
    return nothing
end

"""
    get_global_indices(space) -> Vector{Int}

Global ids corresponding to current `state_vec` order.
"""
get_global_indices(space::StateSpace) = [space.state_to_global[s] for s in space.state_vec]

"""
    expand!(space, model; depth=1)

One- or multi-layer neighborhood expansion using `model.stoichvecs`.
"""
function expand!(space::StateSpace{E,T}, model; depth::Int=1) where {E,T}
    frontier = copy(space.state_vec)
    for _ in 1:depth
        nextf = E[]
        @inbounds for state in frontier
            for ν in model.stoichvecs
                nbr = state + ν
                if nbr ∉ space.active_states
                    add_state!(space, nbr, zero(T))
                    push!(nextf, nbr)
                end
            end
        end
        isempty(nextf) && break
        frontier = nextf
    end
    return space
end

"""
    compress!(space, model, rates, t, prob_quantile; flux_tolerance=1e-6)

Quantile-based candidate pruning with flux-aware protection.
- Candidates: bottom `prob_quantile` fraction by probability.
- Protected: states with `p(x)*Σ_k α_k(x,t) ≥ flux_tolerance * total_flux`.
Returns `(n_removed, total_flux)`.
"""
function compress!(space::StateSpace{E,T}, model, rates, t, prob_quantile; flux_tolerance=1e-6) where {E,T}
    n = length(space.state_vec)
    n_remove = round(Int, prob_quantile * n)
    n_remove == 0 && return (n_removed=0, total_flux=zero(T))

    flux = zeros(T, n)
    @inbounds for i in 1:n
        x = space.state_vec[i]
        # sum of propensities (exit rate)
        exit_rate = sum(prop(x, rates, t) for prop in model.propensities; init=zero(T))
        flux[i] = space.probabilities[i] * exit_rate
    end

    total_flux = sum(flux)
    flux_threshold = total_flux * flux_tolerance

    # bottom by probability
    cand_idx = sortperm(space.probabilities)[1:n_remove]

    to_remove = E[]
    @inbounds for idx in cand_idx
        if flux[idx] < flux_threshold
            push!(to_remove, space.state_vec[idx])
        end
    end

    remove_states!(space, to_remove)
    return (n_removed=length(to_remove), total_flux=total_flux)
end

"""
    renormalize!(space) -> total

Renormalize probabilities in-place. Returns previous total mass.
"""
function renormalize!(space::StateSpace)
    total = sum(space.probabilities)
    if total > 0
        space.probabilities ./= total
    end
    return total
end
export StateSpace, add_state!, remove_states!, expand!, compress!, renormalize!, get_global_indices
