"""
    StateSpace{E, T<:Real}

Tracks active states, their probabilities, and a stable global ID mapping.

- `states[i]` is the i-th active state.
- `probs[i]`  is its probability.
- `index[x]`  gives the compact position of state `x` in `states`.
- `global_id[x]` gives the permanent ID assigned when `x` was first added.
- `id_to_state[id]` maps permanent ID back to state (append-only).
- `next_id` is the next global ID to assign.
"""
mutable struct StateSpace{E, T<:Real}
    states::Vector{E}
    probs::Vector{T}
    index::Dict{E,Int}
    global_id::Dict{E,Int}
    id_to_state::Vector{E}
    next_id::Int
end

function StateSpace{E,T}() where {E,T<:Real}
    StateSpace{E,T}(E[], T[], Dict{E,Int}(), Dict{E,Int}(), E[], 1)
end

Base.length(sp::StateSpace) = length(sp.states)

"""
    add_state!(space, state, prob=0.0)

Add `state` with probability `prob`. No-op if already present. O(1).
"""
function add_state!(sp::StateSpace{E,T}, state::E, prob::T=zero(T)) where {E,T}
    haskey(sp.index, state) && return sp
    push!(sp.states, state)
    push!(sp.probs, prob)
    sp.index[state] = length(sp.states)
    if !haskey(sp.global_id, state)
        sp.global_id[state] = sp.next_id
        push!(sp.id_to_state, state)
        sp.next_id += 1
    end
    sp
end

"""
    remove_states!(space, mask)

Remove states where `mask[i] == true`. Rebuild compact index.
"""
function remove_states!(sp::StateSpace{E,T}, remove_mask::BitVector) where {E,T}
    keep = .!remove_mask
    sp.states = sp.states[keep]
    sp.probs = sp.probs[keep]
    empty!(sp.index)
    for (i, s) in enumerate(sp.states)
        sp.index[s] = i
    end
    sp
end

"""
    expand!(space, model, bc; depth=1)

Add stoichiometric neighbors of all current states that satisfy `bc`.
Probabilities of new states are initialized to zero.
"""
function expand!(sp::StateSpace{E,T}, model, bc; depth::Int=1) where {E,T}
    for _ in 1:depth
        current = copy(sp.states)
        for x in current
            for ν in model.stoichvecs
                y = x + ν
                if !haskey(sp.index, y) && bc(y)
                    add_state!(sp, y)
                end
            end
        end
    end
    sp
end

"""
    compress!(space, model, rates, t, prob_quantile; flux_tolerance=1e-6)

Flux-aware pruning: select the lowest-probability states (up to `prob_quantile`
cumulative mass) as candidates, then only remove those whose flux is below
`flux_tolerance * total_flux`.

Returns `(total_flux, n_removed)`.
"""
function compress!(sp::StateSpace{E,T}, model, rates, t::Real,
                   prob_quantile::Real; flux_tolerance::Real=1e-6) where {E,T}
    n = length(sp)
    n == 0 && return (zero(T), 0)

    # Compute flux for each state: p(x) * sum_k a_k(x)
    flux = Vector{T}(undef, n)
    @inbounds for i in 1:n
        w = zero(T)
        for prop in model.propensities
            w += prop(sp.states[i], rates, t)
        end
        flux[i] = sp.probs[i] * w
    end
    total_flux = sum(flux)
    flux_threshold = total_flux * flux_tolerance

    # Candidate selection: accumulate lowest-probability states up to prob_quantile
    idx = sortperm(sp.probs)
    running_mass = zero(T)
    candidates = Set{Int}()
    for i in idx
        if running_mass + sp.probs[i] > prob_quantile
            break
        end
        running_mass += sp.probs[i]
        push!(candidates, i)
    end

    # Flux protection: only prune candidates with flux below threshold
    remove_mask = falses(n)
    n_removed = 0
    for i in candidates
        if flux[i] < flux_threshold
            remove_mask[i] = true
            n_removed += 1
        end
    end

    remove_states!(sp, remove_mask)
    (total_flux, n_removed)
end

"""
    renormalize!(space)

Scale probabilities so they sum to 1. Returns the sum before normalization.
"""
function renormalize!(sp::StateSpace)
    s = sum(sp.probs)
    if s > 0
        sp.probs ./= s
    end
    s
end

"""
    get_global_ids(space) -> Vector{Int}

Return the global IDs of the current active states (in compact order).
"""
function get_global_ids(sp::StateSpace)
    [sp.global_id[s] for s in sp.states]
end
