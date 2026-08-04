mutable struct FESSNetwork{S}
    states::Vector{S}
    index::Dict{S,Int}
end

FESSNetwork(::Type{S}) where {S} = FESSNetwork(S[], Dict{S,Int}())
FESSNetwork(state::S) where {S} = FESSNetwork([state], Dict(state => 1))

Base.length(network::FESSNetwork) = length(network.states)
Base.isempty(network::FESSNetwork) = isempty(network.states)

function augment!(network::FESSNetwork{S}, state::S) where {S}
    haskey(network.index, state) && return network.index[state]
    push!(network.states, state)
    network.index[state] = length(network)
    length(network)
end

function _remove!(network::FESSNetwork, remove::Int)
    deleteat!(network.states, remove)
    empty!(network.index)
    for (index, state) in enumerate(network.states)
        network.index[state] = index
    end
    network
end

_valid_fess_state(::Any) = true
_valid_fess_state(state::CartesianIndex) = all(value -> value >= 0, Tuple(state))

function _propensities(state, model, rates, time)
    Float64[max(0.0, propensity(state, rates, time))
            for propensity in model.propensities]
end

function _terminal(state, model, rates, time, atol)
    sum(_propensities(state, model, rates, time)) <= atol
end

"""
    augment!(network, entrances, model, rates; depth, time, atol)

Add graph shells around `entrances`. Invalid and globally terminal states remain
outside the network, so transitions to them are treated as exits.
"""
function augment!(
    network::FESSNetwork,
    entrances,
    model,
    rates;
    depth::Int=1,
    time::Real=0.0,
    atol::Real=1.0e-12,
)
    frontier = collect(entrances)
    foreach(state -> augment!(network, state), frontier)
    visited = Set(frontier)

    for _ in 1:depth
        next_frontier = empty(frontier)
        for state in frontier, change in model.stoichvecs
            destination = state + change
            destination in visited && continue
            push!(visited, destination)
            _valid_fess_state(destination) || continue
            _terminal(destination, model, rates, time, atol) && continue
            augment!(network, destination)
            push!(next_frontier, destination)
        end
        frontier = next_frontier
    end
    network
end

mutable struct FESSWorkspace{S}
    network::FESSNetwork{S}
end

FESSWorkspace(::Type{S}) where {S} = FESSWorkspace(FESSNetwork(S))
FESSWorkspace(state::S) where {S} = FESSWorkspace(FESSNetwork(state))

