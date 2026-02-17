"""
    FSPSolution

Stores snapshots from an adaptive FSP solve.

# Fields
- `t::Vector{Float64}`: Time points of saved snapshots.
- `snapshots::Vector{Tuple{Vector{E}, Vector{Float64}}}`: `(states, probs)` at each snapshot.
- `state_space_sizes::Vector{Int}`: State space size at each snapshot.
- `n_species::Int`: Number of species (dimension of CartesianIndex).
"""
struct FSPSolution{E}
    t::Vector{Float64}
    snapshots::Vector{Tuple{Vector{E}, Vector{Float64}}}
    state_space_sizes::Vector{Int}
    n_species::Int
end

Base.length(sol::FSPSolution) = length(sol.t)

"""
    sol[i] -> (states, probs)

Return the i-th snapshot.
"""
Base.getindex(sol::FSPSolution, i::Int) = sol.snapshots[i]

"""
    mean_trajectory(sol) -> Matrix{Float64}

Compute the mean state at each snapshot. Returns an `n_snapshots × n_species` matrix.
"""
function mean_trajectory(sol::FSPSolution)
    ns = sol.n_species
    nt = length(sol)
    result = zeros(nt, ns)
    for i in 1:nt
        states, probs = sol.snapshots[i]
        for (s, p) in zip(states, probs)
            for d in 1:ns
                result[i, d] += s[d] * p
            end
        end
    end
    result
end

"""
    marginal(sol, i, species) -> (values, probs)

Compute the marginal distribution of `species` (1-indexed) at snapshot `i`.
Returns sorted `(values, probs)` vectors.
"""
function marginal(sol::FSPSolution, i::Int, species::Int)
    states, probs = sol.snapshots[i]
    dist = Dict{Int, Float64}()
    for (s, p) in zip(states, probs)
        v = s[species]
        dist[v] = get(dist, v, 0.0) + p
    end
    vals = sort(collect(keys(dist)))
    ps = [dist[v] for v in vals]
    (values=vals, probs=ps)
end
