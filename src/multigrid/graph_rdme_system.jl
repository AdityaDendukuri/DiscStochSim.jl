# ─── Unified channel-based RDME generator ────────────────────────────────────
#
# One generator construction for every spatial example — the RDME stated as a
# plain CME over the extended species set {(species, voxel)}.  It returns a
# `DiscreteStochasticSystem{CartesianIndex{S*K}}` (the same channel-list object
# used by the well-mixed CME solvers), so the unified SA-AMG stepper drives it
# unchanged.  No closure, no Poisson, no mean-field: every channel acts on the
# genuine joint distribution.
#
# State layout: component (v-1)*S + s of the CartesianIndex holds the count of
# species s in voxel v  (v ∈ 1:K, s ∈ 1:S).
#
# Two channel classes, exactly as in the RDME:
#   • local reaction r in voxel v : shift ν_r placed on voxel v's S components,
#                                    propensity a_r(local counts of voxel v)
#   • diffusion hop of species s on edge (i,j): unimolecular transfer
#                                    X_{s,i} → X_{s,j} at rate D_s · n_{s,i}
#
# A `LocalReaction{S}` is the per-voxel template; its `rate` reads the S local
# counts (an NTuple{S,Int}) and returns a propensity.  `local_rxns(v)` may return
# a voxel-dependent list, which is how localized sources (A injected here, B
# injected there) are expressed without breaking the uniform scheme.

"""
    LocalReaction{S}

A single within-voxel reaction template over `S` species.

- `dn`   : net change applied to the voxel's `S` species, `NTuple{S,Int}`
- `rate` : propensity, `NTuple{S,Int} -> Float64`, of the voxel's local counts
"""
struct LocalReaction{S}
    dn   :: NTuple{S,Int}
    rate :: Function
end

"""
    ring_edges(K) -> Vector{Tuple{Int,Int}}

Periodic 1-D chain (a ring) on `K` voxels: edges (1,2),…,(K-1,K),(K,1).
"""
ring_edges(K::Int) = Tuple{Int,Int}[(k, k % K + 1) for k in 1:K]

"""
    chain_edges(K) -> Vector{Tuple{Int,Int}}

Open 1-D chain on `K` voxels (no wrap-around).
"""
chain_edges(K::Int) = Tuple{Int,Int}[(k, k+1) for k in 1:K-1]

@inline _comp(v::Int, s::Int, S::Int) = (v - 1) * S + s

"""
    build_graph_rdme(K, ::Val{S}, edges, local_rxns, D) -> DiscreteStochasticSystem{CartesianIndex{S*K}}

Assemble the joint RDME generator as a channel list.

Arguments
---------
- `K`          : number of voxels
- `Val(S)`     : species per voxel (compile-time)
- `edges`      : voxel adjacency, list of `(i,j)` with `i<j` (e.g. `ring_edges(K)`)
- `local_rxns` : `v -> Vector{LocalReaction{S}}`, the reactions in voxel `v`
                 (pass a constant function for spatially uniform kinetics)
- `D`          : per-species diffusion coefficients, `NTuple{S,Float64}` or vector
                 (`d_s = D_s`, dx = 1); use `0.0` for an immobile species

Returns a `DiscreteStochasticSystem` whose element type is `CartesianIndex{S*K}`.
"""
function build_graph_rdme(K::Int, ::Val{S}, edges::Vector{Tuple{Int,Int}},
                          local_rxns, D) where {S}
    N = S * K
    stoichs = CartesianIndex{N}[]
    props   = Function[]

    eu(c)   = CartesianIndex(ntuple(i -> i == c ? 1 : 0, N))

    # ── local reactions, per voxel ────────────────────────────────────────────
    for v in 1:K
        for rx in local_rxns(v)
            # shift: place dn on voxel v's S components
            shift = CartesianIndex(ntuple(N) do i
                vv, ss = divrem(i - 1, S); vv += 1; ss += 1
                vv == v ? rx.dn[ss] : 0
            end)
            comps = ntuple(s -> _comp(v, s, S), S)
            let rx = rx, comps = comps
                push!(stoichs, shift)
                push!(props, (x, rv, t) -> begin
                    nv = ntuple(s -> Tuple(x)[comps[s]], S)
                    rx.rate(nv)
                end)
            end
        end
    end

    # ── diffusion hops, per edge per mobile species ───────────────────────────
    for (i, j) in edges
        for s in 1:S
            Ds = Float64(D[s])
            Ds > 0 || continue
            ci = _comp(i, s, S); cj = _comp(j, s, S)
            let ci = ci, cj = cj, Ds = Ds
                push!(stoichs, -eu(ci) + eu(cj))               # i → j
                push!(props,  (x, rv, t) -> Ds * Tuple(x)[ci])
                push!(stoichs,  eu(ci) - eu(cj))               # j → i
                push!(props,  (x, rv, t) -> Ds * Tuple(x)[cj])
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, props)
end

# ─── Per-voxel marginal extraction (for plotting; means come OUT, never in) ───

"""
    voxel_species_means(f::UnifiedFSP, K, S) -> Matrix{Float64}  (K × S)

Exact per-voxel, per-species mean ⟨n_{s,v}⟩ = Σ n·p from the realised joint.
"""
function voxel_species_means(f::UnifiedFSP{N}, K::Int, S::Int) where {N}
    @assert N == S * K
    μ = zeros(K, S)
    for (ci, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(ci)
        @inbounds for v in 1:K, s in 1:S
            μ[v, s] += t[_comp(v, s, S)] * p
        end
    end
    μ
end

"""
    voxel_hi_prob(f::UnifiedFSP, K, S, s, n_thresh) -> Vector{Float64}  (length K)

Per-voxel P(n_{s,v} > n_thresh) — the bimodal/front indicator for species `s`.
"""
function voxel_hi_prob(f::UnifiedFSP{N}, K::Int, S::Int, s::Int, n_thresh::Int) where {N}
    ϕ = zeros(K)
    for (ci, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(ci)
        @inbounds for v in 1:K
            t[_comp(v, s, S)] > n_thresh && (ϕ[v] += p)
        end
    end
    ϕ
end

"""
    voxel_species_cov(f::UnifiedFSP, K, S, s) -> Matrix{Float64}  (K × K)

Exact voxel-voxel covariance of species `s`,
`Cov(n_{s,i}, n_{s,j}) = ⟨n_{s,i} n_{s,j}⟩ − ⟨n_{s,i}⟩⟨n_{s,j}⟩`,
read straight from the realised joint. This is the observable a product /
per-voxel-marginal solver structurally cannot produce: under any product
ansatz the off-diagonal is identically zero, so a non-zero off-diagonal here
is direct evidence that the joint distribution (the genuine multigrid target)
is being represented, not a marginal closure.
"""
function voxel_species_cov(f::UnifiedFSP{N}, K::Int, S::Int, s::Int) where {N}
    @assert N == S * K
    comp = ntuple(v -> _comp(v, s, S), K)
    μ  = zeros(K)
    M2 = zeros(K, K)               # ⟨n_i n_j⟩
    for (ci, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(ci)
        @inbounds for i in 1:K
            ni = t[comp[i]]
            μ[i] += ni * p
            for j in i:K
                M2[i, j] += ni * t[comp[j]] * p
            end
        end
    end
    C = zeros(K, K)
    @inbounds for i in 1:K, j in i:K
        c = M2[i, j] - μ[i] * μ[j]
        C[i, j] = c
        C[j, i] = c
    end
    C
end
