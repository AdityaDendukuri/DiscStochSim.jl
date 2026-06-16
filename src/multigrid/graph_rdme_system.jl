# ─── Per-voxel joint observable extraction ───────────────────────────────────
#
# Read exact per-voxel / per-species observables out of a realised joint
# distribution (a `UnifiedFSP` over `CartesianIndex{S*K}` configurations). The
# joint generator itself is assembled by `build_rdme_joint` (unified_rdme.jl).
#
# State layout: component (v-1)*S + s of the CartesianIndex holds the count of
# species s in voxel v  (v ∈ 1:K, s ∈ 1:S).
#
# Means come OUT of the joint, never in: no closure, no Poisson, no mean-field.

@inline _comp(v::Int, s::Int, S::Int) = (v - 1) * S + s

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
