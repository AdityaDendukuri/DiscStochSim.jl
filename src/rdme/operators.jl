"""
    restrict(sp_fine) -> StateSpace{CartesianIndex{K÷2}, T}

Restriction operator  𝓡: fine → coarse.

Marginalizes the fine-level probability distribution over all within-coarse-voxel
splits.  For each fine state x = (n₁, …, n_K), the corresponding coarse state is
x̄ = (n₁+n₂, n₃+n₄, …, n_{K-1}+n_K) and its probability is accumulated.

This is exact: 𝓡 preserves total probability mass.
"""
function restrict(sp_fine::StateSpace{CartesianIndex{K}, T}) where {K, T}
    iseven(K) || error("K = $K must be even for restriction")
    K2 = K ÷ 2
    sp_coarse = StateSpace{CartesianIndex{K2}, T}()

    for i in eachindex(sp_fine.states)
        t_fine       = Tuple(sp_fine.states[i])
        coarse_tuple = ntuple(j -> t_fine[2j-1] + t_fine[2j], Val(K2))
        coarse_state = CartesianIndex(coarse_tuple)
        prob         = sp_fine.probs[i]

        idx = get(sp_coarse.index, coarse_state, 0)
        if idx == 0
            add_state!(sp_coarse, coarse_state, prob)
        else
            sp_coarse.probs[idx] += prob
        end
    end

    sp_coarse
end

# ─── prolongation ──────────────────────────────────────────────────────────────

"""
    prolong(sp_coarse, ::Val{K}; weight_tol=1e-14)
        -> StateSpace{CartesianIndex{K}, T}

Prolongation operator  𝓟: coarse → fine.

For each coarse state x̄ = (n̄₁, …, n̄_{K/2}) with probability p̄, distributes
p̄ across all fine states consistent with x̄ using the binomial stationary
distribution of fast symmetric intra-voxel diffusion:

    P(mⱼ | n̄ⱼ) = C(n̄ⱼ, mⱼ) (1/2)^{n̄ⱼ}

Fine states with weight below `weight_tol` are discarded to keep the
active set manageable.

Note: for large per-voxel counts n̄ⱼ the number of fine splits is ∏ⱼ (n̄ⱼ+1);
use `weight_tol` to prune negligible contributions.
"""
function prolong(sp_coarse::StateSpace{CartesianIndex{K2}, T},
                 ::Val{K};
                 weight_tol::Float64 = 1e-14) where {K2, K, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for one-level prolongation")

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_recurse!(sp_fine, coarse_tup, prob, fine_buf, 1, K2, weight_tol)
    end

    sp_fine
end

"""
Recursive enumeration of all fine splits of the coarse count vector.

At each level `pair_idx` we enumerate m ∈ 0..n̄ⱼ and multiply in the
Binomial(n̄ⱼ, 1/2) weight.  Pruning at `weight_tol` avoids exponential blowup
for high-count states.
"""
function _prolong_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                            coarse_tup::NTuple{K2, Int},
                            weight::T,
                            fine_buf::Vector{Int},
                            pair_idx::Int,
                            n_pairs::Int,
                            weight_tol::Float64) where {K, K2, T}
    if pair_idx > n_pairs
        fine_state = CartesianIndex(NTuple{K, Int}(fine_buf))
        idx = get(sp_fine.index, fine_state, 0)
        if idx == 0
            add_state!(sp_fine, fine_state, weight)
        else
            sp_fine.probs[idx] += weight
        end
        return
    end

    nj        = coarse_tup[pair_idx]
    log_half_n = -nj * log(2.0)    # log( (1/2)^nj )

    for m in 0:nj
        binom_weight = weight * exp(log_half_n) * binomial(nj, m)
        binom_weight > weight_tol || continue
        fine_buf[2*pair_idx - 1] = m
        fine_buf[2*pair_idx]     = nj - m
        _prolong_recurse!(sp_fine, coarse_tup, binom_weight, fine_buf,
                          pair_idx + 1, n_pairs, weight_tol)
    end
end

# ─── lumpability diagnostics ──────────────────────────────────────────────────

"""
    lumpability_ratio(model, fine_grid) -> Float64

Return the lumpability parameter  ε = d_inter / d_intra = (D/(2h)²) / (D/h²) = 1/4

for a uniform grid.  Always 1/4 per level for homogeneous diffusion.
This serves as a sanity check; heterogeneous grids would return spatially varying ε.
"""
lumpability_ratio(model::RDMEModel1D, fine_grid::VoxelGrid) =
    diffusion_rate(model.D, coarsen(fine_grid)) / diffusion_rate(model.D, fine_grid)

"""
    per_molecule_flux(sp, model, grid) -> Vector{Float64}

Compute the per-molecule spatial flux across each interior boundary,
φ_{k,k+1} = Φ_{k,k+1} / max(1, 𝔼[xₖ]).

Used in the two-criterion refinement rule (Definition 6.1 of the paper) to detect
spatial bottlenecks: boundaries with high φ (high throughput per molecule) must be
refined even when the absolute flux Φ is small.
"""
function per_molecule_flux(sp::StateSpace{CartesianIndex{K}, T},
                            model::RDMEModel1D,
                            grid::VoxelGrid) where {K, T}
    K - 1 > 0 || return Float64[]
    d = diffusion_rate(model.D, grid)

    # Expected count per voxel
    mean_count = zeros(Float64, K)
    for i in eachindex(sp.states)
        t = Tuple(sp.states[i])
        p = sp.probs[i]
        for k in 1:K
            mean_count[k] += p * t[k]
        end
    end

    # Absolute flux Φ_{k,k+1} = d · 𝔼[xₖ]  (forward only; symmetric)
    # Per-molecule flux φ_{k,k+1} = Φ / max(1, 𝔼[xₖ])
    phi = Vector{Float64}(undef, K-1)
    for k in 1:K-1
        Phi_k  = d * mean_count[k]
        phi[k] = Phi_k / max(1.0, mean_count[k])
    end
    phi
end
