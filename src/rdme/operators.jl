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
    prolong(sp_coarse, ::Val{K}; weight_tol=1e-14, binom_tol=0.0)
        -> StateSpace{CartesianIndex{K}, T}

Prolongation operator 𝓟: coarse → fine.

Two thresholds control pruning:
- `weight_tol`: skip entire coarse states with probability < weight_tol.
- `binom_tol`: at each voxel-pair, skip Binomial(n̄ⱼ, 1/2) splits where the
  *per-level* Binomial fraction (1/2)^{n̄ⱼ} C(n̄ⱼ, mⱼ) < binom_tol.

Using `binom_tol` (per-pair relative threshold) instead of relying solely on
the accumulated absolute `weight_tol` avoids the K-level compounding problem:
with K/2 pairs, the accumulated weight shrinks by a factor of ~0.375^{K/2} per
centre split, so an absolute threshold silently truncates the recursion mid-way
for any coarse state with small probability.  `binom_tol` ≈ 1e-4 keeps all
Binomial splits contributing at least 1e-4 of the per-pair probability, and
delegates absolute probability cuts to `prune_threshold!` in the solver.
"""
function prolong(sp_coarse::StateSpace{CartesianIndex{K2}, T},
                 ::Val{K};
                 weight_tol::Float64 = 1e-14,
                 binom_tol::Float64  = 0.0) where {K2, K, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for one-level prolongation")

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_recurse!(sp_fine, coarse_tup, prob, fine_buf, 1, K2,
                          weight_tol, binom_tol)
    end

    sp_fine
end

"""
Recursive enumeration of all fine splits of the coarse count vector.

Two pruning thresholds:
- `weight_tol`: prune if accumulated weight (coarse_prob × ∏ Binom fractions) < weight_tol.
- `binom_tol`:  prune if the per-pair Binomial fraction < binom_tol (prevents
  mid-recursion truncation for small-probability coarse states with large K).
"""
function _prolong_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                            coarse_tup::NTuple{K2, Int},
                            weight::T,
                            fine_buf::Vector{Int},
                            pair_idx::Int,
                            n_pairs::Int,
                            weight_tol::Float64,
                            binom_tol::Float64) where {K, K2, T}
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

    nj         = coarse_tup[pair_idx]
    log_half_n = -nj * log(2.0)    # log( (1/2)^nj )

    for m in 0:nj
        # Use log-space to avoid Int64 overflow for large nj
        log_bc     = (sum(log(i) for i in (nj-m+1):nj; init=0.0) -
                      sum(log(i) for i in 1:m;          init=0.0))
        binom_frac = exp(log_half_n + log_bc)
        binom_frac   > binom_tol  || continue               # relative per-pair check
        binom_weight = weight * binom_frac
        binom_weight > weight_tol || continue               # absolute accumulated check
        fine_buf[2*pair_idx - 1] = m
        fine_buf[2*pair_idx]     = nj - m
        _prolong_recurse!(sp_fine, coarse_tup, binom_weight, fine_buf,
                          pair_idx + 1, n_pairs, weight_tol, binom_tol)
    end
end

"""
    prolong_multiplicative(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Multiplicative prolongation: applies a pointwise correction ratio to the pre-smoothed
fine distribution.

    p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)

where:
- p̃^h      = `sp_fine`        (pre-smoothed fine distribution)
- p^{2h}   = `sp_coarse_pre`  (restricted pre-smoothed: 𝓡 p̃^h)
- p̄^{2h}  = `sp_coarse_post` (after coarse solve)

Unlike interpolatory prolongation (Binomial-π), this preserves within-pair asymmetry:
if p̃^h is concentrated at (n_low, n_high), the ratio rescales it without mixing in
the symmetric mode (n_high, n_low) — avoiding the symmetry trap.

Fine states whose coarse denominator is below `prob_tol` are dropped.
"""
function prolong_multiplicative(
    sp_fine::StateSpace{CartesianIndex{K}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{K2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{K2}, T};
    prob_tol::Float64 = 1e-300
) where {K, K2, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for injection prolongation")

    sp_fine_new = StateSpace{CartesianIndex{K}, T}()

    for i in eachindex(sp_fine.states)
        x_fine = sp_fine.states[i]
        p_fine = sp_fine.probs[i]
        p_fine > 0.0 || continue

        # Compute coarse state: x̄_j = x_{2j-1} + x_{2j}
        t_fine   = Tuple(x_fine)
        x_coarse = CartesianIndex(ntuple(j -> t_fine[2j-1] + t_fine[2j], Val(K2)))

        # Look up p^{2h}(x̄) and p̄^{2h}(x̄)
        idx_pre  = get(sp_coarse_pre.index,  x_coarse, 0)
        idx_post = get(sp_coarse_post.index, x_coarse, 0)

        p_pre  = idx_pre  > 0 ? sp_coarse_pre.probs[idx_pre]   : 0.0
        p_post = idx_post > 0 ? sp_coarse_post.probs[idx_post] : 0.0

        # Multiplicative correction; skip if denominator is negligible
        p_new = p_pre > prob_tol ? p_fine * (p_post / p_pre) : 0.0
        p_new > 0.0 && add_state!(sp_fine_new, x_fine, p_new)
    end

    sp_fine_new
end

"""
    prolong_conditional(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Multiplicative prolongation with conditional extension to newly-expanded coarse states.
Uses the pre-smoothed fine conditional p̃^h(x|x̄) as the prolongation kernel for both
covered and uncovered coarse states.

**Covered coarse states** (those in `sp_coarse_pre`):
  Pure injection — p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄).
  Exact and asymmetry-preserving.

**Uncovered coarse states** (added by expand!, not in `sp_coarse_pre`):
  For each new x̄_new, find the nearest covered neighbor x̄_cov (L1 distance 1
  in coarse-voxel coordinates, i.e., differ by ±1 in exactly one pair).
  Carry over the fine conditional from x̄_cov, shifting the second voxel of the
  changed pair to be consistent with nc_new:

      p^h(n₁, nc_new−n₁, ...) ∝ p̃^h(n₁, nc_cov−n₁, ...) / p^{2h}(x̄_cov)

  This avoids the symmetry trap: for a front at (n_low, n_high) with nc_cov=193,
  nc_new=192 gets mass at (n_low, 192−n_low) = (n_low, n_high−1) — NOT the
  symmetric (n_high−1, n_low) that static π would produce.

For uncovered states with no covered neighbor within L1=1, the state is skipped.
"""
function prolong_conditional(
    sp_fine::StateSpace{CartesianIndex{K}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{K2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{K2}, T};
    prob_tol::Float64 = 1e-300
) where {K, K2, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for conditional prolongation")

    # ── Step 1: multiplicative correction for covered coarse states ──────────
    sp_fine_new = prolong_multiplicative(sp_fine, sp_coarse_pre, sp_coarse_post;
                                         prob_tol = prob_tol)

    # ── Step 2: precompute (covered coarse index → fine states + conditional) ─
    # cov_to_fine[i] = [(fine_state, conditional_prob), ...]
    # where conditional_prob = p̃^h(x) / p^{2h}(x̄_cov)
    cov_to_fine = Dict{Int, Vector{Tuple{CartesianIndex{K}, Float64}}}()
    cov_totals  = zeros(Float64, length(sp_coarse_pre))

    for k in eachindex(sp_fine.states)
        x_fine   = sp_fine.states[k]
        p_fine   = sp_fine.probs[k]
        t_fine   = Tuple(x_fine)
        x_coarse = CartesianIndex(ntuple(j -> t_fine[2j-1]+t_fine[2j], Val(K2)))
        idx = get(sp_coarse_pre.index, x_coarse, 0)
        idx == 0 && continue
        cov_totals[idx] += p_fine
        list = get!(cov_to_fine, idx, Tuple{CartesianIndex{K}, Float64}[])
        push!(list, (x_fine, p_fine))
    end
    # Normalize to get conditionals
    for (idx, list) in cov_to_fine
        total = cov_totals[idx]
        total > prob_tol || continue
        for i in eachindex(list)
            (s, p) = list[i]; list[i] = (s, p / total)
        end
    end

    # ── Step 3: conditional transfer for uncovered (expanded) coarse states ───
    n_pre = length(sp_coarse_pre)

    for i in (n_pre+1):length(sp_coarse_post.states)
        x_coarse_new = sp_coarse_post.states[i]
        p_coarse_new = sp_coarse_post.probs[i]
        p_coarse_new > prob_tol || continue
        t_new = Tuple(x_coarse_new)

        # Find covered neighbor at L1 distance 1 (differs by ±1 in one pair)
        cov_idx   = 0
        cov_pair  = 0   # which pair index j (1-based) differs
        for j in 1:K2
            for delta in (1, -1)
                t_try = ntuple(k -> k == j ? t_new[k] - delta : t_new[k], Val(K2))
                all(v -> v ≥ 0, t_try) || continue
                x_try = CartesianIndex(t_try)
                idx   = get(sp_coarse_pre.index, x_try, 0)
                if idx != 0
                    cov_idx  = idx
                    cov_pair = j
                    break
                end
            end
            cov_idx != 0 && break
        end

        (cov_idx != 0 && haskey(cov_to_fine, cov_idx)) || continue

        nc_j_new = t_new[cov_pair]   # new count for pair cov_pair

        for (x_fine, cond_prob) in cov_to_fine[cov_idx]
            cond_prob > prob_tol || continue
            t_fine = Tuple(x_fine)

            # Shift second voxel of pair cov_pair to match nc_j_new;
            # first voxel of the pair (index 2*cov_pair-1) is kept the same.
            n1_j     = t_fine[2*cov_pair - 1]
            n2_j_new = nc_j_new - n1_j
            n2_j_new ≥ 0 || continue   # first voxel exceeds new nc — skip

            t_shifted = ntuple(k -> begin
                if k == 2*cov_pair - 1
                    n1_j
                elseif k == 2*cov_pair
                    n2_j_new
                else
                    t_fine[k]
                end
            end, Val(K))

            x_fine_new = CartesianIndex(t_shifted)
            p_add      = cond_prob * p_coarse_new
            p_add > prob_tol || continue

            idx2 = get(sp_fine_new.index, x_fine_new, 0)
            if idx2 == 0
                add_state!(sp_fine_new, x_fine_new, p_add)
            else
                sp_fine_new.probs[idx2] += p_add
            end
        end
    end

    sp_fine_new
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

"""
    pair_admissibility(sp, model, grid) -> (α, φ_norm)

For each adjacent pair of fine voxels {2j-1, 2j}, compute two diagnostics:

  α_j     = |E[n_{2j-1}] - E[n_{2j}]| / (E[n_{2j-1}] + E[n_{2j}])
            Within-pair asymmetry.  α ≈ 0 → symmetric → Binomial is a
            plausible closure.  α ≈ 1 → one voxel dominates → Binomial wrong.
            This is a first-moment heuristic, not a Binomial validity proof.

  φ_norm_j = φ_j / mean(φ)   (per-molecule flux, normalised)
            Structural-bottleneck indicator.  High φ_norm with low mean count
            signals a boundary that mediates critical transitions even with
            few molecules present (background inactivity vs bottleneck).

Returns two vectors of length K/2.
"""
function pair_admissibility(sp::StateSpace{CartesianIndex{K}, T},
                             model::RDMEModel1D,
                             grid::VoxelGrid) where {K, T}
    @assert iseven(K)

    μ = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p * t[k]; end
    end

    d = diffusion_rate(model.D, grid)
    n_pairs = K ÷ 2

    α      = Vector{Float64}(undef, n_pairs)
    φ_norm = Vector{Float64}(undef, n_pairs)

    φ_raw = Vector{Float64}(undef, K - 1)
    for k in 1:K-1
        φ_raw[k] = d * μ[k] / max(1.0, μ[k])   # ≈ d when μ[k] >> 1
    end
    φ_mean = mean(φ_raw)

    for j in 1:n_pairs
        k1, k2   = 2j-1, 2j
        total    = μ[k1] + μ[k2]
        α[j]     = total > 0.0 ? abs(μ[k1] - μ[k2]) / total : 0.0
        # per-molecule flux at the intra-pair boundary k1→k2
        φ_norm[j] = φ_mean > 0.0 ? φ_raw[k1] / φ_mean : 0.0
    end

    α, φ_norm
end

# ─── partial restriction (adaptive coarsening) ───────────────────────────────

"""
    partial_restrict(sp_fine, mask) -> StateSpace{CartesianIndex{K_mixed}, T}

Partially restrict a fine state space using a pair-level coarsening mask.

For a K-voxel state space with K/2 pairs:
  - mask[j] = true  → pair {2j-1, 2j} is coarsened: nc_j = n_{2j-1} + n_{2j}
  - mask[j] = false → pair kept fine: both n_{2j-1} and n_{2j} pass through

Output dimension K_mixed = K - sum(mask)  (each coarsened pair loses one dim).
Probabilities are accumulated (exact marginalisation over coarsened splits).
"""
function partial_restrict(sp_fine::StateSpace{CartesianIndex{K}, T},
                           mask::BitVector) where {K, T}
    iseven(K) || error("K=$K must be even")
    n_pairs = K ÷ 2
    length(mask) == n_pairs || error("mask length $(length(mask)) ≠ K/2 = $n_pairs")
    KM = K - sum(mask)
    _partial_restrict_val(sp_fine, mask, Val(KM))
end

function _mixed_tuple(t::NTuple{K,Int}, mask::BitVector, ::Val{KM}) where {K, KM}
    n_pairs = K ÷ 2
    out = Vector{Int}(undef, KM)
    pos = 1
    for j in 1:n_pairs
        k1, k2 = 2j-1, 2j
        if mask[j]
            out[pos] = t[k1] + t[k2];  pos += 1
        else
            out[pos] = t[k1];  out[pos+1] = t[k2];  pos += 2
        end
    end
    NTuple{KM,Int}(out)
end

function _partial_restrict_val(sp_fine::StateSpace{CartesianIndex{K}, T},
                                mask::BitVector, ::Val{KM}) where {K, T, KM}
    sp_mixed = StateSpace{CartesianIndex{KM}, T}()
    for i in eachindex(sp_fine.states)
        t     = Tuple(sp_fine.states[i])
        p     = sp_fine.probs[i]
        mt    = _mixed_tuple(t, mask, Val(KM))
        ms    = CartesianIndex(mt)
        idx   = get(sp_mixed.index, ms, 0)
        if idx == 0
            add_state!(sp_mixed, ms, p)
        else
            sp_mixed.probs[idx] += p
        end
    end
    sp_mixed
end

"""
    partial_prolong(sp_mixed, mask, ::Val{K}; weight_tol, binom_tol)
        -> StateSpace{CartesianIndex{K}, T}

Partial prolongation: expand a mixed-dimension state space back to K fine voxels.

  - mask[j] = false → fine pair: both dimensions copied through as-is
  - mask[j] = true  → coarse pair: nc_j split via Binomial(nc_j, 1/2)

`weight_tol` prunes states with accumulated probability below threshold.
`binom_tol`  prunes per-pair Binomial fractions below threshold (avoids
 silent mid-recursion truncation for large nc_j).
"""
function partial_prolong(sp_mixed::StateSpace{CartesianIndex{KM}, T},
                          mask::BitVector,
                          ::Val{K};
                          weight_tol::Float64 = 1e-14,
                          binom_tol::Float64  = 0.0) where {KM, K, T}
    iseven(K) || error("K=$K must be even")
    n_pairs = K ÷ 2
    length(mask) == n_pairs || error("mask length $(length(mask)) ≠ K/2=$n_pairs")
    KM == K - sum(mask) || error("KM=$KM ≠ K - sum(mask) = $(K - sum(mask))")

    # Starting index in the mixed tuple for each pair
    mixed_offsets = Vector{Int}(undef, n_pairs)
    pos = 1
    for j in 1:n_pairs
        mixed_offsets[j] = pos
        pos += mask[j] ? 1 : 2
    end

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_mixed.states)
        prob = sp_mixed.probs[i]
        prob > weight_tol || continue
        t_mixed = Tuple(sp_mixed.states[i])
        _partial_prolong_recurse!(sp_fine, t_mixed, mask, mixed_offsets,
                                   prob, fine_buf, 1, n_pairs, weight_tol, binom_tol)
    end

    sp_fine
end

function _partial_prolong_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                                    t_mixed::NTuple{KM, Int},
                                    mask::BitVector,
                                    mixed_offsets::Vector{Int},
                                    weight::T,
                                    fine_buf::Vector{Int},
                                    pair_idx::Int,
                                    n_pairs::Int,
                                    weight_tol::Float64,
                                    binom_tol::Float64) where {K, KM, T}
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

    mp = mixed_offsets[pair_idx]
    k1 = 2*pair_idx - 1
    k2 = 2*pair_idx

    if !mask[pair_idx]
        # Fine pair: both dims pass through unchanged
        fine_buf[k1] = t_mixed[mp]
        fine_buf[k2] = t_mixed[mp+1]
        _partial_prolong_recurse!(sp_fine, t_mixed, mask, mixed_offsets,
                                   weight, fine_buf, pair_idx+1, n_pairs,
                                   weight_tol, binom_tol)
    else
        # Coarse pair: enumerate Binomial(nc, 1/2) splits
        nj         = t_mixed[mp]
        log_half_n = -nj * log(2.0)
        for m in 0:nj
            log_bc = (sum(log(float(i)) for i in (nj-m+1):nj; init=0.0) -
                      sum(log(float(i)) for i in 1:m;          init=0.0))
            binom_frac   = exp(log_half_n + log_bc)
            binom_frac   > binom_tol  || continue
            binom_weight = weight * binom_frac
            binom_weight > weight_tol || continue
            fine_buf[k1] = m
            fine_buf[k2] = nj - m
            _partial_prolong_recurse!(sp_fine, t_mixed, mask, mixed_offsets,
                                       binom_weight, fine_buf, pair_idx+1, n_pairs,
                                       weight_tol, binom_tol)
        end
    end
end

"""
    select_coarsening_mask(sp, model, grid, prev_mask;
                           α_lo, α_hi, φ_lo, φ_hi) -> BitVector

Hysteretic coarsening-mask update based on pair admissibility diagnostics.

Two independent bottleneck criteria are combined:

  1. Asymmetry  α_eff[j] = max(α_intra[j], β_left[j], β_right[j])
       α_intra[j]  = |E[n_{2j-1}] − E[n_{2j}]| / (E[n_{2j-1}] + E[n_{2j}])
       β_right[j]  = |E[n_{2j}]   − E[n_{2j+1}]| / (E[n_{2j}] + E[n_{2j+1}])  (j < n_pairs)
       β_left[j]   = β_right[j-1]                                                (j > 1)

     Including inter-pair gradients ensures both pairs that touch an active
     pair-boundary front are forced fine, making refinement
     partition-independent.

  2. Per-molecule flux φ_norm[j]: low-count bottleneck guard (unchanged).

Decision rule (hysteresis):
  α_eff[j] > α_hi  OR  φ_norm[j] > φ_hi  →  mask[j] = false  (force fine)
  α_eff[j] < α_lo  AND φ_norm[j] < φ_lo  →  mask[j] = true   (coarsen)
  otherwise                               →  keep prev_mask[j]

`prev_mask` is the mask from the previous time step.  Pass `nothing` to
initialise (treated as all-coarsenable).

Defaults: α_lo=0.10, α_hi=0.25, φ_lo=1.5, φ_hi=3.0.
"""
function select_coarsening_mask(sp::StateSpace{CartesianIndex{K}, T},
                                 model::RDMEModel1D,
                                 grid::VoxelGrid,
                                 prev_mask::Union{BitVector, Nothing} = nothing;
                                 α_lo::Float64 = 0.10,
                                 α_hi::Float64 = 0.25,
                                 φ_lo::Float64 = 1.5,
                                 φ_hi::Float64 = 3.0) where {K, T}
    @assert iseven(K)
    n_pairs = K ÷ 2
    α, φ_norm = pair_admissibility(sp, model, grid)

    # Voxel means (needed for inter-pair boundary gradients)
    μ = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p * t[k]; end
    end

    # β[j] = asymmetry across the boundary between pair j and pair j+1
    # (right boundary of pair j = left boundary of pair j+1)
    β = zeros(n_pairs - 1)
    for j in 1:n_pairs-1
        total = μ[2j] + μ[2j+1]
        β[j] = total > 0.0 ? abs(μ[2j] - μ[2j+1]) / total : 0.0
    end

    mask = isnothing(prev_mask) ? trues(n_pairs) : copy(prev_mask)
    for j in 1:n_pairs
        β_left  = j > 1      ? β[j-1] : 0.0
        β_right = j < n_pairs ? β[j]   : 0.0
        α_eff   = max(α[j], β_left, β_right)

        if α_eff > α_hi || φ_norm[j] > φ_hi
            mask[j] = false   # clearly inadmissible → force fine
        elseif α_eff < α_lo && φ_norm[j] < φ_lo
            mask[j] = true    # clearly admissible  → coarsen
        end
        # else: keep prev_mask[j] (hysteresis band)
    end
    mask
end
