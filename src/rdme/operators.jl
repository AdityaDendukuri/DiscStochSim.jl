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
                           levels::Vector{Int}) where {K, T}
    iseven(K) || error("K=$K must be even")
    n_pairs = K ÷ 2
    length(levels) == n_pairs || error("levels length $(length(levels)) ≠ K/2 = $n_pairs")
    KM = _mixed_dim(levels)
    _partial_restrict_val(sp_fine, levels, Val(KM))
end

function _mixed_tuple(t::NTuple{K,Int}, levels::Vector{Int}, ::Val{KM}) where {K, KM}
    n_pairs = K ÷ 2
    out = Vector{Int}(undef, KM)
    pos = 1
    j = 1
    while j <= n_pairs
        k1, k2 = 2j-1, 2j
        l = levels[j]
        if l == 0
            out[pos] = t[k1]; out[pos+1] = t[k2]; pos += 2; j += 1
        elseif l == 1
            out[pos] = t[k1] + t[k2]; pos += 1; j += 1
        else  # level 2: quad total of 4 voxels
            k3, k4 = 2(j+1)-1, 2(j+1)
            out[pos] = t[k1] + t[k2] + t[k3] + t[k4]; pos += 1; j += 2
        end
    end
    NTuple{KM,Int}(out)
end

function _partial_restrict_val(sp_fine::StateSpace{CartesianIndex{K}, T},
                                levels::Vector{Int}, ::Val{KM}) where {K, T, KM}
    sp_mixed = StateSpace{CartesianIndex{KM}, T}()
    for i in eachindex(sp_fine.states)
        t     = Tuple(sp_fine.states[i])
        p     = sp_fine.probs[i]
        mt    = _mixed_tuple(t, levels, Val(KM))
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
    partial_prolong(sp_mixed, levels, ::Val{K}; weight_tol, binom_tol)
        -> StateSpace{CartesianIndex{K}, T}

Partial prolongation: expand a mixed-dimension state space back to K fine voxels.

  - levels[j] = 0 → fine pair: both dimensions copied through as-is
  - levels[j] = 1 → pair-coarse: nc_j split via Binomial(nc_j, 1/2)
  - levels[j] = 2 → quad-coarse: ñ split via two nested Binomials
                    (ñ → two pair totals → four fine voxels)

`weight_tol` prunes states with accumulated probability below threshold.
`binom_tol`  prunes per-pair Binomial fractions below threshold.
"""
function partial_prolong(sp_mixed::StateSpace{CartesianIndex{KM}, T},
                          levels::Vector{Int},
                          ::Val{K};
                          weight_tol::Float64 = 1e-14,
                          binom_tol::Float64  = 0.0) where {KM, K, T}
    iseven(K) || error("K=$K must be even")
    n_pairs = K ÷ 2
    length(levels) == n_pairs || error("levels length $(length(levels)) ≠ K/2=$n_pairs")
    KM == _mixed_dim(levels) || error("KM=$KM ≠ _mixed_dim(levels) = $(_mixed_dim(levels))")

    mixed_offsets = _mixed_offsets(levels)

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_mixed.states)
        prob = sp_mixed.probs[i]
        prob > weight_tol || continue
        t_mixed = Tuple(sp_mixed.states[i])
        _partial_prolong_recurse!(sp_fine, t_mixed, levels, mixed_offsets,
                                   prob, fine_buf, 1, n_pairs, weight_tol, binom_tol)
    end

    sp_fine
end

function _partial_prolong_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                                    t_mixed::NTuple{KM, Int},
                                    levels::Vector{Int},
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
    l  = levels[pair_idx]

    if l == 0
        # Fine pair: both dims pass through unchanged
        fine_buf[k1] = t_mixed[mp]
        fine_buf[k2] = t_mixed[mp+1]
        _partial_prolong_recurse!(sp_fine, t_mixed, levels, mixed_offsets,
                                   weight, fine_buf, pair_idx+1, n_pairs,
                                   weight_tol, binom_tol)
    elseif l == 1
        # Level-1: Binomial(nc, 1/2) split
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
            _partial_prolong_recurse!(sp_fine, t_mixed, levels, mixed_offsets,
                                       binom_weight, fine_buf, pair_idx+1, n_pairs,
                                       weight_tol, binom_tol)
        end
    else
        # Level-2: quad total → two pair totals → four fine voxels (two nested Binomials)
        nq          = t_mixed[mp]  # quad total (shared by pairs j and j+1)
        k3          = 2*(pair_idx+1) - 1
        k4          = 2*(pair_idx+1)
        log_half_nq = -nq * log(2.0)
        for nb in 0:nq
            # Outer split: Binomial(nq, 1/2) → (nb, nq-nb) pair totals
            log_bc_outer = (sum(log(float(i)) for i in (nq-nb+1):nq; init=0.0) -
                            sum(log(float(i)) for i in 1:nb;          init=0.0))
            bf_outer = exp(log_half_nq + log_bc_outer)
            bf_outer > binom_tol  || continue
            w_outer  = weight * bf_outer
            w_outer  > weight_tol || continue
            nb2      = nq - nb
            log_half_nb  = -nb  * log(2.0)
            log_half_nb2 = -nb2 * log(2.0)
            for m1 in 0:nb
                # Inner split for pair j: Binomial(nb, 1/2)
                log_bc1 = (sum(log(float(i)) for i in (nb-m1+1):nb; init=0.0) -
                           sum(log(float(i)) for i in 1:m1;          init=0.0))
                bf1  = exp(log_half_nb + log_bc1)
                bf1  > binom_tol  || continue
                w1   = w_outer * bf1
                w1   > weight_tol || continue
                fine_buf[k1] = m1; fine_buf[k2] = nb - m1
                for m2 in 0:nb2
                    # Inner split for pair j+1: Binomial(nb2, 1/2)
                    log_bc2 = (sum(log(float(i)) for i in (nb2-m2+1):nb2; init=0.0) -
                               sum(log(float(i)) for i in 1:m2;            init=0.0))
                    bf2  = exp(log_half_nb2 + log_bc2)
                    bf2  > binom_tol  || continue
                    w2   = w1 * bf2
                    w2   > weight_tol || continue
                    fine_buf[k3] = m2; fine_buf[k4] = nb2 - m2
                    _partial_prolong_recurse!(sp_fine, t_mixed, levels, mixed_offsets,
                                               w2, fine_buf, pair_idx+2, n_pairs,
                                               weight_tol, binom_tol)
                end
            end
        end
    end
end

"""
    select_coarsening_mask(sp, model, grid, prev_mask; halo=1) -> Vector{Int}

Compute a level vector using the argmax-α + halo criterion:

  1. Reconstruct voxel means μ from the state.
  2. Compute α_j = |μ[2j-1] − μ[2j]| / (μ[2j-1] + μ[2j]) for all blocks j.
  3. Keep the most-imbalanced block j* and its ±halo neighbors fine (level 0);
     merge everything else (level 1).

`prev_mask` is accepted for interface compatibility but not used.
"""
function select_coarsening_mask(sp::StateSpace{CartesianIndex{K}, T},
                                 model::RDMEModel1D,
                                 grid::VoxelGrid,
                                 prev_levels::Union{Vector{Int}, Nothing} = nothing;
                                 halo::Int = 1) where {K, T}
    @assert iseven(K)
    n_pairs = K ÷ 2

    # Voxel means
    μ = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p * t[k]; end
    end

    # α_j: within-block imbalance
    α = [begin tot = μ[2j-1] + μ[2j]; tot > 0.0 ? abs(μ[2j-1] - μ[2j]) / tot : 0.0 end
         for j in 1:n_pairs]

    # Fine window: most imbalanced block ± halo; everything else merged.
    j_hot   = argmax(α)
    fine_lo = max(1, j_hot - halo)
    fine_hi = min(n_pairs, j_hot + halo)
    [fine_lo <= j <= fine_hi ? 0 : 1 for j in 1:n_pairs]
end

# ─── mixed-state mask selection and adaptation ───────────────────────────────

"""
    _mixed_dim(levels) -> Int

Number of coordinates in the mixed-resolution state tuple for a given level vector.
  level 0 (fine):       2 coords per pair
  level 1 (pair-coarse): 1 coord per pair
  level 2 (quad-coarse): 1 coord per *two* consecutive pairs (paired at j, j+1)
"""
function _mixed_dim(levels::Vector{Int})
    n = length(levels)
    dim = 0
    j = 1
    while j <= n
        l = levels[j]
        if l == 0;     dim += 2; j += 1
        elseif l == 1; dim += 1; j += 1
        else;          dim += 1; j += 2   # level-2: two pairs → one quad coord
        end
    end
    dim
end

function _mixed_offsets(levels::Vector{Int})
    n_pairs = length(levels)
    offsets = Vector{Int}(undef, n_pairs)
    pos = 1
    j = 1
    while j <= n_pairs
        l = levels[j]
        if l == 0
            offsets[j] = pos; pos += 2; j += 1
        elseif l == 1
            offsets[j] = pos; pos += 1; j += 1
        else  # level 2: pairs j and j+1 share one quad coordinate
            offsets[j] = pos; offsets[j+1] = pos; pos += 1; j += 2
        end
    end
    offsets
end

"""
    select_coarsening_mask_mixed(sp_mixed, levels, model, fine_grid; halo) -> Vector{Int}

Compute a new level vector for the mixed-resolution state using a simple
argmax-α + halo criterion:

  1. Reconstruct voxel means μ from the mixed state.
  2. Compute α_j = |μ[2j-1] - μ[2j]| / (μ[2j-1] + μ[2j]) for all blocks j.
  3. Keep the most-imbalanced block j* and its ±halo neighbors fine (level 0);
     merge everything else (level 1).

This is the complete criterion — no β, no flux guard, no hysteresis thresholds.
"""
function select_coarsening_mask_mixed(
    sp_mixed::StateSpace{CartesianIndex{KM}, T},
    levels::Vector{Int},
    model::RDMEModel1D,
    fine_grid::VoxelGrid;
    halo::Int = 1,
) where {KM, T}
    n_pairs = length(levels)
    K       = 2 * n_pairs
    KM == _mixed_dim(levels) || error("KM=$KM inconsistent with levels (expected $(_mixed_dim(levels)))")

    old_off = _mixed_offsets(levels)

    # Reconstruct fine voxel means from mixed state.
    # Level-1: μ[2j-1] = μ[2j] = E[n̄_j]/2
    μ = zeros(K)
    for (i, s) in enumerate(sp_mixed.states)
        t = Tuple(s); p = sp_mixed.probs[i]
        j = 1
        while j <= n_pairs
            op = old_off[j]
            l  = levels[j]
            if l == 0
                μ[2j-1] += p * t[op]; μ[2j] += p * t[op+1]; j += 1
            else  # level 1: both voxels share the total equally
                h = p * t[op] / 2; μ[2j-1] += h; μ[2j] += h; j += 1
            end
        end
    end

    # α_j: within-block imbalance for all blocks.
    α = [begin tot = μ[2j-1] + μ[2j]; tot > 0.0 ? abs(μ[2j-1] - μ[2j]) / tot : 0.0 end
         for j in 1:n_pairs]

    # Fine window: most imbalanced block ± halo; everything else merged.
    j_hot   = argmax(α)
    fine_lo = max(1, j_hot - halo)
    fine_hi = min(n_pairs, j_hot + halo)
    [fine_lo <= j <= fine_hi ? 0 : 1 for j in 1:n_pairs]
end

"""
    adapt_mixed_state(sp_mixed, old_mask, new_mask, ::Val{K}; weight_tol, binom_tol)
        -> StateSpace{CartesianIndex{KM_new}, T}

Transition a mixed-resolution state from `old_mask` to `new_mask` without
materialising the full fine grid.

For each pair:
  - fine → fine (mask unchanged False):     copy both dimensions
  - coarse → coarse (mask unchanged True):  copy the coarse count
  - fine → coarse:                          marginalize (sum the two fine dims)
  - coarse → fine:                          Binomial(n̄, 1/2) expansion

Only pairs that change status trigger any state-space growth.  Pairs that
remain coarse are never expanded.
"""
function adapt_mixed_state(
    sp_mixed::StateSpace{CartesianIndex{KM_old}, T},
    old_levels::Vector{Int},
    new_levels::Vector{Int},
    ::Val{K};
    weight_tol::Float64    = 1e-14,
    binom_tol::Float64     = 1e-6,
    binom_n_sigma::Float64 = Inf,
) where {KM_old, T, K}
    n_pairs = K ÷ 2
    length(old_levels) == n_pairs || error("old_levels length mismatch")
    length(new_levels) == n_pairs || error("new_levels length mismatch")
    KM_new = _mixed_dim(new_levels)
    _adapt_mixed_val(sp_mixed, old_levels, new_levels, Val(K), Val(KM_new);
                     weight_tol, binom_tol, binom_n_sigma)
end

function _adapt_mixed_val(
    sp_mixed::StateSpace{CartesianIndex{KM_old}, T},
    old_levels::Vector{Int},
    new_levels::Vector{Int},
    ::Val{K},
    ::Val{KM_new};
    weight_tol::Float64,
    binom_tol::Float64,
    binom_n_sigma::Float64,
) where {KM_old, T, K, KM_new}
    n_pairs = K ÷ 2
    old_off = _mixed_offsets(old_levels)
    new_off = _mixed_offsets(new_levels)
    sp_new  = StateSpace{CartesianIndex{KM_new}, T}()
    new_buf = zeros(Int, KM_new)

    for i in eachindex(sp_mixed.states)
        prob = sp_mixed.probs[i]
        prob > weight_tol || continue
        t_old = Tuple(sp_mixed.states[i])
        _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                               new_buf, 1, n_pairs, T(prob), weight_tol, binom_tol,
                               binom_n_sigma)
    end
    sp_new
end

# Helper: log-space Binomial coefficient
@inline function _log_binom_frac(n::Int, m::Int)
    exp(-n * log(2.0) +
        sum(log(float(i)) for i in (n-m+1):n; init=0.0) -
        sum(log(float(i)) for i in 1:m;        init=0.0))
end

function _adapt_mixed_recurse!(
    sp_new::StateSpace{CartesianIndex{KM_new}, T},
    t_old::NTuple{KM_old, Int},
    old_levels::Vector{Int},
    new_levels::Vector{Int},
    old_off::Vector{Int},
    new_off::Vector{Int},
    new_buf::Vector{Int},
    pair_idx::Int,
    n_pairs::Int,
    weight::T,
    weight_tol::Float64,
    binom_tol::Float64,
    binom_n_sigma::Float64,
) where {KM_old, KM_new, T}
    if pair_idx > n_pairs
        weight > weight_tol || return
        t_new = CartesianIndex(NTuple{KM_new, Int}(new_buf))
        idx = get(sp_new.index, t_new, 0)
        if idx == 0; add_state!(sp_new, t_new, weight)
        else;        sp_new.probs[idx] += weight; end
        return
    end

    op  = old_off[pair_idx]
    np  = new_off[pair_idx]
    old = old_levels[pair_idx]
    nw  = new_levels[pair_idx]

    # ── level-2 group: process pairs j and j+1 together ─────────────────────
    if old == 2 || nw == 2
        op2 = old == 2 ? op : old_off[pair_idx+1]   # second pair offset (same as op if old==2)
        np2 = nw  == 2 ? np : new_off[pair_idx+1]

        # Determine the quad total and per-pair info for old state
        if old == 2
            nq = t_old[op]   # single quad coord
        elseif old == 1
            nq = t_old[op] + t_old[op2]   # two pair totals
        else  # old == 0: two fine pairs
            nq = t_old[op] + t_old[op+1] + t_old[op2] + t_old[op2+1]
        end

        if nw == 2
            # Target is quad-coarse: just store the quad total
            new_buf[np] = nq
            _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                                   new_buf, pair_idx+2, n_pairs, weight,
                                   weight_tol, binom_tol, binom_n_sigma)
        elseif nw == 1
            # Quad → two pair totals: Binomial(nq, 1/2) split
            for nb in 0:nq
                frac = _log_binom_frac(nq, nb)
                frac > binom_tol || continue
                w2 = weight * frac; w2 > weight_tol || continue
                new_buf[np]  = nb
                new_buf[np2] = nq - nb
                _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                                       new_buf, pair_idx+2, n_pairs, w2,
                                       weight_tol, binom_tol, binom_n_sigma)
            end
        else
            # Quad → four fine voxels: two nested Binomials (same as level-2 prolong)
            for nb in 0:nq
                f1 = _log_binom_frac(nq, nb); f1 > binom_tol || continue
                w1 = weight * f1;               w1 > weight_tol || continue
                nb2 = nq - nb
                for m1 in 0:nb
                    f2 = _log_binom_frac(nb, m1); f2 > binom_tol || continue
                    w2 = w1 * f2;                  w2 > weight_tol || continue
                    new_buf[np] = m1; new_buf[np+1] = nb - m1
                    for m2 in 0:nb2
                        f3 = _log_binom_frac(nb2, m2); f3 > binom_tol || continue
                        w3 = w2 * f3;                    w3 > weight_tol || continue
                        new_buf[np2] = m2; new_buf[np2+1] = nb2 - m2
                        _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                                               new_buf, pair_idx+2, n_pairs, w3,
                                               weight_tol, binom_tol, binom_n_sigma)
                    end
                end
            end
        end
        return
    end

    # ── single-pair transitions (old, nw ∈ {0, 1}) ──────────────────────────
    if old == 0 && nw == 0
        new_buf[np] = t_old[op]; new_buf[np+1] = t_old[op+1]
        _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                               new_buf, pair_idx+1, n_pairs, weight,
                               weight_tol, binom_tol, binom_n_sigma)
    elseif old == 1 && nw == 1
        new_buf[np] = t_old[op]
        _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                               new_buf, pair_idx+1, n_pairs, weight,
                               weight_tol, binom_tol, binom_n_sigma)
    elseif old == 0 && nw == 1
        new_buf[np] = t_old[op] + t_old[op+1]
        _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                               new_buf, pair_idx+1, n_pairs, weight,
                               weight_tol, binom_tol, binom_n_sigma)
    else  # old == 1, nw == 0: Binomial(nc, 1/2) expansion
        nc = t_old[op]
        σ  = sqrt(nc) / 2.0
        m_lo = isinf(binom_n_sigma) ? 0  : max(0,  round(Int, nc/2 - binom_n_sigma*σ))
        m_hi = isinf(binom_n_sigma) ? nc : min(nc, round(Int, nc/2 + binom_n_sigma*σ))
        for m in m_lo:m_hi
            frac = _log_binom_frac(nc, m)
            frac > binom_tol || continue
            w2 = weight * frac; w2 > weight_tol || continue
            new_buf[np] = m; new_buf[np+1] = nc - m
            _adapt_mixed_recurse!(sp_new, t_old, old_levels, new_levels, old_off, new_off,
                                   new_buf, pair_idx+1, n_pairs, w2,
                                   weight_tol, binom_tol, binom_n_sigma)
        end
    end
end

"""
    mixed_coord_bounds(sp, levels; n_max_per_voxel, c_sigma) -> Vector{Int}

Compute per-coordinate upper bounds for the mixed-resolution FSP state space.

Each coordinate k gets an individual bound:

    bound[k] = max(r[k] * n_max_per_voxel,  ceil(μ[k] + c_sigma * σ[k]))

where:
  - `r[k]` is the voxel multiplicity of coordinate k:
       1  for level-0 (fine) coordinates,
       2  for level-1 (pair-coarse) coordinates,
       4  for level-2 (quad-coarse) coordinates.
  - `μ[k]` and `σ[k]` are the marginal mean and standard deviation of
    coordinate k under the current distribution `sp`.
  - `n_max_per_voxel` is the per-fine-voxel floor, so bulk level-1 and
    level-2 coordinates automatically receive proportionally larger floors.
  - `c_sigma` controls the adaptive slack above the mean.

Setting `c_sigma = 0.0` reduces to purely level-scaled floors.
The `n_max_per_voxel` floor ensures the bound never falls below the physical
support even when the distribution is very narrow (e.g., just after IC).
"""
function mixed_coord_bounds(
    sp    ::StateSpace{CartesianIndex{KM}, Float64},
    levels::Vector{Int};
    n_max_per_voxel::Int    = 200,
    c_sigma        ::Float64 = 6.0,
) where {KM}
    off = _mixed_offsets(levels)

    # marginal mean and E[n²] for every mixed coordinate
    μ  = zeros(KM)
    μ2 = zeros(KM)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:KM
            μ[k]  += p * t[k]
            μ2[k] += p * t[k] * t[k]
        end
    end
    σ = sqrt.(max.(μ2 .- μ .* μ, 0.0))

    # voxel multiplicity r[k] for each mixed coordinate
    r = ones(Int, KM)
    j = 1
    while j <= length(levels)
        l  = levels[j]
        op = off[j]
        if l == 0
            r[op] = 1; r[op+1] = 1; j += 1
        elseif l == 1
            r[op] = 2; j += 1
        else   # l == 2: pairs j and j+1 share coord op
            r[op] = 4; j += 2
        end
    end

    [max(r[k] * n_max_per_voxel, ceil(Int, μ[k] + c_sigma * σ[k] + 1))
     for k in 1:KM]
end
