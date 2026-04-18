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
    prolong_injection(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Injection prolongation: multiplicative correction applied to the pre-smoothed fine
distribution.

    p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)

where:
- p̃^h      = `sp_fine`        (pre-smoothed fine distribution)
- p^{2h}   = `sp_coarse_pre`  (restricted pre-smoothed: 𝓡 p̃^h)
- p̄^{2h}  = `sp_coarse_post` (after coarse solve)

Unlike correction prolongation, this preserves asymmetry in the fine level: if p̃^h is
concentrated at (n_low, n_high), the ratio p̄^{2h}(x̄) / p^{2h}(x̄) rescales it without
mixing in the (n_high, n_low) mode — avoiding the symmetry trap.

Fine states whose coarse denominator is below `prob_tol` are dropped.
"""
function prolong_injection(
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
    prolong_fine_conditional(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Fix-3 prolongation: uses the pre-smoothed fine conditional p̃^h(x|x̄) as the
prolongation kernel for both covered and uncovered (expanded) coarse states.

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
function prolong_fine_conditional(
    sp_fine::StateSpace{CartesianIndex{K}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{K2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{K2}, T};
    prob_tol::Float64 = 1e-300
) where {K, K2, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for fine-conditional prolongation")

    # ── Step 1: injection for covered coarse states ───────────────────────────
    sp_fine_new = prolong_injection(sp_fine, sp_coarse_pre, sp_coarse_post;
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
