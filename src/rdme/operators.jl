"""
    PatchPropensity

Describes the per-voxel propensities for a single species in a 2-voxel patch.

Each propensity is a function `n::Int -> Float64` giving the rate of that
channel firing in a voxel holding `n` molecules.  The sign convention for the
change in *total* patch count determines how the channel couples to the QSD
computation:

  - `conserving`   : Δ(n₁+n₂) = 0  (e.g. isomerisation A⇌B)
  - `absorbing_pos`: Δ(n₁+n₂) = +1 (e.g. birth  ∅→A)
  - `absorbing_neg`: Δ(n₁+n₂) = -1 (e.g. death  A→∅, or bimolecular A+A→…)

`absorbing_pos` and `absorbing_neg` channels contribute to the *absorption
rate* that defines the QSD: the system is absorbed whenever a molecule enters
or leaves the patch.  The QSD is the stationary distribution conditioned on
staying in the N-particle subspace.

Fields
------
- `conserving_fwd(n)` : rate of n → n+1 *within a voxel* (N-conserving half-step)
- `conserving_rev(n)` : rate of n → n-1 *within a voxel* (N-conserving half-step)
- `absorbing_pos(n)`  : per-voxel absorption rate that increases total N by 1
- `absorbing_neg(n)`  : per-voxel absorption rate that decreases total N by 1
"""
struct PatchPropensity
    conserving_fwd::Function   # n -> rate  (n → n+1 in this voxel, partner n→n-1)
    conserving_rev::Function   # n -> rate  (n → n-1 in this voxel, partner n→n+1)
    absorbing_pos::Function    # n -> rate  (birth-like, leaves N-subspace upward)
    absorbing_neg::Function    # n -> rate  (death-like, leaves N-subspace downward)
end

PatchPropensity() = PatchPropensity(_ -> 0.0, _ -> 0.0, _ -> 0.0, _ -> 0.0)

"""
    QSDProlongTable

Precomputed quasi-stationary distributions for a 2-voxel patch at each
total count N ∈ {0, …, N_max}.

`table[n1+1, N+1]` = P(n₁ = n1 | n₁+n₂ = N) under the local patch QSD.

Constructed via [`build_qsd_table`](@ref).
"""
struct QSDProlongTable
    table::Matrix{Float64}   # (N_max+1) × (N_max+1)
    N_max::Int
    d::Float64               # internal diffusion rate used to build the table
end

"""
    build_qsd_table(N_max, d; propensity, tol) -> QSDProlongTable

Precompute the quasi-stationary prolongation distribution for every total
count N ∈ {0, …, N_max} in a symmetric 2-voxel patch.

# Arguments
- `N_max`       : maximum total molecule count to tabulate.
- `d`           : internal diffusion rate  D / dx²  between the two voxels.
- `propensity`  : a [`PatchPropensity`](@ref) describing per-voxel reactions.
  Defaults to `PatchPropensity()` (pure diffusion → recovers Binomial).
- `tol`         : entries of the QSD below this threshold are zeroed.

# Algorithm
For each N, build the (N+1)×(N+1) CME generator restricted to the
N-molecule subspace, treating channels that change total count as absorbing.
When the absorption rates are *state-independent* (e.g. symmetric birth-death),
the QSD equals the ordinary stationary distribution of the conserved dynamics
and is found as the null vector of Q_local.  When absorption is state-dependent
(bimolecular, Schlögl), we solve for the leading eigenvector of Q_local - Λ,
where Λ(n₁) is the total absorption rate from state n₁.

# Notes
- Pure diffusion (no propensity) recovers `Binomial(N, 1/2)` exactly.
- For the Schlögl model pass `absorbing_pos(n) = c3` and
  `absorbing_neg(n) = c4*n + c1*n*(n-1)/2 + c2*n*(n-1)*(n-2)/6`.
"""
function build_qsd_table(N_max::Int, d::Float64;
                         propensity::PatchPropensity = PatchPropensity(),
                         tol::Float64 = 1e-14)::QSDProlongTable
    table = zeros(Float64, N_max + 1, N_max + 1)

    for N in 0:N_max
        if N == 0
            table[1, 1] = 1.0
            continue
        end

        # Build (N+1)×(N+1) local generator.
        # Index i corresponds to n₁ = i-1, n₂ = N-(i-1).
        n_states = N + 1
        Q = zeros(n_states, n_states)

        for i in 1:n_states
            n1 = i - 1
            n2 = N - n1

            # ── internal diffusion ─────────────────────────────────────────
            if n1 < N          # hop: voxel 2 → voxel 1  (n₁ increases)
                r = d * n2
                Q[i+1, i] += r;  Q[i, i] -= r
            end
            if n1 > 0          # hop: voxel 1 → voxel 2  (n₁ decreases)
                r = d * n1
                Q[i-1, i] += r;  Q[i, i] -= r
            end

            # ── N-conserving reactions (shift between voxels) ──────────────
            # Forward: molecule "produced" in voxel 1 balanced by "consumed" in voxel 2
            if n1 < N && n2 > 0
                r = propensity.conserving_fwd(n2)   # voxel 2 acts as source
                Q[i+1, i] += r;  Q[i, i] -= r
            end
            if n1 > 0 && n2 < N
                r = propensity.conserving_rev(n1)   # voxel 1 acts as source
                Q[i-1, i] += r;  Q[i, i] -= r
            end

            # ── absorption: channels that change total N ───────────────────
            # These leave the N-subspace; subtract from diagonal only (QSD formulation).
            abs_rate = propensity.absorbing_pos(n1) + propensity.absorbing_pos(n2) +
                       propensity.absorbing_neg(n1) + propensity.absorbing_neg(n2)
            Q[i, i] -= abs_rate
        end

        # ── find leading eigenvector of Q ──────────────────────────────────
        # If absorption is state-independent, Q has a null vector (ordinary ss).
        # If absorption is state-dependent, the leading (least-negative) eigenvector
        # is the QSD.  We use a unified path: find the eigenvector for the eigenvalue
        # closest to zero.
        vals, vecs = eigen(Q)                       # sorted ascending (most negative first)
        idx_lead   = argmax(real.(vals))            # eigenvalue closest to 0
        π          = real.(vecs[:, idx_lead])
        π          = abs.(π)                        # eigenvector sign is arbitrary
        π_sum      = sum(π)
        if π_sum > 0.0
            π ./= π_sum
        else
            π .= 1.0 / n_states                    # degenerate fallback: uniform
        end
        π[π .< tol] .= 0.0
        s = sum(π)
        s > 0.0 && (π ./= s)

        table[1:n_states, N+1] .= π
    end

    QSDProlongTable(table, N_max, d)
end

"""
    prolong_qsd(sp_coarse, ::Val{K}, qt; weight_tol, prob_tol) -> StateSpace

QSD prolongation 𝓟_QSD: coarse → fine using a precomputed [`QSDProlongTable`](@ref).

Replaces the Binomial assumption in [`prolong`](@ref) with the local
quasi-stationary distribution, which is exact in the fast-diffusion limit
(recovers Binomial) and captures reaction-induced asymmetry at arbitrary
Damköhler numbers.

# Arguments
- `sp_coarse`  : coarse-level `StateSpace{CartesianIndex{K2}}`.
- `::Val{K}`   : fine-level dimensionality (K = 2·K2).
- `qt`         : a `QSDProlongTable` built with [`build_qsd_table`](@ref).
- `weight_tol` : skip coarse states with probability below this.
- `prob_tol`   : skip fine states whose QSD weight is below this.
"""
function prolong_qsd(sp_coarse::StateSpace{CartesianIndex{K2}, T},
                     ::Val{K},
                     qt::QSDProlongTable;
                     weight_tol::Float64 = 1e-14,
                     prob_tol::Float64   = 0.0) where {K2, K, T}
    K == 2 * K2 || error("K=$K must equal 2*K2 for one-level QSD prolongation")

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_qsd_recurse!(sp_fine, coarse_tup, prob, fine_buf, 1, K2, qt,
                              weight_tol, prob_tol)
    end

    sp_fine
end

function _prolong_qsd_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                                coarse_tup::NTuple{K2, Int},
                                weight::T,
                                fine_buf::Vector{Int},
                                pair_idx::Int,
                                n_pairs::Int,
                                qt::QSDProlongTable,
                                weight_tol::Float64,
                                prob_tol::Float64) where {K, K2, T}
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

    N = coarse_tup[pair_idx]
    col = N + 1   # column in QSD table

    for n1 in 0:min(N, qt.N_max)
        qsd_w = qt.table[n1 + 1, col]
        qsd_w > prob_tol || continue
        fine_weight = weight * qsd_w
        fine_weight > weight_tol || continue
        fine_buf[2*pair_idx - 1] = n1
        fine_buf[2*pair_idx]     = N - n1
        _prolong_qsd_recurse!(sp_fine, coarse_tup, fine_weight, fine_buf,
                              pair_idx + 1, n_pairs, qt, weight_tol, prob_tol)
    end
end

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

Binomial prolongation 𝓟: coarse → fine via Binomial(n̄ⱼ, 1/2) per pair.

Exact when the within-pair distribution is at diffusion equilibrium (birth-death).
Two pruning thresholds:
- `weight_tol`: skip entire coarse states with probability < weight_tol.
- `binom_tol`: skip per-pair Binomial fractions < binom_tol to avoid mid-recursion
  truncation for small-probability coarse states with large K.
"""
function prolong(sp_coarse::StateSpace{CartesianIndex{K2}, T},
                 ::Val{K};
                 weight_tol::Float64 = 1e-14,
                 binom_tol::Float64  = 0.0,
                 max_states::Int     = 10_000_000) where {K2, K, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for one-level prolongation")

    sp_fine  = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_recurse!(sp_fine, coarse_tup, prob, fine_buf, 1, K2,
                          weight_tol, binom_tol)
        length(sp_fine) > max_states && error("Prolongation exploded: |S_fine| > $max_states")
    end

    sp_fine
end

"""
    prolong_product(sp_fine_pre, sp_coarse; weight_tol)
        → StateSpace{CartesianIndex{K}, Float64}

Product-convolution prolongation: the principled CME analog of polynomial
interpolation in PDE multigrid.

For each pair j = {2j-1, 2j} with coarse total N_j:

  p(n₁, n₂ | N_j) ∝ m_j1(n₁) × m_j2(N_j − n₁)

where m_j1, m_j2 are the marginals of the two fine voxels extracted from the
PRE-SMOOTHED fine state `sp_fine_pre`.

Why this is the right prolongation
-----------------------------------
In PDEs, prolongation interpolates using the SMOOTHNESS of the solution —
polynomials of degree p reproduce exact values at p+1 points.

In the CME, the analog of smoothness is QUASI-STATIONARITY: after the
pre-smoother runs for time τ_pre, the within-pair distribution relaxes toward
its QSD.  At QSD the joint factorises approximately as p(n₁,n₂) ≈ p₁(n₁)⋅p₂(n₂).
The product convolution is therefore exact in the limit τ_pre → ∞.

Contrast with Binomial prolongation (current default):
  Binomial:  p(n₁|N) = Bin(N, 0.5) — assumes EQUAL MEANS within the pair
  Product:   p(n₁|N) ∝ m₁(n₁)⋅m₂(N−n₁) — uses ACTUAL MARGINALS from pre-smooth

For a bistable interface pair with one voxel at n_lo and one at n_hi:
  Binomial gives a unimodal peak at N/2 (the WRONG transition state)
  Product gives a BIMODAL peak at (n_lo, n_hi) and (n_hi, n_lo)  ← correct

Normalization
--------------
For each pair j, the convolution Z_j = Σ_n m_j1(n) ⋅ m_j2(N_j−n) is computed
and divided out so that total probability is exactly preserved.

Arguments
----------
- `sp_fine_pre` : StateSpace at the fine level, POST-pre-smooth (marginals extracted here)
- `sp_coarse`   : StateSpace at the coarse level, post-coarse-solve
- `weight_tol`  : skip coarse states with probability < weight_tol
"""
function prolong_product(
    sp_fine_pre::StateSpace{CartesianIndex{K}, Float64},
    sp_coarse::StateSpace{CartesianIndex{K2}, Float64};
    weight_tol::Float64 = 1e-13,
    max_states::Int     = 10_000_000
) where {K, K2}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for product-convolution prolongation")

    # ── 1. Extract per-voxel marginals from the pre-smoothed fine state ────────
    # Determine support upper bound
    n_max_est = 0
    for ci in sp_fine_pre.states
        for n in Tuple(ci); n > n_max_est && (n_max_est = n); end
    end

    # marginals[k][n+1] = P(voxel k has n molecules) in pre-smoothed distribution
    marginals = [zeros(Float64, n_max_est + 2) for _ in 1:K]
    for (ci, prob) in zip(sp_fine_pre.states, sp_fine_pre.probs)
        t = Tuple(ci)
        for k in 1:K
            n = t[k]; marginals[k][n+1] += prob
        end
    end
    # Normalize (sum might be slightly < 1 after pruning)
    for mk in marginals
        s = sum(mk); s > 0.0 && (mk ./= s)
    end

    # ── 2. Build fine state space via product-convolution over pairs ───────────
    sp_fine = StateSpace{CartesianIndex{K}, Float64}()
    buf     = zeros(Int, K)

    for (cs, cp) in zip(sp_coarse.states, sp_coarse.probs)
        cp > weight_tol || continue
        ct = Tuple(cs)

        # Per-pair normalization constants Z_j = Σ_n m_j1(n)⋅m_j2(N_j−n)
        Z = ntuple(Val(K2)) do j
            N_j = ct[j]; m1 = marginals[2j-1]; m2 = marginals[2j]
            z   = 0.0
            for n1 in 0:min(N_j, length(m1)-1)
                n2 = N_j - n1
                n2 < length(m2) && (z += m1[n1+1] * m2[n2+1])
            end
            z > 0.0 ? z : 1.0
        end

        _prolong_product_recurse!(sp_fine, cp, buf, Val(K), Val(K2),
                                   ct, marginals, Z, weight_tol)
        length(sp_fine) > max_states &&
            error("Product prolongation exploded: |S_fine| > $max_states")
    end

    sp_fine
end

# ── Recursive tensor product over K2 pairs ────────────────────────────────────
# pair_idx goes 1..K2; builds CartesianIndex{K} in buf when all pairs done.
function _prolong_product_recurse!(
    sp_fine::StateSpace{CartesianIndex{K}, Float64},
    weight::Float64,
    buf::Vector{Int},
    ::Val{K}, ::Val{K2},
    ct::NTuple{K2, Int},
    marginals::Vector{Vector{Float64}},
    Z::NTuple{K2, Float64},
    weight_tol::Float64,
    pair_idx::Int = 1
) where {K, K2}
    if pair_idx > K2
        ci  = CartesianIndex(NTuple{K, Int}(buf))
        idx = get(sp_fine.index, ci, 0)
        if idx == 0; add_state!(sp_fine, ci, weight)
        else;        sp_fine.probs[idx] += weight; end
        return
    end

    N_j = ct[pair_idx]; k1 = 2pair_idx - 1; k2 = 2pair_idx
    m1  = marginals[k1]; m2 = marginals[k2]
    Z_j = Z[pair_idx]

    for n1 in 0:min(N_j, length(m1)-1)
        n2 = N_j - n1
        n2 >= length(m2) && continue
        w  = weight * m1[n1+1] * m2[n2+1] / Z_j
        w  > weight_tol || continue
        buf[k1] = n1; buf[k2] = n2
        _prolong_product_recurse!(sp_fine, w, buf, Val(K), Val(K2),
                                   ct, marginals, Z, weight_tol, pair_idx + 1)
    end
end

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
    log_half_n = -nj * log(2.0)

    for m in 0:nj
        log_bc     = (sum(log(i) for i in (nj-m+1):nj; init=0.0) -
                      sum(log(i) for i in 1:m;          init=0.0))
        binom_frac = exp(log_half_n + log_bc)
        binom_frac   > binom_tol  || continue
        binom_weight = weight * binom_frac
        binom_weight > weight_tol || continue
        fine_buf[2*pair_idx - 1] = m
        fine_buf[2*pair_idx]     = nj - m
        _prolong_recurse!(sp_fine, coarse_tup, binom_weight, fine_buf,
                          pair_idx + 1, n_pairs, weight_tol, binom_tol)
    end
end

"""
    prolong_multiplicative(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Multiplicative (injection) prolongation:

    p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)

Preserves within-pair asymmetry: mass at (n_low, n_high) stays there rather than
being symmetrized. Fine states whose coarse denominator is below `prob_tol` are dropped.
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

        t_fine   = Tuple(x_fine)
        x_coarse = CartesianIndex(ntuple(j -> t_fine[2j-1] + t_fine[2j], Val(K2)))

        idx_pre  = get(sp_coarse_pre.index,  x_coarse, 0)
        idx_post = get(sp_coarse_post.index, x_coarse, 0)

        p_pre  = idx_pre  > 0 ? sp_coarse_pre.probs[idx_pre]   : 0.0
        p_post = idx_post > 0 ? sp_coarse_post.probs[idx_post] : 0.0

        p_new = p_pre > prob_tol ? p_fine * (p_post / p_pre) : 0.0
        p_new > 0.0 && add_state!(sp_fine_new, x_fine, p_new)
    end

    sp_fine_new
end

"""
    prolong_conditional(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{K}, T}

Multiplicative prolongation with conditional extension for bistable systems.

**Covered coarse states** (in `sp_coarse_pre`): multiplicative update
    p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)
Asymmetry-preserving: mass at (n_lo, n_hi) stays there.

**Uncovered coarse states** (added by expand!): carry over the fine conditional
from the nearest covered L1=1 neighbor (reaction boundary) or L1=2 neighbor
(coarse diffusion boundary), shifting the changed voxels to match the new pair count.
"""
function prolong_conditional(
    sp_fine::StateSpace{CartesianIndex{K}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{K2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{K2}, T};
    prob_tol::Float64 = 1e-300
) where {K, K2, T}
    K == 2 * K2 || error("K=$K must equal 2*K2=$K2 for conditional prolongation")

    # Step 1: multiplicative correction for covered coarse states
    sp_fine_new = prolong_multiplicative(sp_fine, sp_coarse_pre, sp_coarse_post;
                                         prob_tol = prob_tol)

    # Step 2: build covered-state conditionals p̃^h(x|x̄_cov)
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
        push!(get!(cov_to_fine, idx, Tuple{CartesianIndex{K}, Float64}[]), (x_fine, p_fine))
    end
    for (idx, list) in cov_to_fine
        total = cov_totals[idx]
        total > prob_tol || continue
        for i in eachindex(list)
            (s, p) = list[i]; list[i] = (s, p / total)
        end
    end

    # Step 3: extend to expand!-added coarse states (indices > n_pre).
    # States 1..n_pre are covered (handled by multiplicative step).
    # States n_pre+1..end were added by coarse expand! with p=0 and need
    # fine states initialized from the nearest covered neighbor.
    #
    # Two cases, both principled by homogeneity:
    #
    # L1=1 (one pair changes by ±1 from a reaction):
    #   GAIN (+1): new molecule appears uniformly in the pair → both orientations
    #              (n1+1, n2) and (n1, n2+1) each with weight 1/2.
    #   LOSS (−1): any molecule leaves with probability proportional to its voxel count →
    #              (n1−1, n2) with weight n1/nc, (n1, n2−1) with weight n2/nc.
    #
    # L1=2 (adjacent pairs j and j+1 change by ∓1 — coarse diffusion):
    #   Molecule crossed the specific inter-pair boundary voxel:
    #   - Left-to-right (j loses, j+1 gains): exits RIGHT voxel of j, enters LEFT voxel of j+1.
    #   - Right-to-left (j+1 loses, j gains): exits LEFT voxel of j+1, enters RIGHT voxel of j.
    #   Only fine states with a molecule at the exit voxel contribute.
    n_pre = length(sp_coarse_pre)

    # Helper: accumulate a fine state into sp_fine_new
    function _add!(s::CartesianIndex{K}, w::Float64)
        w > prob_tol || return
        idx2 = get(sp_fine_new.index, s, 0)
        if idx2 == 0; add_state!(sp_fine_new, s, w) else sp_fine_new.probs[idx2] += w end
    end

    for i in (n_pre + 1):length(sp_coarse_post.states)
        p_coarse_new = sp_coarse_post.probs[i]
        p_coarse_new > prob_tol || continue
        x_coarse_new = sp_coarse_post.states[i]
        t_new = Tuple(x_coarse_new)

        # ── Find covered neighbor with fine support ──────────────────────────
        cov_idx   = 0
        cov_pair  = 0
        cov_pair2 = 0
        cov_delta1 = 0
        cov_delta2 = 0

        # L1=1: one pair changes by ±1 (reaction)
        for j in 1:K2
            for delta in (1, -1)
                t_try = ntuple(k -> k == j ? t_new[k] - delta : t_new[k], Val(K2))
                all(v -> v ≥ 0, t_try) || continue
                idx = get(sp_coarse_pre.index, CartesianIndex(t_try), 0)
                if idx != 0 && haskey(cov_to_fine, idx)
                    cov_idx = idx; cov_pair = j; cov_delta1 = delta; break
                end
            end
            cov_idx != 0 && break
        end

        # L1=2: adjacent pairs j and j+1 change by opposite ±1 (coarse diffusion)
        if cov_idx == 0 && K2 >= 2
            found2 = false
            for j in 1:(K2-1)
                found2 && break
                for d in (1, -1)   # d = delta for pair j; pair j+1 changes by -d
                    t_try = ntuple(k -> k == j ? t_new[k] - d :
                                       k == j+1 ? t_new[k] + d : t_new[k], Val(K2))
                    all(v -> v ≥ 0, t_try) || continue
                    idx = get(sp_coarse_pre.index, CartesianIndex(t_try), 0)
                    if idx != 0 && haskey(cov_to_fine, idx)
                        cov_idx = idx; cov_pair = j; cov_pair2 = j+1
                        cov_delta1 = d; cov_delta2 = -d
                        found2 = true; break
                    end
                end
            end
        end

        cov_idx != 0 || continue

        # ── Extend fine states ────────────────────────────────────────────────
        for (x_fine, cond_prob) in cov_to_fine[cov_idx]
            cond_prob > prob_tol || continue
            w0 = cond_prob * p_coarse_new
            t_fine = Tuple(x_fine)

            if cov_pair2 == 0
                # ── L1=1 reaction: pair cov_pair changes by cov_delta1 ──────
                n1 = t_fine[2*cov_pair - 1]
                n2 = t_fine[2*cov_pair]
                nc = n1 + n2

                if cov_delta1 == 1  # new state has one MORE molecule (nc+1 in new)
                    # Molecule gained: equal probability of appearing in either sub-voxel
                    s1 = CartesianIndex(ntuple(k -> k==2*cov_pair-1 ? n1+1 :
                                                    k==2*cov_pair   ? n2   : t_fine[k], Val(K)))
                    s2 = CartesianIndex(ntuple(k -> k==2*cov_pair-1 ? n1   :
                                                    k==2*cov_pair   ? n2+1 : t_fine[k], Val(K)))
                    _add!(s1, w0 * 0.5)
                    _add!(s2, w0 * 0.5)

                else   # cov_delta1 == -1: new state has one FEWER molecule (nc-1 in new)
                    # Molecule lost: probability proportional to local count
                    nc > 0 || continue
                    if n1 > 0
                        s1 = CartesianIndex(ntuple(k -> k==2*cov_pair-1 ? n1-1 :
                                                        k==2*cov_pair   ? n2   : t_fine[k], Val(K)))
                        _add!(s1, w0 * (n1 / nc))
                    end
                    if n2 > 0
                        s2 = CartesianIndex(ntuple(k -> k==2*cov_pair-1 ? n1   :
                                                        k==2*cov_pair   ? n2-1 : t_fine[k], Val(K)))
                        _add!(s2, w0 * (n2 / nc))
                    end
                end

            else
                # ── L1=2 diffusion: pair cov_pair ↔ pair cov_pair+1 ──────────
                # cov_delta1 = how pair cov_pair changed (covered→new): +1 means cov_pair GAINED
                # so diffusion direction: cov_pair2→cov_pair (right-to-left) if cov_delta1=+1
                #                        cov_pair→cov_pair2 (left-to-right)  if cov_delta1=-1
                na1 = t_fine[2*cov_pair - 1];   na2 = t_fine[2*cov_pair]
                nb1 = t_fine[2*cov_pair2 - 1];  nb2 = t_fine[2*cov_pair2]

                if cov_delta1 == -1
                    # LEFT-TO-RIGHT: pair cov_pair loses via its RIGHT voxel (na2)
                    #                pair cov_pair2 gains via its LEFT voxel (nb1)
                    na2 > 0 || continue   # no molecule at exit voxel — skip
                    s = CartesianIndex(ntuple(k ->
                        k == 2*cov_pair-1  ? na1      :
                        k == 2*cov_pair    ? na2 - 1  :
                        k == 2*cov_pair2-1 ? nb1 + 1  :
                        k == 2*cov_pair2   ? nb2       : t_fine[k], Val(K)))
                    _add!(s, w0)
                else
                    # RIGHT-TO-LEFT: pair cov_pair2 loses via its LEFT voxel (nb1)
                    #                pair cov_pair  gains via its RIGHT voxel (na2)
                    nb1 > 0 || continue   # no molecule at exit voxel — skip
                    s = CartesianIndex(ntuple(k ->
                        k == 2*cov_pair-1  ? na1      :
                        k == 2*cov_pair    ? na2 + 1  :
                        k == 2*cov_pair2-1 ? nb1 - 1  :
                        k == 2*cov_pair2   ? nb2       : t_fine[k], Val(K)))
                    _add!(s, w0)
                end
            end
        end
    end

    sp_fine_new
end

# ─── diagnostics ─────────────────────────────────────────────────────────────

"""
    lumpability_ratio(model, fine_grid) -> Float64

ε = d_inter / d_intra = 1/4 for a uniform grid (always, for homogeneous diffusion).
"""
lumpability_ratio(model::RDMEModel1D, fine_grid::VoxelGrid) =
    diffusion_rate(model.D, coarsen(fine_grid)) / diffusion_rate(model.D, fine_grid)

"""
    pair_admissibility(sp, model, grid) -> (α, φ_norm)

For each voxel pair {2j-1, 2j}, compute:
  α_j     = |E[n_{2j-1}] - E[n_{2j}]| / (E[n_{2j-1}] + E[n_{2j}])  (within-pair asymmetry)
  φ_norm_j = normalised per-molecule flux at the intra-pair boundary

α ≈ 0 → Binomial prolongation is appropriate. α ≈ 1 → bistable front, use conditional.
"""
function pair_admissibility(sp::StateSpace{CartesianIndex{K}, T},
                             model::RDMEModel1D,
                             grid::VoxelGrid) where {K, T}
    @assert iseven(K)
    n_pairs = K ÷ 2

    μ = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p * t[k]; end
    end

    d = diffusion_rate(model.D, grid)
    φ_raw  = [d * μ[k] / max(1.0, μ[k]) for k in 1:K-1]
    φ_mean = mean(φ_raw)

    α      = Vector{Float64}(undef, n_pairs)
    φ_norm = Vector{Float64}(undef, n_pairs)
    for j in 1:n_pairs
        k1, k2 = 2j-1, 2j
        total   = μ[k1] + μ[k2]
        α[j]     = total > 0.0 ? abs(μ[k1] - μ[k2]) / total : 0.0
        φ_norm[j] = φ_mean > 0.0 ? φ_raw[k1] / φ_mean : 0.0
    end

    α, φ_norm
end

# ─── 2-species restrict / prolong ────────────────────────────────────────────

"""
    restrict2s(sp_fine) -> StateSpace{CartesianIndex{N÷2}, T}

Restriction for 2-species interleaved states `(nA_1,nB_1,…,nA_K,nB_K)`.
Pairs fine voxels (2j-1, 2j) into coarse voxel j; sums each species independently.
"""
function restrict2s(sp_fine::StateSpace{CartesianIndex{N}, T}) where {N, T}
    N % 4 == 0 || error("N=$N must be divisible by 4 for 2-species restriction")
    N2 = N ÷ 2   # coarse state dimension

    sp_coarse = StateSpace{CartesianIndex{N2}, T}()

    for i in eachindex(sp_fine.states)
        t    = Tuple(sp_fine.states[i])
        prob = sp_fine.probs[i]

        # ci=1,2 → coarse voxel 1; ci=3,4 → coarse voxel 2; etc.
        # Odd ci: A species → t[4j-3]+t[4j-1]  (j = (ci+1)÷2)
        # Even ci: B species → t[4j-2]+t[4j]
        coarse_tuple = ntuple(ci -> begin
            j = (ci + 1) ÷ 2
            isodd(ci) ? t[4j-3] + t[4j-1] : t[4j-2] + t[4j]
        end, Val(N2))

        coarse_state = CartesianIndex(coarse_tuple)
        idx = get(sp_coarse.index, coarse_state, 0)
        if idx == 0; add_state!(sp_coarse, coarse_state, prob)
        else; sp_coarse.probs[idx] += prob; end
    end

    sp_coarse
end

"""
    prolong2s(sp_coarse, ::Val{N}; weight_tol, binom_tol, max_states)
        -> StateSpace{CartesianIndex{N}, T}

Binomial prolongation for 2-species interleaved states.
For each coarse voxel j with counts (nA_bar, nB_bar), draws nA_1 ~ Bin(nA_bar,½)
and nB_1 ~ Bin(nB_bar,½) independently.
"""
function prolong2s(sp_coarse::StateSpace{CartesianIndex{N2}, T},
                   ::Val{N};
                   weight_tol::Float64 = 1e-14,
                   binom_tol::Float64  = 0.0,
                   max_states::Int     = 10_000_000) where {N2, N, T}
    N == 2 * N2 || error("N=$N must equal 2*N2=$N2 for 2-species prolongation")
    N2 % 2 == 0 || error("N2=$N2 must be even")
    K2_vox = N2 ÷ 2

    sp_fine  = StateSpace{CartesianIndex{N}, T}()
    fine_buf = zeros(Int, N)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong2s_recurse!(sp_fine, coarse_tup, prob, fine_buf, 1, K2_vox,
                            weight_tol, binom_tol)
        length(sp_fine) > max_states && error("2-species prolongation exploded: |S| > $max_states")
    end

    sp_fine
end

function _prolong2s_recurse!(sp_fine::StateSpace{CartesianIndex{N}, T},
                              coarse_tup::NTuple{N2, Int},
                              weight::T,
                              fine_buf::Vector{Int},
                              pair_idx::Int,
                              n_pairs::Int,
                              weight_tol::Float64,
                              binom_tol::Float64) where {N, N2, T}
    if pair_idx > n_pairs
        fine_state = CartesianIndex(NTuple{N, Int}(fine_buf))
        idx = get(sp_fine.index, fine_state, 0)
        if idx == 0; add_state!(sp_fine, fine_state, weight)
        else; sp_fine.probs[idx] += weight; end
        return
    end

    nA_bar = coarse_tup[2*pair_idx - 1]
    nB_bar = coarse_tup[2*pair_idx]
    log_half_nA = -nA_bar * log(2.0)
    log_half_nB = -nB_bar * log(2.0)

    for mA in 0:nA_bar
        log_bcA = (sum(log(Float64(i)) for i in (nA_bar-mA+1):nA_bar; init=0.0) -
                   sum(log(Float64(i)) for i in 1:mA; init=0.0))
        fracA = exp(log_half_nA + log_bcA)
        fracA > binom_tol || continue
        wA = weight * fracA
        wA > weight_tol || continue

        for mB in 0:nB_bar
            log_bcB = (sum(log(Float64(i)) for i in (nB_bar-mB+1):nB_bar; init=0.0) -
                       sum(log(Float64(i)) for i in 1:mB; init=0.0))
            fracB = exp(log_half_nB + log_bcB)
            fracB > binom_tol || continue
            w = wA * fracB
            w > weight_tol || continue

            # fine voxel 2j-1: (nA=mA, nB=mB);  fine voxel 2j: (nA=nA_bar-mA, nB=nB_bar-mB)
            fine_buf[4*pair_idx - 3] = mA
            fine_buf[4*pair_idx - 2] = mB
            fine_buf[4*pair_idx - 1] = nA_bar - mA
            fine_buf[4*pair_idx]     = nB_bar - mB
            _prolong2s_recurse!(sp_fine, coarse_tup, w, fine_buf,
                                pair_idx + 1, n_pairs, weight_tol, binom_tol)
        end
    end
end

"""
    prolong_multiplicative_2s(sp_fine, sp_coarse_pre, sp_coarse_post; prob_tol)
        -> StateSpace{CartesianIndex{N}, T}

Multiplicative prolongation for 2-species interleaved states.
"""
function prolong_multiplicative_2s(
    sp_fine::StateSpace{CartesianIndex{N}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{N2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{N2}, T};
    prob_tol::Float64 = 1e-300
) where {N, N2, T}
    N == 2 * N2 || error("N=$N must equal 2*N2=$N2")

    sp_fine_new = StateSpace{CartesianIndex{N}, T}()

    for i in eachindex(sp_fine.states)
        x_fine = sp_fine.states[i]
        p_fine = sp_fine.probs[i]
        p_fine > 0.0 || continue

        t_fine = Tuple(x_fine)
        x_coarse = CartesianIndex(ntuple(ci -> begin
            j = (ci + 1) ÷ 2
            isodd(ci) ? t_fine[4j-3] + t_fine[4j-1] : t_fine[4j-2] + t_fine[4j]
        end, Val(N2)))

        idx_pre  = get(sp_coarse_pre.index,  x_coarse, 0)
        idx_post = get(sp_coarse_post.index, x_coarse, 0)
        p_pre  = idx_pre  > 0 ? sp_coarse_pre.probs[idx_pre]   : 0.0
        p_post = idx_post > 0 ? sp_coarse_post.probs[idx_post] : 0.0

        p_new = p_pre > prob_tol ? p_fine * (p_post / p_pre) : 0.0
        p_new > 0.0 && add_state!(sp_fine_new, x_fine, p_new)
    end

    sp_fine_new
end

# ─── CoarseningMap-based 2-species operators ─────────────────────────────────
# Generic versions that work for any patch size (M=2 for 1D, M=4 for 2D).

"""
    restrict2s(sp_fine, cmap, ::Val{N2}) -> StateSpace{CartesianIndex{N2}, T}

Generic 2-species restriction using CoarseningMap.
For each coarse voxel J: nA_coarse = Σ nA_k, nB_coarse = Σ nB_k over k in patch J.
Works for M=2 (1D pairs) and M=4 (2D 2×2 blocks).
"""
function restrict2s(sp_fine::StateSpace{CartesianIndex{N}, T},
                    cmap::CoarseningMap,
                    ::Val{N2}) where {N, N2, T}
    N2 == 2 * cmap.n_coarse || error("N2=$N2 ≠ 2×n_coarse=$(2*cmap.n_coarse)")
    sp_coarse = StateSpace{CartesianIndex{N2}, T}()

    for i in eachindex(sp_fine.states)
        t    = Tuple(sp_fine.states[i])
        prob = sp_fine.probs[i]
        coarse_tuple = ntuple(ci -> begin
            J     = (ci + 1) ÷ 2
            patch = cmap.coarse_to_fine[J]
            isodd(ci) ? sum(t[2k-1] for k in patch) : sum(t[2k] for k in patch)
        end, Val(N2))
        coarse_state = CartesianIndex(coarse_tuple)
        idx = get(sp_coarse.index, coarse_state, 0)
        if idx == 0; add_state!(sp_coarse, coarse_state, prob)
        else; sp_coarse.probs[idx] += prob; end
    end
    sp_coarse
end

"""
    prolong2s(sp_coarse, cmap, ::Val{N}; kwargs...) -> StateSpace{CartesianIndex{N}, T}

Generic 2-species Multinomial prolongation using CoarseningMap.
For each coarse voxel J with patch size M:
  nA split ~ Multinomial(nA_bar, 1/M,...); same for nB independently.
Supports M=2 (Binomial, 1D) and M=4 (4-way Multinomial, 2D).
"""
function prolong2s(sp_coarse::StateSpace{CartesianIndex{N2}, T},
                   cmap::CoarseningMap,
                   ::Val{N};
                   weight_tol::Float64 = 1e-14,
                   binom_tol::Float64  = 0.0,
                   max_states::Int     = 10_000_000) where {N2, N, T}
    N == 2 * cmap.n_fine   || error("N=$N ≠ 2×n_fine=$(2*cmap.n_fine)")
    N2 == 2 * cmap.n_coarse || error("N2=$N2 ≠ 2×n_coarse=$(2*cmap.n_coarse)")
    sp_fine  = StateSpace{CartesianIndex{N}, T}()
    fine_buf = zeros(Int, N)
    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]
        prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong2s_cmap_core!(sp_fine, coarse_tup, prob, fine_buf, 1, cmap, weight_tol, binom_tol)
        length(sp_fine) > max_states && error("2s prolong exploded: |S| > $max_states")
    end
    sp_fine
end

function _prolong2s_cmap_core!(sp_fine::StateSpace{CartesianIndex{N}, T},
                                coarse_tup::NTuple{N2, Int},
                                weight::T,
                                fine_buf::Vector{Int},
                                J::Int,
                                cmap::CoarseningMap,
                                weight_tol::Float64,
                                binom_tol::Float64) where {N, N2, T}
    if J > cmap.n_coarse
        fine_state = CartesianIndex(NTuple{N, Int}(fine_buf))
        idx = get(sp_fine.index, fine_state, 0)
        if idx == 0; add_state!(sp_fine, fine_state, weight)
        else; sp_fine.probs[idx] += weight; end
        return
    end

    nA_bar = coarse_tup[2J-1]
    nB_bar = coarse_tup[2J]
    patch  = cmap.coarse_to_fine[J]
    M      = cmap.patch_size

    log_A0 = sum(log(Float64(i)) for i in 1:nA_bar; init=0.0) - nA_bar * log(Float64(M))
    log_B0 = sum(log(Float64(i)) for i in 1:nB_bar; init=0.0) - nB_bar * log(Float64(M))

    if M == 2
        k1, k2 = patch[1], patch[2]
        for mA in 0:nA_bar
            lA = log_A0 - (sum(log(Float64(i)) for i in 1:mA; init=0.0) +
                            sum(log(Float64(i)) for i in 1:(nA_bar-mA); init=0.0))
            wA = exp(lA)
            wA > binom_tol || continue
            for mB in 0:nB_bar
                lB = log_B0 - (sum(log(Float64(i)) for i in 1:mB; init=0.0) +
                                sum(log(Float64(i)) for i in 1:(nB_bar-mB); init=0.0))
                w = weight * wA * exp(lB)
                w > weight_tol || continue
                fine_buf[2k1-1]=mA;  fine_buf[2k1]=mB
                fine_buf[2k2-1]=nA_bar-mA;  fine_buf[2k2]=nB_bar-mB
                _prolong2s_cmap_core!(sp_fine, coarse_tup, w, fine_buf, J+1, cmap, weight_tol, binom_tol)
            end
        end

    elseif M == 4
        k1, k2, k3, k4 = patch[1], patch[2], patch[3], patch[4]
        for nA1 in 0:nA_bar
            lA1 = sum(log(Float64(i)) for i in 1:nA1; init=0.0)
            for nA2 in 0:(nA_bar-nA1)
                lA12 = lA1 + sum(log(Float64(i)) for i in 1:nA2; init=0.0)
                for nA3 in 0:(nA_bar-nA1-nA2)
                    nA4  = nA_bar - nA1 - nA2 - nA3
                    lA   = log_A0 - (lA12 + sum(log(Float64(i)) for i in 1:nA3; init=0.0) +
                                             sum(log(Float64(i)) for i in 1:nA4; init=0.0))
                    wA   = exp(lA)
                    wA > binom_tol || continue
                    for nB1 in 0:nB_bar
                        lB1 = sum(log(Float64(i)) for i in 1:nB1; init=0.0)
                        for nB2 in 0:(nB_bar-nB1)
                            lB12 = lB1 + sum(log(Float64(i)) for i in 1:nB2; init=0.0)
                            for nB3 in 0:(nB_bar-nB1-nB2)
                                nB4 = nB_bar - nB1 - nB2 - nB3
                                lB  = log_B0 - (lB12 + sum(log(Float64(i)) for i in 1:nB3; init=0.0) +
                                                        sum(log(Float64(i)) for i in 1:nB4; init=0.0))
                                w   = weight * wA * exp(lB)
                                w > weight_tol || continue
                                fine_buf[2k1-1]=nA1; fine_buf[2k1]=nB1
                                fine_buf[2k2-1]=nA2; fine_buf[2k2]=nB2
                                fine_buf[2k3-1]=nA3; fine_buf[2k3]=nB3
                                fine_buf[2k4-1]=nA4; fine_buf[2k4]=nB4
                                _prolong2s_cmap_core!(sp_fine, coarse_tup, w, fine_buf, J+1, cmap, weight_tol, binom_tol)
                            end
                        end
                    end
                end
            end
        end
    else
        error("Unsupported patch_size=$(cmap.patch_size); only M=2 (1D) and M=4 (2D) supported")
    end
end

"""
    prolong_multiplicative_2s(sp_fine, sp_cpre, sp_cpost, cmap, ::Val{N2}; prob_tol)
        -> StateSpace{CartesianIndex{N}, T}

Multiplicative correction prolongation for 2-species using CoarseningMap.
"""
function prolong_multiplicative_2s(
    sp_fine::StateSpace{CartesianIndex{N}, T},
    sp_coarse_pre::StateSpace{CartesianIndex{N2}, T},
    sp_coarse_post::StateSpace{CartesianIndex{N2}, T},
    cmap::CoarseningMap,
    ::Val{N2};
    prob_tol::Float64 = 1e-300
) where {N, N2, T}
    sp_fine_new = StateSpace{CartesianIndex{N}, T}()
    for i in eachindex(sp_fine.states)
        x_fine = sp_fine.states[i]; p_fine = sp_fine.probs[i]
        p_fine > 0.0 || continue
        t = Tuple(x_fine)
        x_coarse = CartesianIndex(ntuple(ci -> begin
            J     = (ci + 1) ÷ 2
            patch = cmap.coarse_to_fine[J]
            isodd(ci) ? sum(t[2k-1] for k in patch) : sum(t[2k] for k in patch)
        end, Val(N2)))
        idx_pre  = get(sp_coarse_pre.index,  x_coarse, 0)
        idx_post = get(sp_coarse_post.index, x_coarse, 0)
        p_pre    = idx_pre  > 0 ? sp_coarse_pre.probs[idx_pre]   : 0.0
        p_post   = idx_post > 0 ? sp_coarse_post.probs[idx_post] : 0.0
        p_new    = p_pre > prob_tol ? p_fine * (p_post / p_pre) : 0.0
        p_new > 0.0 && add_state!(sp_fine_new, x_fine, p_new)
    end
    sp_fine_new
end

# ─── Adaptive coarsening: marginal reconstruction ─────────────────────────────

"""
    fine_marginals_from_coarse(sp_coarse, cmap, n_max_fine) → Vector{Vector{Float64}}

Compute K_fine per-voxel marginal distributions from a coarse StateSpace,
using uniform Multinomial(n̄_j, 1/Mⱼ, …) splitting within each coarse region j.

Returns `marginals[k][n+1]` = P(nₖ = n) for k = 1…K_fine, n = 0…n_max_fine.

Valid in the fast-diffusion limit: within each merged region the distribution
is near-Multinomial, so marginalising over the Binomial/Multinomial split is
exact for birth-death and a good approximation for Schlögl with fast diffusion.
The key advantage over full prolongation: no joint fine state space is constructed.
Cost is O(K_c × n̄_max × n_max_fine) instead of O(n_max_fine^K_fine).
"""
function fine_marginals_from_coarse(
    sp_coarse::StateSpace{CartesianIndex{K_c}, T},
    cmap::CoarseningMap,
    n_max_fine::Int
) where {K_c, T}
    K_c == cmap.n_coarse || error("K_c=$K_c ≠ cmap.n_coarse=$(cmap.n_coarse)")
    K_f  = cmap.n_fine
    marg = [zeros(Float64, n_max_fine + 1) for _ in 1:K_f]

    for i in eachindex(sp_coarse.states)
        p_c = sp_coarse.probs[i]
        p_c > 0.0 || continue
        tc = Tuple(sp_coarse.states[i])

        for j in 1:K_c
            n̄    = tc[j]
            patch = cmap.coarse_to_fine[j]
            M     = length(patch)
            q     = 1.0 / M

            log_q   = M == 1 ? 0.0  : log(q)
            log_1mq = M == 1 ? -Inf : log(1.0 - q)

            for n in 0:min(n̄, n_max_fine)
                if M == 1
                    bw = n == n̄ ? 1.0 : 0.0
                else
                    # log Bin(n̄, q)(n) = log C(n̄,n) + n log q + (n̄-n) log(1-q)
                    lc = 0.0
                    for ii in 1:n;     lc += log(Float64(n̄ - n + ii)) - log(Float64(ii)); end
                    lp = lc + n * log_q + (n̄ - n) * log_1mq
                    bw = exp(lp)
                end
                bw > 1e-15 || continue
                contrib = p_c * bw
                for k in patch
                    marg[k][n + 1] += contrib
                end
            end
        end
    end
    marg
end
