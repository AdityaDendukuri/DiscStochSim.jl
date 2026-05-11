# ─── Active voxel set ─────────────────────────────────────────────────────────

"""
    ActiveVoxelSet1D

Tracks a contiguous range of "active" voxels in a 1D grid of K_total voxels.
Active voxels k_lo..k_hi are tracked with a full joint FSP; the rest are
mean-field (μA[k], μB[k] only).
"""
struct ActiveVoxelSet1D
    k_lo    :: Int
    k_hi    :: Int
    K_total :: Int
end

n_active(avs::ActiveVoxelSet1D)      = avs.k_hi - avs.k_lo + 1
left_inactive(avs::ActiveVoxelSet1D)  = avs.k_lo > 1          ? avs.k_lo - 1 : 0
right_inactive(avs::ActiveVoxelSet1D) = avs.k_hi < avs.K_total ? avs.k_hi + 1 : 0

# ─── System builder with mean-field boundary coupling ─────────────────────────

"""
    build_active_bottleneck_1d(model, grid_act, left_bc, right_bc, μA, μB, ::Val{N_ACT})

Build the DiscreteStochasticSystem for the active sub-grid of a 1D bottleneck RDME.

- `grid_act`  : VoxelGrid with n_voxels = K_act (active voxels only, same dx as fine grid)
- `left_bc`   : index into μA/μB for the left inactive neighbor (0 = none / domain boundary)
- `right_bc`  : index into μA/μB for the right inactive neighbor (0 = none / domain boundary)
- `μA`, `μB`  : mean counts for ALL K_total voxels, captured by reference

At the boundaries, inactive→active coupling is mean-field:
  gain at active voxel 1 from left  at rate  dA · μA[left_bc]
  loss from active voxel 1 to left  at rate  dA · nA_1
(and symmetrically on the right).  Updating μA/μB in-place keeps the
propensities current without rebuilding the system.
"""
function build_active_bottleneck_1d(
    model::BottleneckModel1D,
    grid_act::VoxelGrid,
    left_bc::Int,
    right_bc::Int,
    μA::Vector{Float64},
    μB::Vector{Float64},
    ::Val{N_ACT};
    source_voxels_full::Union{Nothing, AbstractSet{Int}} = nothing,
    k_lo_full::Int = 1   # full-grid index of active voxel 1 (for source lookup)
) where N_ACT
    K_act = grid_act.n_voxels
    N_ACT == 2 * K_act || error("N_ACT=$N_ACT must equal 2*K_act=$K_act")

    dA = diffusion_rate(model.D_A, grid_act)
    dB = diffusion_rate(model.D_B, grid_act)
    k_prod, k_AB, k_deg = model.k_prod, model.k_AB, model.k_deg

    # Determine which active-voxel indices (1..K_act) are source voxels.
    # source_voxels_full=nothing → all voxels produce (original behavior).
    is_source(k_active::Int) = source_voxels_full === nothing ||
                                (k_lo_full + k_active - 1) in source_voxels_full

    stoichs      = CartesianIndex{N_ACT}[]
    propensities = Function[]
    _eu(i) = CartesianIndex(ntuple(j -> j == i ? 1 : 0, Val(N_ACT)))

    # ── Per-voxel reactions (A→B→∅, optionally ∅→A at source voxels) ────────
    for k in 1:K_act
        let k = k
            iA = 2k - 1; iB = 2k
            eA = _eu(iA); eB = _eu(iB)
            if is_source(k)
                push!(stoichs, eA); push!(propensities, (x, rv, t) -> k_prod)
            end
            push!(stoichs, -eA + eB); push!(propensities, (x, rv, t) -> k_AB * max(0, Tuple(x)[iA]))
            push!(stoichs, -eB);      push!(propensities, (x, rv, t) -> k_deg * max(0, Tuple(x)[iB]))
        end
    end

    # ── Interior diffusion between adjacent active voxels ─────────────────────
    for k in 1:K_act - 1
        let k = k
            iAk  = 2k - 1; iAk1 = 2k + 1
            iBk  = 2k;     iBk1 = 2k + 2
            eAk  = _eu(iAk); eAk1 = _eu(iAk1)
            eBk  = _eu(iBk); eBk1 = _eu(iBk1)
            push!(stoichs, -eAk + eAk1); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAk]))
            push!(stoichs,  eAk - eAk1); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAk1]))
            push!(stoichs, -eBk + eBk1); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBk]))
            push!(stoichs,  eBk - eBk1); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBk1]))
        end
    end

    # ── Left boundary: active voxel 1 ↔ inactive left_bc ─────────────────────
    if left_bc > 0
        let left_bc = left_bc
            eA1 = _eu(1); eB1 = _eu(2)
            push!(stoichs, -eA1); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[1]))
            push!(stoichs, -eB1); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[2]))
            push!(stoichs,  eA1); push!(propensities, (x, rv, t) -> dA * μA[left_bc])
            push!(stoichs,  eB1); push!(propensities, (x, rv, t) -> dB * μB[left_bc])
        end
    end

    # ── Right boundary: active voxel K_act ↔ inactive right_bc ───────────────
    if right_bc > 0
        let right_bc = right_bc, iAK = 2K_act - 1, iBK = 2K_act
            eAK = _eu(iAK); eBK = _eu(iBK)
            push!(stoichs, -eAK); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAK]))
            push!(stoichs, -eBK); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBK]))
            push!(stoichs,  eAK); push!(propensities, (x, rv, t) -> dA * μA[right_bc])
            push!(stoichs,  eBK); push!(propensities, (x, rv, t) -> dB * μB[right_bc])
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N_ACT}}(stoichs, propensities)
end

# ─── Mixed-resolution V-cycle ─────────────────────────────────────────────────

"""
    mixed_vcycle_1d_2s(sp_h, model, hierarchy, left_bc, right_bc, μA, μB, rates, t, dt; kwargs...)

V-cycle multigrid for the 1D 2-species bottleneck RDME restricted to the active sub-grid.

- `sp_h`       : StateSpace{CartesianIndex{N}} for the K_act-voxel active sub-grid (N = 2·K_act)
- `hierarchy`  : VoxelGrid hierarchy for the active sub-grid only
- `left_bc`    : μA/μB index of left inactive neighbor (0 = domain boundary, no coupling)
- `right_bc`   : μA/μB index of right inactive neighbor (0 = no coupling)
- `μA`, `μB`   : means for all K_total voxels (referenced by propensity closures)

At every coarsening level the same left_bc / right_bc are used: the boundary
between active and inactive regions stays at the same spatial location
regardless of how finely the active sub-grid is resolved.
"""
function mixed_vcycle_1d_2s(
    sp_h::StateSpace{CartesianIndex{N}, Float64},
    model::BottleneckModel1D,
    hierarchy::Vector{VoxelGrid},
    left_bc::Int,
    right_bc::Int,
    μA::Vector{Float64},
    μB::Vector{Float64},
    rates, t::Real, dt::Real;
    τ_pre::Real         = 0.0,
    τ_post::Real        = 0.0,
    krylov_m::Int       = 30,
    weight_tol::Float64 = 1e-14,
    binom_tol::Float64  = 0.0,
    expand_coarse::Bool = true,
    coarse_r_depth::Int = 2,
    coarse_d_depth::Int = 1,
    coarse_n_max::Int   = 80,
    max_states::Int     = 1_000_000,
    prune_tol::Float64  = 0.0
) where N
    length(sp_h) <= max_states || error("State space exploded: |S| = $(length(sp_h)) > $max_states")

    N2 = N ÷ 2   # coarse state dim

    # ── Coarsest level: solve directly with full boundary coupling ────────────
    if length(hierarchy) == 1
        grid_c = hierarchy[1]
        sys_c  = build_active_bottleneck_1d(model, grid_c, left_bc, right_bc, μA, μB, Val(N))
        A_c, = build_generator(sp_h, sys_c, rates, t)
        sp_h.probs .= expv(Float64(dt), A_c, sp_h.probs; m = krylov_m)
        return sp_h
    end

    grid_h  = hierarchy[1]
    grid_2h = hierarchy[2]

    # ── 1. Pre-smooth (intra-pair diffusion only) ─────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        intra_sys = build_intra_system_2s(model, grid_h, grid_2h)
        A_intra, = build_generator(sp_h, intra_sys, rates, t)
        sp_s = _copy_sp(sp_h)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_h.probs; m = krylov_m)
        sp_s
    else
        sp_h
    end

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    sp_for_restrict = _filter_positive(sp_smoothed, weight_tol)
    sp_coarse       = restrict2s(sp_for_restrict)
    n_coarse_pre    = length(sp_coarse)
    probs_pre       = copy(sp_coarse.probs)

    # ── 3. Coarsen model (k_prod scales by r=2; diffusion D/dx² automatic) ───
    model_2h = coarsen_model(model, 2.0)

    # ── 4. Coarse expand ──────────────────────────────────────────────────────
    if expand_coarse && (coarse_r_depth > 0 || coarse_d_depth > 0)
        _expand_coarse_rd_2s!(sp_coarse, coarse_n_max, Val(N2);
                              r_depth = coarse_r_depth, d_depth = coarse_d_depth)
        length(sp_coarse) <= max_states ||
            error("Coarse state space exploded: |S_coarse| = $(length(sp_coarse))")
    end

    # ── 5. Recurse (left_bc / right_bc are constant across levels) ───────────
    sp_coarse_post = mixed_vcycle_1d_2s(
        sp_coarse, model_2h, hierarchy[2:end],
        left_bc, right_bc, μA, μB, rates, t, dt;
        τ_pre, τ_post, krylov_m, weight_tol, binom_tol,
        expand_coarse, coarse_r_depth, coarse_d_depth,
        coarse_n_max = 2 * coarse_n_max,
        max_states, prune_tol)

    # ── 6. Prolongate (multiplicative correction + Binomial for new states) ───
    sp_c_pre  = StateSpace{CartesianIndex{N2}, Float64}()
    sp_c_post = StateSpace{CartesianIndex{N2}, Float64}()
    for i in 1:n_coarse_pre
        add_state!(sp_c_pre,  sp_coarse_post.states[i], probs_pre[i])
        add_state!(sp_c_post, sp_coarse_post.states[i], sp_coarse_post.probs[i])
    end

    sp_h_new = prolong_multiplicative_2s(sp_for_restrict, sp_c_pre, sp_c_post;
                                          prob_tol = weight_tol)

    sp_δ_new = StateSpace{CartesianIndex{N2}, Float64}()
    for i in (n_coarse_pre + 1):length(sp_coarse_post.states)
        p = sp_coarse_post.probs[i]
        p > weight_tol && add_state!(sp_δ_new, sp_coarse_post.states[i], p)
    end
    if length(sp_δ_new) > 0
        sp_δf_new = prolong2s(sp_δ_new, Val(N); weight_tol, binom_tol, max_states)
        for i in eachindex(sp_δf_new.states)
            s = sp_δf_new.states[i]; δp = sp_δf_new.probs[i]
            idx = get(sp_h_new.index, s, 0)
            if idx == 0; add_state!(sp_h_new, s, δp) else sp_h_new.probs[idx] += δp end
        end
    end

    for i in eachindex(sp_h_new.probs)
        sp_h_new.probs[i] < weight_tol && (sp_h_new.probs[i] = 0.0)
    end

    # ── 7. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        intra_sys = build_intra_system_2s(model, grid_h, grid_2h)
        A_intra, = build_generator(sp_h_new, intra_sys, rates, t)
        sp_h_new.probs .= expv(Float64(τ_post), A_intra, sp_h_new.probs; m = krylov_m)
    end

    prune_tol > 0.0 && prune_threshold!(sp_h_new, prune_tol)
    sp_h_new
end

# ─── Inactive means ODE ───────────────────────────────────────────────────────

"""
    update_inactive_means_1d!(μA, μB, sp_act, avs, model, dx, dt)

Two-phase mean update for the 1D bottleneck mixed-resolution state:

Phase 1 (exact): recompute E[nA_k], E[nB_k] for active voxels from sp_act.

Phase 2 (exact for linear reactions): advance inactive-voxel means via
explicit-Euler sub-stepping of the mean-field ODEs:
  dμA_k/dt = k_prod - k_AB·μA_k + D_A/dx² · Δμ_A_k
  dμB_k/dt = k_AB·μA_k - k_deg·μB_k + D_B/dx² · Δμ_B_k
where Δ is the 1D discrete Laplacian (no-flux at domain boundaries).
Sub-step count is chosen for explicit-Euler stability.
"""
function update_inactive_means_1d!(
    μA::Vector{Float64}, μB::Vector{Float64},
    sp_act,
    avs::ActiveVoxelSet1D,
    model::BottleneckModel1D,
    dx::Float64,
    dt::Float64;
    source_voxels_full::Union{Nothing, AbstractSet{Int}} = nothing
)
    K      = avs.K_total
    k_lo   = avs.k_lo
    k_hi   = avs.k_hi
    dA     = model.D_A / dx^2
    dB     = model.D_B / dx^2
    k_prod = model.k_prod
    k_AB   = model.k_AB
    k_deg  = model.k_deg

    # Phase 1: active means from FSP
    for k in k_lo:k_hi; μA[k] = 0.0; μB[k] = 0.0; end
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t = Tuple(sp_act.states[idx])
        for (ki, k) in enumerate(k_lo:k_hi)
            μA[k] += p * t[2ki - 1]
            μB[k] += p * t[2ki]
        end
    end

    # Phase 2: inactive means via explicit Euler sub-steps
    max_rate = 2dA + k_AB + k_deg + k_prod / max(1.0, maximum(μA))
    n_sub    = max(1, ceil(Int, max_rate * dt * 2))
    dt_sub   = dt / n_sub

    ΔA = zeros(K); ΔB = zeros(K)
    for _ in 1:n_sub
        fill!(ΔA, 0.0); fill!(ΔB, 0.0)
        for k in 1:K
            k_lo <= k <= k_hi && continue   # active: skip
            is_src = source_voxels_full === nothing || k in source_voxels_full
            ΔA[k] = (is_src ? k_prod : 0.0) - k_AB * μA[k]
            ΔB[k] = k_AB * μA[k] - k_deg * μB[k]
            k > 1 && (ΔA[k] += dA*(μA[k-1] - μA[k]); ΔB[k] += dB*(μB[k-1] - μB[k]))
            k < K && (ΔA[k] += dA*(μA[k+1] - μA[k]); ΔB[k] += dB*(μB[k+1] - μB[k]))
        end
        for k in 1:K
            k_lo <= k <= k_hi && continue
            μA[k] = max(0.0, μA[k] + ΔA[k] * dt_sub)
            μB[k] = max(0.0, μB[k] + ΔB[k] * dt_sub)
        end
    end
end

# ─── Adaptive active set: promote / demote ───────────────────────────────────

@inline function _poisson_pmf(k::Int, λ::Float64)
    λ <= 0.0 && return k == 0 ? 1.0 : 0.0
    exp(-λ + k * log(λ) - sum(log(Float64(i)) for i in 1:k; init = 0.0))
end

# Smallest k such that P(X ≤ k) ≥ prob_thr for X ~ Poisson(λ).
# Used to truncate Poisson sums: for slow diffusion (λ << 1) this is ≈ 1-2,
# preventing the O(n_max²) tensor blow-up in promote_left/right_1d.
function _poisson_upper(λ::Float64, n_max::Int, prob_thr::Float64 = 0.9999)
    λ <= 0.0 && return 0
    cum = exp(-λ); p = cum; k = 0
    while cum < prob_thr && k < n_max
        k += 1; p *= λ / k; cum += p
    end
    k
end

"""
    promote_left_1d(sp_act, λA, λB, n_max, weight_tol) -> StateSpace{CartesianIndex{N+2}}

Extend the active FSP by adding a new voxel to the LEFT.
The new voxel's initial (nA, nB) distribution is Poisson(λA) × Poisson(λB),
tensorially independent of the existing state (mean-field prior).
"""
function promote_left_1d(
    sp_act::StateSpace{CartesianIndex{N}, Float64},
    λA::Float64, λB::Float64,
    n_max::Int,
    weight_tol::Float64 = 1e-14
) where N
    N_NEW  = N + 2
    n_hi_A = _poisson_upper(λA, n_max)
    n_hi_B = _poisson_upper(λB, n_max)
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        for nA in 0:n_hi_A
            pA = _poisson_pmf(nA, λA)
            pA * p > weight_tol || continue
            for nB in 0:n_hi_B
                w = p * pA * _poisson_pmf(nB, λB)
                w > weight_tol || continue
                t_new = NTuple{N_NEW, Int}((nA, nB, t_old...))
                s = CartesianIndex(t_new)
                i = get(sp_new.index, s, 0)
                if i == 0; add_state!(sp_new, s, w) else sp_new.probs[i] += w end
            end
        end
    end
    sp_new
end

"""
    promote_right_1d(sp_act, λA, λB, n_max, weight_tol) -> StateSpace{CartesianIndex{N+2}}

Extend the active FSP by adding a new voxel to the RIGHT.
"""
function promote_right_1d(
    sp_act::StateSpace{CartesianIndex{N}, Float64},
    λA::Float64, λB::Float64,
    n_max::Int,
    weight_tol::Float64 = 1e-14
) where N
    N_NEW  = N + 2
    n_hi_A = _poisson_upper(λA, n_max)
    n_hi_B = _poisson_upper(λB, n_max)
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        for nA in 0:n_hi_A
            pA = _poisson_pmf(nA, λA)
            pA * p > weight_tol || continue
            for nB in 0:n_hi_B
                w = p * pA * _poisson_pmf(nB, λB)
                w > weight_tol || continue
                t_new = NTuple{N_NEW, Int}((t_old..., nA, nB))
                s = CartesianIndex(t_new)
                i = get(sp_new.index, s, 0)
                if i == 0; add_state!(sp_new, s, w) else sp_new.probs[i] += w end
            end
        end
    end
    sp_new
end

"""
    demote_left_1d(sp_act) -> StateSpace{CartesianIndex{N-2}}

Remove the leftmost voxel from the active FSP by marginalizing it out.
The demoted voxel's mean is already held in μA[k_lo] from the last
`update_inactive_means_1d!` call and will be advanced by the ODE from here on.
"""
function demote_left_1d(sp_act::StateSpace{CartesianIndex{N}, Float64}) where N
    N >= 4 || error("Need at least 2 active voxels to demote")
    N_NEW = N - 2
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        t_new = NTuple{N_NEW, Int}(t_old[3:end])
        s = CartesianIndex(t_new)
        i = get(sp_new.index, s, 0)
        if i == 0; add_state!(sp_new, s, p) else sp_new.probs[i] += p end
    end
    sp_new
end

"""
    demote_right_1d(sp_act) -> StateSpace{CartesianIndex{N-2}}

Remove the rightmost voxel from the active FSP by marginalizing it out.
"""
function demote_right_1d(sp_act::StateSpace{CartesianIndex{N}, Float64}) where N
    N >= 4 || error("Need at least 2 active voxels to demote")
    N_NEW = N - 2
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        t_new = NTuple{N_NEW, Int}(t_old[1:end-2])
        s = CartesianIndex(t_new)
        i = get(sp_new.index, s, 0)
        if i == 0; add_state!(sp_new, s, p) else sp_new.probs[i] += p end
    end
    sp_new
end

"""
    adapt_active_1d(sp_act, avs, μA, μB, n_max; kwargs...) -> (new_sp_act, new_avs)

Adaptively expand or contract the active FSP region based on two criteria:

**Promote** (inactive → active):
- Right neighbor: μA[k_hi+1] or μB[k_hi+1] exceeds `promote_thresh`
  (wavefront arriving — need exact joint tracking).
- Left neighbor: same criterion (rare, for backward-propagating fronts).
Promotes in pairs when possible to keep K_act even for V-cycle compatibility.

**Demote** (active → mean-field) — two disjoint cases:

*Ahead-of-front* (right boundary, empty voxels):
  Both μA[k_hi] < `demote_empty_tol` AND μB[k_hi] < `demote_empty_tol`.
  Wavefront has not yet reached this voxel; mean-field is exact (all mass at 0).

*Behind-front* (left boundary, equilibrated voxels):
  Both relative errors vs analytical SS are below `demote_ss_tol`:
    |μA[k_lo]/μA_ss - 1| < `demote_ss_tol`  AND  |μB[k_lo]/μB_ss - 1| < `demote_ss_tol`
  Wavefront has passed; distribution is near-Poisson(μ_ss) → mean-field sufficient.
  Set `μA_ss = μB_ss = Inf` to disable this criterion (use only empty-based demotion).

At least `min_k_act` active voxels (≥ 2) are always retained.
"""
function adapt_active_1d(
    sp_act,
    avs::ActiveVoxelSet1D,
    μA::Vector{Float64}, μB::Vector{Float64},
    n_max::Int;
    μA_ss::Float64           = Inf,
    μB_ss::Float64           = Inf,
    promote_thresh::Float64   = 0.3,
    demote_empty_tol::Float64 = 0.02,
    demote_ss_tol::Float64    = 0.10,
    min_k_act::Int            = 2,
    max_k_act::Int            = 4,    # hard cap — prevents state-space explosion
    weight_tol::Float64       = 1e-14
)
    new_sp  = sp_act
    new_avs = avs

    # ── Demote right first: leading voxels ahead of the wavefront (empty) ────
    # Do this before promoting so we don't grow K_act unnecessarily.
    while n_active(new_avs) > min_k_act
        k_hi = new_avs.k_hi
        μA[k_hi] < demote_empty_tol && μB[k_hi] < demote_empty_tol || break
        new_sp  = demote_right_1d(new_sp)
        new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_hi - 1, new_avs.K_total)
    end

    # ── Demote left: trailing voxels converged to steady state (or empty) ────
    while n_active(new_avs) > min_k_act
        k_lo = new_avs.k_lo
        empty = μA[k_lo] < demote_empty_tol && μB[k_lo] < demote_empty_tol
        at_ss = isfinite(μA_ss) && isfinite(μB_ss) &&
                μA_ss > 0.0 && μB_ss > 0.0 &&
                abs(μA[k_lo] / μA_ss - 1.0) < demote_ss_tol &&
                abs(μB[k_lo] / μB_ss - 1.0) < demote_ss_tol
        (empty || at_ss) || break
        new_sp  = demote_left_1d(new_sp)
        new_avs = ActiveVoxelSet1D(k_lo + 1, new_avs.k_hi, new_avs.K_total)
    end

    # ── Promote right (wavefront arriving) ────────────────────────────────────
    # Promote in pairs when possible so K_act stays even for V-cycle compatibility.
    # Hard-capped at max_k_act; any excess triggers a forced left demotion below.
    if new_avs.k_hi < new_avs.K_total
        k_new = new_avs.k_hi + 1
        if μA[k_new] > promote_thresh || μB[k_new] > promote_thresh
            new_sp  = promote_right_1d(new_sp, μA[k_new], μB[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_new, new_avs.K_total)
            # Grab pair partner if it also warrants promotion (keeps K_act even)
            if new_avs.k_hi < new_avs.K_total
                k_pair = new_avs.k_hi + 1
                if μA[k_pair] > promote_thresh / 2 || μB[k_pair] > promote_thresh / 2
                    new_sp  = promote_right_1d(new_sp, μA[k_pair], μB[k_pair], n_max, weight_tol)
                    new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_pair, new_avs.K_total)
                end
            end
        end
    end

    # ── Promote left (backward front, rare) ──────────────────────────────────
    if new_avs.k_lo > 1
        k_new = new_avs.k_lo - 1
        if μA[k_new] > promote_thresh || μB[k_new] > promote_thresh
            new_sp  = promote_left_1d(new_sp, μA[k_new], μB[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(k_new, new_avs.k_hi, new_avs.K_total)
        end
    end

    # ── Enforce max_k_act cap: force-demote left until within budget ──────────
    # This is a safety valve: if fast diffusion or wide transient regions cause K_act
    # to exceed the cap, we fall back to mean-field for the trailing voxels.
    # For linear reactions this incurs no error in the means; for nonlinear systems
    # there is an approximation, bounded by the width of the non-equilibrium region.
    while n_active(new_avs) > max_k_act
        new_sp  = demote_left_1d(new_sp)
        k_lo    = new_avs.k_lo
        new_avs = ActiveVoxelSet1D(k_lo + 1, new_avs.k_hi, new_avs.K_total)
    end

    new_sp, new_avs
end

# ─── BC for active FSP ────────────────────────────────────────────────────────

"""
    active_bc(n_max) -> Function

Boundary condition for the active-voxel FSP: all molecule counts in [0, n_max].
"""
active_bc(n_max::Int) = (s::CartesianIndex) -> all(c -> 0 ≤ c ≤ n_max, Tuple(s))

# ─── A+B→∅ annihilation model ────────────────────────────────────────────────

"""
    AnnihilationModel1D

1D two-species annihilation RDME:
  ∅ →(k_A) A  [produced at left source voxels]
  ∅ →(k_B) B  [produced at right source voxels]
  A + B →(k_rxn · nA · nB) ∅  [bimolecular annihilation everywhere]
  + independent diffusion DA, DB

State layout: interleaved `(nA_1, nB_1, nA_2, nB_2, …)`.
Steady-state reaction zone width ≈ sqrt(DA / k_rxn) at the domain center.
"""
struct AnnihilationModel1D
    D_A  :: Float64
    D_B  :: Float64
    k_A  :: Float64    # ∅→A rate per source voxel
    k_B  :: Float64    # ∅→B rate per source voxel
    k_rxn :: Float64   # A+B→∅ bimolecular rate constant
end

"""
    build_active_annihilation_1d(model, grid_act, left_bc, right_bc, μA, μB, ::Val{N_ACT};
                                  source_A, source_B, k_lo_full)

System builder for the annihilation model on the active sub-grid.
`source_A` / `source_B` are Sets of FULL-GRID voxel indices where A / B is produced.
"""
function build_active_annihilation_1d(
    model::AnnihilationModel1D,
    grid_act::VoxelGrid,
    left_bc::Int,
    right_bc::Int,
    μA::Vector{Float64},
    μB::Vector{Float64},
    ::Val{N_ACT};
    source_A::AbstractSet{Int} = Set{Int}(),
    source_B::AbstractSet{Int} = Set{Int}(),
    k_lo_full::Int = 1
) where N_ACT
    K_act = grid_act.n_voxels
    N_ACT == 2 * K_act || error("N_ACT=$N_ACT ≠ 2*K_act=$K_act")

    dA   = diffusion_rate(model.D_A, grid_act)
    dB   = diffusion_rate(model.D_B, grid_act)
    k_A  = model.k_A;  k_B = model.k_B;  k_rxn = model.k_rxn

    stoichs      = CartesianIndex{N_ACT}[]
    propensities = Function[]
    _eu(i) = CartesianIndex(ntuple(j -> j == i ? 1 : 0, Val(N_ACT)))

    for k in 1:K_act
        let k = k, kg = k_lo_full + k - 1   # full-grid index
            iA = 2k-1; iB = 2k
            eA = _eu(iA); eB = _eu(iB)
            # Production
            kg in source_A && (push!(stoichs, eA);  push!(propensities, (x,rv,t) -> k_A))
            kg in source_B && (push!(stoichs, eB);  push!(propensities, (x,rv,t) -> k_B))
            # Annihilation A+B→∅
            push!(stoichs, -eA - eB)
            push!(propensities, (x,rv,t) -> k_rxn * max(0, Tuple(x)[iA]) * max(0, Tuple(x)[iB]))
        end
    end

    # Interior diffusion
    for k in 1:K_act-1
        let k = k
            iAk=2k-1; iAk1=2k+1; iBk=2k; iBk1=2k+2
            eAk=_eu(iAk); eAk1=_eu(iAk1); eBk=_eu(iBk); eBk1=_eu(iBk1)
            push!(stoichs,-eAk+eAk1); push!(propensities,(x,rv,t)->dA*max(0,Tuple(x)[iAk]))
            push!(stoichs, eAk-eAk1); push!(propensities,(x,rv,t)->dA*max(0,Tuple(x)[iAk1]))
            push!(stoichs,-eBk+eBk1); push!(propensities,(x,rv,t)->dB*max(0,Tuple(x)[iBk]))
            push!(stoichs, eBk-eBk1); push!(propensities,(x,rv,t)->dB*max(0,Tuple(x)[iBk1]))
        end
    end

    # Left boundary coupling
    if left_bc > 0
        let left_bc=left_bc
            eA1=_eu(1); eB1=_eu(2)
            push!(stoichs,-eA1); push!(propensities,(x,rv,t)->dA*max(0,Tuple(x)[1]))
            push!(stoichs,-eB1); push!(propensities,(x,rv,t)->dB*max(0,Tuple(x)[2]))
            push!(stoichs, eA1); push!(propensities,(x,rv,t)->dA*μA[left_bc])
            push!(stoichs, eB1); push!(propensities,(x,rv,t)->dB*μB[left_bc])
        end
    end

    # Right boundary coupling
    if right_bc > 0
        let right_bc=right_bc, iAK=2K_act-1, iBK=2K_act
            eAK=_eu(iAK); eBK=_eu(iBK)
            push!(stoichs,-eAK); push!(propensities,(x,rv,t)->dA*max(0,Tuple(x)[iAK]))
            push!(stoichs,-eBK); push!(propensities,(x,rv,t)->dB*max(0,Tuple(x)[iBK]))
            push!(stoichs, eAK); push!(propensities,(x,rv,t)->dA*μA[right_bc])
            push!(stoichs, eBK); push!(propensities,(x,rv,t)->dB*μB[right_bc])
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N_ACT}}(stoichs, propensities)
end

"""
    update_inactive_means_annihilation_1d!(μA, μB, sp_act, avs, model, dx, dt;
                                            source_A, source_B)

Phase 1: recompute active means from FSP.
Phase 2: advance inactive means via explicit-Euler sub-steps of
  dμA/dt = (k in source_A ? k_A : 0) - k_rxn·μA·μB + DA·∇²μA
  dμB/dt = (k in source_B ? k_B : 0) - k_rxn·μA·μB + DB·∇²μB

Note: k_rxn·μA·μB is the mean-field closure for the bimolecular reaction.
For inactive voxels dominated by a single species (μA≈0 or μB≈0) this is exact;
near the reaction front the FSP is used instead, so the closure error is confined
to regions handled by exact tracking.
"""
function update_inactive_means_annihilation_1d!(
    μA::Vector{Float64}, μB::Vector{Float64},
    sp_act,
    avs::ActiveVoxelSet1D,
    model::AnnihilationModel1D,
    dx::Float64,
    dt::Float64;
    source_A::AbstractSet{Int} = Set{Int}(),
    source_B::AbstractSet{Int} = Set{Int}()
)
    K    = avs.K_total; k_lo = avs.k_lo; k_hi = avs.k_hi
    dA   = model.D_A / dx^2; dB = model.D_B / dx^2
    k_A  = model.k_A;  k_B  = model.k_B;  k_rxn = model.k_rxn

    # Phase 1: active means from FSP
    for k in k_lo:k_hi; μA[k] = 0.0; μB[k] = 0.0; end
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t = Tuple(sp_act.states[idx])
        for (ki, k) in enumerate(k_lo:k_hi)
            μA[k] += p * t[2ki-1]; μB[k] += p * t[2ki]
        end
    end

    # Phase 2: inactive voxel ODE with explicit Euler
    max_rate = 2*max(dA,dB) + k_rxn * max(maximum(μA), 1.0) + max(k_A, k_B)
    n_sub    = max(1, ceil(Int, max_rate * dt * 2))
    dt_sub   = dt / n_sub

    ΔA = zeros(K); ΔB = zeros(K)
    for _ in 1:n_sub
        fill!(ΔA,0.0); fill!(ΔB,0.0)
        for k in 1:K
            k_lo <= k <= k_hi && continue
            rxn       = k_rxn * μA[k] * μB[k]
            ΔA[k]     = (k in source_A ? k_A : 0.0) - rxn
            ΔB[k]     = (k in source_B ? k_B : 0.0) - rxn
            k>1 && (ΔA[k]+=dA*(μA[k-1]-μA[k]); ΔB[k]+=dB*(μB[k-1]-μB[k]))
            k<K && (ΔA[k]+=dA*(μA[k+1]-μA[k]); ΔB[k]+=dB*(μB[k+1]-μB[k]))
        end
        for k in 1:K
            k_lo <= k <= k_hi && continue
            μA[k] = max(0.0, μA[k] + ΔA[k]*dt_sub)
            μB[k] = max(0.0, μB[k] + ΔB[k]*dt_sub)
        end
    end
end

"""
    adapt_annihilation_1d(sp_act, avs, μA, μB, n_max; kwargs...) -> (new_sp, new_avs)

Adaptive active-set for the annihilation model.

The reaction zone is where BOTH species are present simultaneously.
- Promote right when A arrives at right boundary (μA[k_hi+1] > promote_thresh_A)
- Promote left  when B arrives at left boundary  (μB[k_lo-1] > promote_thresh_B)
- Demote left   when B has not yet reached k_lo  (μB[k_lo] < demote_thresh_B)
  — left side is A-only; mean-field for A is exact (no bimolecular reaction)
- Demote right  when A has not yet reached k_hi  (μA[k_hi] < demote_thresh_A)
  — right side is B-only; same argument
- Hard cap at max_k_act; excess is force-demoted from the species-depleted end.
"""
function adapt_annihilation_1d(
    sp_act,
    avs::ActiveVoxelSet1D,
    μA::Vector{Float64}, μB::Vector{Float64},
    n_max::Int;
    promote_thresh_A::Float64  = 0.1,
    promote_thresh_B::Float64  = 0.1,
    demote_thresh_A::Float64   = 0.02,
    demote_thresh_B::Float64   = 0.02,
    min_k_act::Int             = 2,
    max_k_act::Int             = 4,
    weight_tol::Float64        = 1e-14
)
    new_sp  = sp_act
    new_avs = avs

    # ── Demote right (A-depleted right boundary) ──────────────────────────────
    while n_active(new_avs) > min_k_act
        k_hi = new_avs.k_hi
        μA[k_hi] < demote_thresh_A || break
        new_sp  = demote_right_1d(new_sp)
        new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_hi-1, new_avs.K_total)
    end

    # ── Demote left (B-depleted left boundary) ────────────────────────────────
    while n_active(new_avs) > min_k_act
        k_lo = new_avs.k_lo
        μB[k_lo] < demote_thresh_B || break
        new_sp  = demote_left_1d(new_sp)
        new_avs = ActiveVoxelSet1D(k_lo+1, new_avs.k_hi, new_avs.K_total)
    end

    # ── Promote right (A front arriving from left) ────────────────────────────
    if new_avs.k_hi < new_avs.K_total
        k_new = new_avs.k_hi + 1
        if μA[k_new] > promote_thresh_A
            new_sp  = promote_right_1d(new_sp, μA[k_new], μB[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_new, new_avs.K_total)
        end
    end

    # ── Promote left (B front arriving from right) ────────────────────────────
    if new_avs.k_lo > 1
        k_new = new_avs.k_lo - 1
        if μB[k_new] > promote_thresh_B
            new_sp  = promote_left_1d(new_sp, μA[k_new], μB[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(k_new, new_avs.k_hi, new_avs.K_total)
        end
    end

    # ── Hard cap: force-demote from the less-populated end ───────────────────
    while n_active(new_avs) > max_k_act
        k_lo = new_avs.k_lo; k_hi = new_avs.k_hi
        if μA[k_lo] <= μB[k_hi]   # left end is more A-like (less critical)
            new_sp  = demote_left_1d(new_sp)
            new_avs = ActiveVoxelSet1D(k_lo+1, k_hi, new_avs.K_total)
        else
            new_sp  = demote_right_1d(new_sp)
            new_avs = ActiveVoxelSet1D(k_lo, k_hi-1, new_avs.K_total)
        end
    end

    new_sp, new_avs
end
