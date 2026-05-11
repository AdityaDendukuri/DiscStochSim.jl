# ─── Single-species Schlögl mixed-resolution 1D ──────────────────────────────
#
# Active voxels: K_act voxels tracked with joint single-species FSP.
# Inactive voxels: mean-field ODE  dμ_k/dt = f_sch(μ_k) + D/dx² · Δμ_k.
#
# State layout: CartesianIndex{K_act} where component k = molecule count
# in active voxel k (local index 1..K_act, full-grid index k_lo..k_hi).

"""
    build_active_schlogl_1d(model, grid_act, left_bc, right_bc, μ, ::Val{K_ACT})

Schlögl RDME system for K_act active voxels with mean-field boundary coupling.

- `left_bc`  : full-grid index of left inactive neighbor (0 = domain boundary)
- `right_bc` : full-grid index of right inactive neighbor (0 = domain boundary)
- `μ`        : mean counts for ALL K_total voxels (closured by reference)

Boundary flux between active boundary and inactive neighbor:
  gain into active voxel 1 from left:  d · μ[left_bc]
  loss out of active voxel 1 to left:  d · n_1
(and symmetrically on right).
"""
function build_active_schlogl_1d(
    model::SchloglModel1D,
    grid_act::VoxelGrid,
    left_bc::Int,
    right_bc::Int,
    μ::Vector{Float64},
    ::Val{K_ACT}
) where K_ACT
    K_act = grid_act.n_voxels
    K_ACT == K_act || error("K_ACT=$K_ACT ≠ K_act=$K_act")

    d  = diffusion_rate(model.D, grid_act)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
    _eu(i) = CartesianIndex(ntuple(j -> j == i ? 1 : 0, Val(K_ACT)))

    stoichs      = CartesianIndex{K_ACT}[]
    propensities = Function[]

    for k in 1:K_act
        let k = k
            ek = _eu(k)
            push!(stoichs,  ek); push!(propensities, (x, rv, t) -> c3)
            push!(stoichs, -ek); push!(propensities, (x, rv, t) -> c4 * max(0, Tuple(x)[k]))
            push!(stoichs,  ek); push!(propensities, (x, rv, t) -> begin
                n = Tuple(x)[k]; c1 * max(0, n*(n-1)) / 2 end)
            push!(stoichs, -ek); push!(propensities, (x, rv, t) -> begin
                n = Tuple(x)[k]; c2 * max(0, n*(n-1)*(n-2)) / 6 end)
        end
    end

    for k in 1:K_act-1
        let k = k
            ek = _eu(k); ek1 = _eu(k+1)
            push!(stoichs, -ek + ek1); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[k]))
            push!(stoichs,  ek - ek1); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[k+1]))
        end
    end

    if left_bc > 0
        let left_bc = left_bc
            e1 = _eu(1)
            push!(stoichs, -e1); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[1]))
            push!(stoichs,  e1); push!(propensities, (x, rv, t) -> d * μ[left_bc])
        end
    end

    if right_bc > 0
        let right_bc = right_bc
            eK = _eu(K_act)
            push!(stoichs, -eK); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[K_act]))
            push!(stoichs,  eK); push!(propensities, (x, rv, t) -> d * μ[right_bc])
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_ACT}}(stoichs, propensities)
end

"""
    update_inactive_means_schlogl_1d!(μ, sp_act, avs, model, dx, dt)

Phase 1: recompute active-voxel means from FSP.
Phase 2: advance inactive-voxel means via explicit-Euler sub-steps of
  dμ_k/dt = f_sch(μ_k) + D/dx² · Δμ_k
where f_sch(μ) = c3 + c1·μ(μ-1)/2 - c4·μ - c2·μ(μ-1)(μ-2)/6.

The mean-field closure (E[n^m] ≈ μ^m) is accurate for inactive voxels
near the Schlögl fixed points, where the distribution is narrow (Poisson-like).
"""
function update_inactive_means_schlogl_1d!(
    μ::Vector{Float64},
    sp_act,
    avs::ActiveVoxelSet1D,
    model::SchloglModel1D,
    dx::Float64,
    dt::Float64
)
    K = avs.K_total; k_lo = avs.k_lo; k_hi = avs.k_hi
    d  = model.D / dx^2
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
    f_sch(m) = c3 + c1*m*(m-1)/2 - c4*m - c2*m*(m-1)*(m-2)/6

    # Phase 1: active means from FSP marginals
    for k in k_lo:k_hi; μ[k] = 0.0; end
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t = Tuple(sp_act.states[idx])
        for (ki, k) in enumerate(k_lo:k_hi)
            μ[k] += p * t[ki]
        end
    end

    # Phase 2: inactive voxels — Schlögl mean-field ODE, explicit Euler
    # Stability: |f'(n)| is largest near n_unstable (~70-80 for default params).
    # Conservative bound: 4 × max reaction rate magnitude.
    n_ref = 80.0
    fprime_est = abs(c1*(2n_ref-1)/2 - c2*(3n_ref^2-6n_ref+2)/6 - c4)
    max_rate   = 2d + fprime_est + 1.0
    n_sub = max(1, ceil(Int, max_rate * dt * 2))
    dt_sub = dt / n_sub

    Δ = zeros(K)
    for _ in 1:n_sub
        fill!(Δ, 0.0)
        for k in 1:K
            k_lo <= k <= k_hi && continue
            Δ[k] = f_sch(μ[k])
            k > 1 && (Δ[k] += d * (μ[k-1] - μ[k]))
            k < K && (Δ[k] += d * (μ[k+1] - μ[k]))
        end
        for k in 1:K
            k_lo <= k <= k_hi && continue
            μ[k] = max(0.0, μ[k] + Δ[k] * dt_sub)
        end
    end
end

# ─── Single-species promote / demote ─────────────────────────────────────────

"""
    promote_left_schlogl(sp_act, λ, n_max, weight_tol) -> StateSpace{CartesianIndex{N+1}}

Extend the active FSP by adding a new voxel to the LEFT.
New voxel's initial distribution: Poisson(λ).
"""
function promote_left_schlogl(
    sp_act::StateSpace{CartesianIndex{N}, Float64},
    λ::Float64, n_max::Int,
    weight_tol::Float64 = 1e-14
) where N
    N_NEW  = N + 1
    n_hi_k = _poisson_upper(λ, n_max)
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        for n in 0:n_hi_k
            w = p * _poisson_pmf(n, λ)
            w > weight_tol || continue
            t_new = NTuple{N_NEW, Int}((n, t_old...))
            s = CartesianIndex(t_new)
            i = get(sp_new.index, s, 0)
            if i == 0; add_state!(sp_new, s, w) else sp_new.probs[i] += w end
        end
    end
    sp_new
end

"""
    promote_right_schlogl(sp_act, λ, n_max, weight_tol) -> StateSpace{CartesianIndex{N+1}}

Extend the active FSP by adding a new voxel to the RIGHT.
"""
function promote_right_schlogl(
    sp_act::StateSpace{CartesianIndex{N}, Float64},
    λ::Float64, n_max::Int,
    weight_tol::Float64 = 1e-14
) where N
    N_NEW  = N + 1
    n_hi_k = _poisson_upper(λ, n_max)
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        for n in 0:n_hi_k
            w = p * _poisson_pmf(n, λ)
            w > weight_tol || continue
            t_new = NTuple{N_NEW, Int}((t_old..., n))
            s = CartesianIndex(t_new)
            i = get(sp_new.index, s, 0)
            if i == 0; add_state!(sp_new, s, w) else sp_new.probs[i] += w end
        end
    end
    sp_new
end

"""
    demote_left_schlogl(sp_act) -> StateSpace{CartesianIndex{N-1}}

Remove the leftmost voxel by marginalizing it out.
"""
function demote_left_schlogl(sp_act::StateSpace{CartesianIndex{N}, Float64}) where N
    N >= 2 || error("Need ≥ 2 active voxels to demote")
    N_NEW = N - 1
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        t_new = NTuple{N_NEW, Int}(t_old[2:end])
        s = CartesianIndex(t_new)
        i = get(sp_new.index, s, 0)
        if i == 0; add_state!(sp_new, s, p) else sp_new.probs[i] += p end
    end
    sp_new
end

"""
    demote_right_schlogl(sp_act) -> StateSpace{CartesianIndex{N-1}}

Remove the rightmost voxel by marginalizing it out.
"""
function demote_right_schlogl(sp_act::StateSpace{CartesianIndex{N}, Float64}) where N
    N >= 2 || error("Need ≥ 2 active voxels to demote")
    N_NEW = N - 1
    sp_new = StateSpace{CartesianIndex{N_NEW}, Float64}()
    for idx in eachindex(sp_act.states)
        p = sp_act.probs[idx]; p > 0.0 || continue
        t_old = Tuple(sp_act.states[idx])
        t_new = NTuple{N_NEW, Int}(t_old[1:end-1])
        s = CartesianIndex(t_new)
        i = get(sp_new.index, s, 0)
        if i == 0; add_state!(sp_new, s, p) else sp_new.probs[i] += p end
    end
    sp_new
end

"""
    adapt_schlogl_1d(sp_act, avs, μ, n_max, n_lo, n_hi; kwargs...) -> (new_sp, new_avs)

Adaptive active-set management for the 1D Schlögl mixed-resolution state.

**Promote** (inactive → active):
Inactive neighbor enters the bistable transition zone:
  n_lo + promote_margin < μ[k_neighbor] < n_hi - promote_margin

**Demote** (active → mean-field):
Active boundary voxel has converged to a fixed point:
  μ[k] ≤ n_lo + demote_margin  OR  μ[k] ≥ n_hi - demote_margin

The bistable transition zone is defined by (n_lo+margin, n_hi-margin); voxels
outside this zone have unimodal distributions → mean-field is accurate.
"""
function adapt_schlogl_1d(
    sp_act,
    avs::ActiveVoxelSet1D,
    μ::Vector{Float64},
    n_max::Int,
    n_lo::Int, n_hi::Int;
    promote_margin::Float64  = 15.0,
    demote_margin::Float64   = 10.0,
    min_k_act::Int           = 2,
    max_k_act::Int           = 6,
    weight_tol::Float64      = 1e-14
)
    new_sp  = sp_act
    new_avs = avs

    in_bistable(m) = n_lo + promote_margin < m < n_hi - promote_margin
    at_fixedpt(m)  = m <= n_lo + demote_margin || m >= n_hi - demote_margin

    # ── Demote right: trailing voxel converged to fixed point ─────────────────
    while n_active(new_avs) > min_k_act
        k_hi = new_avs.k_hi
        at_fixedpt(μ[k_hi]) || break
        new_sp  = demote_right_schlogl(new_sp)
        new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_hi - 1, new_avs.K_total)
    end

    # ── Demote left: leading voxel converged to fixed point ───────────────────
    while n_active(new_avs) > min_k_act
        k_lo = new_avs.k_lo
        at_fixedpt(μ[k_lo]) || break
        new_sp  = demote_left_schlogl(new_sp)
        new_avs = ActiveVoxelSet1D(k_lo + 1, new_avs.k_hi, new_avs.K_total)
    end

    # ── Promote right: right inactive neighbor entering bistable zone ──────────
    if new_avs.k_hi < new_avs.K_total
        k_new = new_avs.k_hi + 1
        if in_bistable(μ[k_new])
            new_sp  = promote_right_schlogl(new_sp, μ[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(new_avs.k_lo, k_new, new_avs.K_total)
        end
    end

    # ── Promote left: left inactive neighbor entering bistable zone ───────────
    if new_avs.k_lo > 1
        k_new = new_avs.k_lo - 1
        if in_bistable(μ[k_new])
            new_sp  = promote_left_schlogl(new_sp, μ[k_new], n_max, weight_tol)
            new_avs = ActiveVoxelSet1D(k_new, new_avs.k_hi, new_avs.K_total)
        end
    end

    # ── Enforce max_k_act: force-demote from the less-critical boundary ───────
    while n_active(new_avs) > max_k_act
        k_lo = new_avs.k_lo; k_hi = new_avs.k_hi
        # Demote whichever end is closer to a fixed point
        if at_fixedpt(μ[k_lo]) || (!at_fixedpt(μ[k_hi]) && μ[k_lo] <= n_lo + 2*demote_margin)
            new_sp  = demote_left_schlogl(new_sp)
            new_avs = ActiveVoxelSet1D(k_lo + 1, k_hi, new_avs.K_total)
        else
            new_sp  = demote_right_schlogl(new_sp)
            new_avs = ActiveVoxelSet1D(k_lo, k_hi - 1, new_avs.K_total)
        end
    end

    new_sp, new_avs
end

"""
    schlogl_mixed_bc(n_max) -> Function

FSP truncation boundary condition: all molecule counts in [0, n_max].
"""
schlogl_mixed_bc(n_max::Int) = (s::CartesianIndex) -> all(c -> 0 ≤ c ≤ n_max, Tuple(s))
