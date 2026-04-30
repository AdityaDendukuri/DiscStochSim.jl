"""
    two_level_vcycle(sp_fine, model, fine_grid, coarse_grid, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{K}, Float64}

Two-level multigrid V-cycle for the RDME (Algorithm 6.1 of the paper).

Executes the following steps:

  1. **Pre-smooth**: evolve only intra-coarse-pair diffusion for time τ_pre.
     Brings within-pair distributions close to Binomial(n̄ⱼ, 1/2).

  2. **Restrict**: marginalize fine → coarse  (𝓡 operator).
     Saves p^{2h} = 𝓡 p̃^h (the restricted pre-smoothed distribution).

  3. **Coarse expand** (optional, before solve): add FSP boundary states to the
     coarse state space so probability can flow into them during the solve.

  4. **Coarse solve**: evolve the full coarse-level RDME for time dt → p̃^{2h}.

  5. **Prolongate correction**: compute the coarse correction δ^{2h} = p̃^{2h} - p^{2h}
     and inject it into the fine level:
         p^h_new = p̃^h + 𝓟(δ^{2h})
     (Algorithm 6.1, step 4).  This is the key difference from the operator-splitting
     variant: only the *correction* is prolongated, not the full coarse distribution.
     The fine state space therefore grows only from the small leaked probability, not
     from the entire (large) coarse support.

  6. **Post-smooth**: re-equilibrate within-pair distributions for time τ_post.

Arguments:
- `sp_fine`      : fine-level StateSpace{CartesianIndex{K}, Float64}
- `model`        : RDMEModel1D
- `fine_grid`    : VoxelGrid (finest)
- `coarse_grid`  : VoxelGrid (one level coarser, K/2 voxels)
- `rates`        : rate parameter vector (passed through to propensities)
- `t`            : current time
- `dt`           : time step for the coarse solve

Keyword arguments:
- `τ_pre`             : pre-smoothing time (default 0.0 → skip pre-smooth)
- `τ_post`            : post-smoothing time (default 0.0 → skip post-smooth)
- `krylov_m`          : Krylov subspace dimension for expmv (default 30)
- `weight_tol`        : minimum |δ| for a coarse-state correction to be prolongated
- `binom_tol`         : drop per-pair Binomial fraction < tol (default 0.0 = keep all)
- `expand_coarse`     : expand coarse FSP before coarse solve (default true)
- `coarse_expand_depth` : shell depth for coarse expansion (default 1)
- `coarse_n_max`      : per-voxel FSP truncation for coarse level (default 80)
"""
function two_level_vcycle(sp_fine::StateSpace{CartesianIndex{K}, Float64},
                           model::RDMEModel1D,
                           fine_grid::VoxelGrid,
                           coarse_grid::VoxelGrid,
                           rates, t::Real, dt::Real;
                           τ_pre::Real              = 0.0,
                           τ_post::Real             = 0.0,
                           krylov_m::Int            = 30,
                           weight_tol::Float64      = 1e-14,
                           binom_tol::Float64       = 0.0,
                           expand_coarse::Bool      = true,
                           coarse_expand_depth::Int = 1,
                           coarse_n_max::Int        = 80) where {K}
    @assert fine_grid.n_voxels  == K
    @assert coarse_grid.n_voxels == K ÷ 2
    K2 = K ÷ 2

    # ── 1. Pre-smooth ─────────────────────────────────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        intra_sys = build_intra_system(model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine, intra_sys, rates, t)
        sp_s = _copy_sp(sp_fine)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_fine.probs; m = krylov_m)
        sp_s
    else
        sp_fine
    end

    # ── 2. Restrict: p^{2h} = 𝓡 p̃^h ─────────────────────────────────────────
    sp_coarse = restrict(sp_smoothed)
    n_coarse_pre  = length(sp_coarse)        # number of states before expansion
    probs_pre     = copy(sp_coarse.probs)    # save p^{2h} for correction computation

    # ── 3. Coarse setup + optional expand (BEFORE solve) ─────────────────────
    coarse_sys = build_coarse_system(model, coarse_grid, fine_grid)

    if expand_coarse && coarse_expand_depth > 0
        coarse_bc = state -> rdme_bc(state, coarse_n_max)
        expand!(sp_coarse, coarse_sys, coarse_bc; depth = coarse_expand_depth)
    end

    # ── 4. Coarse solve: p̃^{2h} ──────────────────────────────────────────────
    A_coarse, = build_generator(sp_coarse, coarse_sys, rates, t)
    sp_coarse.probs .= expv(Float64(dt), A_coarse, sp_coarse.probs; m = krylov_m)

    # ── 5. Prolongate correction: p^h_new = p̃^h + 𝓟(p̃^{2h} - p^{2h}) ───────
    #
    # Compute δ^{2h} = p̃^{2h} - p^{2h} for each coarse state.
    # Split into positive and negative parts (prolong is defined for non-negative
    # distributions; we handle the signed correction by prolongating each part
    # separately then combining).
    #
    # New coarse states (added by expand!, indices n_coarse_pre+1 .. end) had
    # p^{2h} = 0, so their correction δ = p̃^{2h} ≥ 0.

    sp_δ_pos = StateSpace{CartesianIndex{K2}, Float64}()   # δ > 0
    sp_δ_neg = StateSpace{CartesianIndex{K2}, Float64}()   # δ < 0 (stored as |δ|)

    for i in eachindex(sp_coarse.states)
        prob_before = i ≤ n_coarse_pre ? probs_pre[i] : 0.0
        δ = sp_coarse.probs[i] - prob_before
        s = sp_coarse.states[i]
        if δ > weight_tol
            add_state!(sp_δ_pos, s, δ)
        elseif δ < -weight_tol
            add_state!(sp_δ_neg, s, -δ)   # store absolute value
        end
    end

    # Prolongate each part (non-negative, so standard prolong applies)
    sp_δ_fine_pos = prolong(sp_δ_pos, Val(K); weight_tol = weight_tol,
                                               binom_tol  = binom_tol)
    sp_δ_fine_neg = prolong(sp_δ_neg, Val(K); weight_tol = weight_tol,
                                               binom_tol  = binom_tol)

    # Inject correction into the pre-smoothed fine distribution
    sp_fine_new = _copy_sp(sp_smoothed)

    for i in eachindex(sp_δ_fine_pos.states)
        s = sp_δ_fine_pos.states[i]; δp = sp_δ_fine_pos.probs[i]
        idx = get(sp_fine_new.index, s, 0)
        if idx == 0
            add_state!(sp_fine_new, s, δp)
        else
            sp_fine_new.probs[idx] += δp
        end
    end

    for i in eachindex(sp_δ_fine_neg.states)
        s = sp_δ_fine_neg.states[i]; δp = sp_δ_fine_neg.probs[i]
        idx = get(sp_fine_new.index, s, 0)
        if idx == 0
            add_state!(sp_fine_new, s, -δp)   # may be negative; pruned externally
        else
            sp_fine_new.probs[idx] -= δp
        end
    end

    # Clip any numerically-negative probabilities to zero
    # (small negative values arise from splitting-error cancellations)
    for i in eachindex(sp_fine_new.probs)
        if sp_fine_new.probs[i] < 0.0
            sp_fine_new.probs[i] = 0.0
        end
    end

    # ── 6. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        intra_sys = build_intra_system(model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine_new, intra_sys, rates, t)
        sp_fine_new.probs .= expv(Float64(τ_post), A_intra, sp_fine_new.probs;
                                   m = krylov_m)
    end

    sp_fine_new
end

# ─── Schlögl V-cycle with dynamic-π prolongation ─────────────────────────────

"""
    two_level_vcycle_schlogl(sp_fine, model, fine_grid, coarse_grid, pi_table,
                              rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{K}, Float64}

Two-level multigrid V-cycle for the Schlögl RDME using dynamic-π prolongation.

Identical structure to `two_level_vcycle` but uses:
  - Pre-smooth:  intra-pair **diffusion only** (reactions are handled by the
    coarse operator).  A proxy `RDMEModel1D(D, 0, 0)` captures only diffusion.
  - Coarse sys:  `build_schlogl_coarse_system_dynamic(model, …, pi_table)`,
    with Galerkin rates precomputed from the dynamic-π table.
  - Prolongation: `prolong_dynamic` — distributes each coarse correction
    according to pi_table rather than Binomial(nc, 1/2).

The `pi_table` should be precomputed once via `compute_dynamic_pi`.

Keyword arguments are the same as `two_level_vcycle`.
"""
function two_level_vcycle_schlogl(sp_fine::StateSpace{CartesianIndex{K}, Float64},
                                   model::SchloglModel1D,
                                   fine_grid::VoxelGrid,
                                   coarse_grid::VoxelGrid,
                                   pi_table::Vector{Vector{Float64}},
                                   rates, t::Real, dt::Real;
                                   τ_pre::Real              = 0.0,
                                   τ_post::Real             = 0.0,
                                   krylov_m::Int            = 30,
                                   weight_tol::Float64      = 1e-14,
                                   binom_tol::Float64       = 0.0,
                                   use_dynamic_pi::Bool     = true,
                                   expand_coarse::Bool      = true,
                                   coarse_expand_depth::Int = 1,
                                   coarse_n_max::Int        = 80) where {K}
    @assert fine_grid.n_voxels   == K
    @assert coarse_grid.n_voxels == K ÷ 2
    K2 = K ÷ 2
    n_max_fine = coarse_n_max ÷ 2   # per-voxel fine truncation

    # ── 1. Pre-smooth: intra-pair diffusion only ──────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        proxy_model = RDMEModel1D(model.D, 0.0, 0.0)   # diffusion only
        intra_sys   = build_intra_system(proxy_model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine, intra_sys, rates, t)
        sp_s = _copy_sp(sp_fine)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_fine.probs; m = krylov_m)
        sp_s
    else
        sp_fine
    end

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    sp_coarse    = restrict(sp_smoothed)
    n_coarse_pre = length(sp_coarse)
    probs_pre    = copy(sp_coarse.probs)

    # ── 3. Coarse system + optional expand ────────────────────────────────────
    coarse_sys = if use_dynamic_pi
        build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
    else
        build_schlogl_coarse_system(model, coarse_grid, fine_grid)
    end
    if expand_coarse && coarse_expand_depth > 0
        coarse_bc = state -> rdme_bc(state, coarse_n_max)
        expand!(sp_coarse, coarse_sys, coarse_bc; depth = coarse_expand_depth)
    end

    # ── 4. Coarse solve ───────────────────────────────────────────────────────
    A_coarse, = build_generator(sp_coarse, coarse_sys, rates, t)
    sp_coarse.probs .= expv(Float64(dt), A_coarse, sp_coarse.probs; m = krylov_m)

    # ── 5. Prolongate correction ───────────────────────────────────────────────
    sp_δ_pos = StateSpace{CartesianIndex{K2}, Float64}()
    sp_δ_neg = StateSpace{CartesianIndex{K2}, Float64}()

    for i in eachindex(sp_coarse.states)
        prob_before = i ≤ n_coarse_pre ? probs_pre[i] : 0.0
        δ = sp_coarse.probs[i] - prob_before
        s = sp_coarse.states[i]
        if δ > weight_tol
            add_state!(sp_δ_pos, s, δ)
        elseif δ < -weight_tol
            add_state!(sp_δ_neg, s, -δ)
        end
    end

    sp_δ_fine_pos = if use_dynamic_pi
        prolong_dynamic(sp_δ_pos, pi_table, n_max_fine; weight_tol = weight_tol)
    else
        prolong(sp_δ_pos, Val(K); weight_tol = weight_tol, binom_tol = binom_tol)
    end
    sp_δ_fine_neg = if use_dynamic_pi
        prolong_dynamic(sp_δ_neg, pi_table, n_max_fine; weight_tol = weight_tol)
    else
        prolong(sp_δ_neg, Val(K); weight_tol = weight_tol, binom_tol = binom_tol)
    end

    # ── 6. Inject correction ──────────────────────────────────────────────────
    sp_fine_new = _copy_sp(sp_smoothed)

    for i in eachindex(sp_δ_fine_pos.states)
        s = sp_δ_fine_pos.states[i]; δp = sp_δ_fine_pos.probs[i]
        idx = get(sp_fine_new.index, s, 0)
        if idx == 0; add_state!(sp_fine_new, s,  δp)
        else;        sp_fine_new.probs[idx] += δp; end
    end
    for i in eachindex(sp_δ_fine_neg.states)
        s = sp_δ_fine_neg.states[i]; δp = sp_δ_fine_neg.probs[i]
        idx = get(sp_fine_new.index, s, 0)
        if idx == 0; add_state!(sp_fine_new, s, -δp)
        else;        sp_fine_new.probs[idx] -= δp; end
    end
    for i in eachindex(sp_fine_new.probs)
        sp_fine_new.probs[i] < 0.0 && (sp_fine_new.probs[i] = 0.0)
    end

    # ── 7. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        proxy_model = RDMEModel1D(model.D, 0.0, 0.0)
        intra_sys   = build_intra_system(proxy_model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine_new, intra_sys, rates, t)
        sp_fine_new.probs .= expv(Float64(τ_post), A_intra, sp_fine_new.probs;
                                   m = krylov_m)
    end

    sp_fine_new
end

# ─── Schlögl V-cycle with injection prolongation ─────────────────────────────

"""
    two_level_vcycle_schlogl_injection(sp_fine, model, fine_grid, coarse_grid,
                                        pi_table, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{K}, Float64}

Two-level V-cycle for the Schlögl RDME using **injection prolongation** (Fix 1).

Instead of prolongating a signed correction δ = p̄^{2h} - p^{2h} via Binomial or
dynamic-π splits, the fine distribution is rescaled multiplicatively:

    p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)

This preserves asymmetry in the fine distribution: a fine state concentrated at
(n_low, n_high) stays there after the coarse update rather than bleeding into the
symmetric (n_high, n_low) mode.

Keyword arguments are the same as `two_level_vcycle_schlogl`.
"""
function two_level_vcycle_schlogl_injection(
    sp_fine::StateSpace{CartesianIndex{K}, Float64},
    model::SchloglModel1D,
    fine_grid::VoxelGrid,
    coarse_grid::VoxelGrid,
    pi_table::Vector{Vector{Float64}},
    rates, t::Real, dt::Real;
    τ_pre::Real              = 0.0,
    τ_post::Real             = 0.0,
    krylov_m::Int            = 30,
    weight_tol::Float64      = 1e-14,
    use_dynamic_pi::Bool     = true,
    expand_coarse::Bool      = true,
    coarse_expand_depth::Int = 1,
    coarse_n_max::Int        = 80
) where {K}
    @assert fine_grid.n_voxels   == K
    @assert coarse_grid.n_voxels == K ÷ 2
    K2 = K ÷ 2

    # ── 1. Pre-smooth: intra-pair diffusion only ──────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        proxy_model = RDMEModel1D(model.D, 0.0, 0.0)
        intra_sys   = build_intra_system(proxy_model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine, intra_sys, rates, t)
        sp_s = _copy_sp(sp_fine)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_fine.probs; m = krylov_m)
        sp_s
    else
        sp_fine
    end

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    sp_coarse     = restrict(sp_smoothed)
    sp_coarse_pre = _copy_sp(sp_coarse)     # save p^{2h} before coarse solve for multiplicative correction

    # ── 3. Coarse system + optional expand ────────────────────────────────────
    coarse_sys = if use_dynamic_pi
        build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
    else
        build_schlogl_coarse_system(model, coarse_grid, fine_grid)
    end
    if expand_coarse && coarse_expand_depth > 0
        coarse_bc = state -> rdme_bc(state, coarse_n_max)
        expand!(sp_coarse, coarse_sys, coarse_bc; depth = coarse_expand_depth)
    end

    # ── 4. Coarse solve ───────────────────────────────────────────────────────
    A_coarse, = build_generator(sp_coarse, coarse_sys, rates, t)
    sp_coarse.probs .= expv(Float64(dt), A_coarse, sp_coarse.probs; m = krylov_m)

    # ── 5. Multiplicative prolongation with conditional extension ─────────────
    #
    # Covered coarse states: exact multiplicative correction (preserves asymmetry).
    # Newly-expanded coarse states: carry over the fine conditional from the
    # nearest covered neighbor, shifted to match the new pair count.
    sp_fine_new = prolong_conditional(sp_smoothed, sp_coarse_pre, sp_coarse;
                                           prob_tol = weight_tol)

    # ── 6. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        proxy_model = RDMEModel1D(model.D, 0.0, 0.0)
        intra_sys   = build_intra_system(proxy_model, fine_grid, coarse_grid)
        A_intra, = build_generator(sp_fine_new, intra_sys, rates, t)
        sp_fine_new.probs .= expv(Float64(τ_post), A_intra, sp_fine_new.probs;
                                   m = krylov_m)
    end

    sp_fine_new
end

# ─── helper: copy a StateSpace (states + probs, fresh index/ids) ─────────────

function _copy_sp(sp::StateSpace{E, T}) where {E, T}
    sp2 = StateSpace{E, T}()
    for i in eachindex(sp.states)
        add_state!(sp2, sp.states[i], sp.probs[i])
    end
    sp2
end

# ─── Adaptive two-level V-cycle ───────────────────────────────────────────────
"""
    two_level_vcycle_adaptive(sp_fine, model, fine_grid, rates, t, dt; kwargs...)

Adaptive V-cycle: selects a coarsening mask, restricts to mixed space, solves
the mixed CME, then prolongs back to the fine space.

Returns `(sp_fine_new, mask)` so the caller can pass `mask` as `prev_mask` on
the next time step (hysteresis).

Keyword arguments:
- `prev_mask`        : mask from previous step (nothing = first step)
- `halo`             : fine-window half-width around most imbalanced block (default 1)
- `krylov_m`         : Krylov subspace dimension for expv
- `weight_tol`       : probability weight cutoff for partial_prolong
- `binom_tol`        : Binomial truncation threshold for partial_prolong
- `expand_mixed`     : whether to expand the mixed state space before solving
- `mixed_expand_depth` : depth for expand!
- `mixed_n_max`      : per-voxel upper bound for mixed boundary condition
"""
function two_level_vcycle_adaptive(
    sp_fine::StateSpace{CartesianIndex{K}, Float64},
    model::SchloglModel1D,
    fine_grid::VoxelGrid,
    rates,
    t::Real,
    dt::Real;
    prev_mask::Union{Vector{Int}, Nothing} = nothing,
    halo::Int     = 1,
    krylov_m::Int = 30,
    weight_tol::Float64 = 1e-14,
    binom_tol::Float64 = 1e-6,
    expand_mixed::Bool = true,
    mixed_expand_depth::Int = 1,
    mixed_n_max::Int = 400,
) where {K}

    # ── 1. Select coarsening levels ───────────────────────────────────────────
    rdme_model = RDMEModel1D(model.D, 0.0, 0.0)   # diffusion-only for admissibility
    levels = select_coarsening_mask(sp_fine, rdme_model, fine_grid, prev_mask;
                                    halo)

    # If nothing is coarsened, skip the mixed machinery and solve directly.
    if !any(l > 0 for l in levels)
        fine_sys = build_schlogl_rdme_system(model, fine_grid)
        A_fine, = build_generator(sp_fine, fine_sys, rates, t)
        sp_out = _copy_sp(sp_fine)
        sp_out.probs .= expv(Float64(dt), A_fine, sp_fine.probs; m = krylov_m)
        return sp_out, levels
    end

    KM = _mixed_dim(levels)

    # ── 2. Partial restriction: fine → mixed ──────────────────────────────────
    sp_mixed = _partial_restrict_val(sp_fine, levels, Val(KM))

    # ── 3. Build mixed-resolution CME system ─────────────────────────────────
    mixed_sys = build_schlogl_mixed_system(model, fine_grid, levels)

    # ── 4. Optionally expand the mixed state space ────────────────────────────
    if expand_mixed
        bc_mixed = s -> all(c -> 0 ≤ c ≤ mixed_n_max, Tuple(s))
        expand!(sp_mixed, mixed_sys, bc_mixed; depth = mixed_expand_depth)
    end

    # ── 5. Solve mixed CME ────────────────────────────────────────────────────
    A_mixed, = build_generator(sp_mixed, mixed_sys, Float64[], Float64(t))
    sp_mixed.probs .= expv(Float64(dt), A_mixed, sp_mixed.probs; m = krylov_m)

    # ── 6. Partial prolongation: mixed → fine ─────────────────────────────────
    sp_fine_new = partial_prolong(sp_mixed, levels, Val(K);
                                   weight_tol, binom_tol)

    return sp_fine_new, levels
end

# ─── Mixed-state cache ────────────────────────────────────────────────────────
"""
    MixedSolverCache

Persistent cache for the per-mask objects that are expensive to rebuild.

Stores, keyed by `BitVector` mask:
  - `mixed_systems`: the mixed-grid CME system (reaction propensities + topology)

Pass a single `MixedSolverCache()` instance across time steps so that operators
for frequently-repeated masks are only built once.
"""
mutable struct MixedSolverCache
    mixed_systems::Dict{Vector{Int}, Any}
    MixedSolverCache() = new(Dict{Vector{Int}, Any}())
end

function _get_mixed_system!(cache::MixedSolverCache, levels, model, grid)
    get!(cache.mixed_systems, levels) do
        build_schlogl_mixed_system(model, grid, levels)
    end
end
_get_mixed_system!(::Nothing, levels, model, grid) = build_schlogl_mixed_system(model, grid, levels)

# ─── Persistent mixed-state solver ───────────────────────────────────────────
"""
    two_level_step_mixed(sp_mixed, mask, model, fine_grid, dt; kwargs...)
        -> (sp_mixed_new, new_mask)

One time step that keeps the mixed-resolution distribution as the primary state.

Unlike `two_level_vcycle_adaptive`, this never materialises the full fine grid:
  1. Conditionally recompute `new_mask` (every `mask_check_interval` steps).
  2. If the mask changed, adapt `sp_mixed` directly: fine→coarse pairs are
     marginalised; coarse→fine pairs are Binomial-expanded one pair at a time,
     optionally bounded to ±`binom_n_sigma` standard deviations.
  3. Optionally expand the mixed state space, build the mixed generator, and
     advance by `dt` via Krylov-based matrix exponentiation.

The full fine distribution is never constructed.  Prolong to the fine grid only
when a fine observable is actually needed (call `partial_prolong` explicitly).

Pass a `MixedSolverCache()` instance via `cache` to reuse mixed-grid operators
across repeated calls with the same mask.

Keyword arguments:
  `cache`                 `MixedSolverCache` or `nothing` (default nothing)
  `mask_check_interval`   recompute mask every N calls; Ref counter is updated (default 1)
  `step_count`            `Ref{Int}` tracking calls so far (default Ref(0))
  `halo`                  fine-window half-width (default 1)
  `krylov_m`              Krylov subspace dimension (default 30)
  `weight_tol`            probability weight pruning threshold (default 1e-14)
  `binom_tol`             Binomial fraction pruning threshold (default 1e-6)
  `binom_n_sigma`         keep ±n_sigma·σ of Binomial support; Inf = full (default Inf)
  `expand_mixed`          whether to expand the mixed state space before solve (default true)
  `mixed_expand_depth`    expansion depth (default 1)
  `mixed_n_max`           uniform per-dimension upper bound; used when `anisotropic_expand=false` (default 400)
  `anisotropic_expand`    if true, use per-coordinate adaptive bounds from `mixed_coord_bounds` (default false)
  `n_max_per_voxel`       per-fine-voxel floor for anisotropic bounds; 0 = use `mixed_n_max` (default 0)
  `c_sigma`               adaptive slack: bound[k] = max(floor[k], μ[k] + c_sigma·σ[k]) (default 6.0)
  `prune_tol`             remove states with prob < prune_tol after solve; 0 = off (default 0)
  `reexpand_depth`        after pruning, re-expand by this depth to track outgoing flux (default 1)
"""
function two_level_step_mixed(
    sp_mixed::StateSpace{CartesianIndex{KM}, Float64},
    levels::Vector{Int},
    model::SchloglModel1D,
    fine_grid::VoxelGrid,
    dt::Real;
    cache::Union{MixedSolverCache, Nothing} = nothing,
    mask_check_interval::Int  = 1,
    step_count::Ref{Int}      = Ref(0),
    krylov_m::Int             = 30,
    weight_tol::Float64       = 1e-14,
    binom_tol::Float64        = 1e-6,
    binom_n_sigma::Float64    = Inf,
    expand_mixed::Bool        = true,
    mixed_expand_depth::Int   = 1,
    mixed_n_max::Int          = 400,
    anisotropic_expand::Bool  = false,
    n_max_per_voxel::Int      = 0,
    c_sigma::Float64          = 6.0,
    prune_tol::Float64        = 0.0,
    reexpand_depth::Int       = 1,
    halo::Int                 = 1,
    window_shift_tol::Float64 = 0.5,
    prev_window_center::Ref{Float64} = Ref(-1.0),
) where {KM}
    K       = 2 * length(levels)
    n_pairs = K ÷ 2
    rdme_m  = RDMEModel1D(model.D, 0.0, 0.0)
    step_count[] += 1

    # 1. Conditionally recompute levels.
    # Cheap center-of-mass pre-check: skip the full diagnostic when the
    # front hasn't moved by more than window_shift_tol blocks.
    new_levels = if step_count[] % mask_check_interval == 0
        skip_full_check = false
        # Always compute pair center-of-mass (cheap O(|S^M|) scan).
        let off_cm = _mixed_offsets(levels)
            μ_pair_cm = zeros(n_pairs)
            for (i, s) in enumerate(sp_mixed.states)
                tt = Tuple(s); p = sp_mixed.probs[i]
                j = 1
                while j <= n_pairs
                    l = levels[j]; op = off_cm[j]
                    if l == 0
                        μ_pair_cm[j] += p * (tt[op] + tt[op+1])
                        j += 1
                    elseif l == 1
                        μ_pair_cm[j] += p * tt[op]
                        j += 1
                    else  # l == 2: pairs j and j+1 share coord op
                        half = p * tt[op] / 2
                        μ_pair_cm[j]   += half
                        μ_pair_cm[j+1] += half
                        j += 2
                    end
                end
            end
            μ_total_cm = sum(μ_pair_cm)
            curr_center = μ_total_cm > 0.0 ?
                sum(j * μ_pair_cm[j] for j in 1:n_pairs) / μ_total_cm :
                (prev_window_center[] >= 0.0 ? prev_window_center[] : 1.0)
            if prev_window_center[] >= 0.0
                skip_full_check = abs(curr_center - prev_window_center[]) < window_shift_tol
            end
            if !skip_full_check
                prev_window_center[] = curr_center
            end
        end

        if skip_full_check
            levels
        else
            select_coarsening_mask_mixed(sp_mixed, levels, rdme_m, fine_grid;
                                          halo)
        end
    else
        levels
    end

    # 2. Adapt mixed state if the levels changed.
    sp_cur = if new_levels == levels
        sp_mixed
    else
        adapt_mixed_state(sp_mixed, levels, new_levels, Val(K);
                          weight_tol, binom_tol, binom_n_sigma)
    end

    # 3. Retrieve (or build) the mixed system for the new levels.
    mixed_sys = _get_mixed_system!(cache, new_levels, model, fine_grid)

    # 4. Optionally expand the mixed state space.
    if expand_mixed
        _base_npv = n_max_per_voxel > 0 ? n_max_per_voxel : mixed_n_max
        bc_mixed = if anisotropic_expand
            let bnds = mixed_coord_bounds(sp_cur, new_levels;
                                          n_max_per_voxel = _base_npv, c_sigma)
                s -> begin t = Tuple(s); all(t[k] ≤ bnds[k] for k in eachindex(t)) end
            end
        else
            s -> all(c -> 0 ≤ c ≤ mixed_n_max, Tuple(s))
        end
        expand!(sp_cur, mixed_sys, bc_mixed; depth = mixed_expand_depth)
    end

    # 5. Solve mixed CME.
    A, = build_generator(sp_cur, mixed_sys, Float64[], 0.0)
    sp_cur.probs .= expv(Float64(dt), A, sp_cur.probs; m = krylov_m)

    # 6. Prune low-probability states (adaptive FSP truncation).
    if prune_tol > 0.0
        keep_idx = findall(p -> p > prune_tol, sp_cur.probs)
        if length(keep_idx) < length(sp_cur.states)
            KM_cur    = length(Tuple(sp_cur.states[1]))
            sp_pruned = StateSpace{CartesianIndex{KM_cur}, Float64}()
            for i in keep_idx
                add_state!(sp_pruned, sp_cur.states[i], sp_cur.probs[i])
            end
            sp_pruned.probs ./= sum(sp_pruned.probs)
            sp_cur = sp_pruned
        end
    end

    # 7. Re-expand around surviving mass to capture outgoing flux next step.
    if prune_tol > 0.0 && reexpand_depth > 0
        _base_npv = n_max_per_voxel > 0 ? n_max_per_voxel : mixed_n_max
        bc_re = if anisotropic_expand
            let bnds = mixed_coord_bounds(sp_cur, new_levels;
                                          n_max_per_voxel = _base_npv, c_sigma)
                s -> begin t = Tuple(s); all(t[k] ≤ bnds[k] for k in eachindex(t)) end
            end
        else
            s -> all(c -> 0 ≤ c ≤ mixed_n_max, Tuple(s))
        end
        expand!(sp_cur, mixed_sys, bc_re; depth = reexpand_depth)
    end

    return sp_cur, new_levels
end
