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

# ─── helper: copy a StateSpace (states + probs, fresh index/ids) ─────────────

function _copy_sp(sp::StateSpace{E, T}) where {E, T}
    sp2 = StateSpace{E, T}()
    for i in eachindex(sp.states)
        add_state!(sp2, sp.states[i], sp.probs[i])
    end
    sp2
end
