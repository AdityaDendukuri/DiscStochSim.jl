"""
    two_level_vcycle(sp_fine, model, fine_grid, coarse_grid, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{K}, Float64}

Two-level multigrid V-cycle for the RDME (Section 5 of the paper).

Executes the following steps:

  1. **Pre-smooth**: evolve only intra-coarse-pair diffusion for time τ_pre.
     Brings within-pair distributions close to Binomial(n̄ⱼ, 1/2).

  2. **Restrict**: marginalize fine → coarse  (𝓡 operator).

  3. **Coarse solve**: evolve the full coarse-level RDME (reactions + inter-pair
     diffusion) for time dt using a single Krylov expmv step.

  4. **Prolongate**: reconstruct fine distribution from coarse solution (𝓟 operator).
     The fine-level distribution is replaced by the prolongated coarse solution.

  5. **Post-smooth**: re-equilibrate within-pair distributions for time τ_post.

This implements the *operator-splitting* interpretation of the V-cycle
(Remark 5.1 of the paper).  The correction-based variant (which additively
updates p̃^h) is left for a future implementation.

Arguments:
- `sp_fine`      : fine-level StateSpace{CartesianIndex{K}, Float64}
- `model`        : RDMEModel1D
- `fine_grid`    : VoxelGrid (finest)
- `coarse_grid`  : VoxelGrid (one level coarser, K/2 voxels)
- `rates`        : rate parameter vector (passed through to propensities)
- `t`            : current time
- `dt`           : time step for the coarse solve

Keyword arguments:
- `τ_pre`        : pre-smoothing time (default 0.0 → skip pre-smooth)
- `τ_post`       : post-smoothing time (default 0.0 → skip post-smooth)
- `krylov_m`     : Krylov subspace dimension for expmv (default 30)
- `weight_tol`   : drop fine states with prolongated weight < tol (default 1e-14)
"""
function two_level_vcycle(sp_fine::StateSpace{CartesianIndex{K}, Float64},
                           model::RDMEModel1D,
                           fine_grid::VoxelGrid,
                           coarse_grid::VoxelGrid,
                           rates, t::Real, dt::Real;
                           τ_pre::Real       = 0.0,
                           τ_post::Real      = 0.0,
                           krylov_m::Int     = 30,
                           weight_tol::Float64 = 1e-14) where {K}
    @assert fine_grid.n_voxels  == K
    @assert coarse_grid.n_voxels == K ÷ 2

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

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    sp_coarse = restrict(sp_smoothed)

    # ── 3. Coarse solve ───────────────────────────────────────────────────────
    coarse_sys = build_coarse_system(model, coarse_grid)
    A_coarse, = build_generator(sp_coarse, coarse_sys, rates, t)
    sp_coarse.probs .= expv(Float64(dt), A_coarse, sp_coarse.probs; m = krylov_m)

    # ── 4. Prolongate ─────────────────────────────────────────────────────────
    sp_fine_new = prolong(sp_coarse, Val(K); weight_tol = weight_tol)

    # ── 5. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        intra_sys = build_intra_system(model, fine_grid, coarse_grid)
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
