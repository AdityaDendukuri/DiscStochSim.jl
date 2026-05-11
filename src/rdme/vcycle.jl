# ─── helpers ─────────────────────────────────────────────────────────────────

# Expand the coarse state space with separate reaction and diffusion BFS depths.
#
# Reaction expansion (r_depth): each state reaches neighbors where any single coarse
# voxel count changes by ±1 (all-stoichiometry-1 reactions of Schlögl/birth-death).
# This is iterated r_depth times so each nc_j can range ±r from any covered state.
#
# Diffusion expansion (d_depth): adjacent pairs (j, j+1) exchange one molecule at a
# time (nc_j ± 1, nc_{j+1} ∓ 1), iterated d_depth times.
#
# Both expansions add new states with p=0; existing states and probabilities are unchanged.
function _expand_coarse_rd!(sp::StateSpace{CartesianIndex{K2}, T},
                             n_max::Int,
                             ::Val{K2};
                             r_depth::Int = 2,
                             d_depth::Int = 1) where {K2, T}
    # BFS over reaction transitions
    for _ in 1:r_depth
        frontier = sp.states[begin:end]   # snapshot current states
        for s in frontier
            t = Tuple(s)
            for j in 1:K2, δ in (-1, 1)
                nj = t[j] + δ
                (nj >= 0 && nj <= n_max) || continue
                s_new = CartesianIndex(ntuple(k -> k == j ? nj : t[k], Val(K2)))
                haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
            end
        end
    end
    # BFS over diffusion transitions (adjacent pairs exchange ±1)
    K2 >= 2 || return
    for _ in 1:d_depth
        frontier = sp.states[begin:end]
        for s in frontier
            t = Tuple(s)
            for j in 1:(K2-1), δ in (-1, 1)
                nj  = t[j]   + δ
                nj1 = t[j+1] - δ
                (nj >= 0 && nj <= n_max && nj1 >= 0 && nj1 <= n_max) || continue
                s_new = CartesianIndex(ntuple(k -> k == j ? nj : k == j+1 ? nj1 : t[k], Val(K2)))
                haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
            end
        end
    end
end

function _copy_sp(sp::StateSpace{E, T}) where {E, T}
    sp2 = StateSpace{E, T}()
    for i in eachindex(sp.states)
        add_state!(sp2, sp.states[i], sp.probs[i])
    end
    sp2
end

# Return a new StateSpace containing only states with prob > tol.
function _filter_positive(sp::StateSpace{E, T}, tol::Float64) where {E, T}
    sp2 = StateSpace{E, T}()
    for i in eachindex(sp.states)
        sp.probs[i] > tol && add_state!(sp2, sp.states[i], sp.probs[i])
    end
    sp2
end

# Build the 2-voxel Schlögl stationary distribution table for dynamic-π prolongation.
function _compute_schlogl_pi_table(model::SchloglModel1D, grid::VoxelGrid, n_max::Int)
    sp2  = StateSpace{CartesianIndex{2}, Float64}()
    sys2 = build_schlogl_rdme_system(model, VoxelGrid(2, grid.dx, grid.level))
    bc2  = s -> rdme_bc(s, n_max)
    add_state!(sp2, CartesianIndex(0, 0), 1.0)
    expand!(sp2, sys2, bc2; depth = n_max)
    A2, = build_generator(sp2, sys2, Float64[], 0.0)
    compute_dynamic_pi(sp2, A2; n_max = n_max)
end

"""
    multi_level_vcycle(sp_fine, model, hierarchy, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{K}, Float64}

N-level recursive multigrid V-cycle for the RDME.

Arguments:
- `sp_fine`      : StateSpace at current level
- `model`        : RDMEModel1D or SchloglModel1D
- `hierarchy`    : Vector{VoxelGrid} from current level down to coarsest
- `rates`        : rate parameter vector
- `t`            : current time
- `dt`           : time step

Keyword arguments:
- `τ_pre`             : pre-smoothing time
- `τ_post`            : post-smoothing time
- `krylov_m`          : Krylov subspace dimension
- `weight_tol`        : correction prolongation weight threshold
- `binom_tol`         : correction prolongation binomial threshold
- `expand_coarse`     : expand coarse FSP before solve/recursion
- `coarse_r_depth`    : reaction expansion depth — each coarse voxel count expands by ±r
- `coarse_d_depth`    : diffusion expansion depth — adjacent pairs exchange ±d molecules
- `coarse_n_max`      : per-voxel FSP truncation (scales by 2 each level)
- `max_states`        : emergency abort if |S| exceeds this (default 10^6)
- `pi_tables`         : reserved for future use; currently unused.
"""
function multi_level_vcycle(sp_h::StateSpace{CartesianIndex{K}, Float64},
                             model::Union{RDMEModel1D, SchloglModel1D},
                             hierarchy::Vector{VoxelGrid},
                             rates, t::Real, dt::Real;
                             τ_pre::Real              = 0.0,
                             τ_post::Real             = 0.0,
                             krylov_m::Int            = 30,
                             weight_tol::Float64      = 1e-14,
                             binom_tol::Float64       = 0.0,
                             expand_coarse::Bool      = true,
                             coarse_r_depth::Int      = 2,
                             coarse_d_depth::Int      = 1,
                             coarse_n_max::Int        = 80,
                             max_states::Int          = 1_000_000,
                             prune_tol::Float64       = 0.0,
                             pi_tables::Vector{Vector{Vector{Float64}}} = Vector{Vector{Vector{Float64}}}()) where {K}
    length(sp_h) <= max_states || error("State space exploded: |S| = $(length(sp_h)) > $max_states")

    # ── Coarsest level: solve directly ───────────────────────────────────────
    if length(hierarchy) == 1
        grid_c = hierarchy[1]
        sys_c  = model isa RDMEModel1D ? build_rdme_system(model, grid_c) :
                                         build_schlogl_rdme_system(model, grid_c)
        A_c, = build_generator(sp_h, sys_c, rates, t)
        sp_h.probs .= expv(Float64(dt), A_c, sp_h.probs; m = krylov_m)
        return sp_h
    end

    grid_h  = hierarchy[1]
    grid_2h = hierarchy[2]
    K2 = grid_2h.n_voxels

    # Extract this level's pi_table; fall back to auto-compute if not provided.
    # pi_tables[2:end] are forwarded to the recursive call so each level uses
    # the table built from its own (coarsened) model.
    pi_table = if model isa SchloglModel1D
        if !isempty(pi_tables) && !isempty(pi_tables[1])
            pi_tables[1]
        else
            _compute_schlogl_pi_table(model, grid_h, coarse_n_max ÷ 2)
        end
    else
        Vector{Vector{Float64}}()
    end
    pi_tables_next = length(pi_tables) >= 2 ? pi_tables[2:end] : Vector{Vector{Vector{Float64}}}()

    # ── 1. Pre-smooth ─────────────────────────────────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        smooth_model = RDMEModel1D(model.D, 0.0, 0.0)
        intra_sys = build_intra_system(smooth_model, grid_h, grid_2h)
        A_intra, = build_generator(sp_h, intra_sys, rates, t)
        sp_s = _copy_sp(sp_h)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_h.probs; m = krylov_m)
        sp_s
    else
        sp_h
    end

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    # Filter positive-prob states: zero-prob expand! boundary states contribute
    # nothing to the coarse distribution and can be dropped before restricting.
    sp_for_restrict = _filter_positive(sp_smoothed, weight_tol)
    sp_coarse = restrict(sp_for_restrict)
    n_coarse_pre = length(sp_coarse)
    probs_pre    = copy(sp_coarse.probs)

    # ── 3. Coarse model (diffusion rate scales as 1/dx²) ─────────────────────
    model_2h = coarsen_model(model, 2.0)

    # ── 4. Coarse expand ─────────────────────────────────────────────────────
    # Controlled expansion: reaction steps change individual coarse counts by ±1
    # per step (stoichiometry for Schlögl/birth-death); diffusion steps exchange ±1
    # between adjacent pairs conserving the pair total.  Keeping these separate
    # avoids the combinatorial explosion from treating all events equally.
    if expand_coarse && (coarse_r_depth > 0 || coarse_d_depth > 0)
        _expand_coarse_rd!(sp_coarse, coarse_n_max, Val(K2);
                           r_depth = coarse_r_depth, d_depth = coarse_d_depth)
        length(sp_coarse) <= max_states ||
            error("State space exploded after coarse expand: |S_coarse| = $(length(sp_coarse)) > $max_states")
    end

    # ── 5. Recurse ────────────────────────────────────────────────────────────
    sp_coarse_post = multi_level_vcycle(sp_coarse, model_2h, hierarchy[2:end], rates, t, dt;
                                         τ_pre, τ_post, krylov_m, weight_tol, binom_tol,
                                         expand_coarse, coarse_r_depth, coarse_d_depth,
                                         coarse_n_max = 2 * coarse_n_max,
                                         max_states, prune_tol,
                                         pi_tables = pi_tables_next)

    # ── 6. Prolongate correction ──────────────────────────────────────────────
    # Part A — Covered coarse states (indices 1..n_coarse_pre):
    #   Multiplicative correction — scales existing fine states by p̄_post/p̄_pre.
    #   Symmetric, no new states, no negatives.
    #
    # Part B — New coarse states (expand!-added, indices > n_coarse_pre):
    #   Birth-death: Binomial initialization (Smith & Grima 2016 — exact in fast
    #                diffusion limit when stable states are near nc/2).
    #   Schlögl:     Conditional extension from nearest covered neighbor — preserves
    #                the bimodal distribution structure; Binomial is wrong here since
    #                the stable states (n_low, n_high) are far from nc/2.

    # Build restricted coarse pre/post (only the n_coarse_pre covered states)
    sp_c_pre  = StateSpace{CartesianIndex{K2}, Float64}()
    sp_c_post = StateSpace{CartesianIndex{K2}, Float64}()
    for i in 1:n_coarse_pre
        add_state!(sp_c_pre,  sp_coarse_post.states[i], probs_pre[i])
        add_state!(sp_c_post, sp_coarse_post.states[i], sp_coarse_post.probs[i])
    end

    if model isa SchloglModel1D
        # Schlögl: prolong_conditional handles both parts in one call —
        # Step 1 applies multiplicative correction for covered states (sp_c_pre),
        # Step 3 applies L1 conditional extension for expand!-added new states.
        # The bimodal stable states (n_low, n_high) are far from nc/2 so Binomial
        # would miss them; conditional extension propagates the actual distribution.
        sp_h_new = prolong_conditional(sp_for_restrict, sp_c_pre, sp_coarse_post;
                                       prob_tol = weight_tol)
    else
        # Birth-death: Part A multiplicative + Part B Binomial for new states.
        # Smith & Grima (2016) Eq. (6-7): Binomial(nc, 1/2) is exact in the fast-
        # diffusion limit for convergent (mass-action) propensities.
        sp_h_new = prolong_multiplicative(sp_for_restrict, sp_c_pre, sp_c_post;
                                          prob_tol = weight_tol)
        sp_δ_new = StateSpace{CartesianIndex{K2}, Float64}()
        for i in (n_coarse_pre + 1):length(sp_coarse_post.states)
            p = sp_coarse_post.probs[i]
            p > weight_tol && add_state!(sp_δ_new, sp_coarse_post.states[i], p)
        end
        if length(sp_δ_new) > 0
            sp_δf_new = prolong(sp_δ_new, Val(K); weight_tol, binom_tol, max_states)
            for i in eachindex(sp_δf_new.states)
                s = sp_δf_new.states[i]; δp = sp_δf_new.probs[i]
                idx = get(sp_h_new.index, s, 0)
                if idx == 0; add_state!(sp_h_new, s, δp) else sp_h_new.probs[idx] += δp end
            end
        end
    end

    for i in eachindex(sp_h_new.probs)
        sp_h_new.probs[i] < weight_tol && (sp_h_new.probs[i] = 0.0)
    end

    # ── 7. Post-smooth ────────────────────────────────────────────────────────
    if τ_post > 0.0
        smooth_model = RDMEModel1D(model.D, 0.0, 0.0)
        intra_sys = build_intra_system(smooth_model, grid_h, grid_2h)
        A_intra, = build_generator(sp_h_new, intra_sys, rates, t)
        sp_h_new.probs .= expv(Float64(τ_post), A_intra, sp_h_new.probs; m = krylov_m)
    end

    # ── 8. Prune ──────────────────────────────────────────────────────────────
    # Remove low-probability states created by prolongation so the returned
    # state space doesn't grow without bound across steps.
    prune_tol > 0.0 && prune_threshold!(sp_h_new, prune_tol)

    sp_h_new
end

# ─── 2-species V-cycle ────────────────────────────────────────────────────────

# Coarse expansion for 2-species interleaved states (nA_1,nB_1,…,nA_K,nB_K).
# Reaction expansion: change any single A or B component by ±1.
# Diffusion expansion: adjacent coarse voxels exchange ±1 of the same species.
function _expand_coarse_rd_2s!(sp::StateSpace{CartesianIndex{N2}, T},
                                n_max::Int,
                                ::Val{N2};
                                r_depth::Int = 2,
                                d_depth::Int = 1) where {N2, T}
    N2 % 2 == 0 || error("N2=$N2 must be even (2-species interleaved)")
    K2_vox = N2 ÷ 2

    for _ in 1:r_depth
        frontier = sp.states[begin:end]
        for s in frontier
            t = Tuple(s)
            for ci in 1:N2, δ in (-1, 1)
                nci = t[ci] + δ
                (nci >= 0 && nci <= n_max) || continue
                s_new = CartesianIndex(ntuple(k -> k == ci ? nci : t[k], Val(N2)))
                haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
            end
        end
    end

    K2_vox >= 2 || return
    for _ in 1:d_depth
        frontier = sp.states[begin:end]
        for s in frontier
            t = Tuple(s)
            for j in 1:(K2_vox-1), δ in (-1, 1)
                # A: coarse voxels j (comp 2j-1) and j+1 (comp 2j+1)
                nAj = t[2j-1] + δ; nAj1 = t[2j+1] - δ
                if nAj >= 0 && nAj <= n_max && nAj1 >= 0 && nAj1 <= n_max
                    s_new = CartesianIndex(ntuple(k -> k==2j-1 ? nAj : k==2j+1 ? nAj1 : t[k], Val(N2)))
                    haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
                end
                # B: coarse voxels j (comp 2j) and j+1 (comp 2j+2)
                nBj = t[2j] + δ; nBj1 = t[2j+2] - δ
                if nBj >= 0 && nBj <= n_max && nBj1 >= 0 && nBj1 <= n_max
                    s_new = CartesianIndex(ntuple(k -> k==2j ? nBj : k==2j+2 ? nBj1 : t[k], Val(N2)))
                    haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
                end
            end
        end
    end
end

"""
    multi_level_vcycle_2s(sp_fine, model, hierarchy, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{N}, Float64}

N-level recursive multigrid V-cycle for the 2-species bottleneck RDME.
State layout: interleaved `(nA_1,nB_1,…,nA_K,nB_K)` with `N = 2K`.

Restriction `restrict2s` sums each species separately over fine voxel pairs.
Prolongation uses independent Binomial(n̄_species, 1/2) per species per pair.
"""
function multi_level_vcycle_2s(sp_h::StateSpace{CartesianIndex{N}, Float64},
                                model::BottleneckModel1D,
                                hierarchy::Vector{VoxelGrid},
                                rates, t::Real, dt::Real;
                                τ_pre::Real          = 0.0,
                                τ_post::Real         = 0.0,
                                krylov_m::Int        = 30,
                                weight_tol::Float64  = 1e-14,
                                binom_tol::Float64   = 0.0,
                                expand_coarse::Bool  = true,
                                coarse_r_depth::Int  = 2,
                                coarse_d_depth::Int  = 1,
                                coarse_n_max::Int    = 80,
                                max_states::Int      = 1_000_000,
                                prune_tol::Float64   = 0.0,
                                patch_qsd::Union{PatchQSD, Nothing} = nothing) where {N}
    length(sp_h) <= max_states || error("State space exploded: |S| = $(length(sp_h)) > $max_states")

    # ── Coarsest level: solve directly ───────────────────────────────────────
    if length(hierarchy) == 1
        grid_c = hierarchy[1]
        sys_c  = build_bottleneck_system(model, grid_c)
        A_c, = build_generator(sp_h, sys_c, rates, t)
        sp_h.probs .= expv(Float64(dt), A_c, sp_h.probs; m = krylov_m)
        return sp_h
    end

    grid_h  = hierarchy[1]
    grid_2h = hierarchy[2]
    N2      = N ÷ 2   # coarse state dimension (still interleaved 2-species)

    # ── 1. Pre-smooth (intra-pair diffusion for both species) ─────────────────
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

    # ── 3. Coarsen model ─────────────────────────────────────────────────────
    model_2h = coarsen_model(model, 2.0)

    # ── 4. Coarse expand ─────────────────────────────────────────────────────
    if expand_coarse && (coarse_r_depth > 0 || coarse_d_depth > 0)
        _expand_coarse_rd_2s!(sp_coarse, coarse_n_max, Val(N2);
                              r_depth = coarse_r_depth, d_depth = coarse_d_depth)
        length(sp_coarse) <= max_states ||
            error("State space exploded after coarse expand: |S_coarse| = $(length(sp_coarse))")
    end

    # ── 5. Recurse ────────────────────────────────────────────────────────────
    sp_coarse_post = multi_level_vcycle_2s(sp_coarse, model_2h, hierarchy[2:end], rates, t, dt;
                                            τ_pre, τ_post, krylov_m, weight_tol, binom_tol,
                                            expand_coarse, coarse_r_depth, coarse_d_depth,
                                            coarse_n_max = 2 * coarse_n_max,
                                            max_states, prune_tol, patch_qsd)

    # ── 6. Prolongate ─────────────────────────────────────────────────────────
    sp_c_pre  = StateSpace{CartesianIndex{N2}, Float64}()
    sp_c_post = StateSpace{CartesianIndex{N2}, Float64}()
    for i in 1:n_coarse_pre
        add_state!(sp_c_pre,  sp_coarse_post.states[i], probs_pre[i])
        add_state!(sp_c_post, sp_coarse_post.states[i], sp_coarse_post.probs[i])
    end

    # Part A: multiplicative correction for covered coarse states (always exact)
    sp_h_new = prolong_multiplicative_2s(sp_for_restrict, sp_c_pre, sp_c_post;
                                          prob_tol = weight_tol)

    # Part B: physics-informed QSD (or Binomial fallback) for expand!-added states
    sp_δ_new = StateSpace{CartesianIndex{N2}, Float64}()
    for i in (n_coarse_pre + 1):length(sp_coarse_post.states)
        p = sp_coarse_post.probs[i]
        p > weight_tol && add_state!(sp_δ_new, sp_coarse_post.states[i], p)
    end
    if length(sp_δ_new) > 0
        sp_δf_new = patch_qsd === nothing ?
            prolong2s(sp_δ_new, Val(N); weight_tol, binom_tol, max_states) :
            prolong_patch_qsd(sp_δ_new, Val(N), patch_qsd; weight_tol, prob_tol = weight_tol)
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

    # ── 8. Prune ──────────────────────────────────────────────────────────────
    prune_tol > 0.0 && prune_threshold!(sp_h_new, prune_tol)

    sp_h_new
end

# ─── 2D V-cycle ──────────────────────────────────────────────────────────────

# Coarse expansion for 2D 2-species interleaved states using CoarseningMap.
# Reaction: change any single A or B component by ±1.
# Diffusion: exchange ±1 of same species between adjacent coarse voxels.
function _expand_coarse_2d_2s!(sp::StateSpace{CartesianIndex{N2}, T},
                                n_max::Int,
                                ::Val{N2},
                                adj::Vector{Tuple{Int,Int}};   # adjacency list of coarse voxel pairs
                                r_depth::Int = 2,
                                d_depth::Int = 1) where {N2, T}
    for _ in 1:r_depth
        frontier = sp.states[begin:end]
        for s in frontier
            t = Tuple(s)
            for ci in 1:N2, δ in (-1, 1)
                nci = t[ci] + δ
                (nci >= 0 && nci <= n_max) || continue
                s_new = CartesianIndex(ntuple(k -> k == ci ? nci : t[k], Val(N2)))
                haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
            end
        end
    end

    isempty(adj) && return
    for _ in 1:d_depth
        frontier = sp.states[begin:end]
        for s in frontier
            t = Tuple(s)
            for (J1, J2) in adj, δ in (-1, 1)
                # A species: coarse voxel J1 (comp 2J1-1) and J2 (comp 2J2-1)
                nA1 = t[2J1-1] + δ; nA2 = t[2J2-1] - δ
                if nA1 >= 0 && nA1 <= n_max && nA2 >= 0 && nA2 <= n_max
                    s_new = CartesianIndex(ntuple(k -> k==2J1-1 ? nA1 : k==2J2-1 ? nA2 : t[k], Val(N2)))
                    haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
                end
                # B species: coarse voxel J1 (comp 2J1) and J2 (comp 2J2)
                nB1 = t[2J1] + δ; nB2 = t[2J2] - δ
                if nB1 >= 0 && nB1 <= n_max && nB2 >= 0 && nB2 <= n_max
                    s_new = CartesianIndex(ntuple(k -> k==2J1 ? nB1 : k==2J2 ? nB2 : t[k], Val(N2)))
                    haskey(sp.index, s_new) || add_state!(sp, s_new, zero(T))
                end
            end
        end
    end
end

# Build the coarse-level adjacency list from Grid2D dimensions.
function _coarse_adjacency(K2_x::Int, K2_y::Int)
    adj = Tuple{Int,Int}[]
    for I in 1:K2_x, J in 1:K2_y
        Jc = (I-1)*K2_y + J
        if J < K2_y; push!(adj, (Jc, Jc+1)); end        # horizontal
        if I < K2_x; push!(adj, (Jc, Jc+K2_y)); end     # vertical
    end
    adj
end

"""
    multi_level_vcycle_2d_2s(sp_fine, model, grids, cmaps, rates, t, dt; kwargs...)
        -> StateSpace{CartesianIndex{N}, Float64}

N-level recursive multigrid V-cycle for the 2D two-species bottleneck RDME.

Arguments:
- `grids`  : Vector{Grid2D}, grids[1]=finest, grids[end]=coarsest
- `cmaps`  : Vector{CoarseningMap}, cmaps[ℓ] maps grids[ℓ] → grids[ℓ+1]
- `source_voxels` : fine-level source voxels for production (passed to each level)
"""
function multi_level_vcycle_2d_2s(
    sp_h::StateSpace{CartesianIndex{N}, Float64},
    model::BottleneckModel1D,
    grids::Vector{Grid2D},
    cmaps::Vector{CoarseningMap},
    rates, t::Real, dt::Real;
    τ_pre::Real          = 0.0,
    τ_post::Real         = 0.0,
    krylov_m::Int        = 30,
    weight_tol::Float64  = 1e-14,
    binom_tol::Float64   = 0.0,
    expand_coarse::Bool  = true,
    coarse_r_depth::Int  = 2,
    coarse_d_depth::Int  = 1,
    coarse_n_max::Int    = 40,
    max_states::Int      = 1_000_000,
    prune_tol::Float64   = 0.0,
    source_voxels        = nothing,
    patch_qsd::Union{PatchQSD, Nothing} = nothing
) where {N}
    length(sp_h) <= max_states || error("State space exploded: |S| = $(length(sp_h)) > $max_states")

    grid_h = grids[1]

    # ── Coarsest level: solve directly ───────────────────────────────────────
    if length(grids) == 1
        # At coarsest level, source_voxels has been propagated through coarsening
        # — just build system with whatever source_voxels covers this level
        sys_c = build_bottleneck_system_2d(model, grid_h; source_voxels=source_voxels)
        A_c, = build_generator(sp_h, sys_c, rates, t)
        sp_h.probs .= expv(Float64(dt), A_c, sp_h.probs; m = krylov_m)
        return sp_h
    end

    grid_2h = grids[2]
    cmap    = cmaps[1]
    N2      = 2 * cmap.n_coarse   # coarse state dim: 2 species × coarse voxels
    N2_val  = Val(N2)

    # ── 1. Pre-smooth ─────────────────────────────────────────────────────────
    sp_smoothed = if τ_pre > 0.0
        intra_sys = build_intra_system_2d(model, grid_h, cmap)
        A_intra, = build_generator(sp_h, intra_sys, rates, t)
        sp_s = _copy_sp(sp_h)
        sp_s.probs .= expv(Float64(τ_pre), A_intra, sp_h.probs; m = krylov_m)
        sp_s
    else
        sp_h
    end

    # ── 2. Restrict ───────────────────────────────────────────────────────────
    sp_for_restrict = _filter_positive(sp_smoothed, weight_tol)
    sp_coarse       = restrict2s(sp_for_restrict, cmap, N2_val)
    n_coarse_pre    = length(sp_coarse)
    probs_pre       = copy(sp_coarse.probs)

    # ── 3. Coarsen model (scale k_prod by patch_size for volume ratio) ────────
    model_2h = coarsen_model(model, Float64(cmap.patch_size))

    # ── 4. Coarse expand ──────────────────────────────────────────────────────
    # Only expand when the NEXT call is the coarsest level (length(grids)==2).
    # At deeper intermediate levels the Multinomial prolong over many coarse
    # voxels creates an exponential state-space explosion — skip it there.
    at_penultimate = (length(grids) == 2)
    if expand_coarse && at_penultimate && (coarse_r_depth > 0 || coarse_d_depth > 0)
        adj = _coarse_adjacency(grid_2h.K_x, grid_2h.K_y)
        _expand_coarse_2d_2s!(sp_coarse, coarse_n_max, N2_val, adj;
                               r_depth = coarse_r_depth, d_depth = coarse_d_depth)
        length(sp_coarse) <= max_states ||
            error("Coarse state space exploded: |S| = $(length(sp_coarse))")
    end

    # ── 5. Recurse ────────────────────────────────────────────────────────────
    # Map source_voxels to coarse level
    source_coarse = source_voxels === nothing ? nothing :
                    unique(cmap.fine_to_coarse[k] for k in source_voxels)
    sp_coarse_post = multi_level_vcycle_2d_2s(
        sp_coarse, model_2h, grids[2:end], cmaps[2:end], rates, t, dt;
        τ_pre, τ_post, krylov_m, weight_tol, binom_tol,
        expand_coarse, coarse_r_depth, coarse_d_depth,
        coarse_n_max = cmap.patch_size * coarse_n_max,
        max_states, prune_tol,
        source_voxels = source_coarse, patch_qsd)

    # ── 6. Prolongate ─────────────────────────────────────────────────────────
    sp_c_pre  = StateSpace{CartesianIndex{N2}, Float64}()
    sp_c_post = StateSpace{CartesianIndex{N2}, Float64}()
    for i in 1:n_coarse_pre
        add_state!(sp_c_pre,  sp_coarse_post.states[i], probs_pre[i])
        add_state!(sp_c_post, sp_coarse_post.states[i], sp_coarse_post.probs[i])
    end

    sp_h_new = prolong_multiplicative_2s(sp_for_restrict, sp_c_pre, sp_c_post, cmap, N2_val;
                                          prob_tol = weight_tol)

    sp_δ_new = StateSpace{CartesianIndex{N2}, Float64}()
    for i in (n_coarse_pre+1):length(sp_coarse_post.states)
        p = sp_coarse_post.probs[i]
        p > weight_tol && add_state!(sp_δ_new, sp_coarse_post.states[i], p)
    end
    if length(sp_δ_new) > 0
        sp_δf_new = patch_qsd === nothing ?
            prolong2s(sp_δ_new, cmap, Val(N); weight_tol, binom_tol, max_states) :
            prolong_patch_qsd(sp_δ_new, cmap, Val(N), patch_qsd; weight_tol, prob_tol = weight_tol)
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
        intra_sys = build_intra_system_2d(model, grid_h, cmap)
        A_intra, = build_generator(sp_h_new, intra_sys, rates, t)
        sp_h_new.probs .= expv(Float64(τ_post), A_intra, sp_h_new.probs; m = krylov_m)
    end

    prune_tol > 0.0 && prune_threshold!(sp_h_new, prune_tol)
    sp_h_new
end
