"""
    RDMEMultigridFSP

Algorithm parameters for the multi-level multigrid RDME FSP solver.

Fields:
- `dt`                 : fixed time step per V-cycle
- `τ_pre`              : pre-smoothing time  (0 → skip, set to log(1/ε)/(2d) for ε-smoothing)
- `τ_post`             : post-smoothing time
- `n_max`              : FSP truncation: max molecules per fine voxel
- `krylov_m`           : Krylov subspace dimension for all expmv calls
- `weight_tol`         : drop prolongated fine states below this weight
- `binom_tol`          : drop per-pair Binomial fraction < binom_tol in prolongation
- `prune_tol`          : post-vcycle: prune fine states with prob < prune_tol
- `expand_coarse`      : expand coarse state space before coarse solve (recommended: true)
- `coarse_r_depth`     : coarse reaction expansion depth (±r in each coarse voxel count)
- `coarse_d_depth`     : coarse diffusion expansion depth (adjacent pairs exchange ±d)
- `expand_depth`       : post-vcycle fine expansion depth (0 = skip; use coarse expansion instead)
- `coarse_only`        : work entirely at coarse level; prolong only at snapshot times (recommended
                         for large K where fine state space is intractable between V-cycles)
- `n_levels`           : number of coarsening levels (1 = K→K/2, 2 = K→K/4, …)
- `coarse_n_total`     : if > 0, also enforce total-molecule FSP boundary (Σnk ≤ coarse_n_total)
                         at the coarsest level in addition to per-voxel.  Cuts corner states like
                         (n_max, n_max, 0, 0) that are physically negligible.  0 = disabled.
- `coarse_n_max_per_voxel` : if > 0, override the default coarsest-level per-voxel bound
                         (which is 2^n_levels * n_max) with this tighter value.
- `max_states`         : emergency abort if |S| exceeds this (default 10^6)
- `save_every`         : save a snapshot every this many steps
"""
Base.@kwdef struct RDMEMultigridFSP
    dt::Float64               = 0.01
    τ_pre::Float64            = 0.0
    τ_post::Float64           = 0.0
    n_max::Int                = 40
    krylov_m::Int             = 30
    prune_tol::Float64        = 1e-12
    expand_coarse::Bool       = true
    coarse_r_depth::Int       = 2
    coarse_d_depth::Int       = 1
    expand_depth::Int         = 0
    coarse_only::Bool              = false
    n_levels::Int                  = 1
    coarse_n_total::Int            = 0
    coarse_n_max_per_voxel::Int    = 0
    max_states::Int                = 1_000_000
    save_every::Int                = 1
    weight_tol::Float64       = 1e-13
    binom_tol::Float64        = 1e-10
end

"""
    RDMESolution

Stores the RDME FSP solution: probability distributions at each saved time point.

Fields:
- `t`           : saved time points
- `snapshots`   : vector of (states, probs) tuples at each saved time
- `K`           : number of fine voxels
- `state_space_sizes` : |S| at each saved time
"""
struct RDMESolution{K}
    t::Vector{Float64}
    snapshots::Vector{Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}}
    state_space_sizes::Vector{Int}
end

"""
    solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)
        -> RDMESolution

Solve the 1D RDME using a multi-level multigrid FSP.

Arguments:
- `model`     : `RDMEModel1D` or `SchloglModel1D`
- `fine_grid` : `VoxelGrid` with `n_voxels = K`
- `u0`        : initial state `CartesianIndex{K}`
- `tspan`     : `(t0, tf)` time interval
- `rates`     : rate parameter vector
- `alg`       : `RDMEMultigridFSP` parameters
"""
function solve_rdme_multigrid(model::Union{RDMEModel1D, SchloglModel1D},
                               fine_grid::VoxelGrid,
                               u0::CartesianIndex{K},
                               tspan::Tuple{<:Real,<:Real},
                               rates,
                               alg::RDMEMultigridFSP) where {K}
    if alg.coarse_only
        return _solve_rdme_coarse_only(model, fine_grid, u0, tspan, rates, alg)
    end

    t0, tf = Float64.(tspan)

    # ── Build grid hierarchy ──────────────────────────────────────────────────
    hierarchy = build_hierarchy(fine_grid, alg.n_levels)


    # ── Build full RDME system (used for expansion) ───────────────────────────
    full_sys = model isa RDMEModel1D ? build_rdme_system(model, fine_grid) :
                                       build_schlogl_rdme_system(model, fine_grid)
    bc = state -> rdme_bc(state, alg.n_max)

    # ── Initialise state space ────────────────────────────────────────────────
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, full_sys, bc; depth = 1)

    t    = t0
    step = 0

    # Storage for snapshots
    saved_t    = Float64[]
    snapshots  = Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}[]
    ss_sizes   = Int[]

    # Save initial snapshot
    _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp)

    # ── Main time-stepping loop ───────────────────────────────────────────────
    while t < tf
        dt_step = min(alg.dt, tf - t)

        sp = multi_level_vcycle(sp, model, hierarchy, rates, t, dt_step;
                                 τ_pre              = alg.τ_pre,
                                 τ_post             = alg.τ_post,
                                 krylov_m           = alg.krylov_m,
                                 weight_tol         = alg.weight_tol,
                                 binom_tol          = alg.binom_tol,
                                 expand_coarse      = alg.expand_coarse,
                                 coarse_r_depth     = alg.coarse_r_depth,
                                 coarse_d_depth     = alg.coarse_d_depth,
                                 coarse_n_max       = 2 * alg.n_max,
                                 max_states         = alg.max_states,
                                 prune_tol          = alg.prune_tol,
)

        t    += dt_step
        step += 1

        # Renormalize and prune negligible states
        renormalize!(sp)
        prune_threshold!(sp, alg.prune_tol)

        # Optionally expand to keep FSP boundary ahead of support
        if alg.expand_depth > 0
            expand!(sp, full_sys, bc; depth = alg.expand_depth)
        end

        # Save snapshot
        if step % alg.save_every == 0 || t ≥ tf
            _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp)
        end
    end

    RDMESolution{K}(saved_t, snapshots, ss_sizes)
end

"""
    _solve_rdme_coarse_only(model, fine_grid, u0, tspan, rates, alg) -> RDMESolution

Coarse-only variant of the multigrid solver.
"""
function _solve_rdme_coarse_only(model::Union{RDMEModel1D, SchloglModel1D},
                                  fine_grid::VoxelGrid,
                                  u0::CartesianIndex{K},
                                  tspan::Tuple{<:Real,<:Real},
                                  rates,
                                  alg::RDMEMultigridFSP) where {K}
    iseven(K) || error("K=$K must be even for multigrid")
    t0, tf = Float64.(tspan)
    n      = alg.n_levels

    # ── Build coarsest grid and system ────────────────────────────────────────
    coarsest_grid = fine_grid
    for _ in 1:n
        coarsest_grid = coarsen(coarsest_grid)
    end
    default_n_max  = (2^n) * alg.n_max
    coarsest_n_max = alg.coarse_n_max_per_voxel > 0 ? alg.coarse_n_max_per_voxel : default_n_max
    
    # Scale model for coarsest level
    coarsest_model = model
    for _ in 1:n
        coarsest_model = coarsen_model(coarsest_model, 2.0)
    end
    
    coarsest_sys = coarsest_model isa RDMEModel1D ? build_rdme_system(coarsest_model, coarsest_grid) :
                                                    build_schlogl_rdme_system(coarsest_model, coarsest_grid)

    # Boundary condition: per-voxel (possibly tightened) + optional total-molecule
    coarsest_bc = if alg.coarse_n_total > 0
        state -> rdme_bc(state, coarsest_n_max) && sum(Tuple(state)) ≤ alg.coarse_n_total
    else
        state -> rdme_bc(state, coarsest_n_max)
    end

    # ── Initialise coarsest state space from fine IC ──────────────────────────
    sp_fine0 = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp_fine0, u0, 1.0)
    sp_c = _multi_restrict(sp_fine0, Val(n))

    t    = t0
    step = 0

    # Storage for fine snapshots
    saved_t   = Float64[]
    snapshots = Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}[]
    ss_sizes  = Int[]

    # ── Save initial snapshot ─────────────────────────────────────────────────
    sp_fine_snap = _multi_prolong(sp_c, Val(K), Val(n);
                                  weight_tol=alg.weight_tol, binom_tol=alg.binom_tol)
    _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp_fine_snap)

    # ── Main loop: work entirely at coarsest level ────────────────────────────
    while t < tf
        dt_step = min(alg.dt, tf - t)

        # Expand coarsest FSP before solve (keep boundary ahead of support)
        if alg.expand_coarse && alg.coarse_r_depth > 0
            expand!(sp_c, coarsest_sys, coarsest_bc; depth = alg.coarse_r_depth)
        end

        # Build generator and evolve
        A_c, = build_generator(sp_c, coarsest_sys, rates, t)
        sp_c.probs .= expv(Float64(dt_step), A_c, sp_c.probs; m = alg.krylov_m)

        t    += dt_step
        step += 1

        # Renormalize and prune coarsest states
        renormalize!(sp_c)
        prune_threshold!(sp_c, alg.prune_tol)

        # Save snapshot: prolong all levels to fine, save, discard fine
        if step % alg.save_every == 0 || t ≥ tf
            sp_fine_snap = _multi_prolong(sp_c, Val(K), Val(n);
                                          weight_tol = alg.weight_tol,
                                          binom_tol  = alg.binom_tol)
            _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp_fine_snap)
        end
    end

    RDMESolution{K}(saved_t, snapshots, ss_sizes)
end

# ─── multi-level restrict / prolong helpers ───────────────────────────────────

_multi_restrict(sp::StateSpace, ::Val{0}) = sp
_multi_restrict(sp::StateSpace, ::Val{N}) where N = _multi_restrict(restrict(sp), Val(N-1))

_multi_prolong(sp, ::Val{K}, ::Val{0}; kw...) where K = sp
function _multi_prolong(sp, ::Val{K}, ::Val{N}; max_states=10_000_000, kw...) where {K, N}
    # Prolong current level (K / 2^N) to next level (K / 2^(N-1))
    sp_next = prolong(sp, Val(K ÷ 2^(N-1)); max_states=max_states, kw...)
    _multi_prolong(sp_next, Val(K), Val(N-1); max_states=max_states, kw...)
end

# ─── internal helpers ─────────────────────────────────────────────────────────

function _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp)
    push!(saved_t,   t)
    push!(snapshots, (copy(sp.states), copy(sp.probs)))
    push!(ss_sizes,  length(sp))
end

# ─── accessors ───────────────────────────────────────────────────────────────

function mean_voxel_counts(sol::RDMESolution{K}) where {K}
    n_t = length(sol.t)
    mu  = zeros(Float64, n_t, K)
    for (ti, (states, probs)) in enumerate(sol.snapshots)
        for (state, p) in zip(states, probs)
            t_state = Tuple(state)
            for k in 1:K
                mu[ti, k] += p * t_state[k]
            end
        end
    end
    mu
end

function marginal_voxel(sol::RDMESolution{K}, ti::Int, k::Int) where {K}
    states, probs = sol.snapshots[ti]
    count_map = Dict{Int, Float64}()
    for (state, p) in zip(states, probs)
        n = Tuple(state)[k]
        count_map[n] = get(count_map, n, 0.0) + p
    end
    ns = sort!(collect(keys(count_map)))
    ps = [count_map[n] for n in ns]
    ns, ps
end
