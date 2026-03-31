"""
    RDMEMultigridFSP

Algorithm parameters for the two-level multigrid RDME FSP solver.

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
- `coarse_expand_depth`: stoichiometric shell depth for coarse expansion
- `expand_depth`       : post-vcycle fine expansion depth (0 = skip; use coarse expansion instead)
- `coarse_only`        : work entirely at coarse level; prolong only at snapshot times (recommended
                         for large K where fine state space is intractable between V-cycles)
- `n_levels`           : number of coarsening levels for coarse_only mode (1 = K→K/2, 2 = K→K/4, …)
- `coarse_n_total`     : if > 0, also enforce total-molecule FSP boundary (Σnk ≤ coarse_n_total)
                         at the coarsest level in addition to per-voxel.  Cuts corner states like
                         (n_max, n_max, 0, 0) that are physically negligible.  0 = disabled.
- `coarse_n_max_per_voxel` : if > 0, override the default coarsest-level per-voxel bound
                         (which is 2^n_levels * n_max) with this tighter value.  Useful when
                         prune_tol is loose enough that counts above this value are always pruned
                         anyway (avoids wasteful expand-then-prune cycling of boundary states).
                         Rule of thumb: set to stationary_mean_per_coarse_voxel + 4*sqrt(mean).
                         0 = disabled (use default 2^n_levels * n_max).
- `save_every`         : save a snapshot every this many steps
"""
Base.@kwdef struct RDMEMultigridFSP
    dt::Float64               = 0.01
    τ_pre::Float64            = 0.0
    τ_post::Float64           = 0.0
    n_max::Int                = 40
    krylov_m::Int             = 30
    weight_tol::Float64       = 1e-14   # prolongation: skip coarse states with prob < weight_tol
    binom_tol::Float64        = 0.0     # prolongation: skip per-pair Binomial fraction < binom_tol
    prune_tol::Float64        = 1e-12   # post-vcycle: prune fine states with prob < prune_tol
    expand_coarse::Bool       = true    # expand coarse FSP before coarse solve (recommended)
    coarse_expand_depth::Int  = 1       # shell depth for coarse expansion
    expand_depth::Int         = 0       # post-vcycle fine expansion depth (0 = skip)
    coarse_only::Bool              = false   # work entirely at coarse level; prolong only at snapshots
    n_levels::Int                  = 1       # coarsening levels for coarse_only (1=K→K/2, 2=K→K/4)
    coarse_n_total::Int            = 0       # additional total-molecule FSP bound (0 = disabled)
    coarse_n_max_per_voxel::Int    = 0       # override per-voxel bound at coarsest level (0 = default)
    save_every::Int                = 1
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

Solve the 1D RDME using a two-level multigrid FSP.

Arguments:
- `model`     : `RDMEModel1D`
- `fine_grid` : `VoxelGrid` with `n_voxels = K` (must be even)
- `u0`        : initial state `CartesianIndex{K}` (all probability starts here)
- `tspan`     : `(t0, tf)` time interval
- `rates`     : rate parameter vector (passed to propensity functions)
- `alg`       : `RDMEMultigridFSP` parameters

The solver:
1. Initialises `StateSpace{CartesianIndex{K}, Float64}` with `u0` and expands
   by one stoichiometric shell (all reachable states in one reaction/diffusion step).
2. Repeatedly applies `two_level_vcycle` with step `dt` until `tf` is reached.
3. After each cycle: renormalizes, prunes states with probability < 1e-12,
   then expands by one shell.
4. Saves a snapshot every `save_every` steps.
"""
function solve_rdme_multigrid(model::RDMEModel1D,
                               fine_grid::VoxelGrid,
                               u0::CartesianIndex{K},
                               tspan::Tuple{<:Real,<:Real},
                               rates,
                               alg::RDMEMultigridFSP) where {K}
    if alg.coarse_only
        return _solve_rdme_coarse_only(model, fine_grid, u0, tspan, rates, alg)
    end

    t0, tf = Float64.(tspan)
    coarse_grid = coarsen(fine_grid)

    # ── Build full RDME system (used for expansion) ───────────────────────────
    full_sys = build_rdme_system(model, fine_grid)
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

        sp = two_level_vcycle(sp, model, fine_grid, coarse_grid, rates, t, dt_step;
                               τ_pre              = alg.τ_pre,
                               τ_post             = alg.τ_post,
                               krylov_m           = alg.krylov_m,
                               weight_tol         = alg.weight_tol,
                               binom_tol          = alg.binom_tol,
                               expand_coarse      = alg.expand_coarse,
                               coarse_expand_depth = alg.coarse_expand_depth,
                               coarse_n_max       = 2 * alg.n_max)

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

Works entirely at the K/2^n_levels-voxel coarsest level: expand → solve → prune.
Prolongation to fine is performed **only at snapshot times** and the fine
state space is discarded immediately after saving; it never accumulates
between time steps.  This keeps memory O(|S_coarse|) regardless of how
many steps are taken, making it tractable for large K where the fine state
space would otherwise explode.

`alg.n_levels` controls depth: 1 = K→K/2 (default), 2 = K→K/4, etc.

The snapshots stored in `RDMESolution{K}` contain the prolongated fine
distributions so all post-processing functions (`mean_voxel_counts`, etc.)
work unchanged.
"""
function _solve_rdme_coarse_only(model::RDMEModel1D,
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
    coarsest_sys   = build_coarse_system(model, coarsest_grid, fine_grid)

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
        if alg.expand_coarse && alg.coarse_expand_depth > 0
            expand!(sp_c, coarsest_sys, coarsest_bc; depth = alg.coarse_expand_depth)
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

_multi_restrict(sp::StateSpace{CartesianIndex{K}, Float64}, ::Val{1}) where K =
    restrict(sp)

_multi_restrict(sp::StateSpace{CartesianIndex{K}, Float64}, ::Val{2}) where K =
    restrict(restrict(sp))

_multi_restrict(sp::StateSpace{CartesianIndex{K}, Float64}, ::Val{3}) where K =
    restrict(restrict(restrict(sp)))

function _multi_prolong(sp, ::Val{K}, ::Val{1}; kw...) where K
    prolong(sp, Val(K); kw...)
end

function _multi_prolong(sp, ::Val{K}, ::Val{2}; kw...) where K
    prolong(prolong(sp, Val(K÷2); kw...), Val(K); kw...)
end

function _multi_prolong(sp, ::Val{K}, ::Val{3}; kw...) where K
    prolong(prolong(prolong(sp, Val(K÷4); kw...), Val(K÷2); kw...), Val(K); kw...)
end

# ─── internal helpers ─────────────────────────────────────────────────────────

function _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp)
    push!(saved_t,   t)
    push!(snapshots, (copy(sp.states), copy(sp.probs)))
    push!(ss_sizes,  length(sp))
end

# ─── accessors ───────────────────────────────────────────────────────────────

"""
    mean_voxel_counts(sol) -> Matrix{Float64}

Return a `(n_snapshots × K)` matrix of mean molecule counts per voxel at each
saved time.
"""
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

"""
    marginal_voxel(sol, ti, k) -> (counts, probs)

Marginal distribution of molecule count in voxel `k` at snapshot index `ti`.
"""
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
