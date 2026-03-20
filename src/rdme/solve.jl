"""
    RDMEMultigridFSP

Algorithm parameters for the two-level multigrid RDME FSP solver.

Fields:
- `dt`          : fixed time step per V-cycle
- `τ_pre`       : pre-smoothing time  (0 → skip, set to log(1/ε)/(2d) for ε-smoothing)
- `τ_post`      : post-smoothing time
- `n_max`       : FSP truncation: max molecules per fine voxel
- `krylov_m`    : Krylov subspace dimension for all expmv calls
- `weight_tol`  : drop prolongated fine states below this weight
- `save_every`  : save a snapshot every this many steps
"""
Base.@kwdef struct RDMEMultigridFSP
    dt::Float64          = 0.01
    τ_pre::Float64       = 0.0
    τ_post::Float64      = 0.0
    n_max::Int           = 40
    krylov_m::Int        = 30
    weight_tol::Float64  = 1e-14
    save_every::Int      = 1
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
                               τ_pre      = alg.τ_pre,
                               τ_post     = alg.τ_post,
                               krylov_m   = alg.krylov_m,
                               weight_tol = alg.weight_tol)

        t    += dt_step
        step += 1

        # Renormalize and prune negligible states
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        # Expand by one shell to keep FSP boundary ahead of support
        expand!(sp, full_sys, bc; depth = 1)

        # Save snapshot
        if step % alg.save_every == 0 || t ≥ tf
            _save_snapshot!(saved_t, snapshots, ss_sizes, t, sp)
        end
    end

    RDMESolution{K}(saved_t, snapshots, ss_sizes)
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
