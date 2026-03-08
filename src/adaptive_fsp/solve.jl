"""
    AdaptiveFSP(; kwargs...)

Algorithm parameters for the adaptive Finite State Projection method.

# Keyword Arguments
- `ε_dt::Float64 = 1.0`: Time-step control numerator.
- `prob_quantile::Float64 = 0.1`: Quantile parameter used by pruning.
- `flux_tolerance::Float64 = 1e-6`: Flux protection threshold (relative to total flux).
- `expansion_depth::Int = 1`: Number of neighbor shells added per step.
- `save_interval::Int = 1000`: Save a snapshot every N iterations.
- `max_iter::Int = typemax(Int)`: Maximum number of iterations.
- `expand_method::Symbol = :stoich`: State space expansion strategy per step.
  `:stoich` adds all stoichiometric neighbors (`expand!`); `:ssa` adds one
  SSA-sampled successor per state (`expand_ssa!`). Use `:ssa` for oscillatory
  systems (e.g. Oregonator) where `:stoich` inflates the state space.
- `dt_max::Float64 = Inf`: Hard upper bound on the time step. Required for
  slow-start systems where `Φ_total → 0` at `t=0` (e.g. the Bottleneck system
  with tiny initial flux), which would otherwise produce a single giant step
  that the depth-1 expansion cannot contain. Set `dt_max ≈ 1/λ_downstream`
  where `λ_downstream` is the rate of the fast downstream reaction.

## Time-step formula

    dt = min(ε_dt / Φ_total,  dt_max,  tf - t)

where `Φ_total = Σᵢ pᵢ·wᵢ` is the total probability-weighted outflux.
With `ε_dt = 1.0` (the default), this gives approximately one expected
reaction per step — the natural unit for CME integration.
"""
Base.@kwdef struct AdaptiveFSP
    ε_dt::Float64 = 1.0
    prob_quantile::Float64 = 0.1
    flux_tolerance::Float64 = 1e-6
    expansion_depth::Int = 1
    save_interval::Int = 1000
    max_iter::Int = typemax(Int)
    expand_method::Symbol = :stoich
    dt_max::Float64 = Inf
    # Optional progress callback: called every save_interval steps with
    # (t_current, t_log, mean_trajectory_matrix).  Use terminal_progress()
    # to get a live UnicodePlots display in the terminal.
    progress_callback = nothing
    # Deprecated: flux_method no longer used; dt = ε_dt/Φ_total always.
    # Retained for backward compatibility with existing call sites.
    flux_method::Symbol = :total
end

"""
    AdaptiveFSPDiagnostics

Per-step diagnostic data from an AdaptiveFSP solve.

- `dt_log`: time step taken at each iteration.
- `t_log`: simulated time after each iteration.
- `size_log`: projection size `|S|` after each iteration.
- `total_iters`: total number of iterations executed.
"""
struct AdaptiveFSPDiagnostics
    dt_log::Vector{Float64}
    t_log::Vector{Float64}
    size_log::Vector{Int}
    total_iters::Int
end

"""
    solve(prob::FSPProblem, alg::AdaptiveFSP) -> (FSPSolution, AdaptiveFSPDiagnostics)

Run the adaptive FSP algorithm:
1. Expand state space (via `alg.expand_method`: `:stoich` = all neighbors, `:ssa` = SSA-sampled).
2. Build the CME generator (incremental).
3. Compute adaptive time step: `dt = ε_dt / Φ_total`.
4. Propagate probability via matrix exponential (`expv`).
5. Flux-aware pruning + renormalization.
6. Repeat until `t >= tf` or `alg.max_iter` iterations reached.

Returns `(solution, diagnostics)` where `diagnostics` contains per-step
`dt_log`, `t_log`, and `size_log`.
"""
function CommonSolve.solve(prob::FSPProblem{E,T}, alg::AdaptiveFSP) where {E,T}
    model = prob.model
    rates = prob.rates
    bc = prob.bc
    t0, tf = prob.tspan
    n_species = length(Tuple(prob.u0))

    # Initialize state space
    sp = StateSpace{E, Float64}()
    add_state!(sp, prob.u0, 1.0)

    # Solution storage
    sol_t = Float64[t0]
    sol_snaps = [(copy(sp.states), copy(sp.probs))]
    sol_sizes = [length(sp)]

    # Per-step diagnostics (lightweight: no state copies)
    dt_log   = Float64[]
    t_log    = Float64[]
    size_log_diag = Int[]

    t = Float64(t0)
    iter = 0

    # Track previous generator for incremental rebuilds
    A_old = spzeros(Float64, 0, 0)
    gids_old = Int[]

    while t < tf && iter < alg.max_iter
        iter += 1

        # 1) Expand state space
        if alg.expand_method === :ssa
            expand_ssa!(sp, model, rates, t, bc; depth=alg.expansion_depth)
        else
            expand!(sp, model, bc; depth=alg.expansion_depth)
        end

        # 2) Build generator (incremental rebuild reusing previous matrix)
        A, in_flow, out_flow = reconstruct_generator(sp, model, rates, t, A_old, gids_old)

        # Cache generator state for next iteration's incremental rebuild
        A_old = A
        gids_old = get_global_ids(sp)

        # 3) Adaptive time step: dt = min(ε_dt/Φ_total, dt_max)
        out_flux = Vector(-out_flow * sp.probs)
        Φ_total  = sum(out_flux)
        dt_flux  = Φ_total > 0 ? alg.ε_dt / Φ_total : (tf - t)
        dt = min(dt_flux, alg.dt_max, tf - t)

        # 4) Propagate probability (subcycle if dt*||A||_1 is too large)
            sp.probs = expv(dt, A, sp.probs)

        # 5) Flux-aware pruning + renormalization
        compress!(sp, model, rates, t, alg.prob_quantile;
                  flux_tolerance=alg.flux_tolerance)
        renormalize!(sp)

        # 6) Advance time
        t += dt

        # 7) Record per-step diagnostics
        push!(dt_log, dt)
        push!(t_log, t)
        push!(size_log_diag, length(sp))

        # 8) Save snapshot
        if iter % alg.save_interval == 0 || t >= tf
            push!(sol_t, t)
            push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
            push!(sol_sizes, length(sp))
            @info "t = $(round(t; sigdigits=4)), |S| = $(length(sp)), dt = $(round(dt; sigdigits=3))"
            if !isnothing(alg.progress_callback)
                partial = FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
                alg.progress_callback(t, sol_t, mean_trajectory(partial))
            end
        end
    end

    # Ensure final snapshot is captured
    if sol_t[end] < t
        push!(sol_t, t)
        push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
        push!(sol_sizes, length(sp))
    end

    sol  = FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
    diag = AdaptiveFSPDiagnostics(dt_log, t_log, size_log_diag, iter)
    (sol, diag)
end
