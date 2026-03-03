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
- `flux_method::Symbol = :total`: Time-step denominator: `:total` uses
  `Φ_total = Σᵢ pᵢ·wᵢ` (probability-weighted sum); `:maximum` uses
  `max_i(pᵢ·wᵢ)` (most active state). Use `:maximum` for stiff systems
  where a single high-rate state would otherwise dominate the sum.
- `expand_method::Symbol = :stoich`: State space expansion strategy per step.
  `:stoich` adds all stoichiometric neighbors (`expand!`); `:ssa` adds one
  SSA-sampled successor per state (`expand_ssa!`). Use `:ssa` for oscillatory
  systems (e.g. Oregonator) where `:stoich` inflates the state space.

## Time-step formula

    dt = min(ε_dt / Φ,  tf - t)

where `Φ = Σᵢ pᵢ·wᵢ` (`:total`) or `Φ = maxᵢ pᵢ·wᵢ` (`:maximum`).
"""
Base.@kwdef struct AdaptiveFSP
    ε_dt::Float64 = 1.0
    prob_quantile::Float64 = 0.1
    flux_tolerance::Float64 = 1e-6
    expansion_depth::Int = 1
    save_interval::Int = 1000
    max_iter::Int = typemax(Int)
    flux_method::Symbol = :total
    expand_method::Symbol = :stoich
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
3. Compute adaptive time step: `dt = ε_dt / Φ` where `Φ` is total or maximum outflux
   depending on `alg.flux_method`.
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

        # 3) Adaptive time step
        out_flux = Vector(-out_flow * sp.probs)
        Φ = alg.flux_method === :maximum ? maximum(out_flux) : sum(out_flux)
        dt = Φ > 0 ? alg.ε_dt / Φ : (tf - t)
        dt = min(dt, tf - t)

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
