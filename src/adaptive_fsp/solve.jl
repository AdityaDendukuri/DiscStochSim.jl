"""
    AdaptiveFSP(; kwargs...)

Algorithm parameters for the adaptive Finite State Projection method.

# Keyword Arguments
- `ε_dt::Float64 = 1.0`: Time-step control parameter.
- `prob_quantile::Float64 = 0.1`: Quantile parameter used by pruning.
- `flux_tolerance::Float64 = 1e-6`: Flux protection threshold (relative to total flux).
- `expansion_depth::Int = 1`: Number of neighbor shells added per step.
- `save_interval::Int = 1000`: Save a snapshot every N iterations.
- `max_iter::Int = typemax(Int)`: Maximum number of iterations.
- `expand_method::Symbol = :stoich`: State space expansion strategy per step.
  `:stoich` adds all stoichiometric neighbors (`expand!`); `:ssa` adds one
  SSA-sampled successor per state (`expand_ssa!`); `:flux` adds all neighbors
  where `α_k(x)/w(x) > expand_threshold` (`expand_flux!`). Use `:ssa` or
  `:flux` for oscillatory systems (e.g. Oregonator) where `:stoich` inflates
  the state space. `:flux` is the deterministic version of `:ssa`.
- `expand_threshold::Float64 = 0.1`: Fractional propensity threshold for
  `:flux` expansion. Reactions with `α_k(x)/w(x) > expand_threshold` add
  their neighbor to J. Ignored by `:stoich` and `:ssa`.
- `dt_max::Float64 = Inf`: Per-depth-unit time step cap: `1/λ_downstream`,
  where `λ_downstream` is the exit rate of the fast downstream states first
  reached by expansion. The actual cap applied is `expansion_depth * dt_max`,
  so increasing expansion depth automatically relaxes the cap proportionally.
  Required for slow-start systems where `Φ_out → 0` at `t=0` (e.g. the
  Bottleneck system), which would otherwise produce a single giant step.
  Set `dt_max = 1/λ_downstream`; leave as `Inf` if the flux step is
  already sufficient (κ = λ_max/Φ_out ≤ 1/ε_dt throughout).

## Time-step formula

    dt = min(ε_dt / Φ,  expansion_depth * dt_max,  tf - t)

where `Φ = Σᵢ pᵢ·wᵢ` is the probability-weighted internal exit rate (compressed
generator diagonal). For stoichiometric expansion all neighbors are in J so Φ ≈ Φ_total
and the bound `‖ε‖₁ ≲ Φ_total·dt ≈ Φ·dt = ε_dt` (Corollary 4.7) holds directly.
For SSA expansion (e.g. Oregonator) Φ < Φ_total, but Φ still tracks instantaneous
system activity and gives empirically good adaptive control.
"""
Base.@kwdef struct AdaptiveFSP
    ε_dt::Float64 = 1.0
    prob_quantile::Float64 = 0.1
    flux_tolerance::Float64 = 1e-6
    expansion_depth::Int = 1
    save_interval::Int = 1000
    max_iter::Int = typemax(Int)
    expand_method::Symbol = :stoich
    expand_threshold::Float64 = 0.1
    dt_max::Float64 = Inf
    # Optional progress callback: called every save_interval steps with
    # (t_current, t_log, mean_trajectory_matrix).  Use terminal_progress()
    # to get a live UnicodePlots display in the terminal.
    progress_callback = nothing
end

"""
    AdaptiveFSPDiagnostics

Per-step diagnostic data from an AdaptiveFSP solve.

- `dt_log`: time step taken at each iteration.
- `t_log`: simulated time after each iteration.
- `size_log`: projection size `|S|` after each iteration.
- `flux_log`: internal flux `Φ = Σᵢ pᵢ·wᵢ` (compressed generator diagonal) at each iteration.
- `total_iters`: total number of iterations executed.
"""
struct AdaptiveFSPDiagnostics
    dt_log::Vector{Float64}
    t_log::Vector{Float64}
    size_log::Vector{Int}
    flux_log::Vector{Float64}
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
    dt_log        = Float64[]
    t_log         = Float64[]
    size_log_diag = Int[]
    flux_log      = Float64[]

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
        elseif alg.expand_method === :flux
            expand_flux!(sp, model, rates, t, bc;
                         depth=alg.expansion_depth, threshold=alg.expand_threshold)
        else
            expand!(sp, model, bc; depth=alg.expansion_depth)
        end

        # 2) Build generator (incremental rebuild reusing previous matrix)
        A, in_flow, out_flow, exit_rate_boundary =
            reconstruct_generator(sp, model, rates, t, A_old, gids_old)

        # Cache generator state for next iteration's incremental rebuild
        A_old = A
        gids_old = get_global_ids(sp)

        # 3) Adaptive time step: dt = min(ε_dt/Φ_total, d·dt_max, tf-t)
        # Compute Φ_total = Σᵢ p[i]·w(xᵢ) directly from the full exit rate w(xᵢ),
        # not from the compressed generator diagonal.
        Φ_total = zero(Float64)
        @inbounds for i in eachindex(sp.probs)
            w = zero(Float64)
            for prop in model.propensities
                w += prop(sp.states[i], rates, t)
            end
            Φ_total += sp.probs[i] * w
        end
        dt_flux = Φ_total > 0 ? alg.ε_dt / Φ_total : (tf - t)
        dt = min(dt_flux, alg.expansion_depth * alg.dt_max, tf - t)

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
        push!(flux_log, Φ_total)

        # 8) Save snapshot
        if iter % alg.save_interval == 0 || t >= tf
            push!(sol_t, t)
            push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
            push!(sol_sizes, length(sp))
            @info "t = $(round(t; sigdigits=4)), |S| = $(length(sp)), dt = $(round(dt; sigdigits=3))"
            if !isnothing(alg.progress_callback)
                partial = FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
                extra = (step_t=t_log, dt_log=dt_log, size_log=size_log_diag, flux_log=flux_log)
                alg.progress_callback(t, sol_t, mean_trajectory(partial), extra)
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
    diag = AdaptiveFSPDiagnostics(dt_log, t_log, size_log_diag, flux_log, iter)
    (sol, diag)
end
