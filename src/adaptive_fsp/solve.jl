"""
    AdaptiveFSP(; kwargs...)

Algorithm parameters for the adaptive Finite State Projection method.

# Keyword Arguments
- `ε_dt::Float64 = 1.0`: Time-step control: `dt = ε_dt / total_flux`.
- `prob_quantile::Float64 = 0.1`: Cumulative probability mass eligible for pruning.
- `flux_tolerance::Float64 = 1e-6`: Flux protection threshold (relative to total flux).
- `expansion_depth::Int = 1`: Number of neighbor shells added per step.
- `save_interval::Int = 1000`: Save a snapshot every N iterations.
"""
Base.@kwdef struct AdaptiveFSP
    ε_dt::Float64 = 1.0
    prob_quantile::Float64 = 0.1
    flux_tolerance::Float64 = 1e-6
    expansion_depth::Int = 1
    save_interval::Int = 1000
end

"""
    solve(prob::FSPProblem, alg::AdaptiveFSP) -> FSPSolution

Run the adaptive FSP algorithm:
1. Expand state space (stoichiometric neighbors).
2. Build/reconstruct the CME generator.
3. Compute adaptive time step from boundary flux.
4. Propagate probability via matrix exponential (`expv`).
5. Flux-aware pruning + renormalization.
6. Repeat until `t >= tf`.
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

    t = Float64(t0)
    iter = 0

    # Generator cache
    A_old = spzeros(Float64, 0, 0)
    gids_old = Int[]

    while t < tf
        iter += 1

        # 1) Expand state space
        expand!(sp, model, bc; depth=alg.expansion_depth)

        # 2) Build/reconstruct generator
        if isempty(gids_old)
            A, in_flow, out_flow = build_generator(sp, model, rates, t)
        else
            A, in_flow, out_flow = reconstruct_generator(sp, model, rates, t, A_old, gids_old)
        end

        # 3) Adaptive time step from total outflux
        out_flux = Vector(-out_flow * sp.probs)
        total_flux = sum(out_flux)
        dt = total_flux > 0 ? alg.ε_dt / total_flux : (tf - t)
        dt = min(dt, tf - t)

        # 4) Propagate probability
        sp.probs = expmv(dt, A, sp.probs)

        # 5) Flux-aware pruning + renormalization
        compress!(sp, model, rates, t, alg.prob_quantile;
                  flux_tolerance=alg.flux_tolerance)
        renormalize!(sp)

        # 6) Cache for next reconstruction
        gids_old = get_global_ids(sp)
        A_old = A

        # 7) Advance time
        t += dt

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

    FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
end
