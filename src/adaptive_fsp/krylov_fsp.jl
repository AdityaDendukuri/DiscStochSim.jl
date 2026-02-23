"""
    KrylovFSP(; kwargs...)

Krylov-FSP-SSA baseline algorithm (Sidje & Vo 2015, reviewed in Dinh & Sidje 2016).

Uses probability-threshold pruning (no flux protection), mass-conservation
accept/reject time stepping, and SSA-driven + reachability expansion.

# Keyword Arguments
- `ε::Float64 = 1e-3`: Total probability mass loss tolerance.
- `τ_init::Float64 = 1.0`: Initial time step.
- `prob_threshold::Float64 = 0.0`: Hard probability threshold for pruning (0 = use ε/|S|).
- `ℓ::Int = 5`: Lookback window length for max-probability tracking.
- `r::Int = 1`: Reachability expansion depth after SSA.
- `save_interval::Int = 1000`: Save a snapshot every N iterations.
- `max_iter::Int = typemax(Int)`: Maximum number of iterations (for diagnostic runs).
"""
Base.@kwdef struct KrylovFSP
    ε::Float64 = 1e-3
    τ_init::Float64 = 1.0
    prob_threshold::Float64 = 0.0
    ℓ::Int = 5
    r::Int = 1
    save_interval::Int = 1000
    max_iter::Int = typemax(Int)
end

"""
    KrylovFSPDiagnostics

Per-iteration diagnostic data from a KrylovFSP solve.
"""
struct KrylovFSPDiagnostics
    τ_log::Vector{Float64}
    mass_log::Vector{Float64}
    time_log::Vector{Float64}
    walltime_log::Vector{Float64}
    size_log::Vector{Int}
    reject_count::Int
    total_iters::Int
end

"""
    solve(prob::FSPProblem, alg::KrylovFSP) -> (FSPSolution, KrylovFSPDiagnostics)

Run the Krylov-FSP-SSA algorithm (Dinh & Sidje 2016, Algorithm [18]).

Returns a tuple of `(solution, diagnostics)` where diagnostics contains
per-iteration τ, mass, timing, and rejection data for failure-mode analysis.
"""
function CommonSolve.solve(prob::FSPProblem{E,T}, alg::KrylovFSP) where {E,T}
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

    # Diagnostics
    τ_log = Float64[]
    mass_log = Float64[]
    time_log = Float64[]
    walltime_log = Float64[]
    size_log = Int[]
    reject_count = 0

    # Lookback window for max-probability tracking
    prob_history = Vector{Dict{E,Float64}}()

    t = Float64(t0)
    τ = Float64(alg.τ_init)
    iter = 0
    wall_start = time()

    while t < tf && iter < alg.max_iter
        iter += 1

        # Determine pruning threshold
        threshold = alg.prob_threshold > 0 ? alg.prob_threshold : alg.ε / max(length(sp), 1)

        # Update lookback window
        current_probs = Dict{E,Float64}(s => sp.probs[sp.index[s]] for s in sp.states)
        push!(prob_history, current_probs)
        if length(prob_history) > alg.ℓ
            popfirst!(prob_history)
        end

        # Compute max probability over lookback window
        max_probs = Dict{E,Float64}()
        for hist in prob_history
            for (s, p) in hist
                max_probs[s] = max(get(max_probs, s, 0.0), p)
            end
        end

        # 1. Threshold prune (with lookback protection)
        prune_threshold!(sp, threshold; max_probs=max_probs)

        # 2. SSA-driven expansion
        expand_ssa!(sp, model, rates, t, bc; depth=1)

        # 3. r-step reachability expansion
        if alg.r > 0
            expand!(sp, model, bc; depth=alg.r)
        end

        # 4. Build generator + propagate
        τ_step = min(τ, tf - t)
        A, _, _ = build_generator(sp, model, rates, t)

        # Attempt matrix exponential; if it fails numerically, reject and halve τ
        local p_new
        expmv_ok = true
        try
            p_new = expmv(τ_step, A, sp.probs)
            if any(isnan, p_new) || any(isinf, p_new)
                expmv_ok = false
            end
        catch
            expmv_ok = false
        end

        if !expmv_ok
            # Treat as rejection
            τ /= 2
            reject_count += 1
            push!(τ_log, τ_step)
            push!(mass_log, NaN)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))
            continue
        end

        # Clamp small negative values from numerical noise
        @inbounds for i in eachindex(p_new)
            if p_new[i] < 0
                p_new[i] = 0.0
            end
        end

        # 5. Accept/reject based on mass conservation
        mass = sum(p_new)
        mass_target = 1.0 - alg.ε

        if mass >= mass_target
            # Accept step
            sp.probs = p_new
            renormalize!(sp)
            t += τ_step

            push!(τ_log, τ_step)
            push!(mass_log, mass)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))

            # Save snapshot
            if iter % alg.save_interval == 0 || t >= tf
                push!(sol_t, t)
                push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
                push!(sol_sizes, length(sp))
                @info "KrylovFSP: t = $(round(t; sigdigits=4)), |S| = $(length(sp)), τ = $(round(τ_step; sigdigits=3)), mass = $(round(mass; sigdigits=6))"
            end
        else
            # Reject: halve τ
            τ /= 2
            reject_count += 1

            push!(τ_log, τ_step)
            push!(mass_log, mass)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))

            if iter % alg.save_interval == 0
                @info "KrylovFSP REJECT: t = $(round(t; sigdigits=4)), τ = $(round(τ; sigdigits=3)), mass = $(round(mass; sigdigits=6))"
            end
        end
    end

    # Ensure final snapshot is captured (avoid duplicate if already saved in loop)
    if length(sol_t) == 1 || sol_t[end] < t - eps(t)
        push!(sol_t, t)
        push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
        push!(sol_sizes, length(sp))
    end

    sol = FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
    diag = KrylovFSPDiagnostics(τ_log, mass_log, time_log, walltime_log, size_log, reject_count, iter)

    (sol, diag)
end
