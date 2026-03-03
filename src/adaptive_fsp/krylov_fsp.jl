"""
    KrylovFSP(; kwargs...)

Faithful implementation of the Adaptive Krylov-FSP-SSA algorithm
(Algorithm 2, Vo & Sidje 2016, WCECS).

Key features:
- Joint probability+derivative pruning: drop state i only if BOTH
  `p_i < droptol` AND `(Ap)_i < droptol_deriv`.
- Time-scaled FSP acceptance: `sum(w) >= 1 - ε*(t+τ)/tf`.
- Adaptive step size: τ grows by `growth_factor` after successful steps,
  halves on FSP rejection.
- droptol auto-reduction: if removing all candidates in L={i: p_i < droptol}
  would violate the mass budget, droptol is reduced by 10× until safe.

# Keyword Arguments
- `ε::Float64 = 1e-6`: FSP tolerance ε_FSP (controls total error at tf).
- `τ_init::Float64 = 1.0`: Initial time step.
- `droptol::Float64 = 1e-10`: Initial probability threshold for pruning.
- `droptol_deriv::Float64 = 1e-16`: Derivative threshold for pruning.
- `Tol::Float64 = 1e-8`: Krylov approximation tolerance (passed to expmv).
- `growth_factor::Float64 = 2.0`: τ multiplier after a successful step.
- `r::Int = 1`: Reachability expansion depth after SSA.
- `save_interval::Int = 1000`: Save a snapshot every N accepted steps.
- `max_iter::Int = typemax(Int)`: Maximum number of iterations.
"""
Base.@kwdef struct KrylovFSP
    ε::Float64 = 1e-6
    τ_init::Float64 = 1.0
    droptol::Float64 = 1e-10
    droptol_deriv::Float64 = 1e-16
    Tol::Float64 = 1e-8
    growth_factor::Float64 = 2.0
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

Run the Adaptive Krylov-FSP-SSA algorithm (Vo & Sidje 2016, Algorithm 2).

Loop structure at each step k:
  1. Build A = A_{J_k} from current projection.
  2. Compute w = exp(τ A) p via Krylov (expmv).
  3. Check time-scaled FSP criterion: sum(w) >= 1 - ε*(t+τ)/tf.
     If fails: halve τ, set iexpand=true, go to 2.
  4. Accept: p ← w, t ← t+τ.
  5. Prune using joint filter: remove i if p_i < droptol AND (Ap)_i < droptol_deriv.
     Auto-reduce droptol if removing all candidates would violate the mass budget.
  6. SSA-guided expansion + r-step reachability for next step.
  7. If step was accepted first try (iexpand=false): τ ← τ * growth_factor.
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

    # Expand once to give Krylov something to work with
    expand!(sp, model, bc; depth=alg.r)

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

    t = Float64(t0)
    τ = Float64(alg.τ_init)
    iter = 0
    accepted_steps = 0
    wall_start = time()

    while t < tf && iter < alg.max_iter
        iter += 1

        # ── Step 1: Build absorbing-form generator for current projection ──
        # absorbing=true: diagonal = total outflow (including to states outside S),
        # so probability leaks out when the projection boundary is hit.
        # This is required for the FSP mass-loss criterion to function correctly.
        A, _, _ = build_generator(sp, model, rates, t; absorbing=true)

        # ── Steps 2–3: Find τ satisfying FSP criterion ────────────────────
        iexpand = false
        τ_step = min(τ, tf - t)
        local p_new
        fsp_ok = false
        n_τ_reduce = 0

        while !fsp_ok && n_τ_reduce < 60
            τ_step = min(τ, tf - t)

            p_new = expv(τ_step, A, sp.probs; tol=alg.Tol)

            # Clamp small negatives from numerical noise
            @inbounds for i in eachindex(p_new)
                if p_new[i] < 0
                    p_new[i] = 0.0
                end
            end

            if any(isnan, p_new) || any(isinf, p_new)
                τ /= 2
                iexpand = true
                reject_count += 1
                n_τ_reduce += 1
                continue
            end

            mass = sum(p_new)
            # Time-scaled FSP criterion (Eq. 8 in Vo & Sidje 2016)
            fsp_threshold = 1.0 - alg.ε * (t + τ_step) / tf

            if mass >= fsp_threshold
                fsp_ok = true
            else
                τ /= 2
                iexpand = true
                reject_count += 1
                n_τ_reduce += 1

                push!(τ_log, τ_step)
                push!(mass_log, mass)
                push!(time_log, t)
                push!(walltime_log, time() - wall_start)
                push!(size_log, length(sp))
            end
        end

        if !fsp_ok
            # Fallback: accept tiny step to avoid infinite loop
            τ_step = min(τ, tf - t)
            p_new = expv(τ_step, A, sp.probs; tol=alg.Tol)
            @inbounds for i in eachindex(p_new)
                if p_new[i] < 0; p_new[i] = 0.0; end
            end
        end

        # ── Step 4: Accept ─────────────────────────────────────────────────
        sp.probs = p_new
        t += τ_step
        accepted_steps += 1

        mass_accepted = sum(sp.probs)

        push!(τ_log, τ_step)
        push!(mass_log, mass_accepted)
        push!(time_log, t)
        push!(walltime_log, time() - wall_start)
        push!(size_log, length(sp))

        # ── Step 5: Prune with joint probability + derivative filter ───────
        # Derivative at current time: dp = A * p  (A is still valid for J_k)
        dp = A * sp.probs

        # Auto-reduce droptol if necessary (Vo & Sidje Section III.C)
        droptol = alg.droptol
        fsp_mass_budget = 1.0 - alg.ε * t / tf
        while droptol > 1e-300
            # Mass of L = all probability-threshold candidates
            mass_L = 0.0
            @inbounds for i in eachindex(sp.probs)
                if sp.probs[i] < droptol
                    mass_L += sp.probs[i]
                end
            end
            if mass_accepted - mass_L >= fsp_mass_budget
                break  # safe to proceed
            end
            droptol /= 10
        end

        # Apply joint filter: remove i iff p_i < droptol AND dp_i < droptol_deriv
        n = length(sp)
        remove_mask = falses(n)
        @inbounds for i in 1:n
            if sp.probs[i] < droptol && dp[i] < alg.droptol_deriv
                remove_mask[i] = true
            end
        end
        remove_states!(sp, remove_mask)

        if accepted_steps % alg.save_interval == 0 || t >= tf
            push!(sol_t, t)
            push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
            push!(sol_sizes, length(sp))
            @info "KrylovFSP: t=$(round(t; sigdigits=4)), |S|=$(length(sp)), " *
                  "τ=$(round(τ_step; sigdigits=3)), mass=$(round(mass_accepted; sigdigits=6)), " *
                  "droptol=$(droptol)"
        end

        # ── Step 6: SSA-guided expansion for next step ─────────────────────
        expand_ssa!(sp, model, rates, t, bc; depth=1)
        if alg.r > 0
            expand!(sp, model, bc; depth=alg.r)
        end

        # ── Step 7: Adapt τ for next step ──────────────────────────────────
        # If accepted first try: grow τ. If had to reduce: keep reduced τ.
        if !iexpand
            τ = min(τ * alg.growth_factor, tf - t)
        end
        # If iexpand=true: τ was already halved; keep as is for next step.
    end

    # Ensure final snapshot
    if length(sol_t) == 1 || sol_t[end] < t - eps(t)
        push!(sol_t, t)
        push!(sol_snaps, (copy(sp.states), copy(sp.probs)))
        push!(sol_sizes, length(sp))
    end

    sol = FSPSolution{E}(sol_t, sol_snaps, sol_sizes, n_species)
    diag = KrylovFSPDiagnostics(τ_log, mass_log, time_log, walltime_log,
                                 size_log, reject_count, iter)
    (sol, diag)
end
