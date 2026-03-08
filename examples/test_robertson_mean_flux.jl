# Test unified time-stepping formula: Δt = ε_dt / Φ̄  where Φ̄ = Φ_total / |J|
# Compare three approaches on Robertson (tf=100 for speed):
#   :maximum  → Δt = ε_dt / max_i(p_i·w_i)          [current Robertson setting]
#   :mean     → Δt = ε_dt / (Φ_total / |J|)  = ε_dt·|J|/Φ_total  [proposed]
#   :total    → Δt = ε_dt / Φ_total                   [current default, very slow]
#
# Expected: :mean ≈ :maximum in dt and accuracy; :total is ~|J|× slower.

using DiscStochSim, Catalyst, CommonSolve, ExponentialUtilities,
      SparseArrays

# Patch in a tiny local solve that supports :mean
function solve_robertson_custom(flux_sym; tf=100.0, ε_dt=0.01, max_iter=typemax(Int))
    rn = @reaction_network begin
        k1, A --> B
        k2, 2B --> C + B
        k3, C + B --> A + C
    end
    rates = [0.04, 3e6, 1.0]
    prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, tf), rates;
                     bounds=(0, 100_001))

    model = prob.model; bc = prob.bc; t0, _tf = prob.tspan
    sp = DiscStochSim.StateSpace{typeof(prob.u0), Float64}()
    DiscStochSim.add_state!(sp, prob.u0, 1.0)

    A_old = spzeros(Float64, 0, 0)
    gids_old = Int[]
    dt_log = Float64[]; t = Float64(t0); iter = 0

    while t < _tf && iter < max_iter
        iter += 1
        DiscStochSim.expand!(sp, model, bc; depth=1)
        A, in_flow, out_flow = DiscStochSim.reconstruct_generator(sp, model, rates, t, A_old, gids_old)
        A_old = A; gids_old = DiscStochSim.get_global_ids(sp)

        out_flux = Vector(-out_flow * sp.probs)
        Φ_total = sum(out_flux)
        Φ = if flux_sym === :maximum
                maximum(out_flux)
            elseif flux_sym === :mean
                Φ_total / length(sp)
            else  # :total
                Φ_total
            end
        dt = Φ > 0 ? ε_dt / Φ : (_tf - t)
        dt = min(dt, _tf - t)
        push!(dt_log, dt)

        sp.probs = expv(dt, A, sp.probs)
        DiscStochSim.compress!(sp, model, rates, t, 0.4; flux_tolerance=1e-9)
        DiscStochSim.renormalize!(sp)
        t += dt
    end

    traj_end = DiscStochSim.mean_trajectory(
        DiscStochSim.FSPSolution{typeof(prob.u0)}(
            [t0, t],
            [(copy(sp.states), copy(sp.probs)), (copy(sp.states), copy(sp.probs))],
            [1, length(sp)],
            3
        )
    )
    (t, iter, dt_log, traj_end[end, :])
end

for sym in [:maximum, :mean]
    println("\n--- flux: :$sym ---")
    t0 = time()
    t_final, iters, dt_log, means = solve_robertson_custom(sym; tf=100.0)
    t_wall = time() - t0
    println("  Wall time:   $(round(t_wall; digits=1))s")
    println("  Total iters: $iters")
    println("  t_final:     $t_final")
    println("  dt range:    [$(round(minimum(dt_log);sigdigits=3)), $(round(maximum(dt_log);sigdigits=3))]")
    println("  Means:       A=$(round(means[1];digits=2)), B=$(round(means[2];sigdigits=3)), C=$(round(means[3];digits=2))")
end
