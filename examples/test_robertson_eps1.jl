# Test unified dt = 1/Phi_total (eps_dt=1.0, :total) on Robertson
# If correct answer matches eps_dt=0.01 :maximum, we have a universal formula.
using DiscStochSim, Catalyst

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
rates = [0.04, 3e6, 1.0]
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 100.0), rates; bounds=(0, 100_001))

# Reference (ODE): at t=100, A≈6142, B≈0.065, C≈3858
println("ODE reference at t=100: A≈6142, B≈0.065, C≈3858\n")

for (label, ε, method) in [
    ("eps=0.01 :maximum (current)",  0.01, :maximum),
    ("eps=1.0  :total  (proposed)",  1.0,  :total),
    ("eps=0.1  :total",              0.1,  :total),
]
    println("--- $label ---")
    t0 = time()
    sol, diag = solve(prob, AdaptiveFSP(ε_dt=ε, prob_quantile=0.4, flux_tolerance=1e-9,
                                        save_interval=100_000, flux_method=method,
                                        expand_method=:stoich))
    t_wall = time() - t0
    traj = mean_trajectory(sol)
    println("  Wall time:   $(round(t_wall; digits=1))s")
    println("  Total iters: $(diag.total_iters)")
    println("  A=$(round(traj[end,1];digits=1)),  B=$(round(traj[end,2];sigdigits=3)),  C=$(round(traj[end,3];digits=1))")
    println("  dt range:    [$(round(minimum(diag.dt_log);sigdigits=3)), $(round(maximum(diag.dt_log);sigdigits=3))]")
    println()
end
