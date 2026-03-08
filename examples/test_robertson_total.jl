# Compare flux_method=:total vs :maximum for Robertson
# Short tf=100 to get results quickly
using DiscStochSim, Catalyst

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

rates = [0.04, 3e6, 1.0]
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 100.0), rates;
                  bounds=(0, 100_001))

for method in [:maximum, :total]
    println("\n--- flux_method=:$method ---")
    t0 = time()
    sol, diag = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
                                        save_interval=1000,
                                        flux_method=method,
                                        expand_method=:stoich))
    t_wall = time() - t0
    traj = mean_trajectory(sol)
    println("  Wall time:   $(round(t_wall; digits=1))s")
    println("  Total iters: $(diag.total_iters)")
    println("  t_final:     $(sol.t[end])")
    println("  Final means: A=$(round(traj[end,1];digits=2)), B=$(round(traj[end,2];sigdigits=3)), C=$(round(traj[end,3];digits=2))")
    println("  |S|_final:   $(sol.state_space_sizes[end])")
    if !isempty(diag.dt_log)
        println("  dt range:    [$(minimum(diag.dt_log)), $(maximum(diag.dt_log))]")
    end
end
