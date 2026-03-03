using DiscStochSim, Catalyst, Random

rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

rates = [1e-6, 1e-1]
prob = FSPProblem(rn, CartesianIndex(1, 0, 0), (0.0, 1e5), rates; bounds=(0, Int(1e11)))

# AdaptiveFSP
println("Running AdaptiveFSP...")
Random.default_rng(1)
t0 = time()
sol_afsp, _ = solve(prob, AdaptiveFSP(ε_dt=1e-3, prob_quantile=0.9, flux_tolerance=1e-6, save_interval=1000,
                                      flux_method=:total, expand_method=:stoich))
t_afsp = time() - t0
traj = mean_trajectory(sol_afsp)
println("  Done: $(length(sol_afsp)) snaps, $(round(t_afsp; digits=1))s")
println("  Final: A=$(round(traj[end,1];sigdigits=4)), B=$(round(traj[end,2];sigdigits=4)), C=$(round(traj[end,3];sigdigits=4))")

# KrylovFSP
println("\nRunning KrylovFSP...")
Random.default_rng(1)
t0 = time()
sol_kfsp, diag = solve(prob, KrylovFSP(ε=1e-4, τ_init=1e3, ℓ=5, r=1, save_interval=100, max_iter=500))
t_kfsp = time() - t0
traj_k = mean_trajectory(sol_kfsp)
println("  Done: $(diag.total_iters) iters, $(diag.reject_count) rejections, $(round(t_kfsp; digits=1))s")
println("  Final: A=$(round(traj_k[end,1];sigdigits=4)), B=$(round(traj_k[end,2];sigdigits=4)), C=$(round(traj_k[end,3];sigdigits=4))")
