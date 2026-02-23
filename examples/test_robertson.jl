using DiscStochSim, Catalyst, Random

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

rates = [0.04, 3e6, 1.0]
tf = 10.0
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, tf), rates; bounds=(0, 100_001))

# AdaptiveFSP
println("Running AdaptiveFSP...")
Random.default_rng(1)
t0 = time()
sol_afsp = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9, save_interval=5000))
t_afsp = time() - t0
traj = mean_trajectory(sol_afsp)
println("  Done: $(length(sol_afsp)) snaps, $(round(t_afsp; digits=1))s, t_final=$(sol_afsp.t[end])")
println("  Final: A=$(round(traj[end,1];digits=1)), B=$(round(traj[end,2];digits=1)), C=$(round(traj[end,3];digits=1))")

# KrylovFSP
println("\nRunning KrylovFSP...")
Random.default_rng(1)
t0 = time()
sol_kfsp, diag = solve(prob, KrylovFSP(ε=1e-2, τ_init=0.1, ℓ=5, r=1, save_interval=50, max_iter=500))
t_kfsp = time() - t0
traj_k = mean_trajectory(sol_kfsp)
println("  Done: $(diag.total_iters) iters, $(diag.reject_count) rejections, $(round(t_kfsp; digits=1))s")
println("  Time reached: $(sol_kfsp.t[end])")
if length(diag.τ_log) > 0
    valid_τ = filter(x -> x > 0, diag.τ_log)
    if length(valid_τ) > 0
        println("  τ range: [$(minimum(valid_τ)), $(maximum(valid_τ))]")
    end
end
