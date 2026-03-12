# Compare :flux vs :ssa expansion on Oregonator (short run)
using DiscStochSim, Catalyst, LinearAlgebra, Random

rn = @reaction_network begin
    @species X(t) Y(t) Z(t); @parameters k1 k2 k3 k4 k5
    k1, Y --> X; k2, X + Y --> 0; k3, X --> 2X + Z; k4, 2X --> 0; k5, Z --> Y
end
y1s,y2s,y3s=500.,1000.,2000.; mu1s,mu2s=2000.,50000.
rates=[mu1s/y2s,mu2s/(y1s*y2s),(mu1s+mu2s)/y1s,2*mu1s/y1s^2,(mu1s+mu2s)/y3s]
prob = FSPProblem(rn, CartesianIndex(500,1000,2000), (0.0,1.0), rates; bounds=(0,50_000))

configs = [
    ("SSA   (1 neighbor/state)",  :ssa,  0.10),
    ("Flux  (thr=0.10, ~3/state)", :flux, 0.10),
    ("Flux  (thr=0.33, ~1/state)", :flux, 0.33),
    ("Flux  (thr=0.40, <1/state)", :flux, 0.40),
]

for (label, method, threshold) in configs
    alg = AdaptiveFSP(ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1.0,
                      expansion_depth=1, save_interval=500,
                      expand_method=method, expand_threshold=threshold,
                      max_iter=2000)
    Random.seed!(42)
    t_start = time()
    sol, diag = solve(prob, alg)
    t_wall = time() - t_start
    traj = mean_trajectory(sol)
    println("=== $label ===")
    println("  iters=$(diag.total_iters), wall=$(round(t_wall;digits=1))s")
    println("  t_reached=$(round(sol.t[end];sigdigits=4))")
    println("  max|S|=$(maximum(sol.state_space_sizes)), final|S|=$(sol.state_space_sizes[end])")
    if !isempty(diag.dt_log)
        println("  dt: min=$(round(minimum(diag.dt_log);sigdigits=3)), max=$(round(maximum(diag.dt_log);sigdigits=3))")
    end
    println("  Final: X=$(round(traj[end,1];digits=1)), Y=$(round(traj[end,2];digits=1)), Z=$(round(traj[end,3];digits=1))")
    println()
end
