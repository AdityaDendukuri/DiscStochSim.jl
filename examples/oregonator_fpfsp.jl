# Oregonator — FP-FSP with total-flux time stepping
# Saves diagnostics to examples/results/oregonator_fpfsp.jls for plotting.
using DiscStochSim, Catalyst, Serialization

rn = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end

y1s, y2s, y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s    = 2000.0, 50000.0
rates = [
    mu1s / y2s,
    mu2s / (y1s * y2s),
    (mu1s + mu2s) / y1s,
    2 * mu1s / y1s^2,
    (mu1s + mu2s) / y3s,
]

tf   = 0.11
prob = FSPProblem(rn, CartesianIndex(500, 1000, 2000), (0.0, tf), rates;
                  bounds=(0, 50_000))

println("Running FP-FSP Oregonator (total flux, tf=$tf)...")
cb = terminal_progress(["X", "Y", "Z"]; tf=tf, title="Oregonator FP-FSP")
t_wall_start = time()
sol, diag = solve(prob, AdaptiveFSP(
    ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1.0,
    expansion_depth=1, save_interval=1000, flux_method=:total, expand_method=:ssa,
    progress_callback=cb,
))
t_wall = time() - t_wall_start

traj = mean_trajectory(sol)
println("Done: $(diag.total_iters) iters, $(round(t_wall; digits=1))s, " *
        "t_reached=$(round(sol.t[end]; sigdigits=4))")
println("Final: X=$(round(traj[end,1]; digits=1)), " *
        "Y=$(round(traj[end,2]; digits=1)), Z=$(round(traj[end,3]; digits=1))")
println("max|S|=$(maximum(sol.state_space_sizes)), final|S|=$(sol.state_space_sizes[end])")
if !isempty(diag.dt_log)
    println("dt range: $(minimum(diag.dt_log)) to $(maximum(diag.dt_log))")
end

results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
out = Dict(
    "dt_log"      => diag.dt_log,
    "t_log"       => diag.t_log,
    "size_log"    => diag.size_log,
    "snap_t"      => sol.t,
    "snap_sizes"  => sol.state_space_sizes,
    "total_iters" => diag.total_iters,
    "wall_time"   => t_wall,
    "mean_X"      => traj[end, 1],
    "mean_Y"      => traj[end, 2],
    "mean_Z"      => traj[end, 3],
    "S_max"       => maximum(sol.state_space_sizes),
    "S_final"     => sol.state_space_sizes[end],
)
serialize(joinpath(results_dir, "oregonator_fpfsp.jls"), out)
println("Saved diagnostics → examples/results/oregonator_fpfsp.jls")
