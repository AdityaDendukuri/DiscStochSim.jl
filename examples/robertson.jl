# Robertson Autocatalytic System (Extreme Stiffness)
# Saves diagnostics to examples/results/robertson_fpfsp.jls for plotting.
using DiscStochSim, Catalyst, Serialization

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

rates = [0.04, 3e6, 1.0]
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates;
                  bounds=(0, 100_001))

cb = terminal_progress(["A", "1000×B", "C"]; tf=1e4, xscale=:log10,
                       transform = m -> (m2 = copy(m); m2[:,2] .*= 10000; m2))
t_wall_start = time()
sol, diag = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
                                    save_interval=1000,
                                    expand_method=:stoich,
                                    progress_callback=cb))
t_wall = time() - t_wall_start

traj = mean_trajectory(sol)
println("Final: A=$(traj[end,1]), B=$(traj[end,2]), C=$(traj[end,3])")
println("Snapshots: $(length(sol)), final |S| = $(sol.state_space_sizes[end])")
println("Wall time: $(round(t_wall; digits=1))s, total_iters=$(diag.total_iters)")
if !isempty(diag.dt_log)
    println("dt range: $(minimum(diag.dt_log)) to $(maximum(diag.dt_log))")
end

# Log-subsample per-step diagnostics: pick ~2000 indices log-spaced in time.
let t_full = diag.t_log, n_save = 2000
    t_grid   = 10 .^ range(log10(t_full[1]), log10(t_full[end]); length=n_save)
    raw_idx  = searchsortedfirst.(Ref(t_full), t_grid)
    global log_idx = unique(clamp.(raw_idx, 1, length(t_full)))
    global dt_log_sub   = diag.dt_log[log_idx]
    global t_log_sub    = t_full[log_idx]
    global size_log_sub = diag.size_log[log_idx]
    global flux_log_sub = diag.flux_log[log_idx]
end
println("Log-subsampled: $(length(t_log_sub)) points from $(diag.total_iters) steps")

# Save diagnostics for plotting
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
out = Dict(
    "dt_log"      => dt_log_sub,
    "t_log"       => t_log_sub,
    "size_log"    => size_log_sub,
    "flux_log"    => flux_log_sub,
    "snap_t"      => sol.t,
    "snap_traj"   => traj,          # mean trajectories at snapshot times
    "snap_sizes"  => sol.state_space_sizes,
    "total_iters" => diag.total_iters,
    "wall_time"   => t_wall,
    "mean_A"      => traj[end, 1],
    "mean_B"      => traj[end, 2],
    "mean_C"      => traj[end, 3],
    "S_final"     => sol.state_space_sizes[end],
)
serialize(joinpath(results_dir, "robertson_fpfsp.jls"), out)
println("Saved diagnostics → examples/results/robertson_fpfsp.jls")
