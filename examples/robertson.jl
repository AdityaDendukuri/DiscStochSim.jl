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

t_wall_start = time()
sol, diag = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
                                    save_interval=1000,
                                    flux_method=:maximum,  # dt = ε_dt / maxᵢ pᵢwᵢ (stiff)
                                    expand_method=:stoich)) # stoichiometric expansion
t_wall = time() - t_wall_start

traj = mean_trajectory(sol)
println("Final: A=$(traj[end,1]), B=$(traj[end,2]), C=$(traj[end,3])")
println("Snapshots: $(length(sol)), final |S| = $(sol.state_space_sizes[end])")
println("Wall time: $(round(t_wall; digits=1))s, total_iters=$(diag.total_iters)")
if !isempty(diag.dt_log)
    println("dt range: $(minimum(diag.dt_log)) to $(maximum(diag.dt_log))")
end

# Save diagnostics for plotting
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
out = Dict(
    "dt_log"      => diag.dt_log,
    "t_log"       => diag.t_log,
    "total_iters" => diag.total_iters,
    "wall_time"   => t_wall,
    "mean_A"      => traj[end, 1],
    "mean_B"      => traj[end, 2],
    "mean_C"      => traj[end, 3],
    "S_final"     => sol.state_space_sizes[end],
)
serialize(joinpath(results_dir, "robertson_fpfsp.jls"), out)
println("Saved diagnostics → examples/results/robertson_fpfsp.jls")
