# Robertson Autocatalytic System (Extreme Stiffness)
using DiscStochSim, Catalyst

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

rates = [0.04, 3e6, 1.0]
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates;
                  bounds=(0, 100_001))
sol = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
                              save_interval=1000))

traj = mean_trajectory(sol)
println("Final: A=$(traj[end,1]), B=$(traj[end,2]), C=$(traj[end,3])")
println("Snapshots: $(length(sol)), final |S| = $(sol.state_space_sizes[end])")
