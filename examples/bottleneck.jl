# Bottleneck System: A→B, B→B+C
# Demonstrates flux-aware pruning preserving bottleneck states
using DiscStochSim, Catalyst

rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

rates = [1e-6, 1e-1]
prob = FSPProblem(rn, CartesianIndex(1, 0, 0), (0.0, 1e6), rates;
                  bounds=(0, Int(1e11)))

# Flux-aware run
sol_flux, _ = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.9, flux_tolerance=1e-6,
                                     save_interval=1000,
                                     flux_method=:total, expand_method=:stoich))
traj_flux = mean_trajectory(sol_flux)

# No-flux run (set flux_tolerance=0 to disable flux protection)
sol_noflux, _ = solve(prob, AdaptiveFSP(ε_dt=0.01, prob_quantile=0.9, flux_tolerance=0.0,
                                       save_interval=1000,
                                       flux_method=:total, expand_method=:stoich))
traj_noflux = mean_trajectory(sol_noflux)

println("Flux-aware final:  A=$(traj_flux[end,1]), B=$(traj_flux[end,2]), C=$(traj_flux[end,3])")
println("No-flux final:     A=$(traj_noflux[end,1]), B=$(traj_noflux[end,2]), C=$(traj_noflux[end,3])")
println("Flux snapshots: $(length(sol_flux)), No-flux snapshots: $(length(sol_noflux))")
