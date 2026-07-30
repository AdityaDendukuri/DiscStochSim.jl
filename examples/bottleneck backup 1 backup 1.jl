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

# dt_max = 1/k2 caps the step at the downstream reaction timescale,
# preventing a single giant step when Phi_total = k1 ≈ 0 at t=0.
dt_max_val = 1.0 / rates[2]   # = 1/k2 = 10.0

# Flux-aware run
cb_flux = terminal_progress(["A", "B", "C"]; tf=1e6, xscale=:log10, title="Bottleneck flux-aware")
sol_flux, _ = solve(prob, AdaptiveFSP(ε_dt=1.0, prob_quantile=0.9, flux_tolerance=1e-6,
                                     save_interval=1000, dt_max=dt_max_val,
                                     expand_method=:stoich, progress_callback=cb_flux))
traj_flux = mean_trajectory(sol_flux)

# No-flux run (set flux_tolerance=0 to disable flux protection)
cb_noflux = terminal_progress(["A", "B", "C"]; tf=1e6, xscale=:log10, title="Bottleneck no-flux")
sol_noflux, _ = solve(prob, AdaptiveFSP(ε_dt=1.0, prob_quantile=0.9, flux_tolerance=0.0,
                                       save_interval=1000, dt_max=dt_max_val,
                                       expand_method=:stoich, progress_callback=cb_noflux))
traj_noflux = mean_trajectory(sol_noflux)

println("Flux-aware final:  A=$(traj_flux[end,1]), B=$(traj_flux[end,2]), C=$(traj_flux[end,3])")
println("No-flux final:     A=$(traj_noflux[end,1]), B=$(traj_noflux[end,2]), C=$(traj_noflux[end,3])")
println("Flux snapshots: $(length(sol_flux)), No-flux snapshots: $(length(sol_noflux))")
