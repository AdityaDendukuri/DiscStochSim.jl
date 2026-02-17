# Oregonator 3-Species Oscillatory System
using DiscStochSim, Catalyst

rn = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end

# Steady-state parametrization
Y1s, Y2s, Y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s = 2000.0, 50000.0
rates = [mu1s/Y2s, mu2s/(Y1s*Y2s), (mu1s+mu2s)/Y1s, 2*mu1s/Y1s^2, (mu1s+mu2s)/Y3s]

prob = FSPProblem(rn, CartesianIndex(500, 1000, 2000), (0.0, 3.0), rates;
                  bounds=(0, 50_000))
sol = solve(prob, AdaptiveFSP(ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1e-4,
                              save_interval=1000))

traj = mean_trajectory(sol)
println("Final: X=$(traj[end,1]), Y=$(traj[end,2]), Z=$(traj[end,3])")
println("Snapshots: $(length(sol)), final |S| = $(sol.state_space_sizes[end])")
