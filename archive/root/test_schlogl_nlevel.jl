using DiscStochSim
using LinearAlgebra

# ── K=4, 1-level V-cycle (4→2), bistable front IC ────────────────────────
K = 4
fine_grid = VoxelGrid(K, 1.0, 0)
model     = SchloglModel1D(0.1)
rates     = Float64[]
tspan     = (0.0, 1.0)

# IC: high state on left pair, low state on right pair
n_high, n_low = 148, 31
u0_front = CartesianIndex(n_high, n_high, n_low, n_low)

alg = RDMEMultigridFSP(
    dt                = 0.1,
    n_levels          = 1,
    n_max             = 250,
    save_every        = 2,
    max_states        = 200_000,
    expand_depth      = 0,          # let prolong_conditional drive state space growth
    prune_tol         = 1e-10,
    coarse_expand_depth = 2,
    τ_pre             = 0.5,
)

println("════ n-level V-cycle: K=4, 1 level (4→2), bistable-front IC ════")
println("  u0 = (148,148,31,31)  [high|low]")
println()

sol = solve_rdme_multigrid(model, fine_grid, u0_front, tspan, rates, alg)

println("State space sizes: ", sol.state_space_sizes)
println()
mu = mean_voxel_counts(sol)
for (ti, t) in enumerate(sol.t)
    println("  t=$(round(t,digits=1)):  μ=$(round.(mu[ti,:], digits=1))  |S|=$(sol.state_space_sizes[ti])")
end
println()
println("(expect: pair 1 diffusing toward Schlögl low state, pair 2 reacting)")
