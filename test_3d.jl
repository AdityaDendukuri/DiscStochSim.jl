using DiscStochSim, Printf
model = SchloglModel1D(0.1)
n_low, n_uns, n_high = schlogl_fixed_points(model)
println("n_low=$n_low  n_uns=$n_uns  n_high=$n_high")

# Try a 3x3x3 grid (27 voxels) — too large for joint FSP
# Use 2x2x2 = 8 voxels for joint FSP
g3 = grid_3d(2,2,2)
println("grid_3d(2,2,2): $(g3.n_voxels) voxels, $(length(g3.edges)) edges")
println("Degrees: ", degrees(g3))

# Supercritical IC: 4 voxels at n_high, 4 at n_low
# In 2x2x2 grid: voxels 1-4 = "front face", 5-8 = "back face"
u0_sup = CartesianIndex(n_high, n_high, n_high, n_high, n_low, n_low, n_low, n_low)
# Subcritical IC: 4 at n_uns, 4 at n_low
u0_sub = CartesianIndex(n_uns, n_uns, n_uns, n_uns, n_low, n_low, n_low, n_low)

alg = RDMEMultigridFSP(dt=1.0, n_levels=1, n_max=220, save_every=1,
    max_states=300_000, expand_depth=0, coarse_r_depth=1, coarse_d_depth=1,
    prune_tol=1e-5, weight_tol=1e-5, τ_pre=0.3)

println("\n--- Supercritical IC (4 voxels at n_high) ---")
@time sol_sup = solve_rdme_graph(model, g3, u0_sup, (0.0, 30.0), alg)
mu_sup = mean_voxel_counts(sol_sup)
@printf("t=30: E[n] = %s\n", join(round.(mu_sup[end,:], digits=1), ", "))
@printf("|S| range: %d → %d\n", minimum(sol_sup.state_space_sizes), maximum(sol_sup.state_space_sizes))
cp_sup = config_probs(sol_sup, n_uns, length(sol_sup.t))
println("Config probs at t=30: ", [(m, round(cp_sup[m+1], digits=3)) for m in 0:g3.n_voxels if cp_sup[m+1]>0.01])

println("\n--- Subcritical IC (4 voxels at n_uns) ---")
@time sol_sub = solve_rdme_graph(model, g3, u0_sub, (0.0, 30.0), alg)
mu_sub = mean_voxel_counts(sol_sub)
@printf("t=30: E[n] = %s\n", join(round.(mu_sub[end,:], digits=1), ", "))
cp_sub = config_probs(sol_sub, n_uns, length(sol_sub.t))
println("Config probs at t=30: ", [(m, round(cp_sub[m+1], digits=3)) for m in 0:g3.n_voxels if cp_sub[m+1]>0.01])
