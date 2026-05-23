using DiscStochSim, Printf
model = SchloglModel1D(1.0)
n_low, n_uns, n_high = schlogl_fixed_points(model)
K = 6; fine_grid = VoxelGrid(K, 1.0, 0)
u0 = CartesianIndex(n_low, n_low, n_high, n_high, n_low, n_low)
println("IC = $u0   n_low=$n_low  n_uns=$n_uns  n_high=$n_high")

function run_droplet(u0_ic, D_val, t_end_val; label="")
    m_mdl = SchloglModel1D(D_val)
    alg = RDMEMultigridFSP(dt=1.0, n_levels=1, n_max=220, save_every=5,
                            max_states=400_000, expand_depth=0,
                            coarse_r_depth=1, coarse_d_depth=1,
                            prune_tol=1e-5, weight_tol=1e-5, τ_pre=0.3)
    sol = solve_rdme_multigrid(m_mdl, fine_grid, u0_ic, (0.0, t_end_val), Float64[], alg)
    mu  = mean_voxel_counts(sol)
    states_f, probs_f = sol.snapshots[end]
    p_m = zeros(K+1)
    for (s,pr) in zip(states_f, probs_f)
        m_hi = count(n -> n > n_uns, Tuple(s)); p_m[m_hi+1] += pr
    end
    @printf("%s D=%.3f t=%.0f: E=%.0f,%.0f,%.0f,%.0f,%.0f,%.0f  |S|_max=%d\n",
        label, D_val, t_end_val, mu[end,:]..., maximum(sol.state_space_sizes))
    @printf("  P(m hi): %s\n", join([@sprintf("m=%d:%.3f",mi,p_m[mi+1]) for mi in 0:K], "  "))
    sol, mu
end

println("\n=== 2-voxel droplet (stable pair): varying D ===")
u0_2vox = CartesianIndex(n_low, n_low, n_high, n_high, n_low, n_low)
for D_val in [0.01, 0.1, 0.5, 1.0, 2.0]
    run_droplet(u0_2vox, D_val, 100.0; label="2vox")
end

println("\n=== 4-voxel droplet: should always expand ===")
u0_4vox = CartesianIndex(n_low, n_high, n_high, n_high, n_high, n_low)
run_droplet(u0_4vox, 0.1, 100.0; label="4vox")

println("\n=== Unstable IC (all at n_uns): spontaneous bifurcation ===")
u0_uns = CartesianIndex(n_uns, n_uns, n_uns, n_uns, n_uns, n_uns)
run_droplet(u0_uns, 0.1, 50.0; label="uns")
