# Quick feasibility test: K=16 V-cycle with larger D
using DiscStochSim
using ExponentialUtilities
using Printf

D = 0.02; K = 16; h = 1.0/K; n_max = 180
model     = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)   # K=8
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=250)
n_lo, _, n_hi = fp
@printf("D=%.3f  K=%d  h=%.4f  d=D/h²=%.2f\n", D, K, h, D/h^2)
@printf("n_lo=%d  n_hi=%d\n", n_lo, n_hi)

# IC: front at pair 5 (voxels 9-10), so 8 lo voxels then 8 hi voxels
u0 = CartesianIndex(ntuple(k -> k <= 8 ? n_lo : n_hi, K))
@printf("IC: v1-v8=%d, v9-v16=%d  (front at pair 5)\n", n_lo, n_hi)

# Compute dynamic-pi for fine voxel size h
println("\nComputing dynamic-pi table (h=$(h))...")
pair_sys = build_schlogl_rdme_system(model, VoxelGrid(2, h, 0))
sp_pair  = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair, CartesianIndex(n1,n2), 0.0)
end
sp_pair.probs[sp_pair.index[CartesianIndex(n_lo,n_hi)]] = 1.0
A_p, = build_generator(sp_pair, pair_sys, rates, 0.0)
pi_table = compute_dynamic_pi(sp_pair, A_p; n_max=n_max)
println("  Done.")

# Run 5 V-cycle steps and check state space size / timing
sp = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp, u0, 1.0)
dt = 0.2; t = 0.0

println("\nRunning K=$K V-cycle steps...")
for step in 1:15
    global sp, t
    t1 = time()
    sp = two_level_vcycle_schlogl_injection(
        sp, model, fine_grid, coarse_grid, pi_table, rates, t, dt;
        coarse_n_max=2*n_max)
    t += dt
    renormalize!(sp); prune_threshold!(sp, 1e-7)
    mu = zeros(K)
    for (i,s) in enumerate(sp.states); tv=Tuple(s); p=sp.probs[i]
        for k in 1:K; mu[k] += p*tv[k]; end; end
    @printf("  step %d  t=%.1f  |S|=%d  %.2fs  v8=%.1f v9=%.1f\n",
            step, t, length(sp), time()-t1, mu[8], mu[9])
end
println("\nFeasibility check done.")
