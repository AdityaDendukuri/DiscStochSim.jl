# Schlögl bistability: all voxels start at unstable fixed point n_un
# Watch probability band SPLIT into two bands (n_lo and n_hi) over time
# This directly shows bistability in the RDME

using DiscStochSim
using ExponentialUtilities
using Printf

D = 0.1; K = 4; h = 1.0/K; n_max = 200
c3    = 15.0
model = SchloglModel1D(D, 0.028, 3.2e-4, c3, 1.0)
fine_grid = VoxelGrid(K, h, 0)
sys       = build_schlogl_rdme_system(model, fine_grid)
bc        = s -> rdme_bc(s, n_max)
rates     = Float64[]

fp = schlogl_fixed_points(model; n_max=300)
n_lo, n_un, n_hi = fp[1], fp[2], fp[3]
@printf("c3=%.1f  n_lo=%d  n_un=%d  n_hi=%d  d=%.1f\n",
        c3, n_lo, n_un, n_hi, D/h^2)

# IC: ALL voxels at unstable fixed point
# Stochastic fluctuations will tip each voxel to n_lo or n_hi
u0 = CartesianIndex(n_un, n_un, n_un, n_un)
@printf("IC: all voxels at n_un=%d\n\n", n_un)

sp = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp, u0, 1.0)

dt = 0.1; T = 5.0; n_steps = round(Int, T/dt)
t_cur = 0.0

function vmeans(states, probs, K)
    mu = zeros(K)
    for (i,s) in enumerate(states)
        tv=Tuple(s); p=probs[i]
        for k in 1:K; mu[k]+=p*tv[k]; end
    end
    mu
end

println("  t     |S|    v1      v2      v3      v4")
@printf("  %.1f  %5d  %6.1f  %6.1f  %6.1f  %6.1f\n",
        0.0, 1, Float64(n_un), Float64(n_un), Float64(n_un), Float64(n_un))

for step in 1:n_steps
    global t_cur, sp
    dts = min(dt, T-t_cur)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t_cur)
    sp.probs .= max.(expv(dts, A, sp.probs; m=40), 0.0)
    t_cur += dts
    renormalize!(sp); prune_threshold!(sp, 1e-6)

    if mod(step, 5) == 0
        mu = vmeans(sp.states, sp.probs, K)
        @printf("  %.1f  %5d  %6.1f  %6.1f  %6.1f  %6.1f\n",
                t_cur, length(sp), mu...)
        length(sp) > 100_000 && (println("  too large, stopping"); break)
    end
end
println("Done.")
