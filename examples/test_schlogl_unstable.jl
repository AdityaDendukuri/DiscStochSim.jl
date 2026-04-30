# K=8 Schlögl, all voxels at n_un: watch bistability develop spatially
# Use K=4 coarse solver, prolong at snapshots

using DiscStochSim
using ExponentialUtilities
using Printf

D = 0.005; K = 8; h = 1.0/K; n_max = 200
model = SchloglModel1D(D)
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)   # K=4
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=300)
n_lo, n_un, n_hi = fp[1], fp[2], fp[3]
@printf("n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

# IC: v1=n_un (unstable), v2=n_hi, v3..v8=n_lo
# Pair 1 = (n_un, n_hi): nc1 = n_un+n_hi = 234 (above coarse unstable ~144 -> tips to hi)
# Pairs 2-4 = (n_lo,n_lo): nc = 62, stable
# Watch v1 get pulled from n_un to n_hi by its neighbor -> bistability developing
u0_fine = CartesianIndex(n_un, n_hi, n_lo, n_lo, n_lo, n_lo, n_lo, n_lo)
nc1 = n_un + n_hi  # = 234
u0_c = CartesianIndex(nc1, 2*n_lo, 2*n_lo, 2*n_lo)
@printf("Fine IC: v1=%d (n_un), v2=%d (n_hi), v3..v8=%d (n_lo)\n", n_un, n_hi, n_lo)
@printf("Coarse IC: nc = [%d, %d, %d, %d]\n\n", Tuple(u0_c)...)

sys_c = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
bc_c  = s -> rdme_bc(s, 2*n_max)

sp_c = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_c, u0_c, 1.0)

dt = 0.5; T_check = [0.5, 1.0, 2.0, 5.0, 10.0]
n_steps_total = round(Int, maximum(T_check)/dt)
t_cur = 0.0

function coarse_means(sp)
    mu = zeros(4)
    for (i,s) in enumerate(sp.states); tv=Tuple(s); p=sp.probs[i]
        for j in 1:4; mu[j]+=p*tv[j]; end; end
    mu
end

println("  t      |S_c|   nc1    nc2    nc3    nc4   (nc should split to ~$(2*n_lo) or ~$(2*n_hi))")

for step in 1:n_steps_total
    global t_cur, sp_c
    dts = min(dt, maximum(T_check)-t_cur)
    expand!(sp_c, sys_c, bc_c; depth=1)
    A, = build_generator(sp_c, sys_c, rates, t_cur)
    sp_c.probs .= max.(expv(dts, A, sp_c.probs; m=40), 0.0)
    t_cur += dts
    renormalize!(sp_c); prune_threshold!(sp_c, 1e-7)

    if any(abs(t_cur - t) < dt/2 for t in T_check)
        mu = coarse_means(sp_c)
        @printf("  %.1f   %5d   %.1f   %.1f   %.1f   %.1f\n",
                t_cur, length(sp_c), mu[1], mu[2], mu[3], mu[4])
        length(sp_c) > 100_000 && (println("  too large"); break)
    end
end
println("Done.")
