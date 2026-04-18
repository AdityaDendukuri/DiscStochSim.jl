"""
Birth-Death-Diffusion: Step 2 — full adaptive time loop, wall time comparison.

Same K=4 / K=2 setup as step 1, but now both solvers run an adaptive
expand → solve → prune loop over [0, T_END].

At each snapshot we prolong the coarse distribution and compute TV vs the
fine FSP.  We report:
  - |S| at each snapshot for both solvers
  - TV(coarse+prolong, fine) at each snapshot
  - Wall time for each solver
"""

using DiscStochSim
using ExponentialUtilities
using Distributions
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

K      = 4
D      = 1.0
k_b    = 2.0
k_d    = 1.0
h      = 1.0 / K
N_MAX  = 20
T_END  = 2.0
dt     = 0.1          # time step per iteration
T_SNAP = [0.5, 1.0, 1.5, 2.0]

model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

d_fine   = D / h^2
d_coarse = D / (2h)^2

println("="^60)
println("Birth-Death-Diffusion  K=$K  D=$D  k_b=$k_b  k_d=$k_d")
@printf("d_fine=%.1f  d_coarse=%.1f  dt=%.2f  T=%.1f\n", d_fine, d_coarse, dt, T_END)
println("IC: product-Poisson gradient μ=[4,2,1,1]")
println("="^60)

# ─── systems ──────────────────────────────────────────────────────────────────

fine_sys   = build_rdme_system(model, fine_grid)
coarse_sys = build_coarse_system(model, coarse_grid, fine_grid)
bc_fine    = state -> rdme_bc(state, N_MAX)
bc_coarse  = state -> rdme_bc(state, 2*N_MAX)

# Smooth IC: spatial gradient Pois(4)×Pois(2)×Pois(1)×Pois(1)
# Restricts cleanly to coarse; Binomial prolong is a good approximation
# because the fine distribution is already spread across many states.
μ_IC = [4.0, 2.0, 1.0, 1.0]   # voxel means at t=0

# ─── analytical stationary distribution ──────────────────────────────────────
# Linear system → stationary = ∏_k Pois(μ*_k), where μ* solves the steady-state
# discrete diffusion-reaction equation with reflecting BCs.

function stationary_means(K, d, k_b, k_d)
    A = zeros(K, K)
    b = fill(-k_b, K)
    for k in 1:K
        A[k,k] = -(k_d + (k==1 || k==K ? d : 2d))
        k > 1 && (A[k,k-1] = d)
        k < K && (A[k,k+1] = d)
    end
    A \ b
end

μ_star = stationary_means(K, d_fine, k_b, k_d)
@printf("Analytical stationary means: %s\n\n", string(round.(μ_star, digits=3)))

function tv_vs_analytical(sp, μ)
    tv = 0.0
    for (i, s) in enumerate(sp.states)
        t   = Tuple(s)
        p_a = prod(pdf(Poisson(μ[k]), t[k]) for k in 1:length(μ))
        tv += abs(sp.probs[i] - p_a)
    end
    tv += abs(1.0 - sum(sp.probs))   # mass outside truncated support
    tv / 2
end

# ─── TV helper ────────────────────────────────────────────────────────────────

function tv(sp_a, sp_b)
    total = 0.0
    for (i, s) in enumerate(sp_b.states)
        pa = get(sp_a.index, s, 0)
        pb = sp_b.probs[i]
        total += abs((pa != 0 ? sp_a.probs[pa] : 0.0) - pb)
    end
    for (i, s) in enumerate(sp_a.states)
        get(sp_b.index, s, 0) != 0 && continue
        total += abs(sp_a.probs[i])
    end
    total / 2
end

# ─── Run: Fine FSP ────────────────────────────────────────────────────────────

# Build fine IC from product-Poisson
sp_fine = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:N_MAX, n2 in 0:N_MAX, n3 in 0:N_MAX, n4 in 0:N_MAX
    p = prod(pdf(Poisson(μ_IC[k]), (n1,n2,n3,n4)[k]) for k in 1:K)
    p > 1e-14 && add_state!(sp_fine, CartesianIndex(n1,n2,n3,n4), p)
end
renormalize!(sp_fine)

snaps_fine = Dict{Float64, StateSpace{CartesianIndex{K}, Float64}}()

t_fine = @elapsed let sp = sp_fine, t_cur = 0.0
    while t_cur < T_END - 1e-10
        dt_step = min(dt, T_END - t_cur)
        expand!(sp, fine_sys, bc_fine; depth=1)
        A, = build_generator(sp, fine_sys, rates, t_cur)
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=30), 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        t_snap = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        if abs(t_cur - t_snap) < dt/2 && !haskey(snaps_fine, t_snap)
            sp2 = StateSpace{CartesianIndex{K}, Float64}()
            for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
            snaps_fine[t_snap] = sp2
        end
    end
    global sp_fine = sp
end

# ─── Run: Coarse FSP ──────────────────────────────────────────────────────────

# Coarse IC: restrict the fine IC exactly
sp_coarse = restrict(sp_fine)

snaps_coarse = Dict{Float64, StateSpace{CartesianIndex{2}, Float64}}()

t_coarse = @elapsed let sp = sp_coarse, t_cur = 0.0
    while t_cur < T_END - 1e-10
        dt_step = min(dt, T_END - t_cur)
        expand!(sp, coarse_sys, bc_coarse; depth=1)
        A, = build_generator(sp, coarse_sys, rates, t_cur)
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=30), 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        t_snap = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        if abs(t_cur - t_snap) < dt/2 && !haskey(snaps_coarse, t_snap)
            sp2 = StateSpace{CartesianIndex{2}, Float64}()
            for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
            snaps_coarse[t_snap] = sp2
        end
    end
    global sp_coarse = sp
end

# ─── Compare ──────────────────────────────────────────────────────────────────

println()
println("  t    |S_fine|  |S_coarse|  reduc   TV(fine,π*)  TV(coarse+prol,fine)  TV(coarse+prol,π*)")
println("  ────────────────────────────────────────────────────────────────────────────────────────")

for t in T_SNAP
    sf = snaps_fine[t];   nf = length(sf)
    sc = snaps_coarse[t]; nc = length(sc)
    sp_prol = prolong(sc, Val(K); weight_tol=1e-14, binom_tol=1e-4)
    renormalize!(sp_prol)
    tv_fc   = tv(sp_prol, sf)
    tv_f_pi = tv_vs_analytical(sf, μ_star)
    tv_c_pi = tv_vs_analytical(sp_prol, μ_star)
    @printf("  %4.1f  %7d    %7d    %4.0fx    %.5f     %.5f               %.5f\n",
            t, nf, nc, nf/nc, tv_f_pi, tv_fc, tv_c_pi)
end

println()
@printf("  Fine   FSP wall time: %.3f s   final |S|=%d\n", t_fine,   length(sp_fine))
@printf("  Coarse FSP wall time: %.3f s   final |S|=%d\n", t_coarse, length(sp_coarse))
@printf("  Speedup: %.1fx\n", t_fine / t_coarse)
println("\nDone.")
