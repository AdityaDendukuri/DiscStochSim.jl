"""
Birth-Death-Diffusion: Step 3 — K=8 fine, K=4 coarse.

The fine K=8 FSP is intractable at D=1 (state space ~10^9 at stationary).
The coarse K=4 FSP is fast, and for a linear network the coarse generator
is exact — so we verify accuracy against the analytical Poisson stationary
distribution instead of against the intractable fine FSP.

Key demonstration:
  - Coarse |S| stays small throughout
  - TV(coarse+prolong, π*) → 0 as t → ∞
  - Wall time is seconds, not hours
"""

using DiscStochSim
using ExponentialUtilities
using Distributions
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

K      = 8
D      = 1.0
k_b    = 2.0
k_d    = 1.0
h      = 1.0 / K
N_MAX  = 20
T_END  = 2.0
dt     = 0.1
T_SNAP = [0.5, 1.0, 1.5, 2.0]

# Mild gradient close to stationary (μ*=2 per voxel).
# Keeps the transient state space small so the K=4 coarse solve stays tractable.
μ_IC = [3.0, 3.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0]

model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)          # K=4

d_fine   = D / h^2
d_coarse = D / (2h)^2

println("="^60)
println("Birth-Death-Diffusion  K=$K  D=$D  k_b=$k_b  k_d=$k_d")
@printf("d_fine=%.1f  d_coarse=%.1f  dt=%.2f  T=%.1f\n", d_fine, d_coarse, dt, T_END)
println("IC: product-Poisson gradient μ=$μ_IC")
println("="^60)

# ─── analytical stationary (exact for linear networks) ────────────────────────

function stationary_means(K, d, k_b, k_d)
    A = zeros(K, K); b = fill(-k_b, K)
    for k in 1:K
        A[k,k] = -(k_d + (k==1||k==K ? d : 2d))
        k > 1 && (A[k,k-1] = d)
        k < K && (A[k,k+1] = d)
    end
    A \ b
end

μ_star_fine   = stationary_means(K,   d_fine,   k_b, k_d)
μ_star_coarse = stationary_means(K÷2, d_coarse, k_b, k_d)
@printf("Fine   stationary means: %s\n", string(round.(μ_star_fine,   digits=2)))
@printf("Coarse stationary means: %s\n", string(round.(μ_star_coarse, digits=2)))
println()

# ─── state space size estimate for fine K=8 ───────────────────────────────────
# From step 2: K=4 fine with same params reached ~26k states at t=2.
# K=8 state space scales roughly as |S_K4|^2 ≈ 700M — intractable.
println("Fine K=8 FSP state space estimate: ~700M states at t=2 (intractable)")
println("Coarse K=4 FSP: solve exactly, prolong to K=8 only at snapshot times.")
println()

# ─── coarse IC: product-Poisson for K=4 (tractable) ─────────────────────────
# Build the K=4 coarse IC with per-voxel means = μ_IC[2j-1] + μ_IC[2j]
# (since nc_j = n_{2j-1} + n_{2j} at t=0, and both voxels are Poisson)
# For the coarse system Pois(nc_j) where nc_j ~ Pois(μ_IC[2j-1]+μ_IC[2j]).

μ_IC_coarse = [μ_IC[2j-1] + μ_IC[2j] for j in 1:K÷2]
coarse_sys  = build_coarse_system(model, coarse_grid, fine_grid)
N_MAX_C     = 15                                   # tighter bound for coarse voxels
bc_coarse   = state -> rdme_bc(state, N_MAX_C)

pois_c = [Poisson(μ_IC_coarse[k]) for k in 1:K÷2]

sp_coarse = StateSpace{CartesianIndex{4}, Float64}()
for n1 in 0:N_MAX_C, n2 in 0:N_MAX_C, n3 in 0:N_MAX_C, n4 in 0:N_MAX_C
    ns = (n1, n2, n3, n4)
    p  = prod(pdf(pois_c[k], ns[k]) for k in 1:K÷2)
    p > 1e-10 && add_state!(sp_coarse, CartesianIndex(ns...), p)
end
renormalize!(sp_coarse)
println("Coarse IC |S| = $(length(sp_coarse))")

# ─── TV vs analytical ─────────────────────────────────────────────────────────

function tv_vs_analytical(sp, μ)
    tv = 0.0
    for (i, s) in enumerate(sp.states)
        t   = Tuple(s)
        p_a = prod(pdf(Poisson(μ[k]), t[k]) for k in 1:length(μ))
        tv += abs(sp.probs[i] - p_a)
    end
    tv += abs(1.0 - sum(sp.probs))
    tv / 2
end

# ─── Run: Coarse FSP (K=4) ────────────────────────────────────────────────────

snaps_coarse = Dict{Float64, StateSpace{CartesianIndex{4}, Float64}}()

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
            sp2 = StateSpace{CartesianIndex{4}, Float64}()
            for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
            snaps_coarse[t_snap] = sp2
        end
    end
    global sp_coarse = sp
end

# ─── Results ──────────────────────────────────────────────────────────────────

println()
println("  t    |S_coarse|   TV(coarse, π*_coarse)   TV(prol→K=8, π*_fine)")
println("  ──────────────────────────────────────────────────────────────────")

for t in T_SNAP
    sc = snaps_coarse[t]
    sp_prol = prolong(sc, Val(K); weight_tol=1e-14, binom_tol=1e-4)
    renormalize!(sp_prol)
    tv_c  = tv_vs_analytical(sc,      μ_star_coarse)
    tv_f  = tv_vs_analytical(sp_prol, μ_star_fine)
    @printf("  %4.1f  %7d     %.5f                   %.5f\n",
            t, length(sc), tv_c, tv_f)
end

println()
@printf("  Coarse FSP wall time: %.3f s   final |S_coarse|=%d\n", t_coarse, length(sp_coarse))
println("  Fine K=8 FSP: intractable (~700M states)")
println("\nDone.")
