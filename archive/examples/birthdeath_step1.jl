"""
Birth-Death-Diffusion: Step 1 — single coarse solve + prolong vs fine FSP.

K=4 fine voxels, K=2 coarse.  Linear network so coarse generator is exact.
We compare:
  - Fine FSP (K=4 direct)
  - Coarse solve + Binomial prolong (K=2 → K=4)

No V-cycle machinery — just restrict / solve / prolong manually.
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
h      = 1.0 / K        # voxel size = 0.25
N_MAX  = 20             # per-voxel truncation
T_EVAL = [0.5, 1.0, 2.0]

model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

d_fine   = D / h^2           # fine diffusion rate
d_coarse = D / (2h)^2        # coarse diffusion rate
println("="^60)
println("Birth-Death-Diffusion  K=$K  D=$D  k_b=$k_b  k_d=$k_d")
@printf("d_fine=%.2f   d_coarse=%.2f   ratio=%.2f\n", d_fine, d_coarse, d_fine/d_coarse)
println("="^60)

# ─── analytical stationary means ──────────────────────────────────────────────
# Linear system → stationary = product of Poissons with mean μ*
# μ* solves:  d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0  (reflecting BCs)

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

μ_fine   = stationary_means(K,   d_fine,   k_b, k_d)
μ_coarse = stationary_means(K÷2, d_coarse, k_b, k_d)
@printf("Fine   stationary means: %s\n", string(round.(μ_fine,   digits=3)))
@printf("Coarse stationary means: %s\n", string(round.(μ_coarse, digits=3)))

# ─── build fine FSP (full K=4 state space) ────────────────────────────────────

fine_sys = build_rdme_system(model, fine_grid)
bc_fine  = state -> rdme_bc(state, N_MAX)

sp_fine = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:N_MAX, n2 in 0:N_MAX, n3 in 0:N_MAX, n4 in 0:N_MAX
    add_state!(sp_fine, CartesianIndex(n1,n2,n3,n4), 0.0)
end
# IC: molecules concentrated in voxel 1
u0 = CartesianIndex(5, 0, 0, 0)
sp_fine.probs[sp_fine.index[u0]] = 1.0

A_fine, = build_generator(sp_fine, fine_sys, rates, 0.0)
println("\nFine   state space |S| = $(length(sp_fine))")

# ─── build coarse FSP (K=2) ───────────────────────────────────────────────────

coarse_sys = build_coarse_system(model, coarse_grid, fine_grid)
bc_coarse  = state -> rdme_bc(state, 2*N_MAX)

sp_coarse = StateSpace{CartesianIndex{2}, Float64}()
for nc1 in 0:2*N_MAX, nc2 in 0:2*N_MAX
    add_state!(sp_coarse, CartesianIndex(nc1, nc2), 0.0)
end
# IC: restrict fine IC to coarse  (nc1 = n1+n2 = 5+0, nc2 = n3+n4 = 0+0)
u0_coarse = CartesianIndex(5, 0)
sp_coarse.probs[sp_coarse.index[u0_coarse]] = 1.0

A_coarse, = build_generator(sp_coarse, coarse_sys, rates, 0.0)
println("Coarse state space |S| = $(length(sp_coarse))")
@printf("State space reduction:   %.1fx\n\n", length(sp_fine)/length(sp_coarse))

# ─── TV helper ────────────────────────────────────────────────────────────────

function tv_vs_fine(sp_approx, sp_fine_ref, p_fine_ref)
    tv = 0.0
    for (i, s) in enumerate(sp_fine_ref.states)
        p_ref = p_fine_ref[i]
        idx   = get(sp_approx.index, s, 0)
        p_app = idx != 0 ? sp_approx.probs[idx] : 0.0
        tv   += abs(p_ref - p_app)
    end
    # mass outside fine reference support
    for (i, s) in enumerate(sp_approx.states)
        get(sp_fine_ref.index, s, 0) != 0 && continue
        tv += abs(sp_approx.probs[i])
    end
    tv / 2
end

function voxel_means(sp, K)
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

# ─── evolve and compare ───────────────────────────────────────────────────────

println("  t    Method           ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩   TV vs fine FSP")
println("  ─────────────────────────────────────────────────────────────────")

let p_fine = copy(sp_fine.probs), p_coarse = copy(sp_coarse.probs), t_prev = 0.0
for t in T_EVAL
    dt = t - t_prev

    # Fine FSP
    p_fine = max.(expv(dt, A_fine, p_fine; m=30), 0.0)
    p_fine ./= sum(p_fine)
    sp_fine.probs .= p_fine
    mu_f = voxel_means(sp_fine, K)

    # Coarse solve
    p_coarse = max.(expv(dt, A_coarse, p_coarse; m=30), 0.0)
    p_coarse ./= sum(p_coarse)
    sp_coarse.probs .= p_coarse

    # Prolong coarse → fine (Binomial)
    sp_prol = prolong(sp_coarse, Val(K); weight_tol=1e-14, binom_tol=1e-4)
    renormalize!(sp_prol)
    mu_c = voxel_means(sp_prol, K)
    tv   = tv_vs_fine(sp_prol, sp_fine, p_fine)

    @printf("  %4.1f  Fine FSP (K=4)    %5.2f %5.2f %5.2f %5.2f   —\n",
            t, mu_f[1], mu_f[2], mu_f[3], mu_f[4])
    @printf("  %4s  Coarse+Prolong   %5.2f %5.2f %5.2f %5.2f   %.5f   (|S_prol|=%d)\n\n",
            "", mu_c[1], mu_c[2], mu_c[3], mu_c[4], tv, length(sp_prol))

    t_prev = t
end
end  # let

println("Fine   |S| = $(length(sp_fine))  (per-voxel N_MAX=$N_MAX, K=$K)")
println("Coarse |S| = $(length(sp_coarse))  (per-voxel 2*N_MAX=$(2*N_MAX), K=$(K÷2))")
@printf("Reduction:   %.1fx\n", length(sp_fine)/length(sp_coarse))
println("\nDone.")
