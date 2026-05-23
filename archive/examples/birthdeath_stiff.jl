"""
Birth-Death-Diffusion: stiff regime benchmark.

Shows that the V-cycle takes time steps governed by the slow reaction /
inter-voxel-diffusion dynamics, not by the fast intra-pair diffusion that
constrains explicit fine-grid solvers.

Model: ∅ →(k_b) A →(k_d) ∅, diffusion D.  K=8 fine voxels, h=1/8.

Stiffness numbers (D=1):
  d_intra = D/h²     = 64   (fast, within coarse pair; marginalized by 𝓡)
  d_inter = D/(2h)²  = 16   (slow, between coarse voxels; stays in coarse solve)
  ε       = 1/4              (lumpability ratio)
  Explicit CFL  dt_CFL ≈ 1 / (2·d_intra·n_max) ← fine-grid stability limit

V-cycle uses dt = O(1/k_rxn) — reaction timescale — orders of magnitude larger.
"""

using DiscStochSim
using ExponentialUtilities
using Distributions
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

const D    = 1.0
const K    = 8
const h    = 1.0 / K
const k_b  = 2.0
const k_d  = 1.0
const n_max = 30

model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

d_intra = D / h^2           # 64  per molecule
d_inter = D / (2h)^2        # 16  per molecule
ε       = lumpability_ratio(model, fine_grid)
cfl     = 1.0 / (2 * d_intra * n_max)

t_max   = 5.0
dt_vcyc = 0.5          # reaction-timescale dt for V-cycle
dt_fine = 0.5          # same dt for fine FSP
n_steps = round(Int, t_max / dt_vcyc)

# Point IC near μ* = k_b/k_d = 2 per voxel
u0 = CartesianIndex(ntuple(_ -> round(Int, k_b/k_d), K))

println("=" ^ 65)
println("Birth-Death-Diffusion  K=$K  D=$D  k_b=$k_b  k_d=$k_d")
println("=" ^ 65)
@printf("  d_intra (fast, within-pair)   = %5.1f  per molecule\n", d_intra)
@printf("  d_inter (slow, between pairs) = %5.1f  per molecule\n", d_inter)
@printf("  Lumpability ratio  ε          = %.4f\n", ε)
@printf("  Stiffness  d_intra / k_d      = %.0f×\n", d_intra / k_d)
@printf("  Explicit CFL limit (fine)     ≈ %.2e\n", cfl)
@printf("  V-cycle step size             =  %.2f  (%.0f× larger than CFL)\n\n",
        dt_vcyc, dt_vcyc / cfl)

# ─── analytical stationary means ─────────────────────────────────────────────

function stationary_means()
    d = d_intra
    A = zeros(K, K)
    b = fill(-k_b, K)
    for k in 1:K
        A[k,k] = -(k_d + (k == 1 || k == K ? d : 2d))
        k > 1 && (A[k, k-1] = d)
        k < K && (A[k, k+1] = d)
    end
    A \ b
end

μ_star = stationary_means()

function tv_vs_poisson(sp, μ)
    tv = 0.0
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        p_true = prod(pdf(Poisson(μ[k]), t[k]) for k in 1:K)
        tv += abs(sp.probs[i] - p_true)
    end
    tv += abs(1.0 - sum(sp.probs))
    tv / 2
end

# ─── V-cycle (large dt, reaction timescale) ───────────────────────────────────

println("── V-cycle  (dt=$(dt_vcyc), $(round(Int, dt_vcyc/cfl))× CFL) ──────────────────")

function run_vcycle()
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    hierarchy = build_hierarchy(fine_grid, 1)
    t_cur = 0.0
    for _ in 1:n_steps
        dt_step = min(dt_vcyc, t_max - t_cur)
        sp = multi_level_vcycle(sp, model, hierarchy, rates, t_cur, dt_step;
                                coarse_n_max  = 2 * n_max,
                                expand_coarse = true)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
    end
    sp
end

t_vc  = @elapsed sp_vc = run_vcycle()
tv_vc = tv_vs_poisson(sp_vc, μ_star)
@printf("  Wall time: %.2fs   |S|=%d   TV(π*) = %.5f   steps=%d\n\n",
        t_vc, length(sp_vc), tv_vc, n_steps)

# ─── Fine adaptive FSP (same dt — stiff expv must handle fast modes) ──────────

println("── Fine adaptive FSP  (dt=$(dt_fine), same steps) ───────────────────")

fine_sys = build_rdme_system(model, fine_grid)
fine_bc  = state -> rdme_bc(state, n_max)

function run_fine()
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    t_cur = 0.0
    for _ in 1:n_steps
        dt_step = min(dt_fine, t_max - t_cur)
        expand!(sp, fine_sys, fine_bc; depth=1)
        A, = build_generator(sp, fine_sys, rates, t_cur)
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=50), 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
    end
    sp
end

t_fine_wall = @elapsed sp_fine = run_fine()
tv_fine = tv_vs_poisson(sp_fine, μ_star)
@printf("  Wall time: %.2fs   |S|=%d   TV(π*) = %.5f   steps=%d\n\n",
        t_fine_wall, length(sp_fine), tv_fine, n_steps)

# ─── summary ──────────────────────────────────────────────────────────────────

println("=" ^ 65)
@printf("  Explicit CFL limit:  dt_CFL ≈ %.2e\n", cfl)
@printf("  V-cycle step used:   dt     = %.2f   (%.0f× larger)\n",
        dt_vcyc, dt_vcyc / cfl)
println()
println("  Method        Wall time   |S|      TV(π*)   Steps")
println("  ──────────────────────────────────────────────────────")
@printf("  V-cycle       %6.2fs    %5d    %.5f    %d\n",
        t_vc, length(sp_vc), tv_vc, n_steps)
@printf("  Fine FSP      %6.2fs    %5d    %.5f    %d\n",
        t_fine_wall, length(sp_fine), tv_fine, n_steps)
@printf("  Speedup: %.1f×\n", t_fine_wall / t_vc)
println()
@printf("  Fine generator spectral radius ≈ %.0f  (2·d_intra·n_max)\n",
        2 * d_intra * n_max)
@printf("  Coarse generator spectral radius ≈ %.0f  (2·d_inter·2n_max)\n",
        2 * d_inter * 2n_max)
println("Done.")
