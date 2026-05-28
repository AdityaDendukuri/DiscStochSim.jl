"""
RDME Multigrid FSP — Two-Species Bottleneck Benchmark

Model (K voxels, 2 species A and B):
  ∅  →^(k_prod) A     (zero-order, per voxel)
  A  →^(k_AB)   B     (first-order — the bottleneck, slow)
  B  →^(k_deg)  ∅     (first-order)
  A diffuses at D_A,  B diffuses at D_B

Stiffness: D_A, D_B >> k_AB (fast spatial mixing, slow conversion).
The V-cycle marginalizes over within-pair within-species distributions (each ~Binomial)
and solves the slow k_AB dynamics on the coarse grid.

Analytical steady state (linear network):
  μ_A = k_prod / k_AB  per voxel
  μ_B = k_prod / k_deg per voxel
  P = ∏_k Pois(μ_A)(nA_k) · Pois(μ_B)(nB_k)

Run from project root:
  julia --project examples/bottleneck_rdme.jl
"""

using DiscStochSim
using ExponentialUtilities
using Distributions
using Printf

# ─── Parameters ───────────────────────────────────────────────────────────────

const K      = 4        # fine voxels (must be even; start small)
const D_A    = 5.0      # fast A diffusion  (stiffness = d_intra / k_AB ≈ 3200)
const D_B    = 5.0      # fast B diffusion
const k_prod = 0.2      # ∅ → A per voxel
const k_AB   = 0.1      # A → B  (bottleneck: slow)
const k_deg  = 1.0      # B → ∅
const T_END  = 30.0
const N_MAX  = 10       # max molecules per species per fine voxel

# Analytical steady-state means
const μ_A = k_prod / k_AB   # = 10
const μ_B = k_prod / k_deg  # = 1

# Initial condition: near steady-state for A, no B
const u0 = CartesianIndex(ntuple(i -> isodd(i) ? round(Int, μ_A) : 0, Val(2K)))

model     = BottleneckModel1D(D_A, D_B, k_prod, k_AB, k_deg)
fine_grid = VoxelGrid(K, 1.0 / K, 0)
rates     = Float64[]

# Stiffness diagnostics
h       = fine_grid.dx
d_intra = D_A / h^2
d_inter = D_A / (2h)^2
cfl     = 1.0 / (2 * d_intra * N_MAX)

println("=" ^ 65)
println("Bottleneck RDME  K=$K  D_A=$D_A  D_B=$D_B  k_AB=$k_AB")
println("=" ^ 65)
@printf("  μ_A (steady state) = %.1f  per voxel\n", μ_A)
@printf("  μ_B (steady state) = %.1f  per voxel\n", μ_B)
@printf("  d_intra (fast)  = %g  per molecule\n", d_intra)
@printf("  d_inter (slow)  = %g  per molecule\n", d_inter)
@printf("  Explicit CFL    ≈ %.2e\n", cfl)
println()

# ─── TV distance from analytical steady state ─────────────────────────────────

function tv_analytical(sp::StateSpace{CartesianIndex{N}, Float64}) where {N}
    tv = 0.0
    for (i, state) in enumerate(sp.states)
        t = Tuple(state)
        p_true = 1.0
        for k in 1:K
            nA = t[2k-1]; nB = t[2k]
            p_true *= pdf(Poisson(μ_A), nA) * pdf(Poisson(μ_B), nB)
        end
        tv += abs(sp.probs[i] - p_true)
    end
    tv += abs(1.0 - sum(sp.probs))
    tv / 2
end

# ─── Multigrid FSP ────────────────────────────────────────────────────────────

println("── Multigrid FSP (K=$K → K=$(K÷2), 1 V-cycle level) ───────────────")

hierarchy = build_hierarchy(fine_grid, 1)

dt_mg   = 1.0
n_steps = round(Int, T_END / dt_mg)

function run_multigrid()
    sp = StateSpace{CartesianIndex{2K}, Float64}()
    add_state!(sp, u0, 1.0)
    t_cur = 0.0
    for _ in 1:n_steps
        dt_step = min(dt_mg, T_END - t_cur)
        sp = multi_level_vcycle_2s(sp, model, hierarchy, rates, t_cur, dt_step;
                                    τ_pre          = 0.0,
                                    coarse_n_max   = 2 * N_MAX,
                                    coarse_r_depth = 2,
                                    coarse_d_depth = 1,
                                    expand_coarse  = true,
                                    binom_tol      = 1e-8,
                                    weight_tol     = 1e-12,
                                    prune_tol      = 1e-11)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-11)
    end
    sp
end

t_mg  = @elapsed sp_mg = run_multigrid()
tv_mg = tv_analytical(sp_mg)
@printf("  Wall time:  %.2fs   |S| = %d   TV(π*) = %.5f   steps = %d\n\n",
        t_mg, length(sp_mg), tv_mg, n_steps)

# ─── Fine FSP (reference) ─────────────────────────────────────────────────────

println("── Fine FSP (no multigrid, reference) ────────────────────────────")

fine_sys = build_bottleneck_system(model, fine_grid)
bc       = state -> bottleneck_bc(state, N_MAX)

function run_fine()
    sp = StateSpace{CartesianIndex{2K}, Float64}()
    add_state!(sp, u0, 1.0)
    t_cur = 0.0
    for _ in 1:n_steps
        dt_step = min(dt_mg, T_END - t_cur)
        expand!(sp, fine_sys, bc; depth = 1)
        A, = build_generator(sp, fine_sys, rates, t_cur)
        sp.probs .= max.(expv(dt_step, A, sp.probs; m = 40), 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-11)
    end
    sp
end

t_fine  = @elapsed sp_fine = run_fine()
tv_fine = tv_analytical(sp_fine)
@printf("  Wall time:  %.2fs   |S| = %d   TV(π*) = %.5f   steps = %d\n\n",
        t_fine, length(sp_fine), tv_fine, n_steps)

# ─── Summary ──────────────────────────────────────────────────────────────────

println("=" ^ 65)
@printf("  Explicit CFL limit:  dt_CFL ≈ %.2e  (fine grid)\n", cfl)
@printf("  V-cycle step used:   dt     = %.2f   (%.0f× larger)\n",
        dt_mg, dt_mg / cfl)
println()
println("  Method        Wall time   |S|     TV(π*)   Steps")
println("  ──────────────────────────────────────────────────────")
@printf("  Multigrid     %6.2fs    %5d   %.5f   %d\n",
        t_mg, length(sp_mg), tv_mg, n_steps)
@printf("  Fine FSP      %6.2fs    %5d   %.5f   %d\n",
        t_fine, length(sp_fine), tv_fine, n_steps)
@printf("  Speedup: %.1f×\n", t_fine / t_mg)
println()

# ─── Per-voxel mean counts ────────────────────────────────────────────────────

function mean_counts_2s(sp::StateSpace{CartesianIndex{N}, Float64}) where {N}
    μA = zeros(K); μB = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K
            μA[k] += p * t[2k-1]
            μB[k] += p * t[2k]
        end
    end
    μA, μB
end

μA_mg, μB_mg   = mean_counts_2s(sp_mg)
μA_ref, μB_ref = mean_counts_2s(sp_fine)

println("  Final mean A per voxel:")
@printf("    Multigrid: %s\n", string(round.(μA_mg,  digits=2)))
@printf("    Fine FSP:  %s\n", string(round.(μA_ref, digits=2)))
@printf("    Exact:     %.2f (uniform)\n", μ_A)
println()
println("  Final mean B per voxel:")
@printf("    Multigrid: %s\n", string(round.(μB_mg,  digits=2)))
@printf("    Fine FSP:  %s\n", string(round.(μB_ref, digits=2)))
@printf("    Exact:     %.2f (uniform)\n", μ_B)
println("Done.")
