"""
RDME Multigrid FSP — Benchmark 1: 1D Birth-Death-Diffusion

Model (K voxels, 1 species A):
  ∅  →  A   rate k_b  (per voxel, zero-order)
  A  →  ∅   rate k_d  (first-order, per molecule)
  A_k ⇌ A_{k±1}  rate d = D/h²  (diffusion)

This is a linear (zero + first-order) network, so:
  1. The stationary distribution is a product of independent Poissons
     π(x) = ∏_k Pois(μ_k), where μ solves the discrete diffusion-reaction
     equation  d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0.
  2. The coarse-level generator is EXACT (Proposition 4.1 of the paper).
  3. The V-cycle introduces only splitting errors → vanish as τ_pre → 0.

We compare:
  - Multigrid FSP (2-level): K=8 fine voxels, K=4 coarse voxels
  - Fine-level FSP (no MG): full K=8 voxel solve
  - Analytical stationary distribution (Poisson)

Run from the project root:
  julia --project examples/rdme_birthdeath_1d.jl
"""

using DiscStochSim
using LinearAlgebra
using SparseArrays
using Distributions   # for Poisson

# ─── Model parameters ─────────────────────────────────────────────────────────

const K     = 8       # fine voxels (must be even)
const L     = 1       # number of coarsening levels  (gives K/2 = 4 coarse voxels)
const D     = 1.0     # diffusion coefficient
const k_b   = 2.0     # birth rate per voxel
const k_d   = 1.0     # death rate per molecule
const T_END = 5.0     # simulation end time
const N_MAX = 30      # FSP molecule-count truncation per voxel

# Initial condition: all molecules in voxel 1
const u0 = CartesianIndex(ntuple(k -> k == 1 ? 5 : 0, K))

# ─── Analytical stationary distribution ───────────────────────────────────────

"""
Solve the steady-state diffusion-reaction equation for mean counts μ:
  d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0
with reflecting boundary conditions μ_0 = μ_1, μ_{K+1} = μ_K.

Returns: mean count per voxel at steady state.
"""
function stationary_means(K, D, dx, k_b, k_d)
    d = D / dx^2
    # Steady-state: d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0
    # Rearranged as A*μ = -k_b (note sign: rhs is -k_b, not +k_b)
    A = zeros(K, K)
    b = fill(-k_b, K)
    for k in 1:K
        A[k, k] = -(k_d + (k == 1 || k == K ? d : 2d))
        if k > 1; A[k, k-1] = d; end
        if k < K; A[k, k+1] = d; end
    end
    A \ b
end

"""
Analytical stationary TV distance from the numerical solution.
For linear networks the stationary distribution is ∏_k Pois(μ_k).
"""
function analytical_tv(sp::StateSpace{CartesianIndex{K}, Float64},
                        mu::Vector{Float64}) where K
    tv = 0.0
    # Iterate over all states in the active set
    for (i, state) in enumerate(sp.states)
        t   = Tuple(state)
        p_analytical = prod(pdf(Poisson(mu[k]), t[k]) for k in 1:K)
        tv += abs(sp.probs[i] - p_analytical)
    end
    # Also accumulate probability mass outside active set
    mass_in = sum(sp.probs)
    tv += abs(1.0 - mass_in)   # probability on truncated states
    tv / 2
end

# ─── Setup ────────────────────────────────────────────────────────────────────

fine_grid   = VoxelGrid(K, 1.0 / K, 0)
coarse_grid = coarsen(fine_grid)
model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]   # all rates baked into propensities via closure

println("=" ^ 60)
println("RDME Multigrid FSP — Birth-Death-Diffusion Benchmark")
println("K=$K fine voxels, D=$D, k_b=$k_b, k_d=$k_d, T=$T_END")
println("Lumpability ratio ε = ", lumpability_ratio(model, fine_grid))
println("=" ^ 60)

# ─── Analytical stationary means ──────────────────────────────────────────────

μ_star = stationary_means(K, D, fine_grid.dx, k_b, k_d)
println("\nAnalytical stationary mean counts per voxel:")
println("  μ* = ", round.(μ_star, digits=3))

# ─── Run: Multigrid FSP ───────────────────────────────────────────────────────

println("\n── Multigrid FSP (2-level) ──────────────────────────────")
alg_mg = RDMEMultigridFSP(
    dt         = 0.05,
    τ_pre      = 0.01,   # small pre-smoothing: ≈ 1 intra-pair diffusion event
    τ_post     = 0.01,
    n_max      = N_MAX,
    krylov_m   = 20,
    save_every = 10,
)

t_mg = @elapsed sol_mg = solve_rdme_multigrid(model, fine_grid, u0,
                                               (0.0, T_END), rates, alg_mg)

println("  Wall time:           $(round(t_mg, digits=2)) s")
println("  Steps taken:         $(length(sol_mg.t) - 1)")
println("  Final |S| (fine):    $(sol_mg.state_space_sizes[end])")

# Final mean counts
mu_mg = mean_voxel_counts(sol_mg)
println("  Final mean counts:   ", round.(mu_mg[end, :], digits=3))

# ─── Run: Fine FSP (reference, no multigrid) ──────────────────────────────────

println("\n── Fine-level FSP (no multigrid, reference) ─────────────")

full_sys = build_rdme_system(model, fine_grid)
bc       = state -> rdme_bc(state, N_MAX)

sp_ref = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp_ref, u0, 1.0)
expand!(sp_ref, full_sys, bc; depth = 1)

t_ref_total = 0.0
step_ref    = 0
t_now       = 0.0
dt_ref      = 0.05

t_fine = @elapsed begin
    while t_now < T_END
        dt_step = min(dt_ref, T_END - t_now)
        A_fine, = build_generator(sp_ref, full_sys, rates, t_now)
        sp_ref.probs .= expv(dt_step, A_fine, sp_ref.probs; m=20)
        t_now  += dt_step
        step_ref += 1
        renormalize!(sp_ref)
        prune_threshold!(sp_ref, 1e-12)
        expand!(sp_ref, full_sys, bc; depth = 1)
    end
end

println("  Wall time:           $(round(t_fine, digits=2)) s")
println("  Steps taken:         $step_ref")
println("  Final |S| (fine):    $(length(sp_ref))")

mu_ref = zeros(K)
for (state, p) in zip(sp_ref.states, sp_ref.probs)
    t = Tuple(state)
    for k in 1:K; mu_ref[k] += p * t[k]; end
end
println("  Final mean counts:   ", round.(mu_ref, digits=3))

# ─── Comparison ───────────────────────────────────────────────────────────────

println("\n── Accuracy comparison at t = $T_END ────────────────────")

# TV distance from analytical stationary distribution
tv_mg  = analytical_tv(StateSpace{CartesianIndex{K}, Float64}() |>
                        sp -> (for (s,p) in zip(sol_mg.snapshots[end]...)
                                   add_state!(sp,s,p); end; sp),
                        μ_star)

# Build sp from final MG snapshot for TV calculation
sp_mg_final = StateSpace{CartesianIndex{K}, Float64}()
for (s, p) in zip(sol_mg.snapshots[end][1], sol_mg.snapshots[end][2])
    add_state!(sp_mg_final, s, p)
end
tv_mg  = analytical_tv(sp_mg_final, μ_star)
tv_ref = analytical_tv(sp_ref,      μ_star)

println("  TV(multigrid, analytical):  $(round(tv_mg, digits=5))")
println("  TV(fine FSP,  analytical):  $(round(tv_ref, digits=5))")
println("  TV(multigrid, fine FSP):    ",
        round(abs(tv_mg - tv_ref), digits=5))

println("\n── Lumpability diagnostics ──────────────────────────────")
phi = per_molecule_flux(sp_ref, model, fine_grid)
println("  Per-molecule flux φ_{k,k+1}:")
for (k, φ) in enumerate(phi)
    println("    k=$(k)→$(k+1): φ = $(round(φ, digits=4))")
end
println("  (Uniform → no spatial bottlenecks → all pairs lumpable)")

println("\nDone.")
