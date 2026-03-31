"""
Multigrid FSP — Systematic Debugging Experiments

Run from project root:
  julia --project examples/multigrid_experiments.jl

Experiments:
  EXP 1: Restrict/prolong roundtrip (mass conservation, marginal preservation)
  EXP 2: Coarse system stationary means vs fine system
         → Checks the birth-rate coarsening correctness
  EXP 3: Single V-cycle accuracy on K=2 (minimal case, analytic reference)
  EXP 4: Long-time convergence of multigrid vs fine FSP (K=4, K=8)
  EXP 5: Error vs dt for multigrid (how splitting error scales)
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Distributions   # Poisson
using Printf

# ────────────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────────────

"""Build a StateSpace from explicit (state → prob) pairs."""
function make_sp(::Val{K}, pairs) where K
    sp = StateSpace{CartesianIndex{K}, Float64}()
    for (t, p) in pairs
        add_state!(sp, CartesianIndex(t), p)
    end
    sp
end

"""Total variation distance between two StateSpaces (same fine dimension)."""
function tv_distance(sp1::StateSpace{CartesianIndex{K},Float64},
                     sp2::StateSpace{CartesianIndex{K},Float64}) where K
    all_states = union(Set(sp1.states), Set(sp2.states))
    tv = 0.0
    for s in all_states
        p1 = get(sp1.index, s, 0) > 0 ? sp1.probs[sp1.index[s]] : 0.0
        p2 = get(sp2.index, s, 0) > 0 ? sp2.probs[sp2.index[s]] : 0.0
        tv += abs(p1 - p2)
    end
    tv / 2
end

"""Mean per-voxel counts from a StateSpace."""
function voxel_means(sp::StateSpace{CartesianIndex{K},Float64}) where K
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

"""
Analytical stationary means for 1D birth-death-diffusion.

Steady-state equation for mean μ_k (per voxel):
  d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0
with reflecting BCs: μ_0 = μ_1, μ_{K+1} = μ_K.

Returns positive values (e.g. k_b/k_d for uniform diffusion).
"""
function stationary_means(K, D, dx, k_b, k_d)
    d = D / dx^2
    # A*μ = b  where A[k,k] = -(k_d + 2d) etc., b = -k_b  →  μ = A\b
    A = zeros(K, K)
    b = fill(-k_b, K)   # sign: steady-state rearranges to A*μ = -k_b
    for k in 1:K
        A[k, k] = -(k_d + (k == 1 || k == K ? d : 2d))
        if k > 1; A[k, k-1] = d; end
        if k < K; A[k, k+1] = d; end
    end
    A \ b
end

"""TV from product-Poisson stationary distribution."""
function analytical_tv(sp::StateSpace{CartesianIndex{K},Float64}, mu) where K
    tv = 0.0
    for i in eachindex(sp.states)
        t = Tuple(sp.states[i])
        p_an = prod(pdf(Poisson(mu[k]), t[k]) for k in 1:K)
        tv += abs(sp.probs[i] - p_an)
    end
    tv += abs(1.0 - sum(sp.probs))
    tv / 2
end

"""Run a reference fine-FSP solve and return the final StateSpace."""
function fine_fsp_solve(model, fine_grid, u0, tf; dt=0.05, n_max=30, krylov_m=20)
    K = fine_grid.n_voxels
    full_sys = build_rdme_system(model, fine_grid)
    bc = s -> rdme_bc(s, n_max)
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, full_sys, bc; depth=1)
    t = 0.0
    rates = Float64[]
    while t < tf
        dt_step = min(dt, tf - t)
        A, = build_generator(sp, full_sys, rates, t)
        sp.probs .= expv(dt_step, A, sp.probs; m=krylov_m)
        t += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
        expand!(sp, full_sys, bc; depth=1)
    end
    sp
end

divider(s) = println("\n", "═"^60, "\n  $s\n", "═"^60)
section(s)  = println("\n── $s ─────────────────────────────────")

# ════════════════════════════════════════════════════════════════════════════
# EXP 1: Restrict / Prolong roundtrip
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 1: Restrict / Prolong Roundtrip (K=4)")

let K = 4, K2 = 2
    # Construct a known distribution over K=4 voxels.
    # Use a product-Poisson(λ) for clarity.
    λ = 2.0
    max_n = 8
    states_probs = Pair{NTuple{K,Int},Float64}[]
    for n1 in 0:max_n, n2 in 0:max_n, n3 in 0:max_n, n4 in 0:max_n
        p = prod(pdf(Poisson(λ), n) for n in (n1,n2,n3,n4))
        p > 1e-10 || continue
        push!(states_probs, (n1,n2,n3,n4) => p)
    end
    sp_fine = make_sp(Val(K), states_probs)
    renormalize!(sp_fine)

    println("  Fine states: $(length(sp_fine)),  mass = $(round(sum(sp_fine.probs), digits=6))")
    println("  Fine means:  ", round.(voxel_means(sp_fine), digits=4))

    # Restrict
    sp_c = restrict(sp_fine)
    println("  Coarse states: $(length(sp_c)),  mass = $(round(sum(sp_c.probs), digits=6))")
    mu_c = voxel_means(sp_c)
    println("  Coarse means:  ", round.(mu_c, digits=4),
            "  (expected ≈ [$(2λ), $(2λ)])")

    # Prolong back
    sp_fine2 = prolong(sp_c, Val(K))
    println("  Prolonged fine states: $(length(sp_fine2)),  mass = $(round(sum(sp_fine2.probs), digits=6))")
    println("  Prolonged fine means:  ", round.(voxel_means(sp_fine2), digits=4))

    tv = tv_distance(sp_fine, sp_fine2)
    println("  TV(original, restrict→prolong) = $(round(tv, digits=6))")

    # Check marginals for voxel-pair (1,2): should be Binomial(n1+n2, 1/2)
    # if the original was already a product-Poisson (which is not Binomial-splittable in general)
    println("\n  Checking voxel-pair (1,2) intra-pair distribution after prolong:")
    marginal_check = Dict{Int,Float64}()  # n1 distribution given pair-sum
    for i in eachindex(sp_c.states)
        n_pair = Tuple(sp_c.states[i])[1]  # coarse voxel 1 = n1+n2
        p_c    = sp_c.probs[i]
        for m in 0:n_pair
            w = p_c * exp(-n_pair * log(2.0)) * binomial(n_pair, m)
            marginal_check[m] = get(marginal_check, m, 0.0) + w
        end
    end
    # Compare to true marginal of n1
    marginal_true = Dict{Int,Float64}()
    for i in eachindex(sp_fine.states)
        n1 = Tuple(sp_fine.states[i])[1]
        marginal_true[n1] = get(marginal_true, n1, 0.0) + sp_fine.probs[i]
    end
    tv_marginal = sum(abs(get(marginal_check, n, 0.0) - get(marginal_true, n, 0.0))
                      for n in union(keys(marginal_check), keys(marginal_true))) / 2
    println("    TV(prolong marginal n1, true marginal n1) = $(round(tv_marginal, digits=6))")
    println("    (Non-zero because product-Poisson ≠ product of Binomials in general)")
    println("    This error is small when intra-pair diffusion is fast (pre-smooth equilibrates)")
end

# ════════════════════════════════════════════════════════════════════════════
# EXP 2: Coarse system stationary means
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 2: Coarse System Stationary Means (Birth Rate Bug Check)")

let K = 4, D = 1.0, k_b = 2.0, k_d = 1.0, n_max = 30
    fine_grid   = VoxelGrid(K, 1.0/K, 0)
    coarse_grid = coarsen(fine_grid)
    model       = RDMEModel1D(D, k_b, k_d)
    rates       = Float64[]

    # Fine stationary means (analytical)
    μ_fine = stationary_means(K, D, fine_grid.dx, k_b, k_d)
    println("  Fine stationary means μ* = ", round.(μ_fine, digits=4))
    println("  Total: $(round(sum(μ_fine), digits=4))")

    # Coarse: each coarse voxel = 2 fine voxels, expected mean = sum of 2 fine voxels
    K2 = coarse_grid.n_voxels
    μ_coarse_expected = [μ_fine[2j-1] + μ_fine[2j] for j in 1:K2]
    println("\n  Expected coarse means (sum of pairs): ", round.(μ_coarse_expected, digits=4))
    println("  Total: $(round(sum(μ_coarse_expected), digits=4))")

    # Coarse stationary means from the current coarse system
    μ_coarse_analytical = stationary_means(K2, D, coarse_grid.dx, k_b, k_d)
    println("\n  Coarse system stationary means (current code, k_b=$k_b per coarse voxel):")
    println("  ", round.(μ_coarse_analytical, digits=4))
    println("  Total: $(round(sum(μ_coarse_analytical), digits=4))")
    println("  Error vs expected: ", round.(μ_coarse_analytical .- μ_coarse_expected, digits=4))

    # What if birth rate is scaled by coarsening ratio?
    coarsen_ratio = coarse_grid.dx / fine_grid.dx   # = 2
    μ_coarse_corrected = stationary_means(K2, D, coarse_grid.dx, k_b * coarsen_ratio, k_d)
    println("\n  Coarse system stationary means (corrected: k_b × $(coarsen_ratio) = $(k_b*coarsen_ratio)):")
    println("  ", round.(μ_coarse_corrected, digits=4))
    println("  Total: $(round(sum(μ_coarse_corrected), digits=4))")
    println("  Error vs expected: ", round.(μ_coarse_corrected .- μ_coarse_expected, digits=4))

    println("\n  *** VERDICT ***")
    err_current   = norm(μ_coarse_analytical .- μ_coarse_expected)
    err_corrected = norm(μ_coarse_corrected  .- μ_coarse_expected)
    println("  ‖error‖ current code:   $(@sprintf("%.6f", err_current))")
    println("  ‖error‖ corrected code: $(@sprintf("%.6f", err_corrected))")
    if err_corrected < err_current
        println("  → Birth rate scaling DOES improve coarse system accuracy!")
    end
end

# ════════════════════════════════════════════════════════════════════════════
# EXP 3: Single V-cycle on K=2 (minimal analytical case)
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 3: Single V-cycle — K=2 Minimal Case")

let K = 2, D = 1.0, k_b = 1.0, k_d = 0.5, n_max = 20, dt = 0.1
    fine_grid   = VoxelGrid(K, 1.0/K, 0)    # dx = 0.5, d = D/0.25 = 4
    coarse_grid = coarsen(fine_grid)         # K2=1, dx=1.0, d=D/1=1
    model       = RDMEModel1D(D, k_b, k_d)
    rates       = Float64[]
    u0          = CartesianIndex(2, 0)       # start with 2 molecules in voxel 1

    println("  K=$K, D=$D, k_b=$k_b, k_d=$k_d, dt=$dt")
    println("  Fine diffusion rate d = $(D/fine_grid.dx^2)")
    println("  Coarse diffusion rate d = $(D/coarse_grid.dx^2) (K2=1, no diffusion)")

    full_sys = build_rdme_system(model, fine_grid)
    bc = s -> rdme_bc(s, n_max)

    sp0 = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp0, u0, 1.0)
    expand!(sp0, full_sys, bc; depth=2)

    println("  Initial |S| = $(length(sp0))")

    # Fine FSP reference
    sp_fine_ref = deepcopy(sp0)
    A_fine, = build_generator(sp_fine_ref, full_sys, rates, 0.0)
    sp_fine_ref.probs .= expv(dt, A_fine, sp_fine_ref.probs; m=20)
    renormalize!(sp_fine_ref)
    println("  Fine FSP result means: ", round.(voxel_means(sp_fine_ref), digits=4))

    # V-cycle (no pre/post smoothing)
    section("V-cycle: τ_pre=0 (no smoothing)")
    sp_mg0 = deepcopy(sp0)
    sp_mg0 = two_level_vcycle(sp_mg0, model, fine_grid, coarse_grid, rates, 0.0, dt;
                               τ_pre=0.0, τ_post=0.0, krylov_m=20)
    renormalize!(sp_mg0)
    println("  MG result means (no smooth): ", round.(voxel_means(sp_mg0), digits=4))
    println("  TV(fine, MG no-smooth) = $(round(tv_distance(sp_fine_ref, sp_mg0), digits=6))")

    # V-cycle with moderate pre/post smoothing
    section("V-cycle: τ_pre=0.05, τ_post=0.05")
    sp_mg1 = deepcopy(sp0)
    sp_mg1 = two_level_vcycle(sp_mg1, model, fine_grid, coarse_grid, rates, 0.0, dt;
                               τ_pre=0.05, τ_post=0.05, krylov_m=20)
    renormalize!(sp_mg1)
    println("  MG result means (smooth): ", round.(voxel_means(sp_mg1), digits=4))
    println("  TV(fine, MG smooth) = $(round(tv_distance(sp_fine_ref, sp_mg1), digits=6))")

    # V-cycle with long pre-smoothing (should equilibrate pairs well)
    section("V-cycle: τ_pre=1.0 (long pre-smooth)")
    sp_mg2 = deepcopy(sp0)
    sp_mg2 = two_level_vcycle(sp_mg2, model, fine_grid, coarse_grid, rates, 0.0, dt;
                               τ_pre=1.0, τ_post=1.0, krylov_m=20)
    renormalize!(sp_mg2)
    println("  MG result means (long smooth): ", round.(voxel_means(sp_mg2), digits=4))
    println("  TV(fine, MG long-smooth) = $(round(tv_distance(sp_fine_ref, sp_mg2), digits=6))")
end

# ════════════════════════════════════════════════════════════════════════════
# EXP 4: Long-time convergence K=4 — multigrid vs fine FSP
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 4: Long-time Convergence, K=4")

let K = 4, D = 1.0, k_b = 2.0, k_d = 1.0, tf = 5.0, n_max = 20
    fine_grid   = VoxelGrid(K, 1.0/K, 0)
    coarse_grid = coarsen(fine_grid)
    model       = RDMEModel1D(D, k_b, k_d)
    rates       = Float64[]
    u0          = CartesianIndex(ntuple(k -> k == 1 ? 4 : 0, K))

    μ_star = stationary_means(K, D, fine_grid.dx, k_b, k_d)
    println("  Analytical stationary means: ", round.(μ_star, digits=4))

    # Fine-FSP reference
    section("Fine FSP (reference)")
    t_fine = @elapsed sp_ref = fine_fsp_solve(model, fine_grid, u0, tf;
                                               dt=0.05, n_max=n_max)
    println("  Wall: $(round(t_fine, digits=2))s,  |S|=$(length(sp_ref))")
    println("  Means: ", round.(voxel_means(sp_ref), digits=4))
    tv_ref = analytical_tv(sp_ref, μ_star)
    println("  TV(fine FSP, analytical) = $(round(tv_ref, digits=6))")

    # Multigrid FSP — current code (k_b unchanged in coarse)
    for (dt, τ_pre, τ_post) in [(0.05, 0.0, 0.0), (0.05, 0.01, 0.01), (0.05, 0.1, 0.1)]
        section("Multigrid FSP: dt=$dt, τ_pre=$τ_pre, τ_post=$τ_post")
        alg = RDMEMultigridFSP(dt=dt, τ_pre=τ_pre, τ_post=τ_post,
                                n_max=n_max, krylov_m=20, save_every=10)
        t_mg = @elapsed sol_mg = solve_rdme_multigrid(model, fine_grid, u0,
                                                       (0.0, tf), rates, alg)
        # Reconstruct final StateSpace from snapshot
        sp_mg = StateSpace{CartesianIndex{K},Float64}()
        for (s,p) in zip(sol_mg.snapshots[end]...)
            add_state!(sp_mg, s, p)
        end
        println("  Wall: $(round(t_mg, digits=2))s,  |S|=$(sol_mg.state_space_sizes[end])")
        println("  Means: ", round.(voxel_means(sp_mg), digits=4))
        tv_mg = analytical_tv(sp_mg, μ_star)
        tv_vs_ref = tv_distance(sp_mg, sp_ref)
        println("  TV(MG, analytical)  = $(round(tv_mg, digits=6))")
        println("  TV(MG, fine FSP)    = $(round(tv_vs_ref, digits=6))")
    end
end

# ════════════════════════════════════════════════════════════════════════════
# EXP 5: Error vs dt (splitting error analysis)
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 5: Error vs dt — K=4, single time point t=0.5")

let K = 4, D = 1.0, k_b = 2.0, k_d = 1.0, tf = 0.5, n_max = 20
    fine_grid   = VoxelGrid(K, 1.0/K, 0)
    model       = RDMEModel1D(D, k_b, k_d)
    rates       = Float64[]
    u0          = CartesianIndex(ntuple(k -> k == 1 ? 4 : 0, K))

    # Fine reference (tight dt)
    sp_ref = fine_fsp_solve(model, fine_grid, u0, tf; dt=0.01, n_max=n_max)
    println("  Fine FSP reference means: ", round.(voxel_means(sp_ref), digits=4))

    println("\n  dt        TV(MG, fine)   MG means (voxel 1)")
    println("  ─────────────────────────────────────────────")
    for dt in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01]
        alg = RDMEMultigridFSP(dt=dt, τ_pre=0.01, τ_post=0.01,
                                n_max=n_max, krylov_m=20, save_every=999)
        sol = solve_rdme_multigrid(model, fine_grid, u0, (0.0, tf), rates, alg)
        sp_mg = StateSpace{CartesianIndex{K},Float64}()
        for (s,p) in zip(sol.snapshots[end]...)
            add_state!(sp_mg, s, p)
        end
        tv = tv_distance(sp_mg, sp_ref)
        mu1 = voxel_means(sp_mg)[1]
        @printf("  dt=%.3f   TV=%.6f    μ₁=%.4f\n", dt, tv, mu1)
    end
end

# ════════════════════════════════════════════════════════════════════════════
# EXP 6: K=8 multigrid-only (fine FSP is intractable at K=8 due to 8D state space)
# Compare multigrid solutions against the analytical stationary distribution.
# ════════════════════════════════════════════════════════════════════════════
divider("EXP 6: K=8 — Multigrid Only (fine FSP intractable)")

let K = 8, D = 1.0, k_b = 2.0, k_d = 1.0, tf = 5.0, n_max = 20
    fine_grid = VoxelGrid(K, 1.0/K, 0)
    model     = RDMEModel1D(D, k_b, k_d)
    rates     = Float64[]
    u0        = CartesianIndex(ntuple(k -> k == 1 ? 5 : 0, K))

    μ_star = stationary_means(K, D, fine_grid.dx, k_b, k_d)
    println("  Analytical stationary means: ", round.(μ_star, digits=3))
    println("  Lumpability ratio ε = ", lumpability_ratio(model, fine_grid))
    println("  (Fine FSP omitted — K=8 state space is too large for fine-level reference)")

    configs = [
        (dt=0.05, τ_pre=0.0,  τ_post=0.0,  label="no smooth"),
        (dt=0.05, τ_pre=0.01, τ_post=0.01, label="τ=0.01"),
        (dt=0.05, τ_pre=0.05, τ_post=0.05, label="τ=0.05"),
    ]

    println("\n  Config         Wall   |S|     TV(MG,analyt)   MG means[1:4]")
    println("  ──────────────────────────────────────────────────────────────")
    for cfg in configs
        alg = RDMEMultigridFSP(dt=cfg.dt, τ_pre=cfg.τ_pre, τ_post=cfg.τ_post,
                                n_max=n_max, krylov_m=20, save_every=10)
        t_mg = @elapsed sol_mg = solve_rdme_multigrid(model, fine_grid, u0,
                                                       (0.0, tf), rates, alg)
        sp_mg = StateSpace{CartesianIndex{K},Float64}()
        for (s,p) in zip(sol_mg.snapshots[end]...)
            add_state!(sp_mg, s, p)
        end
        tv_an = analytical_tv(sp_mg, μ_star)
        mu    = voxel_means(sp_mg)
        @printf("  %-14s %4.1fs   %5d   %s         %s\n",
                cfg.label, t_mg, sol_mg.state_space_sizes[end],
                @sprintf("%.5f", tv_an),
                string(round.(mu[1:4], digits=3)))
    end
end

println("\n", "═"^60)
println("All experiments complete.")
println("═"^60)
