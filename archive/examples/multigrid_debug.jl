"""
Multigrid FSP — Detailed Step-by-Step Diagnostic

Instruments every step of the V-cycle: state space size, total probability mass,
per-voxel means, and coarse-level statistics.  Prints a line per V-cycle step
so we can see exactly where the simulation goes wrong.

Run from project root:
  julia --project examples/multigrid_debug.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Distributions
using Printf

# ────────────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────────────

vmeans(sp::StateSpace{CartesianIndex{K},Float64}) where K = begin
    mu = zeros(K)
    for (i,s) in enumerate(sp.states)
        for k in 1:K; mu[k] += sp.probs[i]*Tuple(s)[k]; end
    end
    mu
end

function stationary_means(K, D, dx, k_b, k_d)
    d = D/dx^2; A = zeros(K,K); b = fill(-k_b,K)
    for k in 1:K
        A[k,k] = -(k_d + (k==1||k==K ? d : 2d))
        if k>1; A[k,k-1]=d; end; if k<K; A[k,k+1]=d; end
    end
    A\b
end

# ────────────────────────────────────────────────────────────────────────────
# Instrumented V-cycle
# ────────────────────────────────────────────────────────────────────────────

function vcycle_debug(sp_fine::StateSpace{CartesianIndex{K},Float64},
                      model, fine_grid, coarse_grid, rates, t, dt;
                      τ_pre=0.0, τ_post=0.0, krylov_m=20,
                      weight_tol=1e-14, binom_tol=0.0,
                      expand_coarse=false, coarse_n_max=80, coarse_expand_depth=1,
                      step_label="") where K

    function fmt(sp, label)
        mu = vmeans(sp)
        @printf("    %-22s |S|=%6d  mass=%.6f  μ_tot=%.3f  μ[1]=%.3f  μ[end]=%.3f\n",
                label, length(sp), sum(sp.probs), sum(mu), mu[1], mu[end])
        flush(stdout)
    end

    @printf("  Step %-5s t=%.3f → %.3f\n", step_label, t, t+dt); flush(stdout)
    fmt(sp_fine, "fine_in")

    # 1. Pre-smooth
    sp_s = if τ_pre > 0
        intra = build_intra_system(model, fine_grid, coarse_grid)
        A_i, = build_generator(sp_fine, intra, rates, t)
        sp2 = DiscStochSim._copy_sp(sp_fine)
        sp2.probs .= expv(Float64(τ_pre), A_i, sp_fine.probs; m=krylov_m)
        fmt(sp2, "after_presmooth")
        sp2
    else
        sp_fine
    end

    # 2. Restrict
    sp_c = restrict(sp_s)
    mu_c = vmeans(sp_c)
    @printf("    %-22s |S|=%6d  mass=%.6f  μ_tot=%.3f  μ_c[1]=%.3f\n",
            "after_restrict", length(sp_c), sum(sp_c.probs), sum(mu_c), mu_c[1])
    flush(stdout)

    # 3. Coarse solve (expand BEFORE solve so prob can flow into new states)
    coarse_sys = build_coarse_system(model, coarse_grid, fine_grid)

    # 3a. Coarse expand BEFORE solve (Fix 2)
    if expand_coarse && coarse_expand_depth > 0
        coarse_bc = state -> rdme_bc(state, coarse_n_max)
        expand!(sp_c, coarse_sys, coarse_bc; depth=coarse_expand_depth)
        mu_cx = vmeans(sp_c)
        @printf("    %-22s |S|=%6d  mass=%.6f  μ_tot=%.3f  μ_c[1]=%.3f\n",
                "after_coarse_expand", length(sp_c), sum(sp_c.probs), sum(mu_cx), mu_cx[1])
        flush(stdout)
    end

    A_c, = build_generator(sp_c, coarse_sys, rates, t)
    sp_c.probs .= expv(Float64(dt), A_c, sp_c.probs; m=krylov_m)
    mu_c2 = vmeans(sp_c)
    @printf("    %-22s |S|=%6d  mass=%.6f  μ_tot=%.3f  μ_c[1]=%.3f\n",
            "after_coarse_solve", length(sp_c), sum(sp_c.probs), sum(mu_c2), mu_c2[1])
    flush(stdout)

    # 4. Prolong
    sp_fn = prolong(sp_c, Val(K); weight_tol=weight_tol, binom_tol=binom_tol)
    fmt(sp_fn, "after_prolong")

    # 5. Post-smooth
    if τ_post > 0
        intra = build_intra_system(model, fine_grid, coarse_grid)
        A_i, = build_generator(sp_fn, intra, rates, t)
        sp_fn.probs .= expv(Float64(τ_post), A_i, sp_fn.probs; m=krylov_m)
        fmt(sp_fn, "after_postsmooth")
    end

    sp_fn
end

# ────────────────────────────────────────────────────────────────────────────
# Instrumented solver loop
# ────────────────────────────────────────────────────────────────────────────

function solve_debug(model, fine_grid, u0::CartesianIndex{K}, tspan, rates;
                     dt=0.2, τ_pre=0.0, τ_post=0.0, n_max=20,
                     krylov_m=20, weight_tol=1e-14, binom_tol=0.0,
                     prune_tol=1e-12, expand_depth=0,
                     expand_coarse=true, coarse_expand_depth=1,
                     max_steps=25) where K

    coarse_grid = coarsen(fine_grid)
    full_sys    = build_rdme_system(model, fine_grid)
    bc          = s -> rdme_bc(s, n_max)

    sp = StateSpace{CartesianIndex{K},Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, full_sys, bc; depth=1)

    t0, tf = Float64.(tspan)
    t = t0; step = 0

    println("  Initial state space: |S|=$(length(sp))  mass=$(sum(sp.probs))")
    println("  Initial means: ", round.(vmeans(sp), digits=3))
    println()
    flush(stdout)

    while t < tf && step < max_steps
        dt_s = min(dt, tf - t)
        sp = vcycle_debug(sp, model, fine_grid, coarse_grid, rates, t, dt_s;
                          τ_pre=τ_pre, τ_post=τ_post, krylov_m=krylov_m,
                          weight_tol=weight_tol, binom_tol=binom_tol,
                          expand_coarse=expand_coarse,
                          coarse_n_max=2*n_max, coarse_expand_depth=coarse_expand_depth,
                          step_label="$(step+1)")

        renormalize!(sp)
        prune_threshold!(sp, prune_tol)
        if expand_depth > 0
            expand!(sp, full_sys, bc; depth=expand_depth)
        end

        mu = vmeans(sp)
        @printf("    %-22s |S|=%6d  mass=%.6f  μ_tot=%.3f  μ[1]=%.3f\n",
                "after_prune+expand", length(sp), sum(sp.probs), sum(mu), mu[1])
        println()
        flush(stdout)

        t += dt_s; step += 1
    end

    sp
end

# ════════════════════════════════════════════════════════════════════════════
# RUN 1: K=4, well-behaved case (sanity check everything is correct)
# ════════════════════════════════════════════════════════════════════════════
println("="^65)
println("RUN 1: K=4  (sanity check, should converge to μ*=[2,2,2,2])")
println("="^65)
K=4; D=1.0; k_b=2.0; k_d=1.0
fg=VoxelGrid(K,1.0/K,0); model=RDMEModel1D(D,k_b,k_d)
u0=CartesianIndex(4,0,0,0); rates=Float64[]
mu_star = stationary_means(K,D,fg.dx,k_b,k_d)
println("Stationary means: ", round.(mu_star,digits=3))
println()
flush(stdout)

sp_k4 = solve_debug(model, fg, u0, (0.0,5.0), rates;
                    dt=0.5, τ_pre=0.01, τ_post=0.01, n_max=20,
                    krylov_m=20, weight_tol=1e-14, binom_tol=0.0,
                    prune_tol=1e-12, expand_depth=1, max_steps=10)

println("FINAL K=4: means=", round.(vmeans(sp_k4), digits=4))
println()

# ════════════════════════════════════════════════════════════════════════════
# RUN 2: K=8, expand_depth=0 — show why it fails
# ════════════════════════════════════════════════════════════════════════════
println("="^65)
println("RUN 2: K=8, expand_depth=0  (should FAIL — state space frozen)")
println("="^65)
K=8; fg8=VoxelGrid(K,1.0/K,0); model8=RDMEModel1D(D,k_b,k_d)
u0_8=CartesianIndex(5,0,0,0,0,0,0,0)
mu_star8=stationary_means(K,D,fg8.dx,k_b,k_d)
println("Stationary means: ", round.(mu_star8,digits=3))
println()
flush(stdout)

sp_k8_noexp = solve_debug(model8, fg8, u0_8, (0.0,2.0), rates;
                          dt=0.5, τ_pre=0.01, τ_post=0.01, n_max=12,
                          krylov_m=15, weight_tol=1e-14, binom_tol=0.0,
                          prune_tol=1e-12, expand_depth=0, max_steps=4)

println("FINAL K=8 no-expand: means=", round.(vmeans(sp_k8_noexp),digits=4))
println()

# ════════════════════════════════════════════════════════════════════════════
# RUN 3: K=8, expand_depth=1, binom_tol=1e-3 — track state space growth
# ════════════════════════════════════════════════════════════════════════════
println("="^65)
println("RUN 3: K=8, expand_depth=1, binom_tol=1e-3, prune_tol=1e-8")
println("="^65)
println("(Only 5 steps to see growth pattern before explosion)")
println()
flush(stdout)

sp_k8_grow = solve_debug(model8, fg8, u0_8, (0.0,2.5), rates;
                         dt=0.5, τ_pre=0.01, τ_post=0.01, n_max=12,
                         krylov_m=15, weight_tol=1e-14, binom_tol=1e-3,
                         prune_tol=1e-8, expand_depth=1, expand_coarse=false, max_steps=5)

println("FINAL K=8 with fine-expand: means=", round.(vmeans(sp_k8_grow),digits=4))
println()

# ════════════════════════════════════════════════════════════════════════════
# RUN 4: K=8, FIX 2 — coarse expansion only, fine expand_depth=0
# ════════════════════════════════════════════════════════════════════════════
println("="^65)
println("RUN 4: K=8, FIX 2 — coarse expand only (expand_coarse=true, expand_depth=0)")
println("="^65)
println("Goal: controlled state-space growth, correct spatial spread, t=0→5")
println()
flush(stdout)

# RUN 4a: coarse-only, original params (dt=0.5, krylov_m=20, binom_tol=1e-3)
sp_k8_fix2 = solve_debug(model8, fg8, u0_8, (0.0,5.0), rates;
                          dt=0.5, τ_pre=0.01, τ_post=0.01, n_max=12,
                          krylov_m=20, weight_tol=1e-14, binom_tol=1e-3,
                          prune_tol=1e-10, expand_depth=0,
                          expand_coarse=true, coarse_expand_depth=1,
                          max_steps=10)

println("FINAL K=8 Fix2a: means=", round.(vmeans(sp_k8_fix2),digits=4))
println("Expected (stationary): ", round.(mu_star8, digits=4))
println()

# ════════════════════════════════════════════════════════════════════════════
# RUN 5: K=8, Fix 2b — smaller dt, larger krylov_m, tighter binom_tol
# ════════════════════════════════════════════════════════════════════════════
println("="^65)
println("RUN 5: K=8, Fix 2b — dt=0.1, krylov_m=35, binom_tol=0.05")
println("="^65)
println("Goal: accurate Krylov, controlled fine state space, t=0→2 (20 steps)")
println()
flush(stdout)

sp_k8_fix2b = solve_debug(model8, fg8, u0_8, (0.0,2.0), rates;
                           dt=0.1, τ_pre=0.0, τ_post=0.0, n_max=12,
                           krylov_m=35, weight_tol=1e-14, binom_tol=0.05,
                           prune_tol=1e-8, expand_depth=0,
                           expand_coarse=true, coarse_expand_depth=1,
                           max_steps=20)

println("FINAL K=8 Fix2b: means=", round.(vmeans(sp_k8_fix2b),digits=4))
println("Expected at t=2 (partial convergence): all voxels should be spreading")
println()
println("="^65)
println("DONE")
println("="^65)
