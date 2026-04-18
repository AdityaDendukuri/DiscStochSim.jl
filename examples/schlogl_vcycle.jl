"""
Schlögl V-cycle Benchmark: K=4 with Dynamic-π Prolongation

Demonstrates that the two-level V-cycle with dynamic-π prolongation correctly
tracks the spatial front for K=4 Schlögl RDME.

Model: K=4 voxels, Schlögl bistable chemistry + diffusion  (D=0.005)
IC:    (n_low, n_low, n_high, n_high) — front at the coarse boundary (between
       coarse voxels 1={fine 1,2} and 2={fine 3,4}).

Methods compared:
  1. Coarse-only Binomial-π   (K=2, no V-cycle)
  2. Coarse-only Dynamic-π    (K=2, no V-cycle)
  3. V-cycle Binomial-π       (K=4→K=2, τ_pre=τ_post=0)
  4. V-cycle Dynamic-π        (K=4→K=2, dynamic-π correction)
  5. V-cycle Fix-3            (K=4→K=2, fine-conditional prolongation)

Since K=4 direct FSP (n_max=200: 201^4 ~ 1.6B states) is intractable, we
use the true K=4 sparse adaptive FSP as ground truth (builds state space
adaptively from IC, reaching ~22K states at t=2).

Run: julia --project examples/schlogl_vcycle.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D    = 0.005    # diffusion coefficient; set to 0.00125 for d_fine=0.02 (K=2 MSA matching)
K    = 4
h    = 1.0 / K        # fine voxel size = 0.25
d    = D / h^2        # fine diffusion rate (=0.08 for D=0.005, =0.02 for D=0.00125)

model     = SchloglModel1D(D)
rates     = Float64[]
n_max     = 180       # per-voxel FSP truncation (covers n_high=162 + ~18 margin)
fine_grid = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)          # K=2

fp = schlogl_fixed_points(model)
n_low, n_uns, n_high = fp

println("="^70)
println("Schlögl V-cycle Benchmark  K=$K voxels  D=$D  d=$(round(d,digits=3))")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low, n_low, n_high, n_high) — front at coarse boundary (symmetric pairs)")
println("="^70)

# ─── initial condition ────────────────────────────────────────────────────────

u0 = CartesianIndex(n_low, n_low, n_high, n_high)   # front at coarse boundary
t_solve  = [0.5, 1.0, 2.0]
t_max    = maximum(t_solve)
dt_c     = 0.2       # coarse step
n_steps  = round(Int, t_max / dt_c)

# ─── pair-level ground truth: two independent K=2 FSPs ───────────────────────
# Pair 1: (n_low, n_low),  Pair 2: (n_high, n_high)
# Each pair evolves as an isolated 2-voxel Schlögl system.
# This is a good test when both pairs are SYMMETRIC (nc = 2*n_low or 2*n_high):
# Binom and dynamic-π both give the correct conditional (peaked at fixed point),
# so TV measures inter-pair coupling error, not prolongation error.

println("\n── Pair-level ground truth (2 × K=2 FSP, symmetric IC) ──")

pair_grid = VoxelGrid(2, h, 0)    # same voxel size as fine K=4
pair_sys  = build_schlogl_rdme_system(model, pair_grid)

# Pair 1: starts at (n_low, n_low)
sp_pair1 = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair1, CartesianIndex(n1, n2), 0.0)
end
sp_pair1.probs[sp_pair1.index[CartesianIndex(n_low, n_low)]] = 1.0
A_pair1, = build_generator(sp_pair1, pair_sys, rates, 0.0)

# Pair 2: starts at (n_high, n_high)
sp_pair2 = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair2, CartesianIndex(n1, n2), 0.0)
end
sp_pair2.probs[sp_pair2.index[CartesianIndex(n_high, n_high)]] = 1.0
A_pair2, = build_generator(sp_pair2, pair_sys, rates, 0.0)

t0 = time()
probs_gt = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
let p1 = copy(sp_pair1.probs), p2 = copy(sp_pair2.probs), t_prev = 0.0
    for t in t_solve
        p1 = max.(expv(t - t_prev, A_pair1, p1; m=50), 0.0); p1 ./= sum(p1)
        p2 = max.(expv(t - t_prev, A_pair2, p2; m=50), 0.0); p2 ./= sum(p2)
        probs_gt[t] = (copy(p1), copy(p2))
        t_prev = t
    end
end
t_gt = time() - t0
@printf("  Wall time: %.2fs  (2 × %d states)\n", t_gt, (n_max+1)^2)

# Extract voxel means from pair FSPs
function pair_means(sp, probs)
    mu1 = sum(Tuple(sp.states[i])[1] * probs[i] for i in eachindex(probs))
    mu2 = sum(Tuple(sp.states[i])[2] * probs[i] for i in eachindex(probs))
    mu1, mu2
end

println("  t      ⟨n1⟩  ⟨n2⟩  (pair1)    ⟨n3⟩  ⟨n4⟩  (pair2)")
println("  ─────────────────────────────────────────────────────────")
for t in t_solve
    p1, p2 = probs_gt[t]
    m1a, m1b = pair_means(sp_pair1, p1)
    m2a, m2b = pair_means(sp_pair2, p2)
    @printf("  %4.1f   %5.1f %5.1f              %5.1f %5.1f\n",
            t, m1a, m1b, m2a, m2b)
end

# ─── compute dynamic-π table (shared by all methods) ─────────────────────────

println("\n── Computing dynamic-π table (K=2 pair stationary) ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_pair1, A_pair1; n_max=n_max)
t_pi = time() - t0
@printf("  Stationary solve + slice extraction: %.2fs\n", t_pi)

# Spot-check: dynamic-π peaks at nc=2*n_low and nc=2*n_high (symmetric pairs)
nc_lo = 2 * n_low; nc_hi = 2 * n_high
pi_lo = pi_table[nc_lo + 1]; pi_hi = pi_table[nc_hi + 1]
@printf("  Dynamic-π peak at nc=%d: n1=%d (π=%.4f)  [Binom peak: n1=%d]\n",
        nc_lo, argmax(pi_lo)-1, maximum(pi_lo), nc_lo÷2)
@printf("  Dynamic-π peak at nc=%d: n1=%d (π=%.4f)  [Binom peak: n1=%d]\n",
        nc_hi, argmax(pi_hi)-1, maximum(pi_hi), nc_hi÷2)

# ─── helper: compute voxel means from fine K=4 state space ───────────────────

function voxel_means_k4(states, probs)
    mu = zeros(4)
    for (i, s) in enumerate(states)
        t = Tuple(s)
        for k in 1:4; mu[k] += probs[i] * t[k]; end
    end
    mu
end

# TV distance between a K=4 fine distribution and the pair-level ground truth
# p_gt = (probs_pair1, probs_pair2) evaluated at time t
function tv_vs_pairs(fine_states, fine_probs, sp_p1, pp1, sp_p2, pp2, n_max)
    # Build combined K=4 ground truth by independence assumption:
    # p_gt(n1,n2,n3,n4) = p1(n1,n2) * p2(n3,n4)
    tv = 0.0
    for (i, s) in enumerate(fine_states)
        t4 = Tuple(s)
        n1, n2, n3, n4 = t4
        i1 = get(sp_p1.index, CartesianIndex(n1, n2), 0)
        i2 = get(sp_p2.index, CartesianIndex(n3, n4), 0)
        p_gt = (i1 != 0 ? pp1[i1] : 0.0) * (i2 != 0 ? pp2[i2] : 0.0)
        tv += abs(fine_probs[i] - p_gt)
    end
    tv / 2
end

# ─── true K=4 sparse adaptive FSP (ground truth) ─────────────────────────────

println("\n── True K=4 sparse adaptive FSP (ground truth) ──")
println("  Adaptive: expand from IC, prune small states each step.")

fine_sys_k4 = build_schlogl_rdme_system(model, fine_grid)
fine_bc_k4  = state -> rdme_bc(state, n_max)

sp_gt = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_gt, u0, 1.0)

probs_gt_k4 = Dict{Float64, Tuple{Vector{CartesianIndex{4}}, Vector{Float64}}}()
t0_gt = time()
let sp = sp_gt, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        expand!(sp, fine_sys_k4, fine_bc_k4; depth=1)
        A_gt, = build_generator(sp, fine_sys_k4, rates, t_cur)
        sp.probs .= expv(dt_step, A_gt, sp.probs; m=50)
        sp.probs  = max.(sp.probs, 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-14)

        if step % 2 == 0
            @printf("  step %2d/%d  |S|=%d\r", step, n_steps, length(sp))
            flush(stdout)
        end

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(probs_gt_k4, t_snap)
            probs_gt_k4[t_snap] = (copy(sp.states), copy(sp.probs))
        end
    end
    global sp_gt = sp
end
t_gt_k4 = time() - t0_gt
@printf("\n  Done in %.1fs  final |S|=%d\n", t_gt_k4, length(sp_gt))

# Voxel means from true K=4 FSP
function voxel_means_from_snap(states, probs)
    mu = zeros(4)
    for (i, s) in enumerate(states)
        t = Tuple(s)
        for k in 1:4; mu[k] += probs[i] * t[k]; end
    end
    mu
end

println("  t     ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩")
println("  ────────────────────────────────")
for t in t_solve
    if haskey(probs_gt_k4, t)
        st, pr = probs_gt_k4[t]
        mu = voxel_means_from_snap(st, pr)
        @printf("  %4.1f  %5.1f %5.1f %5.1f %5.1f\n", t, mu[1], mu[2], mu[3], mu[4])
    end
end

# ─── run all five methods ─────────────────────────────────────────────────────

# Common setup: coarse initial condition
function init_coarse()
    sp = StateSpace{CartesianIndex{2}, Float64}()
    nc1 = n_low + n_low    # = 62, symmetric pair 1
    nc2 = n_high + n_high  # = 324, symmetric pair 2
    add_state!(sp, CartesianIndex(nc1, nc2), 1.0)
    sp
end

coarse_bc = state -> rdme_bc(state, 2*n_max)

# ── Method 1: Coarse-only Binomial-π ─────────────────────────────────────────
coarse_sys_binom = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
sp_c_binom = init_coarse()
expand!(sp_c_binom, coarse_sys_binom, coarse_bc; depth=1)
result_binom = Dict{Float64, Vector{Float64}}()

t0_b = time()
let t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        expand!(sp_c_binom, coarse_sys_binom, coarse_bc; depth=1)
        A_cb, = build_generator(sp_c_binom, coarse_sys_binom, rates, t_cur)
        sp_c_binom.probs .= expv(dt_step, A_cb, sp_c_binom.probs; m=30)
        t_cur += dt_step
        renormalize!(sp_c_binom)
        prune_threshold!(sp_c_binom, 1e-12)

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            if !haskey(result_binom, t_snap)
                # Prolong to fine with Binomial and record means
                sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
                for (i, s) in enumerate(sp_c_binom.states)
                    add_state!(sp_tmp, s, sp_c_binom.probs[i])
                end
                sp_prol = prolong(sp_tmp, Val(4); binom_tol=1e-6)
                renormalize!(sp_prol)
                result_binom[t_snap] = voxel_means_k4(sp_prol.states, sp_prol.probs)
            end
        end
    end
end
t_binom = time() - t0_b

# ── Method 2: Coarse-only Dynamic-π ──────────────────────────────────────────
coarse_sys_dyn = build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
sp_c_dyn = init_coarse()
expand!(sp_c_dyn, coarse_sys_dyn, coarse_bc; depth=1)
result_dyn = Dict{Float64, Vector{Float64}}()

t0_d = time()
let t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        expand!(sp_c_dyn, coarse_sys_dyn, coarse_bc; depth=1)
        A_cd, = build_generator(sp_c_dyn, coarse_sys_dyn, rates, t_cur)
        sp_c_dyn.probs .= expv(dt_step, A_cd, sp_c_dyn.probs; m=30)
        t_cur += dt_step
        renormalize!(sp_c_dyn)
        prune_threshold!(sp_c_dyn, 1e-12)

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            if !haskey(result_dyn, t_snap)
                sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
                for (i, s) in enumerate(sp_c_dyn.states)
                    add_state!(sp_tmp, s, sp_c_dyn.probs[i])
                end
                sp_prol = prolong_dynamic(sp_tmp, pi_table, n_max)
                renormalize!(sp_prol)
                result_dyn[t_snap] = voxel_means_k4(sp_prol.states, sp_prol.probs)
            end
        end
    end
end
t_dyn = time() - t0_d

# ── Method 3: V-cycle Binomial-π ─────────────────────────────────────────────
# Uses two_level_vcycle_schlogl with use_dynamic_pi=false (Binomial correction).
# No explicit fine-level expand! — new fine states come only from the coarse correction.
τ_pre_vc = 1.0 / (d * n_low + 1.0)   # ≈ 1/(d·n̄) — one intra-pair relaxation

sp_vc_binom = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_vc_binom, u0, 1.0)
result_vc_binom = Dict{Float64, Vector{Float64}}()

t0_vb = time()
let sp = sp_vc_binom, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl(sp, model,
                                       fine_grid, coarse_grid, pi_table,
                                       rates, t_cur, dt_step;
                                       τ_pre=τ_pre_vc, krylov_m=30,
                                       use_dynamic_pi=false, weight_tol=1e-8,
                                       coarse_n_max=2*n_max,
                                       expand_coarse=true, coarse_expand_depth=1)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        if step % 10 == 0
            @printf("  [V-cyc Binom] step %d/%d  fine states: %d\n",
                    step, n_steps, length(sp.states))
            flush(stdout)
        end

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            !haskey(result_vc_binom, t_snap) &&
                (result_vc_binom[t_snap] = voxel_means_k4(sp.states, sp.probs))
        end
    end
    global sp_vc_binom = sp
end
t_vc_binom = time() - t0_vb
@printf("  V-cycle Binom: %.1fs  final fine states: %d\n", t_vc_binom, length(sp_vc_binom.states))

# ── Method 4: V-cycle Dynamic-π ──────────────────────────────────────────────
sp_vc_dyn = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_vc_dyn, u0, 1.0)
result_vc_dyn = Dict{Float64, Vector{Float64}}()

t0_vd = time()
let sp = sp_vc_dyn, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl(sp, model,
                                       fine_grid, coarse_grid, pi_table,
                                       rates, t_cur, dt_step;
                                       τ_pre=τ_pre_vc, krylov_m=30,
                                       use_dynamic_pi=true, weight_tol=1e-8,
                                       coarse_n_max=2*n_max,
                                       expand_coarse=true, coarse_expand_depth=1)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        if step % 10 == 0
            @printf("  [V-cyc Dyn-π] step %d/%d  fine states: %d\n",
                    step, n_steps, length(sp.states))
            flush(stdout)
        end

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            !haskey(result_vc_dyn, t_snap) &&
                (result_vc_dyn[t_snap] = voxel_means_k4(sp.states, sp.probs))
        end
    end
    global sp_vc_dyn = sp
end
t_vc_dyn = time() - t0_vd
@printf("  V-cycle Dyn-π: %.1fs  final fine states: %d\n", t_vc_dyn, length(sp_vc_dyn.states))

# ── Method 5: V-cycle Fix-3 (fine-conditional prolongation) ──────────────────
sp_vc_fix3 = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_vc_fix3, u0, 1.0)
result_vc_fix3 = Dict{Float64, Vector{Float64}}()

t0_vf = time()
let sp = sp_vc_fix3, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl_injection(sp, model,
                                                fine_grid, coarse_grid, pi_table,
                                                rates, t_cur, dt_step;
                                                τ_pre=τ_pre_vc, krylov_m=30,
                                                use_dynamic_pi=true, weight_tol=1e-8,
                                                coarse_n_max=2*n_max,
                                                expand_coarse=true, coarse_expand_depth=1)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        if step % 10 == 0
            @printf("  [V-cyc Fix-3] step %d/%d  fine states: %d\n",
                    step, n_steps, length(sp.states))
            flush(stdout)
        end

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            !haskey(result_vc_fix3, t_snap) &&
                (result_vc_fix3[t_snap] = voxel_means_k4(sp.states, sp.probs))
        end
    end
    global sp_vc_fix3 = sp
end
t_vc_fix3 = time() - t0_vf
@printf("  V-cycle Fix-3: %.1fs  final fine states: %d\n", t_vc_fix3, length(sp_vc_fix3.states))

# ─── Print comparison table ───────────────────────────────────────────────────

println("\n── Method comparison: voxel means ──")
println("  t     Method               ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩    TV vs K4-GT   |S|")
println("  ─────────────────────────────────────────────────────────────────────────")

for t in t_solve
    !haskey(probs_gt_k4, t) && continue
    st_gt, pr_gt = probs_gt_k4[t]
    mu_gt = voxel_means_from_snap(st_gt, pr_gt)

    # Build lookup for K=4 GT TV
    gt_index = Dict{CartesianIndex{4}, Float64}()
    for (i, s) in enumerate(st_gt); gt_index[s] = pr_gt[i]; end

    function tv_vs_k4gt(states, probs)
        overlap = 0.0
        for (i, s) in enumerate(states)
            p_gt = get(gt_index, s, 0.0)
            overlap += min(probs[i], p_gt)
        end
        1.0 - overlap
    end

    @printf("  %4.1f  %-20s %5.1f %5.1f %5.1f %5.1f   —            %d\n",
            t, "GT (K=4 sparse FSP)", mu_gt[1], mu_gt[2], mu_gt[3], mu_gt[4], length(st_gt))

    # Coarse Binomial
    sp_tmp_b = StateSpace{CartesianIndex{2}, Float64}()
    for (i, s) in enumerate(sp_c_binom.states); add_state!(sp_tmp_b, s, sp_c_binom.probs[i]); end
    sp_prol_b = prolong(sp_tmp_b, Val(4); binom_tol=1e-6); renormalize!(sp_prol_b)
    mu_cb = get(result_binom, t, fill(NaN, 4))
    tv_cb = tv_vs_k4gt(sp_prol_b.states, sp_prol_b.probs)
    @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f       %d\n",
            "", "Coarse Binom-π", mu_cb[1], mu_cb[2], mu_cb[3], mu_cb[4], tv_cb, length(sp_c_binom))

    # Coarse Dynamic-π
    sp_tmp_d = StateSpace{CartesianIndex{2}, Float64}()
    for (i, s) in enumerate(sp_c_dyn.states); add_state!(sp_tmp_d, s, sp_c_dyn.probs[i]); end
    sp_prol_d = prolong_dynamic(sp_tmp_d, pi_table, n_max); renormalize!(sp_prol_d)
    mu_cd = get(result_dyn, t, fill(NaN, 4))
    tv_cd = tv_vs_k4gt(sp_prol_d.states, sp_prol_d.probs)
    @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f       %d\n",
            "", "Coarse Dyn-π", mu_cd[1], mu_cd[2], mu_cd[3], mu_cd[4], tv_cd, length(sp_c_dyn))

    # V-cycle Binomial
    mu_vb = get(result_vc_binom, t, fill(NaN, 4))
    tv_vb = tv_vs_k4gt(sp_vc_binom.states, sp_vc_binom.probs)
    @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f       %d\n",
            "", "V-cyc Binom-π", mu_vb[1], mu_vb[2], mu_vb[3], mu_vb[4], tv_vb, length(sp_vc_binom))

    # V-cycle Dynamic-π
    mu_vd = get(result_vc_dyn, t, fill(NaN, 4))
    tv_vd = tv_vs_k4gt(sp_vc_dyn.states, sp_vc_dyn.probs)
    @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f       %d\n",
            "", "V-cyc Dyn-π", mu_vd[1], mu_vd[2], mu_vd[3], mu_vd[4], tv_vd, length(sp_vc_dyn))

    # V-cycle Fix-3
    mu_vf = get(result_vc_fix3, t, fill(NaN, 4))
    tv_vf = tv_vs_k4gt(sp_vc_fix3.states, sp_vc_fix3.probs)
    @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f       %d\n\n",
            "", "V-cyc Fix-3", mu_vf[1], mu_vf[2], mu_vf[3], mu_vf[4], tv_vf, length(sp_vc_fix3))
end

# ─── Summary ─────────────────────────────────────────────────────────────────

println()
println("── Summary ──────────────────────────────────────────────────────────────")
@printf("  True K=4 sparse FSP (GT):   %.2fs  final |S|=%d\n", t_gt_k4, length(sp_gt))
@printf("  Pair-level GT (2×K=2 FSP):  %.2fs  (approximate)\n", t_gt)
@printf("  Dynamic-π compute:          %.2fs\n", t_pi)
@printf("  Coarse-only Binomial-π:     %.2fs\n", t_binom)
@printf("  Coarse-only Dynamic-π:      %.2fs\n", t_dyn)
@printf("  V-cycle Binomial-π:         %.2fs  (τ_pre=%.3f)  final |S|=%d\n",
        t_vc_binom, τ_pre_vc, length(sp_vc_binom))
@printf("  V-cycle Dynamic-π:          %.2fs  (τ_pre=%.3f)  final |S|=%d\n",
        t_vc_dyn, τ_pre_vc, length(sp_vc_dyn))
@printf("  V-cycle Fix-3:              %.2fs  (τ_pre=%.3f)  final |S|=%d\n",
        t_vc_fix3, τ_pre_vc, length(sp_vc_fix3))
println()
println("="^70)
println("DONE")
println("="^70)
