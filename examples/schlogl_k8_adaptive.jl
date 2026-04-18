"""
Schlögl K=8 Adaptive V-cycle Benchmark

IC: front placed WITHIN pair 2 — voxels (n_low, n_low, n_low, n_high, n_high, n_high, n_high, n_high)

  pair_admissibility detects within-pair asymmetry α[j] = |μ1-μ2|/(μ1+μ2).
  With this IC:
    pair 1: (n_low,  n_low)  → α ≈ 0           → coarsenable
    pair 2: (n_low,  n_high) → α ≈ 0.68 > α_hi → kept fine  ← front pair
    pair 3: (n_high, n_high) → α ≈ 0           → coarsenable
    pair 4: (n_high, n_high) → α ≈ 0           → coarsenable

  Expected adaptive mask: [true, false, true, true]  →  KM = 5

Pair-level GT: 4 independent K=2 FSPs (p_lo_lo, p_lo_hi, p_hi_hi, p_hi_hi).

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/schlogl_k8_adaptive.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf
using Statistics: mean

flush_print(args...)   = (print(args...);   flush(stdout))
flush_println(args...) = (println(args...); flush(stdout))

# ─── parameters ───────────────────────────────────────────────────────────────

D   = 0.005
K   = 8
h   = 1.0 / K
model       = SchloglModel1D(D)
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)
rates       = Float64[]

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_max = n_high + 20   # margin above high fixed point

# IC: front within pair 2  →  (n_low, n_low, n_low, n_high, ...)
u0 = CartesianIndex(n_low, n_low, n_low, n_high, n_high, n_high, n_high, n_high)

t_solve = [0.5, 1.0, 2.0]
t_max   = maximum(t_solve)
dt_c    = 0.2
n_steps = round(Int, t_max / dt_c)

flush_println("="^72)
flush_println("Schlögl K=$K Adaptive V-cycle  D=$D  h=$h  d=$(round(D/h^2, digits=3))")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d  n_max=%d\n",
        n_low, n_uns, n_high, n_max); flush(stdout)
flush_println("IC: (n_low×3, n_high×5) — front within pair 2")
flush_println("Expected mask after step 1: [true, false, true, true]  (pair 2 fine)")
flush_println("="^72)

# ─── pair-level ground truth ──────────────────────────────────────────────────

flush_println("\n── Pair-level GT (4 × K=2 FSP) ──")
pair_grid = VoxelGrid(2, h, 0)
pair_sys  = build_schlogl_rdme_system(model, pair_grid)

function make_pair_sp(n1, n2; depth=3)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    add_state!(sp, CartesianIndex(n1, n2), 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max, Tuple(s))
    expand!(sp, pair_sys, bc; depth=depth)
    sp
end

sp_ll = make_pair_sp(n_low,  n_low)    # pair 1
sp_lh = make_pair_sp(n_low,  n_high)   # pair 2 (front)
sp_hh = make_pair_sp(n_high, n_high)   # pairs 3,4

@printf("  |S_ll|=%d  |S_lh|=%d  |S_hh|=%d\n",
        length(sp_ll), length(sp_lh), length(sp_hh)); flush(stdout)

A_ll, = build_generator(sp_ll, pair_sys, rates, 0.0)
A_lh, = build_generator(sp_lh, pair_sys, rates, 0.0)
A_hh, = build_generator(sp_hh, pair_sys, rates, 0.0)

probs_gt = Dict{Float64, NTuple{3, Vector{Float64}}}()
t0 = time()
let p_ll = copy(sp_ll.probs), p_lh = copy(sp_lh.probs),
    p_hh = copy(sp_hh.probs), t_prev = 0.0
    for t in t_solve
        p_ll = max.(expv(t-t_prev, A_ll, p_ll; m=50), 0.0); p_ll ./= sum(p_ll)
        p_lh = max.(expv(t-t_prev, A_lh, p_lh; m=50), 0.0); p_lh ./= sum(p_lh)
        p_hh = max.(expv(t-t_prev, A_hh, p_hh; m=50), 0.0); p_hh ./= sum(p_hh)
        probs_gt[t] = (copy(p_ll), copy(p_lh), copy(p_hh))
        t_prev = t
    end
end
t_gt = time() - t0
@printf("  Done in %.2fs\n", t_gt); flush(stdout)

function pmean(sp, p, dim)
    sum(Tuple(sp.states[i])[dim] * p[i] for i in eachindex(p))
end

function voxel_means_8(sp)
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

# TV relative to independent-pair GT
# layout: pair1=(1,2)→ll, pair2=(3,4)→lh, pair3=(5,6)→hh, pair4=(7,8)→hh
function tv_vs_pairs(sp_k8, t)
    p_ll, p_lh, p_hh = probs_gt[t]
    overlap = 0.0
    for (i, s) in enumerate(sp_k8.states)
        x = Tuple(s)
        f(sp, p, a, b) = begin
            idx = get(sp.index, CartesianIndex(a,b), 0)
            idx > 0 ? p[idx] : 0.0
        end
        p_gt = f(sp_ll, p_ll, x[1], x[2]) *
               f(sp_lh, p_lh, x[3], x[4]) *
               f(sp_hh, p_hh, x[5], x[6]) *
               f(sp_hh, p_hh, x[7], x[8])
        overlap += min(sp_k8.probs[i], p_gt)
    end
    1.0 - overlap
end

# ─── dynamic-π table  ─────────────────────────────────────────────────────────
flush_println("\n── Computing dynamic-π table ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_ll, A_ll; n_max=n_max)
@printf("  Done in %.2fs\n", time() - t0); flush(stdout)

# ─── local helpers ────────────────────────────────────────────────────────────
rates_fine = build_schlogl_rdme_system(model, fine_grid)

function copy_sp(sp::StateSpace{E,T}) where {E,T}
    sp2 = StateSpace{E,T}()
    for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
    sp2
end

function make_fine_sp(u0; depth=1)
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max, Tuple(s))
    expand!(sp, rates_fine, bc; depth=depth)
    sp
end

# ─── Method 1: Coarse-only Dyn-π ─────────────────────────────────────────────
flush_println("\n── Method 1: Coarse-only Dyn-π (K=4) ──")
coarse_sys_dyn = build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
coarse_bc      = s -> rdme_bc(s, 2*n_max)

# Coarse IC: pair sums — pair2 gets n_low+n_high
sp_c = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_c, CartesianIndex(2n_low, n_low+n_high, 2n_high, 2n_high), 1.0)
expand!(sp_c, coarse_sys_dyn, coarse_bc; depth=1)

result_coarse = Dict{Float64, Vector{Float64}}()
snap_coarse   = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()
t0_c = time()
let t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        expand!(sp_c, coarse_sys_dyn, coarse_bc; depth=1)
        A_c, = build_generator(sp_c, coarse_sys_dyn, rates, t_cur)
        sp_c.probs .= expv(dt_step, A_c, sp_c.probs; m=30)
        t_cur += dt_step
        renormalize!(sp_c)
        prune_threshold!(sp_c, 1e-12)

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_coarse, t_snap)
            sp_prol = prolong_dynamic(sp_c, pi_table, n_max)
            renormalize!(sp_prol)
            result_coarse[t_snap] = voxel_means_8(sp_prol)
            snap_coarse[t_snap]   = sp_prol
            @printf("  t=%.1f  |S_c|=%d  |S_prol|=%d\n",
                    t_snap, length(sp_c), length(sp_prol)); flush(stdout)
        end
    end
end
t_coarse = time() - t0_c
@printf("  Done in %.1fs  final |S_c|=%d\n", t_coarse, length(sp_c)); flush(stdout)

# ─── Method 2: Non-adaptive V-cycle Mult. prolong. ───────────────────────────
flush_println("\n── Method 2: V-cycle Mult. prolong. (K=8→K=4, non-adaptive) ──")
sp_vc = make_fine_sp(u0; depth=1)
result_vc = Dict{Float64, Vector{Float64}}()
snap_vc   = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()

t0_v = time()
let sp = sp_vc, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                                pi_table, rates, t_cur, dt_step;
                                                τ_pre=0.0, τ_post=0.0, krylov_m=30,
                                                use_dynamic_pi=true, weight_tol=1e-8,
                                                coarse_n_max=2*n_max,
                                                expand_coarse=true, coarse_expand_depth=1)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_vc, t_snap)
            result_vc[t_snap] = voxel_means_8(sp)
            snap_vc[t_snap]   = copy_sp(sp)
            @printf("  t=%.1f  |S|=%d\n", t_snap, length(sp)); flush(stdout)
        end
    end
    global sp_vc = sp
end
t_vc = time() - t0_v
@printf("  Done in %.1fs  final |S|=%d\n", t_vc, length(sp_vc)); flush(stdout)

# ─── Method 3: Adaptive V-cycle ──────────────────────────────────────────────
flush_println("\n── Method 3: Adaptive V-cycle (K=8, mask per step) ──")
flush_println("  Expected: mask[2]=false (front pair), mask[1,3,4]=true (bulk)")
sp_ad = make_fine_sp(u0; depth=1)
result_ad    = Dict{Float64, Vector{Float64}}()
snap_ad      = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()
mask_history = Dict{Float64, BitVector}()

t0_a = time()
let sp = sp_ad, prev_mask = nothing, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp, prev_mask = two_level_vcycle_adaptive(sp, model, fine_grid, rates, t_cur, dt_step;
                                                   prev_mask,
                                                   α_lo=0.10, α_hi=0.25,
                                                   φ_lo=1.5,  φ_hi=3.0,
                                                   krylov_m=30,
                                                   expand_mixed=true,
                                                   mixed_expand_depth=1,
                                                   mixed_n_max=2*n_max,
                                                   binom_tol=1e-4)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        KM = K - sum(prev_mask)
        @printf("  step %2d/%d  t=%.1f  |S|=%6d  mask=%s  KM=%d\n",
                step, n_steps, t_cur, length(sp.states), string(prev_mask), KM)
        flush(stdout)

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_ad, t_snap)
            result_ad[t_snap]    = voxel_means_8(sp)
            snap_ad[t_snap]      = copy_sp(sp)
            mask_history[t_snap] = copy(prev_mask)
        end
    end
    global sp_ad   = sp
    global final_mask = prev_mask
end
t_ad = time() - t0_a
@printf("  Done in %.1fs  final |S|=%d  mask=%s\n",
        t_ad, length(sp_ad), string(final_mask)); flush(stdout)

# ─── Results table ────────────────────────────────────────────────────────────
flush_println()
flush_println("── Results ──────────────────────────────────────────────────────────────────")
flush_println("  t    Method                      ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩  ⟨n5..8⟩  TV   |S|")
flush_println("  ─────────────────────────────────────────────────────────────────────────────")

for t in t_solve
    p_ll, p_lh, p_hh = probs_gt[t]
    @printf("  %.1f  %-28s  %5.1f %5.1f %5.1f %5.1f   %5.1f   —      —\n",
            t, "Pair-GT (4×K=2)",
            pmean(sp_ll, p_ll, 1), pmean(sp_ll, p_ll, 2),
            pmean(sp_lh, p_lh, 1), pmean(sp_lh, p_lh, 2),
            pmean(sp_hh, p_hh, 1))

    mu_c  = get(result_coarse, t, fill(NaN, K))
    sp_cs = get(snap_coarse,  t, nothing)
    tv_c  = sp_cs !== nothing ? tv_vs_pairs(sp_cs, t) : NaN
    @printf("  %.1f  %-28s  %5.1f %5.1f %5.1f %5.1f   %5.1f  %.4f  %d\n",
            t, "Coarse Dyn-π (K=4)",
            mu_c[1], mu_c[2], mu_c[3], mu_c[4], mean(mu_c[5:8]), tv_c, length(sp_c))

    mu_v  = get(result_vc, t, fill(NaN, K))
    sp_vs = get(snap_vc,   t, nothing)
    tv_v  = sp_vs !== nothing ? tv_vs_pairs(sp_vs, t) : NaN
    @printf("  %.1f  %-28s  %5.1f %5.1f %5.1f %5.1f   %5.1f  %.4f  %d\n",
            t, "V-cyc Mult. prolong.",
            mu_v[1], mu_v[2], mu_v[3], mu_v[4], mean(mu_v[5:8]), tv_v, length(sp_vc))

    mu_a  = get(result_ad, t, fill(NaN, K))
    sp_as = get(snap_ad,   t, nothing)
    tv_a  = sp_as !== nothing ? tv_vs_pairs(sp_as, t) : NaN
    mk    = get(mask_history, t, nothing)
    @printf("  %.1f  %-28s  %5.1f %5.1f %5.1f %5.1f   %5.1f  %.4f  %d  mask=%s\n\n",
            t, "Adaptive V-cycle",
            mu_a[1], mu_a[2], mu_a[3], mu_a[4], mean(mu_a[5:8]),
            tv_a, length(sp_ad), mk !== nothing ? string(mk) : "?")
end; flush(stdout)

flush_println("── Timing ───────────────────────────────────────────────────────────────────")
@printf("  Pair-level GT:             %.2fs\n", t_gt)
@printf("  Coarse Dyn-π (K=4):        %.1fs  |S|=%d\n", t_coarse, length(sp_c))
@printf("  V-cyc Mult. prolong.:      %.1fs  |S|=%d\n", t_vc, length(sp_vc))
@printf("  Adaptive V-cycle:          %.1fs  |S|=%d\n", t_ad, length(sp_ad))
flush(stdout)
flush_println("\n" * "="^72)
flush_println("DONE")
flush_println("="^72)
