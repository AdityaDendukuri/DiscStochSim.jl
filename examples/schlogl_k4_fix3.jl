"""
Schlögl K=4 Fix-3 Benchmark: Asymmetric Within-Pair Front

Tests Fix-3 (fine-conditional prolongation) against V-cycle Dynamic-π for a case
that triggers the symmetry trap: IC = (n_low, n_high, n_low, n_high).

Each pair has nc_j = n_low + n_high = 193. The static dynamic-π assigns 50/50 weight
between (n_low, n_high) and (n_high, n_low) → symmetry trap → ⟨n1⟩ → nc/2 ≈ 96.
Fix-3 carries the fine conditional from the covered state to new states, keeping
the asymmetry alive.

Note: τ_pre=0 (no pre-smooth) to avoid the pre-smoother itself introducing symmetry
at nc=193 (the intra-pair stationary for nc=193 is bimodal at n_low and n_high).

Ground truth: 2 × independent K=2 FSPs, each from IC (n_low, n_high).

Run: JULIA_PKG_PRECOMPILE_AUTO=0 julia --project examples/schlogl_k4_fix3.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D    = 0.005
K    = 4
h    = 1.0 / K
model      = SchloglModel1D(D)
rates      = Float64[]
n_max      = 180
fine_grid  = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model)
n_low, n_uns, n_high = fp

u0 = CartesianIndex(n_low, n_high, n_low, n_high)  # within-pair fronts

t_solve  = [0.5, 1.0, 2.0]
t_max    = maximum(t_solve)
dt_c     = 0.2
n_steps  = round(Int, t_max / dt_c)

println("="^70)
println("Schlögl K=$K Fix-3 Benchmark — Asymmetric Within-Pair IC")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low, n_high, n_low, n_high) — both pairs are spatial fronts")
println("Each coarse pair: nc = n_low+n_high = $(n_low+n_high) (symmetry trap nc)")
println("τ_pre = 0 (pre-smooth disabled to avoid within-pair symmetrization)")
println("="^70)

# ─── pair-level ground truth: 2 × K=2 FSP from (n_low, n_high) ──────────────

println("\n── Pair-level ground truth (2 × K=2 FSP from front IC) ──")

pair_grid = VoxelGrid(2, h, 0)
pair_sys  = build_schlogl_rdme_system(model, pair_grid)

sp_pair = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair, CartesianIndex(n1, n2), 0.0)
end
sp_pair.probs[sp_pair.index[CartesianIndex(n_low, n_high)]] = 1.0

A_pair, = build_generator(sp_pair, pair_sys, rates, 0.0)

t0 = time()
probs_gt = Dict{Float64, Vector{Float64}}()
let p = copy(sp_pair.probs), t_prev = 0.0
    for t in t_solve
        p = max.(expv(t - t_prev, A_pair, p; m=50), 0.0); p ./= sum(p)
        probs_gt[t] = copy(p)
        t_prev = t
    end
end
@printf("  Wall time: %.2fs  (%d states)\n", time()-t0, (n_max+1)^2)

println("  t     ⟨n1⟩  ⟨n2⟩  (each pair independently)")
println("  ──────────────────────────────────────────")
for t in t_solve
    p = probs_gt[t]
    mu1 = sum(Tuple(sp_pair.states[i])[1] * p[i] for i in eachindex(p))
    mu2 = sum(Tuple(sp_pair.states[i])[2] * p[i] for i in eachindex(p))
    @printf("  %4.1f   %5.1f %5.1f\n", t, mu1, mu2)
end

# ─── compute dynamic-π (for dynamic-π V-cycle and Fix-3) ─────────────────────

println("\n── Computing dynamic-π table ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_pair, A_pair; n_max=n_max)
@printf("  Done in %.2fs\n", time()-t0)

# ─── helpers ──────────────────────────────────────────────────────────────────

function voxel_means(sp)
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

# TV vs independent-pair ground truth at time t
function tv_fine(sp, t)
    p_gt_vec = probs_gt[t]
    overlap  = 0.0
    for (i, s) in enumerate(sp.states)
        t4 = Tuple(s)
        i1 = get(sp_pair.index, CartesianIndex(t4[1], t4[2]), 0)
        i2 = get(sp_pair.index, CartesianIndex(t4[3], t4[4]), 0)
        p_gt = (i1 != 0 ? p_gt_vec[i1] : 0.0) * (i2 != 0 ? p_gt_vec[i2] : 0.0)
        overlap += min(sp.probs[i], p_gt)
    end
    1.0 - overlap
end

# ─── Method 1: V-cycle Dynamic-π (existing correction) ───────────────────────

println("\n── Method 1: V-cycle Dynamic-π correction (τ_pre=0) ──")

sp_vd = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_vd, u0, 1.0)
result_vd = Dict{Float64, NamedTuple}()

t0 = time()
let sp = sp_vd, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid, pi_table,
                                      rates, t_cur, dt_step;
                                      τ_pre=0.0, τ_post=0.0, krylov_m=30,
                                      use_dynamic_pi=true, weight_tol=1e-8,
                                      coarse_n_max=2*n_max,
                                      expand_coarse=true, coarse_expand_depth=1)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        if step % 2 == 0
            @printf("  step %2d/%d  |S|=%d\r", step, n_steps, length(sp.states))
            flush(stdout)
        end

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_vd, t_snap)
            mu = voxel_means(sp)
            tv = tv_fine(sp, t_snap)
            result_vd[t_snap] = (mu=mu, tv=tv, n=length(sp))
        end
    end
    global sp_vd = sp
end
t_vd = time() - t0
@printf("\n  Done in %.1fs  final |S|=%d\n", t_vd, length(sp_vd))

# ─── Method 2: V-cycle Fix-3 (fine-conditional prolongation) ─────────────────

println("\n── Method 2: V-cycle Fix-3 fine-conditional prolongation (τ_pre=0) ──")

sp_f3 = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_f3, u0, 1.0)
result_f3 = Dict{Float64, NamedTuple}()

t0 = time()
let sp = sp_f3, t_cur = 0.0
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

        if step % 2 == 0
            @printf("  step %2d/%d  |S|=%d\r", step, n_steps, length(sp.states))
            flush(stdout)
        end

        t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_f3, t_snap)
            mu = voxel_means(sp)
            tv = tv_fine(sp, t_snap)
            result_f3[t_snap] = (mu=mu, tv=tv, n=length(sp))
        end
    end
    global sp_f3 = sp
end
t_f3 = time() - t0
@printf("\n  Done in %.1fs  final |S|=%d\n", t_f3, length(sp_f3))

# ─── Results ─────────────────────────────────────────────────────────────────

println()
println("── Results ──────────────────────────────────────────────────────────────")
println("  Pair means — GT: ⟨n1⟩≈⟨n3⟩ (low voxels stay near n_low=31);")
println("  symmetry trap flips them toward nc/2≈96.")
println()
println("  t     Method               ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩    TV vs GT   |S|")
println("  ───────────────────────────────────────────────────────────────────────")

for t in t_solve
    p_gt = probs_gt[t]
    mu1_gt = sum(Tuple(sp_pair.states[i])[1] * p_gt[i] for i in eachindex(p_gt))
    mu2_gt = sum(Tuple(sp_pair.states[i])[2] * p_gt[i] for i in eachindex(p_gt))
    @printf("  %4.1f  %-20s %5.1f %5.1f %5.1f %5.1f   —\n",
            t, "GT (each pair)", mu1_gt, mu2_gt, mu1_gt, mu2_gt)

    vd = get(result_vd, t, nothing)
    f3 = get(result_f3, t, nothing)

    if !isnothing(vd)
        @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f    %d\n",
                "", "Dyn-π correction",
                vd.mu[1], vd.mu[2], vd.mu[3], vd.mu[4], vd.tv, vd.n)
    end
    if !isnothing(f3)
        @printf("  %4s  %-20s %5.1f %5.1f %5.1f %5.1f   %.4f    %d\n\n",
                "", "Fix-3 fine-cond.",
                f3.mu[1], f3.mu[2], f3.mu[3], f3.mu[4], f3.tv, f3.n)
    end
end

println("── Timing ───────────────────────────────────────────────────────────────")
@printf("  V-cycle Dyn-π:  %.1fs  final |S|=%d\n", t_vd, length(sp_vd))
@printf("  V-cycle Fix-3:  %.1fs  final |S|=%d\n", t_f3, length(sp_f3))
println()
println("  Key: ⟨n1⟩ and ⟨n3⟩ should stay near n_low=31 (low voxels of each front).")
println("  Fix-3 preserves the within-pair asymmetry; Dyn-π correction symmetrizes.")
println("="^70)
println("DONE")
println("="^70)
