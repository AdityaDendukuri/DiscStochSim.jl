"""
Schlögl Step 2 — K=4 asymmetric IC, front within a coarse pair.

IC: (n_low, n_high, n_high, n_high)
Coarse pairs: {1,2} → nc_1 = n_low+n_high = 193  (asymmetric — the front)
              {3,4} → nc_2 = n_high+n_high = 324  (symmetric — homogeneous)

Ground truth: K=4 sparse adaptive FSP.

Methods compared:
  1. K=4 sparse adaptive FSP   (ground truth)
  2. V-cycle Binomial-π        (K=4→K=2, symmetry trap on pair 1)
  3. V-cycle Fix-3             (K=4→K=2, preserves asymmetry)
"""

using DiscStochSim
using ExponentialUtilities
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D     = 0.005
K     = 4
h     = 1.0 / K
n_max = 180

model       = SchloglModel1D(D)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)       # K=2

fp = schlogl_fixed_points(model)
n_low, n_uns, n_high = fp

u0 = CartesianIndex(n_low, n_high, n_high, n_high)   # front in pair 1

dt     = 0.2
T_SNAP = [0.5, 1.0, 2.0]
t_max  = maximum(T_SNAP)
n_steps = round(Int, t_max / dt)
cnmax   = 2 * n_max

println("="^60)
println("Schlögl K=$K  D=$D")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low, n_high, n_high, n_high) — front in coarse pair 1")
@printf("nc_1=%d (asymmetric)  nc_2=%d (symmetric)\n", n_low+n_high, 2*n_high)
println("="^60)

# ─── pre-compute dynamic-π (from K=2 pair at fine voxel size) ────────────────

println("\n── Pre-computing dynamic-π ──")
pair_grid = VoxelGrid(2, h, 0)
pair_sys  = build_schlogl_rdme_system(model, pair_grid)
sp_pair   = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair, CartesianIndex(n1,n2), 0.0)
end
sp_pair.probs[sp_pair.index[CartesianIndex(n_low, n_low)]] = 1.0
A_pair, = build_generator(sp_pair, pair_sys, rates, 0.0)
pi_table = compute_dynamic_pi(sp_pair, A_pair; n_max=n_max)
println("  Done.")

# ─── helpers ──────────────────────────────────────────────────────────────────

function voxel_means_k4(states, probs)
    mu = zeros(4)
    for (i, s) in enumerate(states)
        t = Tuple(s)
        for k in 1:4; mu[k] += probs[i] * t[k]; end
    end
    mu
end

function tv_vs_gt(sp_approx, gt_index)
    overlap = 0.0
    for (i, s) in enumerate(sp_approx.states)
        overlap += min(sp_approx.probs[i], get(gt_index, s, 0.0))
    end
    1.0 - overlap
end

# ─── ground truth: K=4 sparse adaptive FSP ───────────────────────────────────

println("\n── Ground truth: K=4 sparse adaptive FSP ──")
fine_sys_k4 = build_schlogl_rdme_system(model, fine_grid)
fine_bc_k4  = state -> rdme_bc(state, n_max)

sp_gt = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_gt, u0, 1.0)

gt_snaps = Dict{Float64, Tuple{Vector{CartesianIndex{4}}, Vector{Float64}}}()

t_gt = @elapsed let sp = sp_gt, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt, t_max - t_cur)
        expand!(sp, fine_sys_k4, fine_bc_k4; depth=1)
        A, = build_generator(sp, fine_sys_k4, rates, t_cur)
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=50), 0.0)
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-14)

        t_snap = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        if abs(t_cur - t_snap) < dt/2 && !haskey(gt_snaps, t_snap)
            gt_snaps[t_snap] = (copy(sp.states), copy(sp.probs))
        end
    end
    global sp_gt = sp
end
@printf("  %.2fs  final |S|=%d\n", t_gt, length(sp_gt))

# ─── V-cycle variants ─────────────────────────────────────────────────────────

function run_vcycle(mode::Symbol)
    sp = StateSpace{CartesianIndex{4}, Float64}()
    add_state!(sp, u0, 1.0)
    snaps = Dict{Float64, Tuple{Vector{CartesianIndex{4}}, Vector{Float64}}}()
    t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt, t_max - t_cur)
        sp = if mode == :binom
            two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid,
                                     pi_table, rates, t_cur, dt_step;
                                     use_dynamic_pi=false, coarse_n_max=cnmax)
        else  # :fix3
            two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                               pi_table, rates, t_cur, dt_step;
                                               coarse_n_max=cnmax)
        end
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)

        t_snap = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        if abs(t_cur - t_snap) < dt/2 && !haskey(snaps, t_snap)
            snaps[t_snap] = (copy(sp.states), copy(sp.probs))
        end
    end
    snaps, sp
end

println("\n── V-cycle Binomial-π ──")
t_binom = @elapsed (snaps_binom, sp_binom) = run_vcycle(:binom)
@printf("  %.2fs  final |S|=%d\n", t_binom, length(sp_binom))

println("\n── V-cycle Fix-3 ──")
t_fix3 = @elapsed (snaps_fix3, sp_fix3) = run_vcycle(:fix3)
@printf("  %.2fs  final |S|=%d\n", t_fix3, length(sp_fix3))

# ─── comparison table ─────────────────────────────────────────────────────────

println()
println("  t    Method          ⟨n1⟩  ⟨n2⟩  ⟨n3⟩  ⟨n4⟩   TV vs K4-GT   |S|")
println("  ─────────────────────────────────────────────────────────────────────")

for t in T_SNAP
    sts_gt, prs_gt = gt_snaps[t]
    mu_gt = voxel_means_k4(sts_gt, prs_gt)
    gt_idx = Dict(s => prs_gt[i] for (i,s) in enumerate(sts_gt))

    @printf("  %4.1f  %-15s %5.1f %5.1f %5.1f %5.1f  —             %d\n",
            t, "GT (K=4 FSP)", mu_gt[1], mu_gt[2], mu_gt[3], mu_gt[4], length(sts_gt))

    for (label, snaps) in [("V-cyc Binom-π", snaps_binom), ("V-cyc Fix-3", snaps_fix3)]
        sts, prs = snaps[t]
        mu = voxel_means_k4(sts, prs)
        sp_tmp = StateSpace{CartesianIndex{4}, Float64}()
        for i in eachindex(sts); add_state!(sp_tmp, sts[i], prs[i]); end
        tv = tv_vs_gt(sp_tmp, gt_idx)
        @printf("  %4s  %-15s %5.1f %5.1f %5.1f %5.1f  %.4f        %d\n",
                "", label, mu[1], mu[2], mu[3], mu[4], tv, length(sts))
    end
    println()
end

@printf("  GT: %.2fs   Binom: %.2fs   Fix-3: %.2fs\n", t_gt, t_binom, t_fix3)
println("\nDone.")
