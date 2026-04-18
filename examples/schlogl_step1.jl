"""
Schlögl Step 1 — K=2 spatial front, symmetry trap and Fix-3.

IC: (n_low, n_high) — one voxel at low fixed point, one at high.
Coarse variable nc = n1 + n2 = 193 cannot distinguish this from (n_high, n_low).

Ground truth: direct K=2 FSP (tractable: ~200^2 = 40k states).

Methods:
  1. Direct FSP                     (ground truth)
  2. V-cycle Binomial-π             (symmetry trap: TV ≈ 0.5 forever)
  3. V-cycle Fix-3 (injection)      (tracks the asymmetric front)

Each method runs the same number of steps with the same dt.
"""

using DiscStochSim
using ExponentialUtilities
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D     = 0.005
K     = 2
h     = 1.0 / K
n_max = 200

model       = SchloglModel1D(D)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model; n_max=n_max)
n_low, n_uns, n_high = fp

println("="^60)
println("Schlögl K=$K  D=$D")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low=$n_low, n_high=$n_high)  nc = $(n_low+n_high)")
println("="^60)

# ─── ground truth: direct K=2 FSP ─────────────────────────────────────────────

sp_fsp = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_fsp, CartesianIndex(n1, n2), 0.0)
end
sp_fsp.probs[sp_fsp.index[CartesianIndex(n_low, n_high)]] = 1.0

fine_sys  = build_schlogl_rdme_system(model, fine_grid)
A_fsp, = build_generator(sp_fsp, fine_sys, rates, 0.0)

println("\n── Pre-computing dynamic-π table ──")
pi_table = compute_dynamic_pi(sp_fsp, A_fsp; n_max=n_max)

# ─── helpers ──────────────────────────────────────────────────────────────────

function tv_vs_fsp(sp_approx, p_fsp)
    tv = 0.0
    for (i, s) in enumerate(sp_fsp.states)
        idx = get(sp_approx.index, s, 0)
        tv += abs(p_fsp[i] - (idx != 0 ? sp_approx.probs[idx] : 0.0))
    end
    for (i, s) in enumerate(sp_approx.states)
        get(sp_fsp.index, s, 0) != 0 && continue
        tv += abs(sp_approx.probs[i])
    end
    tv / 2
end

function voxel_means(sp)
    mu1 = sum(Tuple(sp.states[i])[1] * sp.probs[i] for i in eachindex(sp.probs))
    mu2 = sum(Tuple(sp.states[i])[2] * sp.probs[i] for i in eachindex(sp.probs))
    mu1, mu2
end

# ─── run FSP ground truth ─────────────────────────────────────────────────────

dt      = 0.2
T_SNAP  = [0.5, 1.0, 2.0, 5.0]
n_steps = round(Int, maximum(T_SNAP) / dt)
cnmax   = 2 * n_max

fsp_snaps = Dict{Float64, Vector{Float64}}()
let p = copy(sp_fsp.probs), t_prev = 0.0
    for t in T_SNAP
        p = max.(expv(t - t_prev, A_fsp, p; m=50), 0.0)
        p ./= sum(p)
        fsp_snaps[t] = copy(p)
        t_prev = t
    end
end
println("Direct FSP: $(length(sp_fsp)) states (pre-allocated)")

# ─── run V-cycle variants ─────────────────────────────────────────────────────

function run_vcycle(mode::Symbol, T_end, dt)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    add_state!(sp, CartesianIndex(n_low, n_high), 1.0)
    t_cur   = 0.0
    n_steps = round(Int, T_end / dt)
    snaps   = Dict{Float64, Tuple{Vector{CartesianIndex{2}}, Vector{Float64}}}()

    for _ in 1:n_steps
        sp = if mode == :binom
            two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid,
                                     pi_table, rates, t_cur, dt;
                                     use_dynamic_pi=false, coarse_n_max=cnmax)
        else  # :fix3
            two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                               pi_table, rates, t_cur, dt;
                                               coarse_n_max=cnmax)
        end
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
        t_cur += dt

        t_snap = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        if abs(t_cur - t_snap) < dt/2 && !haskey(snaps, t_snap)
            snaps[t_snap] = (copy(sp.states), copy(sp.probs))
        end
    end
    snaps
end

println("\n── Running V-cycle Binomial-π ──")
t_binom = @elapsed snaps_binom = run_vcycle(:binom, maximum(T_SNAP), dt)

println("── Running V-cycle Fix-3 ──")
t_fix3  = @elapsed snaps_fix3  = run_vcycle(:fix3,  maximum(T_SNAP), dt)

# ─── comparison table ─────────────────────────────────────────────────────────

println()
println("  t     Method           ⟨n1⟩    ⟨n2⟩    TV vs FSP   |S|")
println("  ──────────────────────────────────────────────────────────")

for t in T_SNAP
    p_gt = fsp_snaps[t]
    mu1_gt = sum(Tuple(sp_fsp.states[i])[1] * p_gt[i] for i in eachindex(p_gt))
    mu2_gt = sum(Tuple(sp_fsp.states[i])[2] * p_gt[i] for i in eachindex(p_gt))
    @printf("  %4.1f  %-16s %6.1f  %6.1f  —           %d\n",
            t, "Direct FSP", mu1_gt, mu2_gt, length(sp_fsp))

    # Binomial
    if haskey(snaps_binom, t)
        sts, prs = snaps_binom[t]
        sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
        for i in eachindex(sts); add_state!(sp_tmp, sts[i], prs[i]); end
        mu1, mu2 = voxel_means(sp_tmp)
        tv = tv_vs_fsp(sp_tmp, p_gt)
        @printf("  %4s  %-16s %6.1f  %6.1f  %.4f      %d\n",
                "", "V-cyc Binom-π", mu1, mu2, tv, length(sts))
    end

    # Fix-3
    if haskey(snaps_fix3, t)
        sts, prs = snaps_fix3[t]
        sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
        for i in eachindex(sts); add_state!(sp_tmp, sts[i], prs[i]); end
        mu1, mu2 = voxel_means(sp_tmp)
        tv = tv_vs_fsp(sp_tmp, p_gt)
        @printf("  %4s  %-16s %6.1f  %6.1f  %.4f      %d\n\n",
                "", "V-cyc Fix-3", mu1, mu2, tv, length(sts))
    end
end

@printf("  Timings — Binom: %.2fs   Fix-3: %.2fs\n", t_binom, t_fix3)
println("\nDone.")
