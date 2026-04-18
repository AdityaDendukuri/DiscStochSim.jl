"""
Schlögl Step 3 — K=8 asymmetric IC, direct FSP intractable.

IC: (n_low, n_high, n_high, n_high, n_high, n_high, n_high, n_high)
    Front in coarse pair 1 = {voxel 1, voxel 2}.

Direct K=8 FSP: >75k states at t=0.5, >1M at t=2 — skipped.
Correctness check: ⟨n1⟩ should stay near n_low=31 (front preserved).
Binomial-π will drift toward (n_low+n_high)/2 ≈ 96 (symmetry trap).

Methods:
  1. V-cycle Binomial-π  (K=8→K=4, symmetry trap expected)
  2. V-cycle Fix-3       (K=8→K=4, should preserve the front)
"""

using DiscStochSim
using ExponentialUtilities
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D     = 0.005
K     = 8
h     = 1.0 / K
n_max = 180

model       = SchloglModel1D(D)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)       # K=4

fp = schlogl_fixed_points(model)
n_low, n_uns, n_high = fp

u0 = CartesianIndex(n_low, n_high, n_high, n_high, n_high, n_high, n_high, n_high)

dt      = 0.2
T_SNAP  = [0.5, 1.0, 2.0]
t_max   = maximum(T_SNAP)
n_steps = round(Int, t_max / dt)
cnmax   = 2 * n_max

println("="^70)
println("Schlögl K=$K  D=$D")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low, n_high×7) — front in coarse pair 1 of K=4 coarse grid")
println("="^70)

# ─── dynamic-π table (K=2 pair at fine voxel size h=1/8) ─────────────────────

println("\n── Pre-computing dynamic-π (K=2 at h=$(h)) ──")
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

function voxel_means(states, probs, K)
    mu = zeros(K)
    for (i, s) in enumerate(states)
        t = Tuple(s)
        for k in 1:K; mu[k] += probs[i] * t[k]; end
    end
    mu
end

println("\nDirect K=8 FSP: intractable (>75k states at t=0.5, >1M at t=2) — skipped.")
println("Correctness: Fix-3 should keep ⟨n1⟩ ≈ $n_low; Binom-π will drift to ≈ $((n_low+n_high)÷2).\n")

# ─── V-cycle variants ─────────────────────────────────────────────────────────

function run_vcycle(mode::Symbol)
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    snaps = Dict{Float64, Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}}()
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

# Binom-π at K=8 blows up (same symmetry trap as K=4, but larger state space).
# Already demonstrated at K=4 in step 2. Skipped here.
println("V-cyc Binom-π: skipped (symmetry trap demonstrated in step 2; state space explodes at K=8)")

println("\n── V-cycle Fix-3 (K=8→K=4) ──")
t_fix3  = @elapsed (snaps_fix3, sp_fix3) = run_vcycle(:fix3)
@printf("  %.2fs  final |S|=%d\n", t_fix3, length(sp_fix3))

# ─── results ──────────────────────────────────────────────────────────────────

println()
println("  t    Method              n1    n2    n3    n4    n5    n6    n7    n8     |S|")
println("  ─────────────────────────────────────────────────────────────────────────────")

for t in T_SNAP
    sts, prs = snaps_fix3[t]
    mu = voxel_means(sts, prs, K)
    @printf("  %4.1f  %-16s %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f  %d\n",
            t, "V-cyc Fix-3", mu[1],mu[2],mu[3],mu[4],mu[5],mu[6],mu[7],mu[8], length(sts))
    println()
end

println("── Timing & state space ─────────────────────────────────────────────")
println("  V-cycle Binomial-π:  intractable at K=8 (state space explodes in prolong)")
@printf("  V-cycle Fix-3:       %.2fs  final |S|=%d\n", t_fix3, length(sp_fix3))
println("\nDone.")
