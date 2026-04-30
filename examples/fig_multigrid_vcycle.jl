"""
fig_multigrid_vcycle.jl  --  V-cycle benchmark, research style

Tests the two-level V-cycle (two_level_vcycle_schlogl_injection) against:
  1. K=4 direct FSP  (ground truth)
  2. K=4 coarse-only at K=2  (fast, but TV=1 at front)
  3. K=4 V-cycle (K=4->K=2->K=4 with multiplicative prolongation)

Then scales to K=8:
  4. K=8 V-cycle (K=8->K=4->K=8)  -- is it tractable?

At each step: |S|, wall time, voxel means, TV vs direct FSP.
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf

const n_max  = 180
const dt     = 0.2
const T      = 2.0
const n_steps = round(Int, T/dt)
const T_SNAP  = [0.5, 1.0, 2.0]

D     = 0.005
model = SchloglModel1D(D)
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=250)
n_lo, n_un, n_hi = fp
@printf("Fixed points: n_lo=%d  n_un=%d  n_hi=%d\n\n", n_lo, n_un, n_hi)

# ── helpers ────────────────────────────────────────────────────────────────────
function vmeans(states, probs, K)
    mu = zeros(K)
    for (i,s) in enumerate(states)
        tv = Tuple(s); p = probs[i]
        for k in 1:K; mu[k] += p*tv[k]; end
    end
    mu
end

function tv_dist(sa, pa, sb, pb)
    idxb = Dict(sb[i] => pb[i] for i in eachindex(sb))
    idxa = Dict(sa[i] => true  for i in eachindex(sa))
    tv   = 0.0
    for (i,s) in enumerate(sa); tv += abs(pa[i] - get(idxb, s, 0.0)); end
    for (i,s) in enumerate(sb); !haskey(idxa, s) && (tv += pb[i]); end
    tv / 2
end

function fsp_step!(sp, sys, bc, t, dt_s; prune=1e-10)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t)
    sp.probs .= max.(expv(dt_s, A, sp.probs; m=40), 0.0)
    renormalize!(sp); prune_threshold!(sp, prune)
end

println("="^60)

# ══════════════════════════════════════════════════════════════
# K=4 DIRECT FSP  (ground truth)
# ══════════════════════════════════════════════════════════════
println("\n── K=4 direct FSP (ground truth) ──")
K4 = 4; h4 = 1.0/K4
g4   = VoxelGrid(K4, h4, 0)
sys4 = build_schlogl_rdme_system(model, g4)
bc4  = s -> rdme_bc(s, n_max)
u0_4 = CartesianIndex(n_lo, n_hi, n_hi, n_hi)

sp4  = StateSpace{CartesianIndex{K4}, Float64}()
add_state!(sp4, u0_4, 1.0)
snaps4 = Dict{Float64, Tuple{Vector,Vector}}()
t_cur  = 0.0

t_direct = @elapsed for step in 1:n_steps
    global t_cur
    dt_s = min(dt, T-t_cur)
    fsp_step!(sp4, sys4, bc4, t_cur, dt_s)
    t_cur += dt_s
    snap_t = T_SNAP[argmin(abs.(T_SNAP.-t_cur))]
    abs(t_cur-snap_t)<dt/2 && !haskey(snaps4,snap_t) &&
        (snaps4[snap_t] = (copy(sp4.states), copy(sp4.probs)))
end
@printf("  %.2fs   |S|=%d   means: %s\n", t_direct, length(sp4),
        join([@sprintf("%.1f",x) for x in vmeans(sp4.states,sp4.probs,K4)], " "))

# ══════════════════════════════════════════════════════════════
# K=4 COARSE-ONLY at K=2  (cheap baseline)
# ══════════════════════════════════════════════════════════════
println("\n── K=4 coarse-only (K=2 solve) ──")
g4c   = coarsen(g4)
sys4c = build_schlogl_coarse_system(model, g4c, g4)
bc4c  = s -> rdme_bc(s, 2*n_max)
u0_4c = CartesianIndex(n_lo+n_hi, 2*n_hi)

sp4c = StateSpace{CartesianIndex{2}, Float64}()
add_state!(sp4c, u0_4c, 1.0)
t_cur = 0.0

t_coarse = @elapsed for step in 1:n_steps
    global t_cur
    dt_s = min(dt, T-t_cur)
    fsp_step!(sp4c, sys4c, bc4c, t_cur, dt_s)
    t_cur += dt_s
end
mu4c = vmeans(sp4c.states, sp4c.probs, 2)
@printf("  %.2fs   |S_coarse|=%d   coarse means: nc1=%.1f nc2=%.1f\n",
        t_coarse, length(sp4c), mu4c[1], mu4c[2])
@printf("  (Binomial split -> voxel means would be: %.1f %.1f %.1f %.1f)\n",
        mu4c[1]/2, mu4c[1]/2, mu4c[2]/2, mu4c[2]/2)

# ══════════════════════════════════════════════════════════════
# K=4 V-CYCLE  (K=4->K=2, multiplicative prolongation)
# ══════════════════════════════════════════════════════════════
println("\n── K=4 V-cycle (K=4->K=2, multiplicative prolongation) ──")

# Need dynamic-pi table for the fine voxel size h4
println("  Computing dynamic-pi table...")
pair_sys4 = build_schlogl_rdme_system(model, VoxelGrid(2, h4, 0))
sp_pair4  = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair4, CartesianIndex(n1,n2), 0.0)
end
sp_pair4.probs[sp_pair4.index[CartesianIndex(n_lo,n_hi)]] = 1.0
A_pair4, = build_generator(sp_pair4, pair_sys4, rates, 0.0)
pi4 = compute_dynamic_pi(sp_pair4, A_pair4; n_max=n_max)
println("  Done.")

sp_vc4 = StateSpace{CartesianIndex{K4}, Float64}()
add_state!(sp_vc4, u0_4, 1.0)
snaps_vc4 = Dict{Float64, Tuple{Vector,Vector}}()
t_cur = 0.0

t_vc4 = @elapsed for step in 1:n_steps
    global t_cur, sp_vc4
    dt_s = min(dt, T-t_cur)
    sp_vc4 = two_level_vcycle_schlogl_injection(
        sp_vc4, model, g4, g4c, pi4, rates, t_cur, dt_s;
        coarse_n_max=2*n_max)
    t_cur += dt_s
    renormalize!(sp_vc4); prune_threshold!(sp_vc4, 1e-10)
    snap_t = T_SNAP[argmin(abs.(T_SNAP.-t_cur))]
    abs(t_cur-snap_t)<dt/2 && !haskey(snaps_vc4,snap_t) &&
        (snaps_vc4[snap_t] = (copy(sp_vc4.states), copy(sp_vc4.probs)))
end
mu_vc4 = vmeans(sp_vc4.states, sp_vc4.probs, K4)
@printf("  %.2fs   |S|=%d   means: %s\n", t_vc4, length(sp_vc4),
        join([@sprintf("%.1f",x) for x in mu_vc4], " "))

println("\n  TV vs direct at snapshots:")
for ts in T_SNAP
    if haskey(snaps_vc4,ts) && haskey(snaps4,ts)
        tv = tv_dist(snaps_vc4[ts]..., snaps4[ts]...)
        sts,prs = snaps_vc4[ts]
        mu = vmeans(sts,prs,K4)
        @printf("    t=%.1f  TV=%.4f  |S|=%d  means: %s\n",
                ts, tv, length(sts),
                join([@sprintf("%.1f",x) for x in mu], " "))
    end
end

# ══════════════════════════════════════════════════════════════
# K=8 V-CYCLE  (K=8->K=4, multiplicative prolongation)
# ══════════════════════════════════════════════════════════════
println("\n── K=8 V-cycle (K=8->K=4, multiplicative prolongation) ──")

K8 = 8; h8 = 1.0/K8
g8  = VoxelGrid(K8, h8, 0)
g8c = coarsen(g8)   # K=4
u0_8 = CartesianIndex(n_lo, n_lo, n_lo, n_hi, n_hi, n_hi, n_hi, n_hi)

println("  Computing dynamic-pi table for h=$(h8)...")
pair_sys8 = build_schlogl_rdme_system(model, VoxelGrid(2, h8, 0))
sp_pair8  = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair8, CartesianIndex(n1,n2), 0.0)
end
sp_pair8.probs[sp_pair8.index[CartesianIndex(n_lo,n_hi)]] = 1.0
A_pair8, = build_generator(sp_pair8, pair_sys8, rates, 0.0)
pi8 = compute_dynamic_pi(sp_pair8, A_pair8; n_max=n_max)
println("  Done.")

sp_vc8 = StateSpace{CartesianIndex{K8}, Float64}()
add_state!(sp_vc8, u0_8, 1.0)
snaps_vc8 = Dict{Float64, Tuple{Vector,Vector}}()
t_cur = 0.0

t_vc8 = @elapsed for step in 1:n_steps
    global t_cur, sp_vc8
    dt_s = min(dt, T-t_cur)
    sp_vc8 = two_level_vcycle_schlogl_injection(
        sp_vc8, model, g8, g8c, pi8, rates, t_cur, dt_s;
        coarse_n_max=2*n_max)
    t_cur += dt_s
    renormalize!(sp_vc8); prune_threshold!(sp_vc8, 1e-10)
    snap_t = T_SNAP[argmin(abs.(T_SNAP.-t_cur))]
    abs(t_cur-snap_t)<dt/2 && !haskey(snaps_vc8,snap_t) &&
        (snaps_vc8[snap_t] = (copy(sp_vc8.states), copy(sp_vc8.probs)))
    if step % 5 == 0
        @printf("  step %d/%.0f  t=%.1f  |S|=%d\n",
                step, n_steps, t_cur, length(sp_vc8))
    end
end
mu_vc8 = vmeans(sp_vc8.states, sp_vc8.probs, K8)
@printf("  %.2fs   |S|=%d   means: %s\n", t_vc8, length(sp_vc8),
        join([@sprintf("%.1f",x) for x in mu_vc8], " "))

println("\n  Snapshots:")
for ts in T_SNAP
    if haskey(snaps_vc8, ts)
        sts,prs = snaps_vc8[ts]
        mu = vmeans(sts,prs,K8)
        @printf("    t=%.1f  |S|=%d  means: %s\n", ts, length(sts),
                join([@sprintf("%.1f",x) for x in mu], " "))
    end
end

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════
println("\n" * "="^60)
println("SUMMARY  (t=2.0)")
println("="^60)
@printf("  K=4 direct FSP        |S|=%6d   %.2fs\n", length(sp4), t_direct)
@printf("  K=4 coarse (K=2)      |S|=%6d   %.2fs   (cheap, TV=1)\n", length(sp4c), t_coarse)
@printf("  K=4 V-cycle           |S|=%6d   %.2fs\n", length(sp_vc4), t_vc4)
@printf("  K=8 V-cycle           |S|=%6d   %.2fs\n", length(sp_vc8), t_vc8)
println("\nDone.")
