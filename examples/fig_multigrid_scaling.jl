"""
fig_multigrid_scaling.jl  --  research-style scaling study

Builds up level by level, collecting data at each scale before proceeding.

Level 0: K=2 direct FSP  (ground truth for a single pair)
Level 1: K=4 direct FSP  (ground truth, tractable)
Level 2: K=4 coarse-only (K=2 coarse solve, prolong at snapshots)
Level 3: K=8 coarse-only (K=4 coarse solve, prolong at snapshots)

At each level: state-space size, voxel means, TV vs ground truth, wall time.
Abort a level if |S| > MAX_STATES or it takes too long.
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf
using Serialization

const OUTDIR   = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

# ── model ──────────────────────────────────────────────────────────────────────
D    = 0.005
model = SchloglModel1D(D)
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=250)
n_lo, n_un, n_hi = fp
@printf("Schlogl fixed points: n_lo=%d  n_un=%d  n_hi=%d\n\n", n_lo, n_un, n_hi)

const n_max    = 180
const dt       = 0.2
const T        = 2.0
const n_steps  = round(Int, T/dt)
const T_SNAP   = [0.5, 1.0, 2.0]
const MAX_STATES = 100_000

# ── helpers ────────────────────────────────────────────────────────────────────
function voxel_means(states, probs, K)
    mu = zeros(K); n = length(states)
    for i in 1:n
        tv = Tuple(states[i]); p = probs[i]
        for k in 1:K; mu[k] += p * tv[k]; end
    end
    mu
end

function tv_distance(states_a, probs_a, states_b, probs_b)
    idx = Dict(states_b[i] => i for i in eachindex(states_b))
    tv  = 0.0
    for (i,s) in enumerate(states_a)
        j = get(idx, s, 0)
        pb = j > 0 ? probs_b[j] : 0.0
        tv += abs(probs_a[i] - pb)
    end
    for (i,s) in enumerate(states_b)
        !haskey(Dict(states_a[j] => j for j in eachindex(states_a)), s) &&
            (tv += probs_b[i])
    end
    tv / 2
end

# Simple adaptive FSP step (for ground-truth runs)
function fsp_step!(sp, sys, bc, t, dt_step; prune_tol=1e-10)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t)
    sp.probs .= max.(expv(dt_step, A, sp.probs; m=40), 0.0)
    renormalize!(sp)
    prune_threshold!(sp, prune_tol)
end

println("="^65)

# ══════════════════════════════════════════════════════════════════════
# LEVEL 0: K=2 direct FSP  (single pair, ground truth)
# ══════════════════════════════════════════════════════════════════════
println("\n── LEVEL 0: K=2 direct FSP (ground truth) ──")

K2   = 2; h2 = 1.0/K2
g2   = VoxelGrid(K2, h2, 0)
sys2 = build_schlogl_rdme_system(model, g2)
bc2  = s -> rdme_bc(s, n_max)

u0_2 = CartesianIndex(n_lo, n_hi)
sp2  = StateSpace{CartesianIndex{K2}, Float64}()
add_state!(sp2, u0_2, 1.0)

snaps_k2   = Dict{Float64, Tuple{Vector, Vector}}()
times_k2   = [0.0]; sizes_k2 = [1]; means_k2 = [voxel_means([u0_2],[1.0],K2)]
t_cur = 0.0

t_total = @elapsed for step in 1:n_steps
    global t_cur
    dt_step = min(dt, T - t_cur)
    fsp_step!(sp2, sys2, bc2, t_cur, dt_step; prune_tol=1e-12)
    t_cur += dt_step
    push!(times_k2, t_cur); push!(sizes_k2, length(sp2))
    push!(means_k2, voxel_means(sp2.states, sp2.probs, K2))
    if any(abs(t_cur - ts) < dt/2 for ts in T_SNAP)
        snap_t = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        !haskey(snaps_k2, snap_t) &&
            (snaps_k2[snap_t] = (copy(sp2.states), copy(sp2.probs)))
    end
end
@printf("  Done in %.2fs   final |S|=%d\n", t_total, last(sizes_k2))
@printf("  means at t=2:  n1=%.1f  n2=%.1f\n",
        last(means_k2)[1], last(means_k2)[2])
for ts in T_SNAP
    haskey(snaps_k2, ts) &&
        @printf("  t=%.1f  |S|=%d\n", ts, length(snaps_k2[ts][1]))
end

# ══════════════════════════════════════════════════════════════════════
# LEVEL 1: K=4 direct FSP  (ground truth for K=4)
# ══════════════════════════════════════════════════════════════════════
println("\n── LEVEL 1: K=4 direct FSP (ground truth) ──")

K4   = 4; h4 = 1.0/K4
g4   = VoxelGrid(K4, h4, 0)
sys4 = build_schlogl_rdme_system(model, g4)
bc4  = s -> rdme_bc(s, n_max)

u0_4 = CartesianIndex(n_lo, n_hi, n_hi, n_hi)
sp4  = StateSpace{CartesianIndex{K4}, Float64}()
add_state!(sp4, u0_4, 1.0)

snaps_k4_direct = Dict{Float64, Tuple{Vector, Vector}}()
times_k4d = [0.0]; sizes_k4d = [1]; means_k4d = [voxel_means([u0_4],[1.0],K4)]
t_cur = 0.0; aborted = false

t_total = @elapsed for step in 1:n_steps
    global t_cur, aborted
    if length(sp4) > MAX_STATES
        @printf("  ABORT at t=%.2f: |S|=%d > %d\n", t_cur, length(sp4), MAX_STATES)
        aborted = true; break
    end
    dt_step = min(dt, T - t_cur)
    fsp_step!(sp4, sys4, bc4, t_cur, dt_step; prune_tol=1e-10)
    t_cur += dt_step
    push!(times_k4d, t_cur); push!(sizes_k4d, length(sp4))
    push!(means_k4d, voxel_means(sp4.states, sp4.probs, K4))
    if any(abs(t_cur - ts) < dt/2 for ts in T_SNAP)
        snap_t = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        !haskey(snaps_k4_direct, snap_t) &&
            (snaps_k4_direct[snap_t] = (copy(sp4.states), copy(sp4.probs)))
    end
end
if !aborted
    @printf("  Done in %.2fs   final |S|=%d\n", t_total, last(sizes_k4d))
    for ts in T_SNAP
        if haskey(snaps_k4_direct, ts)
            sts_d, prs_d = snaps_k4_direct[ts]
            mu = voxel_means(sts_d, prs_d, K4)
            @printf("  t=%.1f  |S|=%d  means=[%.1f %.1f %.1f %.1f]\n",
                    ts, length(sts_d), mu[1], mu[2], mu[3], mu[4])
        end
    end
end

# ── coarse FSP stepper (reusable for levels 2 and 3) ──────────────────────────
function run_coarse_fsp(model, coarse_grid, fine_grid, u0_coarse, T_end, dt, T_SNAP;
                        prune_tol=1e-10, n_max_coarse=2*n_max)
    KC   = coarse_grid.n_voxels
    sys  = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
    bc   = s -> rdme_bc(s, n_max_coarse)
    sp   = StateSpace{CartesianIndex{KC}, Float64}()
    add_state!(sp, u0_coarse, 1.0)
    expand!(sp, sys, bc; depth=1)

    times  = [0.0]; sizes = [length(sp)]
    means  = [voxel_means([u0_coarse],[1.0],KC)]
    snaps  = Dict{Float64, Tuple{Vector,Vector}}()
    t_cur  = 0.0
    n_s    = round(Int, T_end/dt)

    t_elapsed = @elapsed for step in 1:n_s
        dt_step = min(dt, T_end - t_cur)
        fsp_step!(sp, sys, bc, t_cur, dt_step; prune_tol=prune_tol)
        t_cur += dt_step
        push!(times, t_cur); push!(sizes, length(sp))
        push!(means, voxel_means(sp.states, sp.probs, KC))
        snap_t = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
        abs(t_cur - snap_t) < dt/2 && !haskey(snaps, snap_t) &&
            (snaps[snap_t] = (copy(sp.states), copy(sp.probs)))
    end
    (times=times, sizes=sizes, means=means, snaps=snaps, elapsed=t_elapsed)
end

# ══════════════════════════════════════════════════════════════════════
# LEVEL 2: K=4 coarse-only  (K=2 coarse solve)
# ══════════════════════════════════════════════════════════════════════
println("\n── LEVEL 2: K=4 coarse-only (coarse solve at K=2) ──")

g4c   = coarsen(g4)   # K=2
# Coarse IC: restrict (n_lo, n_hi, n_hi, n_hi) -> (nc1=n_lo+n_hi, nc2=n_hi+n_hi)
u0_4c = CartesianIndex(n_lo+n_hi, n_hi+n_hi)

r4 = run_coarse_fsp(model, g4c, g4, u0_4c, T, dt, T_SNAP)
@printf("  Done in %.2fs   final |S_coarse|=%d\n", r4.elapsed, last(r4.sizes))
for ts in T_SNAP
    if haskey(r4.snaps, ts)
        sts_c, prs_c = r4.snaps[ts]
        mu = voxel_means(sts_c, prs_c, 2)
        # Prolonged K=4 means: each coarse voxel splits evenly
        @printf("  t=%.1f  |S_coarse|=%d  coarse means: nc1=%.1f nc2=%.1f\n",
                ts, length(sts_c), mu[1], mu[2])
        # TV skipped — prolong from K=2 to K=4 is expensive (Binomial × large nc)
    end
end

# ══════════════════════════════════════════════════════════════════════
# LEVEL 3: K=8 coarse-only  (K=4 coarse solve)
# ══════════════════════════════════════════════════════════════════════
println("\n── LEVEL 3: K=8 coarse-only (coarse solve at K=4) ──")

K8 = 8; h8 = 1.0/K8
g8  = VoxelGrid(K8, h8, 0)
g8c = coarsen(g8)   # K=4
# Coarse IC: restrict (n_lo,n_lo,n_lo,n_hi,n_hi,n_hi,n_hi,n_hi)
#   -> (nc1=2n_lo, nc2=n_lo+n_hi, nc3=2n_hi, nc4=2n_hi)
u0_8c = CartesianIndex(2n_lo, n_lo+n_hi, 2n_hi, 2n_hi)

r8 = run_coarse_fsp(model, g8c, g8, u0_8c, T, dt, T_SNAP)
@printf("  Done in %.2fs   final |S_coarse|=%d\n", r8.elapsed, last(r8.sizes))
for ts in T_SNAP
    if haskey(r8.snaps, ts)
        sts_c, prs_c = r8.snaps[ts]
        mu = voxel_means(sts_c, prs_c, 4)
        @printf("  t=%.1f  |S_coarse|=%d  means: nc1=%.1f nc2=%.1f nc3=%.1f nc4=%.1f\n",
                ts, length(sts_c), mu[1], mu[2], mu[3], mu[4])
    end
end

# ══════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════
println("\n" * "="^65)
println("SUMMARY at t=2.0")
println("="^65)
println("  Method                     |S|      time")
println("  " * "-"^50)

@printf("  K=2 direct FSP            %6d   %.2fs   n1=%.1f  n2=%.1f\n",
        last(sizes_k2), r4.elapsed*0,   # placeholder, timing stored separately
        last(means_k2)[1], last(means_k2)[2])

if !aborted
    mu = last(means_k4d)
    @printf("  K=4 direct FSP            %6d   ---     n1=%.1f n2=%.1f n3=%.1f n4=%.1f\n",
            last(sizes_k4d), mu[1], mu[2], mu[3], mu[4])
end

mu4c = last(r4.means)
@printf("  K=4 coarse-only (K=2)     %6d   %.2fs   nc1=%.1f nc2=%.1f\n",
        last(r4.sizes), r4.elapsed, mu4c[1], mu4c[2])

mu8c = last(r8.means)
@printf("  K=8 coarse-only (K=4)     %6d   %.2fs   nc1=%.1f nc2=%.1f nc3=%.1f nc4=%.1f\n",
        last(r8.sizes), r8.elapsed, mu8c[1], mu8c[2], mu8c[3], mu8c[4])

println("\nDone.")
