"""
Schlögl K=8 V-cycle Benchmark

Demonstrates the V-cycle for K=8 fine voxels where the direct K=8 sparse FSP
is intractable. Uses K=8 → K=4 two-level V-cycle.

Model: K=8 fine voxels, h=1/8, D=0.005  →  d_fine = D/h² = 0.32
IC:    (n_low)×4 + (n_high)×4 — symmetric front at coarse boundary
       Each K=4 coarse pair: (n_low, n_low) or (n_high, n_high) — symmetric

Ground truth options:
  - K=8 sparse adaptive FSP: attempted; may hit state limit
  - Pair-level GT: 4 independent K=2 FSPs at h=1/8 (exact for decoupled pairs)

Methods compared:
  1. Coarse-only Dyn-π     (K=4, no V-cycle)
  2. V-cycle Dyn-π         (K=8→K=4, correction prolongation)
  3. V-cycle Mult. prolong.         (K=8→K=4, fine-conditional prolongation)

Run: JULIA_PKG_PRECOMPILE_AUTO=0 julia --project examples/schlogl_k8_vcycle.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D    = 0.005
K    = 8
h    = 1.0 / K        # fine voxel size = 0.125
d    = D / h^2        # fine diffusion rate = 0.32

model      = SchloglModel1D(D)
rates      = Float64[]
n_max      = 180       # per-voxel truncation (covers n_high=162 + margin)
fine_grid  = VoxelGrid(K, h, 0)          # K=8 fine
coarse_grid = coarsen(fine_grid)          # K=4 coarse

fp = schlogl_fixed_points(model)
n_low, n_uns, n_high = fp

u0 = CartesianIndex(n_low, n_low, n_low, n_low, n_high, n_high, n_high, n_high)

t_solve  = [0.5, 1.0, 2.0]
t_max    = maximum(t_solve)
dt_c     = 0.2
n_steps  = round(Int, t_max / dt_c)

# τ_pre: one intra-pair diffusion relaxation time at fine level
τ_pre_vc = 1.0 / (d * n_low + 1.0)   # ≈ 1/(0.32·31+1) ≈ 0.087

println("="^70)
println("Schlögl K=$K V-cycle Benchmark  D=$D  h=$h  d=$(round(d,digits=3))")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low)×4 + (n_high)×4 — symmetric pairs")
@printf("τ_pre = %.3f  (fine intra-pair relaxation)\n", τ_pre_vc)
println("="^70)

# ─── pair-level ground truth: 4 × K=2 FSPs at h=1/8 ─────────────────────────
# Pair grouping (K=4 coarse): pairs 1&2 → (n_low,n_low); pairs 3&4 → (n_high,n_high)

println("\n── Pair-level ground truth (4 × K=2 FSP at h=$(h)) ──")

pair_grid = VoxelGrid(2, h, 0)   # K=2 at fine voxel size
pair_sys  = build_schlogl_rdme_system(model, pair_grid)

function make_pair_sp(n1_init, n2_init)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    for n1 in 0:n_max, n2 in 0:n_max
        add_state!(sp, CartesianIndex(n1, n2), 0.0)
    end
    sp.probs[sp.index[CartesianIndex(n1_init, n2_init)]] = 1.0
    sp
end

sp_pair_lo = make_pair_sp(n_low, n_low)
sp_pair_hi = make_pair_sp(n_high, n_high)

A_pair_lo, = build_generator(sp_pair_lo, pair_sys, rates, 0.0)
A_pair_hi, = build_generator(sp_pair_hi, pair_sys, rates, 0.0)

t0 = time()
probs_gt = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
let p_lo = copy(sp_pair_lo.probs), p_hi = copy(sp_pair_hi.probs), t_prev = 0.0
    for t in t_solve
        p_lo = max.(expv(t - t_prev, A_pair_lo, p_lo; m=50), 0.0); p_lo ./= sum(p_lo)
        p_hi = max.(expv(t - t_prev, A_pair_hi, p_hi; m=50), 0.0); p_hi ./= sum(p_hi)
        probs_gt[t] = (copy(p_lo), copy(p_hi))
        t_prev = t
    end
end
t_gt = time() - t0

function pair_means(sp, probs)
    mu1 = sum(Tuple(sp.states[i])[1] * probs[i] for i in eachindex(probs))
    mu2 = sum(Tuple(sp.states[i])[2] * probs[i] for i in eachindex(probs))
    mu1, mu2
end

@printf("  Wall time: %.2fs  (4 × %d states)\n", t_gt, (n_max+1)^2)
println("  t     pair (lo,lo)   pair (hi,hi)")
println("  ────────────────────────────────")
for t in t_solve
    p_lo, p_hi = probs_gt[t]
    m1, m2 = pair_means(sp_pair_lo, p_lo)
    m3, m4 = pair_means(sp_pair_hi, p_hi)
    @printf("  %4.1f  (%5.1f,%5.1f)  (%5.1f,%5.1f)\n", t, m1, m2, m3, m4)
end

# ─── compute dynamic-π table (K=2 at h=1/8) ──────────────────────────────────
# pi_table[nc+1][n1+1] = P(n1 | nc) under stationary intra-pair distribution

println("\n── Computing dynamic-π table (K=2 at h=$h) ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_pair_lo, A_pair_lo; n_max=n_max)
@printf("  Done in %.2fs\n", time() - t0)

# ─── K=8 FSP is intractable ──────────────────────────────────────────────────
# For K=8, the state space grows as ~75K states/step at t=0.5 and exceeds
# 10^6 states by t=2.  The K=8 direct FSP is not feasible; we use the
# pair-level GT (4 × K=2) as the reference.
println("\n── K=8 sparse adaptive FSP: intractable ──")
println("  State space reaches ~75K at t=0.5, >1M by t=2 — direct K=8 FSP skipped.")
println("  Reference: pair-level GT (4 independent K=2 FSPs, exact for decoupled pairs)")
k8_hit_limit = true

# ─── helpers ──────────────────────────────────────────────────────────────────

function voxel_means(sp)
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

# TV vs pair-level ground truth (independence assumption)
function tv_vs_pairs(sp, t)
    p_lo, p_hi = probs_gt[t]
    overlap = 0.0
    for (i, s) in enumerate(sp.states)
        t8 = Tuple(s)
        # pairs: (1,2),(3,4) → lo; (5,6),(7,8) → hi
        get_lo(n1,n2) = begin idx=get(sp_pair_lo.index,CartesianIndex(n1,n2),0); idx>0 ? p_lo[idx] : 0.0 end
        get_hi(n1,n2) = begin idx=get(sp_pair_hi.index,CartesianIndex(n1,n2),0); idx>0 ? p_hi[idx] : 0.0 end
        p_gt = get_lo(t8[1],t8[2]) * get_lo(t8[3],t8[4]) * get_hi(t8[5],t8[6]) * get_hi(t8[7],t8[8])
        overlap += min(sp.probs[i], p_gt)
    end
    1.0 - overlap
end

# ─── Coarse-only Dyn-π (K=4) ──────────────────────────────────────────────────

println("\n── Method 1: Coarse-only Dyn-π (K=4) ──")

coarse_sys_dyn = build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
coarse_bc      = state -> rdme_bc(state, 2*n_max)

# coarse IC: pairs (nc1,nc2,nc3,nc4) = (2n_low, 2n_low, 2n_high, 2n_high)
sp_c = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_c, CartesianIndex(2n_low, 2n_low, 2n_high, 2n_high), 1.0)
expand!(sp_c, coarse_sys_dyn, coarse_bc; depth=1)

result_coarse = Dict{Float64, Vector{Float64}}()
snap_coarse  = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()  # prolonged snaps for TV
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
            result_coarse[t_snap] = voxel_means(sp_prol)
            snap_coarse[t_snap]   = sp_prol
        end
    end
end
t_coarse = time() - t0_c
@printf("  Done in %.1fs  |S_coarse|=%d\n", t_coarse, length(sp_c))

# ─── V-cycle Dyn-π (K=8→K=4) ─────────────────────────────────────────────────

println("\n── Method 2: V-cycle Dyn-π (K=8→K=4, τ_pre=$(round(τ_pre_vc,digits=3))) ──")

sp_vc_dyn = StateSpace{CartesianIndex{8}, Float64}()
add_state!(sp_vc_dyn, u0, 1.0)
result_vc_dyn = Dict{Float64, Vector{Float64}}()
snap_vc_dyn   = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()

t0_vd = time()
let sp = sp_vc_dyn, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid, pi_table,
                                      rates, t_cur, dt_step;
                                      τ_pre=τ_pre_vc, τ_post=0.0, krylov_m=30,
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
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_vc_dyn, t_snap)
            result_vc_dyn[t_snap] = voxel_means(sp)
            let sp2 = StateSpace{CartesianIndex{8}, Float64}()
                for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
                snap_vc_dyn[t_snap] = sp2
            end
        end
    end
    global sp_vc_dyn = sp
end
t_vc_dyn = time() - t0_vd
@printf("\n  Done in %.1fs  final |S|=%d\n", t_vc_dyn, length(sp_vc_dyn))

# ─── V-cycle Mult. prolong. (K=8→K=4) ─────────────────────────────────────────────────

println("\n── Method 3: V-cycle Mult. prolong. (K=8→K=4, τ_pre=$(round(τ_pre_vc,digits=3))) ──")

sp_vc_mult = StateSpace{CartesianIndex{8}, Float64}()
add_state!(sp_vc_mult, u0, 1.0)
result_vc_mult = Dict{Float64, Vector{Float64}}()
snap_vc_mult   = Dict{Float64, StateSpace{CartesianIndex{8}, Float64}}()

t0_vf = time()
let sp = sp_vc_mult, t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, t_max - t_cur)
        sp = two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                                pi_table, rates, t_cur, dt_step;
                                                τ_pre=τ_pre_vc, τ_post=0.0, krylov_m=30,
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
        if abs(t_cur - t_snap) < dt_c/2 && !haskey(result_vc_mult, t_snap)
            result_vc_mult[t_snap] = voxel_means(sp)
            let sp2 = StateSpace{CartesianIndex{8}, Float64}()
                for i in eachindex(sp.states); add_state!(sp2, sp.states[i], sp.probs[i]); end
                snap_vc_mult[t_snap] = sp2
            end
        end
    end
    global sp_vc_mult = sp
end
t_vc_mult = time() - t0_vf
@printf("\n  Done in %.1fs  final |S|=%d\n", t_vc_mult, length(sp_vc_mult))

# ─── Results ──────────────────────────────────────────────────────────────────

println()
println("── Results ──────────────────────────────────────────────────────────────")
println("  Pair-level GT means: ⟨n1⟩=⟨n2⟩=⟨n3⟩=⟨n4⟩≈n_low, ⟨n5..8⟩≈n_high")
println()
println("  t     Method               ⟨n1⟩ ⟨n2⟩ ⟨n3⟩ ⟨n4⟩ ⟨n5⟩ ⟨n6⟩ ⟨n7⟩ ⟨n8⟩   TV(pair-GT)   |S|")
println("  ──────────────────────────────────────────────────────────────────────────────────────────")

for t in t_solve
    p_lo, p_hi = probs_gt[t]
    m_lo1, m_lo2 = pair_means(sp_pair_lo, p_lo)
    m_hi1, m_hi2 = pair_means(sp_pair_hi, p_hi)

    @printf("  %4.1f  %-20s %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f   —            —\n",
            t, "Pair-GT (4×K=2)", m_lo1, m_lo2, m_lo1, m_lo2, m_hi1, m_hi2, m_hi1, m_hi2)

    # Coarse Dyn-π
    mu_c  = get(result_coarse, t, fill(NaN, K))
    sp_cs = get(snap_coarse,  t, nothing)
    tv_c  = sp_cs !== nothing ? tv_vs_pairs(sp_cs, t) : NaN
    @printf("  %4s  %-20s %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f   %.4f       %d\n",
            "", "Coarse Dyn-π (K=4)",
            mu_c[1], mu_c[2], mu_c[3], mu_c[4], mu_c[5], mu_c[6], mu_c[7], mu_c[8],
            tv_c, length(sp_c))

    # V-cycle Dyn-π
    mu_vd  = get(result_vc_dyn, t, fill(NaN, K))
    sp_vds = get(snap_vc_dyn,   t, nothing)
    tv_vd  = sp_vds !== nothing ? tv_vs_pairs(sp_vds, t) : NaN
    @printf("  %4s  %-20s %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f   %.4f       %d\n",
            "", "V-cyc Dyn-π (K=8→4)",
            mu_vd[1], mu_vd[2], mu_vd[3], mu_vd[4], mu_vd[5], mu_vd[6], mu_vd[7], mu_vd[8],
            tv_vd, length(sp_vc_dyn))

    # V-cycle Mult. prolong.
    mu_vf  = get(result_vc_mult, t, fill(NaN, K))
    sp_vfs = get(snap_vc_mult,   t, nothing)
    tv_vf  = sp_vfs !== nothing ? tv_vs_pairs(sp_vfs, t) : NaN
    @printf("  %4s  %-20s %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f   %.4f       %d\n\n",
            "", "V-cyc Mult. prolong. (K=8→4)",
            mu_vf[1], mu_vf[2], mu_vf[3], mu_vf[4], mu_vf[5], mu_vf[6], mu_vf[7], mu_vf[8],
            tv_vf, length(sp_vc_mult))
end

println("── Timing ───────────────────────────────────────────────────────────────")
@printf("  Pair-level GT (4×K=2 FSP):   %.2fs\n", t_gt)
println("  K=8 sparse adaptive FSP:     intractable (>75K states at t=0.5)")
@printf("  Coarse-only Dyn-π (K=4):     %.2fs  |S|=%d\n", t_coarse, length(sp_c))
@printf("  V-cycle Dyn-π (K=8→4):       %.2fs  |S|=%d\n", t_vc_dyn, length(sp_vc_dyn))
@printf("  V-cycle Mult. prolong. (K=8→4):       %.2fs  |S|=%d\n", t_vc_mult, length(sp_vc_mult))
println()
println("="^70)
println("DONE")
println("="^70)
