"""
FP-FSP for RDME: fine-grid adaptive FSP vs multigrid coarse-only

Model: 1D birth-death-diffusion on K voxels
  ∅  →^{k_b} A,   A →^{k_d} ∅,   A_k ⇌ A_{k±1}  rate d = D/h²

This script benchmarks three approaches on K=4 and K=8:
  1. Direct FSP: enumerate all states up to n_max; build generator once; evolve.
     (Ground truth, tractable only for small K.)
  2. FP-FSP (fine grid): adaptive dt = ε_dt/Φ; expand one shell per step;
     flux-aware compress!; work entirely on the K-voxel fine state space.
  3. Multigrid coarse-only: work at K/2-voxel coarse level; fixed dt;
     prolongate to fine only at snapshot times.

Key result: for diffusion-dominated RDMEs, Φ ∝ d = D/h² grows quadratically
with resolution. At K=8, d=64 forces dt~5×10⁻⁴ and ~10 000 iterations to
reach t=5, making fine-grid FP-FSP infeasible. Coarse-graining removes this
stiffness: d_coarse = D/h_c² is 4× smaller at each level, restoring tractability.

FP-FSP remains the right tool for reaction-dominated systems with spatial
bottlenecks (e.g. bistable Schlögl fronts) where reaction rates — not
diffusion — dominate Φ and coarse-graining introduces significant error.

Run: julia --project examples/rdme_fpfsp.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf
using Statistics: median

# ─── FP-FSP RDME solver ───────────────────────────────────────────────────────

"""
Adaptive FSP on the fine RDME grid.
Loop: expand → build_generator → dt = ε_dt/Φ → expv → compress! → renormalize!
"""
function solve_fpfsp_rdme(model::RDMEModel1D, grid::VoxelGrid,
                           u0::CartesianIndex{K},
                           tspan::Tuple{<:Real,<:Real},
                           rates;
                           n_max::Int              = 15,
                           ε_dt::Float64           = 1.0,
                           prob_quantile::Float64  = 0.1,
                           flux_tolerance::Float64 = 1e-6,
                           t_save::Vector{Float64} = Float64[],
                           krylov_m::Int           = 30,
                           max_iter::Int           = typemax(Int)) where {K}
    sys = build_rdme_system(model, grid)
    bc  = s -> rdme_bc(s, n_max)

    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, sys, bc; depth = 1)

    t, tf = Float64.(tspan)
    iter  = 0

    snaps    = Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}[]
    saved_t  = Float64[]
    ss_sizes = Int[]
    t_log    = Float64[]
    dt_log   = Float64[]
    sz_log   = Int[]

    save_ptr = 1
    function _maybe_save()
        while save_ptr <= length(t_save) && t >= t_save[save_ptr]
            push!(saved_t,  t)
            push!(snaps,   (copy(sp.states), copy(sp.probs)))
            push!(ss_sizes, length(sp))
            save_ptr += 1
        end
    end
    _maybe_save()   # t = 0

    while t < tf && iter < max_iter
        iter += 1
        expand!(sp, sys, bc; depth = 1)
        A, = build_generator(sp, sys, rates, t)
        Φ  = dot(-diag(A), sp.probs)
        dt = Φ > 0 ? min(ε_dt / Φ, tf - t) : (tf - t)
        sp.probs .= expv(dt, A, sp.probs; m = krylov_m)
        t += dt
        compress!(sp, sys, rates, t, prob_quantile; flux_tolerance = flux_tolerance)
        renormalize!(sp)
        push!(t_log, t); push!(dt_log, dt); push!(sz_log, length(sp))
        _maybe_save()
    end

    (saved_t, snaps, ss_sizes, iter, t_log, dt_log, sz_log)
end

# ─── helpers ──────────────────────────────────────────────────────────────────

function snapshot_means(snap::Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}) where K
    states, probs = snap
    mu = zeros(K)
    for (s, p) in zip(states, probs)
        tup = Tuple(s); for k in 1:K; mu[k] += p * tup[k]; end
    end
    mu
end

function tv_snapshots(snap_a, snap_b)
    sa, pa = snap_a;  sb, pb = snap_b
    ia = Dict(s => i for (i,s) in enumerate(sa))
    ib = Dict(s => i for (i,s) in enumerate(sb))
    tv = 0.0
    for s in union(keys(ia), keys(ib))
        pav = haskey(ia,s) ? pa[ia[s]] : 0.0
        pbv = haskey(ib,s) ? pb[ib[s]] : 0.0
        tv += abs(pav - pbv)
    end
    0.5 * tv
end

# ─── multigrid coarse-only ────────────────────────────────────────────────────

function coarse_only_means(model, fine_grid::VoxelGrid,
                            u0::CartesianIndex{K},
                            tspan, rates, alg) where K
    K2 = K ÷ 2
    coarse_grid = coarsen(fine_grid)
    coarse_sys  = build_coarse_system(model, coarse_grid, fine_grid)
    n_max_c  = alg.coarse_n_max_per_voxel > 0 ? alg.coarse_n_max_per_voxel : 2*alg.n_max
    coarse_bc = alg.coarse_n_total > 0 ?
        (s -> rdme_bc(s, n_max_c) && sum(Tuple(s)) ≤ alg.coarse_n_total) :
        (s -> rdme_bc(s, n_max_c))

    sp0 = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp0, u0, 1.0)
    sp_c = restrict(sp0)
    t_cur, tf = Float64.(tspan)
    step = 0
    saved_t = Float64[]; cmeans = Vector{Vector{Float64}}(); ss_sizes = Int[]

    function save_coarse!()
        mc = zeros(K2)
        for (s, p) in zip(sp_c.states, sp_c.probs)
            tup = Tuple(s); for k in 1:K2; mc[k] += p * tup[k]; end
        end
        push!(saved_t, t_cur); push!(cmeans, mc); push!(ss_sizes, length(sp_c))
    end
    save_coarse!()

    while t_cur < tf
        dt_step = min(alg.dt, tf - t_cur)
        alg.expand_coarse && alg.coarse_expand_depth > 0 &&
            expand!(sp_c, coarse_sys, coarse_bc; depth = alg.coarse_expand_depth)
        A_c, = build_generator(sp_c, coarse_sys, rates, t_cur)
        sp_c.probs .= expv(Float64(dt_step), A_c, sp_c.probs; m = alg.krylov_m)
        t_cur += dt_step; step += 1
        renormalize!(sp_c); prune_threshold!(sp_c, alg.prune_tol)
        (step % alg.save_every == 0 || t_cur ≥ tf) && save_coarse!()
    end
    (saved_t, cmeans, ss_sizes)
end

# ─── direct FSP (K=4 ground truth) ───────────────────────────────────────────

function full_state_space(K::Int, n_max::Int)
    sp = StateSpace{CartesianIndex{K}, Float64}()
    function recurse!(tup, k)
        k > K && (add_state!(sp, CartesianIndex(NTuple{K,Int}(tup)), 0.0); return)
        for n in 0:n_max; tup[k] = n; recurse!(tup, k+1); end
    end
    recurse!(zeros(Int, K), 1); sp
end

function solve_direct_fsp(model, grid, u0::CartesianIndex{K}, t_save;
                           n_max=12, krylov_m=60) where K
    sp = full_state_space(K, n_max)
    sp.probs[sp.index[u0]] = 1.0
    sys = build_rdme_system(model, grid)
    A,  = build_generator(sp, sys, Float64[], 0.0)
    snaps = Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}[]
    t_prev = 0.0
    for t in t_save
        dt = t - t_prev
        if dt > 0
            sp.probs .= max.(expv(dt, A, sp.probs; m=krylov_m), 0.0)
            sp.probs ./= sum(sp.probs)
        end
        push!(snaps, (copy(sp.states), copy(sp.probs))); t_prev = t
    end
    snaps
end

# ═════════════════════════════════════════════════════════════════════════════
D = 1.0;  k_b = 2.0;  k_d = 1.0;  μ★ = k_b / k_d
rates = Float64[]
t_save = [0.5, 1.0, 2.0, 3.0, 5.0]

# ═════════════════════════════════════════════════════════════════════════════
# K = 4
# ═════════════════════════════════════════════════════════════════════════════
K4 = 4;  g4 = VoxelGrid(K4, 1.0/K4, 0);  u0_4 = CartesianIndex(4, 0, 0, 0)
n_max4 = 12;  d4 = D / (1.0/K4)^2

println("="^65)
println("K=4  d=D/h²=$(d4)  k_b=$(k_b)  k_d=$(k_d)  μ*=$(μ★)  Σμ*=$(K4*μ★)")
println("="^65)

# ── Direct FSP ────────────────────────────────────────────────────────────────
println("\n── Direct FSP (ground truth, $((n_max4+1)^K4) states) ──")
t0 = time()
snaps_gt = solve_direct_fsp(RDMEModel1D(D,k_b,k_d), g4, u0_4, t_save; n_max=n_max4)
t_gt = time() - t0
@printf("  Wall time: %.2fs\n", t_gt)
for (t, snap) in zip(t_save, snaps_gt)
    mu = snapshot_means(snap)
    @printf("  t=%4.1f  Σμ=%.3f  means=%s\n", t, sum(mu), string(round.(mu, digits=3)))
end

# ── FP-FSP fine K=4 ──────────────────────────────────────────────────────────
println("\n── FP-FSP (fine K=$(K4), adaptive dt) ──")
t0 = time()
saved4, snaps4, sizes4, iters4, tlog4, dtlog4, szlog4 =
    solve_fpfsp_rdme(RDMEModel1D(D,k_b,k_d), g4, u0_4, (0.0, maximum(t_save)), rates;
                     n_max=n_max4, ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1e-6,
                     t_save=t_save, krylov_m=30)
t_fp4 = time() - t0
@printf("  Wall time: %.2fs,  iters=%d,  dt in [%.3e, %.3e]\n",
        t_fp4, iters4, minimum(dtlog4), maximum(dtlog4))
for (t, snap, sz) in zip(saved4, snaps4, sizes4)
    mu  = snapshot_means(snap)
    idx = argmin(abs.(t_save .- t))
    tv  = tv_snapshots(snap, snaps_gt[idx])
    @printf("  t=%4.1f  Σμ=%.3f  |S|=%5d  TV(gt)=%.2e\n", t, sum(mu), sz, tv)
end

# ── Multigrid coarse-only K=4 → K=2 ─────────────────────────────────────────
println("\n── Multigrid coarse-only (K=$(K4÷2)) ──")
alg4 = RDMEMultigridFSP(dt=0.1, n_max=n_max4, krylov_m=40,
                         prune_tol=1e-12, expand_coarse=true, coarse_expand_depth=1,
                         coarse_only=true, save_every=1)
t0 = time()
st4, cm4, ss4 = coarse_only_means(RDMEModel1D(D,k_b,k_d), g4, u0_4,
                                   (0.0, maximum(t_save)), rates, alg4)
t_mg4 = time() - t0
@printf("  Wall time: %.2fs\n", t_mg4)
for t_want in t_save
    idx = argmin(abs.(st4 .- t_want))
    mc  = cm4[idx]
    mf  = vcat([[mc[c]/2, mc[c]/2] for c in 1:K4÷2]...)
    @printf("  t=%4.1f  Σμ=%.3f  |S_c|=%3d  means=%s\n",
            st4[idx], sum(mc), ss4[idx], string(round.(mf, digits=3)))
end

println("\n── K=4 summary ──────────────────────────────────────────")
@printf("  Direct FSP:      %5.2fs  (%d states, 1 generator build)\n", t_gt, (n_max4+1)^K4)
@printf("  FP-FSP fine:     %5.2fs  (%d iters, |S|~%d, TV<%.0e)\n",
        t_fp4, iters4, round(Int,median(szlog4)), maximum(map(i->tv_snapshots(snaps4[i],snaps_gt[i]),1:length(snaps4))))
@printf("  Multigrid K=%d:   %5.2fs  (%d coarse states)\n",
        K4÷2, t_mg4, maximum(ss4))
@printf("\n  Diffusion stiffness: d=%g → Φ~%.0f → dt~%.3f\n",
        d4, 1/median(dtlog4), median(dtlog4))
println("  → FP-FSP is correct but $(round(Int,t_fp4/t_gt))× slower than direct FSP,")
println("    $(round(Int,t_fp4/t_mg4))× slower than multigrid, due to diffusion-forced small dt.")

# ═════════════════════════════════════════════════════════════════════════════
# K = 8
# ═════════════════════════════════════════════════════════════════════════════
println()
println("="^65)
K8 = 8;  g8 = VoxelGrid(K8, 1.0/K8, 0);  u0_8 = CartesianIndex(5,0,0,0,0,0,0,0)
n_max8 = 12;  d8 = D / (1.0/K8)^2;  d8c = D / (1.0/(K8÷2))^2

println("K=8  d=D/h²=$(d8)  (coarse d_c=$(d8c))")
println("Fine state space $(n_max8+1)^8 ≈ $(round(Int,(n_max8+1)^8/1e6))M — direct FSP intractable")
println("="^65)

# ── FP-FSP K=8 analysis (not run — shown analytically) ───────────────────────
println("\n── FP-FSP fine K=8 (analytical estimate) ──")
# At stationarity: N_total ~ K8*μ★ = 16 molecules; w(x) ~ (k_d + 2d8) per mol ~ 129
# Φ ~ N_total * (k_d + 2*d8) ~ 16 * 129 ~ 2064; dt ~ ε_dt/Φ ~ 5e-4
Φ_est8 = K8 * μ★ * (k_d + 2*d8)
dt_est8 = 1.0 / Φ_est8
iters_est8 = round(Int, maximum(t_save) / dt_est8)
@printf("  Φ ≈ K·μ*·(k_d+2d) = %d·%.1f·%.0f ≈ %.0f\n",
        K8, μ★, k_d+2*d8, Φ_est8)
@printf("  dt ≈ ε_dt/Φ ≈ %.2e  →  ~%d iterations to t=%.1f\n",
        dt_est8, iters_est8, maximum(t_save))
@printf("  Compare K=4: Φ~%.0f, dt~%.3f, ~%d iterations\n",
        1/median(dtlog4), median(dtlog4), iters4)
@printf("  Scaling: d ∝ 1/h² → dt ∝ h² → iters ∝ K²  (%.0f× more at K=8 vs K=4)\n",
        (d8/d4)^1)
println("  Fine-grid FP-FSP on K=8 is dominated by diffusion stiffness.")
println("  Stoichiometric expansion of CartesianIndex{8} states also")
println("  generates O(K·n_max) new candidates per state per iteration.")

# ── Multigrid coarse-only K=8 → K=4 ─────────────────────────────────────────
println("\n── Multigrid coarse-only (K=$(K8÷2)) ──")
alg8 = RDMEMultigridFSP(
    dt=0.02, n_max=12, krylov_m=30,
    prune_tol=1e-7, expand_coarse=true, coarse_expand_depth=1,
    coarse_only=true, n_levels=1,
    coarse_n_total=32, coarse_n_max_per_voxel=14,
    save_every=25,
)
t0 = time()
st8, cm8, ss8 = coarse_only_means(RDMEModel1D(D,k_b,k_d), g8, u0_8,
                                   (0.0, maximum(t_save)), rates, alg8)
t_mg8 = time() - t0
@printf("  Coarse d_c=%.0f (%.0fx smaller → %.0fx larger dt feasible)\n",
        d8c, d8/d8c, d8/d8c)
@printf("  Wall time: %.2fs\n", t_mg8)
for t_want in t_save
    idx = argmin(abs.(st8 .- t_want))
    mc  = cm8[idx]
    mf  = vcat([[mc[c]/2, mc[c]/2] for c in 1:K8÷2]...)
    @printf("  t=%4.1f  Σμ=%.3f  |S_c|=%5d  means=%s\n",
            st8[idx], sum(mc), ss8[idx], string(round.(mf, digits=3)))
end

println("\n── K=8 summary ──────────────────────────────────────────")
t_fp8_est = t_fp4 / iters4 * iters_est8 * (d8/d4)   # scale by |S| growth too
@printf("  FP-FSP fine K=8:  intractable  (est. ~%.0fs, ~%d iters)\n",
        t_fp8_est, iters_est8)
@printf("  Multigrid K=4:    %.2fs actual\n", t_mg8)
@printf("  Multigrid speedup vs FP-FSP estimate: ~%.0fx\n", t_fp8_est / t_mg8)
println()
println("Note: FP-FSP becomes competitive when reaction flux (not diffusion)")
println("dominates Φ, e.g. Schlögl bistable fronts where spatial inhomogeneity")
println("creates transit states between attractors. There, coarse-graining")
println("introduces O(1) error while FP-FSP resolves the fine-grained front.")
println()
println("="^65); println("DONE"); println("="^65)
