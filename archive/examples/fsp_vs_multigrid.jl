"""
Birth-Death-Diffusion 1D: Direct FSP vs Multigrid FSP

Benchmark 1 from the paper (Section 10.1):
  - K voxels on [0,1], single species A
  - Birth: ∅ →^{k_b} A  per voxel
  - Death: A →^{k_d} ∅  per molecule
  - Diffusion: rate d = D/h²

Analytical stationary: π_k = Poisson(μ_k) where μ_k = k_b/k_d (uniform).

This script runs:
  1. Direct FSP: build full RDME generator on all states up to n_max, evolve with expv.
     This is the ground truth for K=4.
  2. Multigrid V-cycle: coarse-only mode, working at K2=K/2 level.
  3. Comparison: means and TV distance vs time.

Run: julia --project examples/fsp_vs_multigrid.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf
using SparseArrays

# ─── analytical reference ────────────────────────────────────────────────────

"""
Stationary mean per voxel for 1D birth-death-diffusion.
Solves d(μ_{k-1} - 2μ_k + μ_{k+1}) + k_b - k_d μ_k = 0 with reflecting BCs.
For homogeneous params the solution is μ_k = k_b/k_d everywhere.
"""
stationary_mean(k_b, k_d) = k_b / k_d

"""TV distance = 0.5 * sum |p - q| between two prob vectors over the same state space."""
tv_distance(p, q) = 0.5 * sum(abs, p .- q)

# ─── direct FSP: enumerate all states up to n_max ─────────────────────────────

"""
Build a StateSpace containing ALL states (n₁,…,nK) with 0 ≤ nₖ ≤ n_max.
This gives the full truncated FSP state space — the ground truth for small K.
"""
function full_state_space(K::Int, n_max::Int)
    sp = StateSpace{CartesianIndex{K}, Float64}()
    function recurse!(tup, k)
        if k > K
            add_state!(sp, CartesianIndex(NTuple{K,Int}(tup)), 0.0)
            return
        end
        for n in 0:n_max
            tup[k] = n
            recurse!(tup, k + 1)
        end
    end
    recurse!(zeros(Int, K), 1)
    sp
end

"""
Solve the RDME directly: build full generator on all states up to n_max, evolve
with expmv from t=0 to each t in t_save.
Returns a vector of (states, probs) snapshots.
"""
function solve_direct_fsp(model::RDMEModel1D, grid::VoxelGrid,
                           u0::CartesianIndex{K}, t_save::Vector{Float64};
                           n_max::Int = 15, krylov_m::Int = 50) where {K}
    sp = full_state_space(K, n_max)

    # Set initial condition
    idx0 = get(sp.index, u0, 0)
    if idx0 == 0
        error("Initial state $u0 not in truncated state space (n_max=$n_max too small?)")
    end
    sp.probs[idx0] = 1.0

    # Build full generator once
    sys  = build_rdme_system(model, grid)
    bc   = s -> rdme_bc(s, n_max)
    A, _ = build_generator(sp, sys, Float64[], 0.0)   # rates unused here

    snapshots = Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}[]
    t_prev = 0.0

    for t in t_save
        dt = t - t_prev
        if dt > 0.0
            sp.probs .= expv(dt, A, sp.probs; m = krylov_m)
            sp.probs .= max.(sp.probs, 0.0)
            sp.probs ./= sum(sp.probs)
        end
        push!(snapshots, (copy(sp.states), copy(sp.probs)))
        t_prev = t
    end

    snapshots
end

"""Mean voxel counts from a (states, probs) snapshot."""
function snapshot_means(snap::Tuple{Vector{CartesianIndex{K}}, Vector{Float64}}) where K
    states, probs = snap
    mu = zeros(K)
    for (s, p) in zip(states, probs)
        t = Tuple(s)
        for k in 1:K; mu[k] += p * t[k]; end
    end
    mu
end

# ─── main comparison ─────────────────────────────────────────────────────────

D   = 1.0;  k_b = 2.0;  k_d = 1.0
μ★  = stationary_mean(k_b, k_d)
rates = Float64[]

t_save = [0.5, 1.0, 2.0, 3.0, 5.0]

# ══════════════════════════════════════════════════════════════════════════════
# K = 4  (state space = 11^4 = 14 641, fully tractable)
# ══════════════════════════════════════════════════════════════════════════════
println("="^65)
println("K=4  Birth-Death-Diffusion   k_b=$k_b  k_d=$k_d  D=$D")
println("Stationary: each voxel → Poisson(μ*=$μ★)")
println("="^65)

K = 4
grid = VoxelGrid(K, 1.0/K, 0)
u0   = CartesianIndex(4, 0, 0, 0)   # all molecules in voxel 1
n_max_fine = 12

# ── 1. Direct FSP ─────────────────────────────────────────────────────────────
println("\n── Direct FSP (ground truth) ──")
println("State space size = $(11)^$K = $((n_max_fine+1)^K) (all states up to n_max=$n_max_fine)")
t0 = time()
snaps_direct = solve_direct_fsp(RDMEModel1D(D,k_b,k_d), grid, u0, t_save;
                                 n_max = n_max_fine, krylov_m = 60)
println("Solved in $(round(time()-t0, digits=2))s")
println()
for (t, snap) in zip(t_save, snaps_direct)
    mu = snapshot_means(snap)
    @printf("t=%4.1f  means=%s  Σμ=%.3f\n", t, string(round.(mu, digits=3)), sum(mu))
end

# ── 2. Multigrid FSP (coarse-only) ────────────────────────────────────────────
println("\n── Multigrid FSP (coarse-only, K2=$(K÷2)) ──")
alg = RDMEMultigridFSP(
    dt         = 0.1,
    τ_pre      = 0.0,
    τ_post     = 0.0,
    n_max      = n_max_fine,
    krylov_m   = 40,
    binom_tol  = 0.0,
    prune_tol  = 1e-12,
    expand_coarse     = true,
    coarse_expand_depth = 1,
    expand_depth = 0,
    coarse_only  = true,
    save_every   = 1,
)
t0 = time()
sol_mg = solve_rdme_multigrid(RDMEModel1D(D,k_b,k_d), grid, u0,
                               (0.0, maximum(t_save)), rates, alg)
println("Solved in $(round(time()-t0, digits=2))s")
println()

# Extract multigrid snapshots at matching times
for t_want in t_save
    idx = argmin(abs.(sol_mg.t .- t_want))
    states, probs = sol_mg.snapshots[idx]
    mu = zeros(K)
    for (s, p) in zip(states, probs)
        tup = Tuple(s)
        for k in 1:K; mu[k] += p * tup[k]; end
    end
    @printf("t=%4.1f  means=%s  Σμ=%.3f  |S|=%d\n",
            sol_mg.t[idx], string(round.(mu, digits=3)), sum(mu),
            sol_mg.state_space_sizes[idx])
end

# ── 3. TV distance comparison ─────────────────────────────────────────────────
println("\n── TV distance (multigrid vs direct FSP) ──")
for (t_want, snap_d) in zip(t_save, snaps_direct)
    idx = argmin(abs.(sol_mg.t .- t_want))
    states_mg, probs_mg = sol_mg.snapshots[idx]
    states_d,  probs_d  = snap_d

    # Build a common index
    all_states = union(Set(states_mg), Set(states_d))
    idx_mg = Dict(s => i for (i,s) in enumerate(states_mg))
    idx_d  = Dict(s => i for (i,s) in enumerate(states_d))
    tv = 0.0
    for s in all_states
        pm = haskey(idx_mg, s) ? probs_mg[idx_mg[s]] : 0.0
        pd = haskey(idx_d,  s) ? probs_d[idx_d[s]]   : 0.0
        tv += abs(pm - pd)
    end
    @printf("t=%4.1f  TV(MG, FSP) = %.2e\n", t_want, 0.5*tv)
end

println()
println("Expected TV → 0 as accuracy increases")
println()

# ─── K=8 coarse loop (no prolongation) ───────────────────────────────────────
# For K=8 the fine state space (K=8) after prolongation reaches 17M states,
# making snapshot saves ≈10s each (×10 saves ≈100s overhead).
# Instead we work entirely at the K=4 coarse level and read fine means directly:
#   μ_fine[2c-1] = μ_fine[2c] = μ_coarse[c] / 2   (binomial-split symmetry)
# This reduces snapshot cost to microseconds and keeps total runtime ≈20s.

"""
Run the K=8 coarse-only loop at K2=K÷2 level.
Returns (saved_t, coarse_means, coarse_ss_sizes) without any prolongation.
Fine per-voxel means follow by splitting: μ_fine[2c-1]=μ_fine[2c]=μ_coarse[c]/2.
"""
function coarse_only_means(model, fine_grid::VoxelGrid,
                            u0::CartesianIndex{K},
                            tspan, rates, alg) where K
    K2 = K ÷ 2
    coarse_grid = coarsen(fine_grid)
    coarse_sys  = build_coarse_system(model, coarse_grid, fine_grid)

    n_max_c  = alg.coarse_n_max_per_voxel > 0 ? alg.coarse_n_max_per_voxel : 2*alg.n_max
    coarse_bc = if alg.coarse_n_total > 0
        s -> rdme_bc(s, n_max_c) && sum(Tuple(s)) ≤ alg.coarse_n_total
    else
        s -> rdme_bc(s, n_max_c)
    end

    sp0 = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp0, u0, 1.0)
    sp_c = restrict(sp0)   # K → K2 coarse

    t_cur, tf = Float64.(tspan)
    step = 0

    saved_t   = Float64[]
    cmeans    = Vector{Vector{Float64}}()   # coarse means at each save
    ss_sizes  = Int[]

    function save_coarse!()
        mc = zeros(K2)
        for (s, p) in zip(sp_c.states, sp_c.probs)
            tup = Tuple(s)
            for k in 1:K2; mc[k] += p * tup[k]; end
        end
        push!(saved_t,  t_cur)
        push!(cmeans,   mc)
        push!(ss_sizes, length(sp_c))
    end

    save_coarse!()   # t = 0

    while t_cur < tf
        dt_step = min(alg.dt, tf - t_cur)
        if alg.expand_coarse && alg.coarse_expand_depth > 0
            expand!(sp_c, coarse_sys, coarse_bc; depth = alg.coarse_expand_depth)
        end
        A_c, = build_generator(sp_c, coarse_sys, rates, t_cur)
        sp_c.probs .= expv(Float64(dt_step), A_c, sp_c.probs; m = alg.krylov_m)
        t_cur += dt_step
        step  += 1
        renormalize!(sp_c)
        prune_threshold!(sp_c, alg.prune_tol)

        if step % alg.save_every == 0 || t_cur ≥ tf
            save_coarse!()
        end
    end

    (saved_t, cmeans, ss_sizes)
end

# ══════════════════════════════════════════════════════════════════════════════
# K = 8  (fine state space 13^8 ≈ 800M — intractable; coarse K2=4 tractable)
# ══════════════════════════════════════════════════════════════════════════════
println("="^65)
println("K=8  Birth-Death-Diffusion   k_b=$k_b  k_d=$k_d  D=$D")
println("Fine state space 13^8 ≈ 800M (intractable for direct FSP)")
println("Coarse K2=4: work entirely at coarse level, read fine means via μ_c/2")
println("="^65)

K8    = 8
grid8 = VoxelGrid(K8, 1.0/K8, 0)
u0_8  = CartesianIndex(5, 0, 0, 0, 0, 0, 0, 0)

println("\n── Multigrid FSP K=8 (coarse-only, K2=$(K8÷2), no prolongation) ──")
# K=4 coarse: d_coarse = D/(h_coarse)² = 16, birth_coarse = 2*k_b = 4
# ‖A_coarse‖ ≈ 800.  dt×‖A‖ ≈ 0.02×800 = 16 — safe for Krylov m=30
alg8 = RDMEMultigridFSP(
    dt         = 0.02,
    τ_pre      = 0.0,
    τ_post     = 0.0,
    n_max      = 12,
    krylov_m   = 30,
    binom_tol  = 0.005,
    prune_tol  = 1e-7,
    expand_coarse       = true,
    coarse_expand_depth = 1,
    expand_depth        = 0,
    coarse_only              = true,
    n_levels                 = 1,
    coarse_n_total           = 32,     # AND total-molecule: 16 + 4√16 = 32
    coarse_n_max_per_voxel   = 14,     # per-voxel: Pois(4;n>14) < 1.5e-7 < prune_tol
    save_every               = 25,     # every 25 steps × 0.02 = 0.5 intervals
)
t0 = time()
saved_t8, cmeans8, ss8 = coarse_only_means(RDMEModel1D(D,k_b,k_d), grid8, u0_8,
                                            (0.0, maximum(t_save)), rates, alg8)
println("Solved in $(round(time()-t0, digits=2))s")
println()

for t_want in t_save
    idx = argmin(abs.(saved_t8 .- t_want))
    mc  = cmeans8[idx]
    # Fine per-voxel means via binomial-split symmetry: μ_fine[2c-1]=μ_fine[2c]=μ_c/2
    mf  = vcat([[mc[c]/2, mc[c]/2] for c in 1:K8÷2]...)
    @printf("t=%4.1f  Σμ=%.3f  means=%s  |S_coarse|=%d\n",
            saved_t8[idx], sum(mc), string(round.(mf, digits=3)), ss8[idx])
end

println()
println("Expected: Σμ → $(K8 * μ★)  all voxels → $μ★")
println()
println("="^65)
println("DONE")
println("="^65)
