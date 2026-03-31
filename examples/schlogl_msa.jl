"""
Schlögl MSA (Markov State Aggregation) Benchmark

Demonstrates the failure modes of MSA aggregation for bistable spatial RDME
and the partial fix from dynamic-π prolongation.

Model: K=2 voxels, Schlögl bistable chemistry + diffusion  (D=0.005)
  Per-voxel reactions: ∅→X (c3), X→∅ (c4), 2X→3X (c1), 3X→2X (c2)
  Diffusion: X_1 ⇌ X_2 at rate d = D/h²
  Stable states: n_low ≈ 31, n_high ≈ 162

Setup: IC at (n_low, n_high) — a spatial front between the two stable states.
With D=0.005 (d=0.02), the reaction restoring force (≈1.8/molecule near n_low)
dominates diffusion flux (d·n_high ≈ 3.2), so the front is metastable.

Four methods:
  1. Direct FSP (ground truth): enumerate all (n1,n2) up to n_max.
  2. Coarse MSA (Binomial-π): nc=n1+n2; prolong via Binomial(nc,1/2).
     Completely wrong: unimodal at n1=nc/2=97, TV≈1.
  3. Dynamic-π (stationary): solve A^T π=0 for the 2-voxel stationary,
     extract conditional slices p(n1|nc). Bimodal at n1≈31 and n1≈162 — but
     assigns equal 50/50 weight, while the transient distribution is 100% at
     (31,162) initially. TV≈0.8 — much better than Binomial but not perfect.
  4. [Key insight] Perfect prolongation requires time-dependent π: the coarse
     variable nc=n1+n2 loses the directionality of the front. For K=4→K=2
     multigrid, spatial structure is preserved and dynamic-π per pair is exact.

Run: julia --project examples/schlogl_msa.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D   = 0.005  # diffusion coefficient — slow enough that front is metastable in stationary
K   = 2      # fine voxels
h   = 1.0/K  # voxel size
d   = D/h^2  # diffusion jump rate

model = SchloglModel1D(D)   # default: c1=0.028, c2=3.2e-4, c3=19.5, c4=1.0
rates = Float64[]            # all rates baked into closures
n_max = 200                  # per-voxel FSP truncation

fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)   # K2 = 1 voxel

# ─── find stable states ───────────────────────────────────────────────────────

println("="^60)
println("Schlögl MSA Benchmark  K=$K voxels  D=$D  d=$(round(d,digits=3))")
println("="^60)

fp = schlogl_fixed_points(model; n_max=n_max)
n_low, n_uns, n_high = fp
@printf("\nSingle-voxel fixed points:\n")
@printf("  n_low  ≈ %d  (stable)\n", n_low)
@printf("  n_uns  ≈ %d  (unstable)\n", n_uns)
@printf("  n_high ≈ %d  (stable)\n", n_high)

# ─── direct FSP (ground truth) ───────────────────────────────────────────────

println("\n── Direct FSP (ground truth, $(n_max+1)^2 = $((n_max+1)^2) states) ──")

# Build full state space: all (n1,n2) with 0 ≤ n1,n2 ≤ n_max
sp_fine = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_fine, CartesianIndex(n1, n2), 0.0)
end

# Initial condition: low–high front
u0 = CartesianIndex(n_low, n_high)
sp_fine.probs[sp_fine.index[u0]] = 1.0

fine_sys = build_schlogl_rdme_system(model, fine_grid)
t0 = time()
A_fine, = build_generator(sp_fine, fine_sys, rates, 0.0)
t_gen = time() - t0
@printf("  Generator build: %.2fs  (|S|=%d, nnz=%d)\n",
        t_gen, length(sp_fine), nnz(A_fine))

# Evolve to several times
t_solve = [0.5, 1.0, 2.0, 5.0]
probs_fine = Dict{Float64, Vector{Float64}}()
t0 = time()
let p = copy(sp_fine.probs), t_prev = 0.0
    for t in t_solve
        p = max.(expv(t - t_prev, A_fine, p; m=50), 0.0)
        p ./= sum(p)
        probs_fine[t] = copy(p)
        t_prev = t
    end
end
t_fine_total = time() - t0
@printf("  Solve to t=%.1f: %.2fs\n", maximum(t_solve), t_fine_total)

# Report marginal means
for t in t_solve
    p = probs_fine[t]
    mu1 = sum(Tuple(sp_fine.states[i])[1] * p[i] for i in eachindex(p))
    mu2 = sum(Tuple(sp_fine.states[i])[2] * p[i] for i in eachindex(p))
    @printf("  t=%4.1f  ⟨n1⟩=%.1f  ⟨n2⟩=%.1f  Σ=%.3f\n", t, mu1, mu2, sum(p))
end

# ─── Coarse MSA (Binomial π) ─────────────────────────────────────────────────

println("\n── Coarse MSA (K2=1 voxel, Binomial-π prolongation) ──")

coarse_sys = build_schlogl_coarse_system(model, coarse_grid, fine_grid)

# Initial coarse state: nc = n_low + n_high
nc0 = n_low + n_high
sp_coarse = StateSpace{CartesianIndex{1}, Float64}()
add_state!(sp_coarse, CartesianIndex(nc0), 1.0)

# Expand coarse state space
coarse_bc = s -> rdme_bc(s, 2*n_max)
expand!(sp_coarse, coarse_sys, coarse_bc; depth=1)

# Evolve the coarse system with fixed dt steps
dt_c    = 0.05
n_steps = round(Int, maximum(t_solve) / dt_c)
probs_coarse = Dict{Float64, Tuple{Vector{CartesianIndex{1}}, Vector{Float64}}}()
t_cur = 0.0
t0 = time()
let t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, maximum(t_solve) - t_cur)
        expand!(sp_coarse, coarse_sys, coarse_bc; depth=1)
        A_c, = build_generator(sp_coarse, coarse_sys, rates, t_cur)
        sp_coarse.probs .= expv(dt_step, A_c, sp_coarse.probs; m=30)
        t_cur += dt_step
        renormalize!(sp_coarse)
        prune_threshold!(sp_coarse, 1e-12)

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            if !haskey(probs_coarse, t_snap)
                probs_coarse[t_snap] = (copy(sp_coarse.states), copy(sp_coarse.probs))
            end
        end
    end
end
t_coarse_total = time() - t0
@printf("  Solve to t=%.1f: %.2fs  (max|S_c|=up to %d coarse states)\n",
        maximum(t_solve), t_coarse_total, 2*n_max+1)

# Report coarse marginal means
for t in t_solve
    if haskey(probs_coarse, t)
        states_c, probs_c = probs_coarse[t]
        mu_c = sum(Tuple(states_c[i])[1] * probs_c[i] for i in eachindex(probs_c))
        # Split: ⟨n1⟩ = ⟨n2⟩ = ⟨nc⟩/2 under Binomial(nc, 1/2)
        @printf("  t=%4.1f  ⟨nc⟩=%.1f  ⟨n1⟩=⟨n2⟩=%.1f (Binomial split)  |S_c|=%d\n",
                t, mu_c, mu_c/2, length(states_c))
    end
end

# ─── Comparison: conditional distribution at the front ───────────────────────

println("\n── Conditional distribution: p(n1 | nc = $(n_low+n_high)) at t=0.5 ──")
println("  (nc = n1+n2 = n_low+n_high = $(n_low+n_high) — the spatial front)")

t_diag = 0.5
if haskey(probs_fine, t_diag)
    p_fine = probs_fine[t_diag]
    nc_target = n_low + n_high

    # Collect p(n1, n2) for all states with n1+n2 = nc_target
    cond = Dict{Int,Float64}()
    for (i, s) in enumerate(sp_fine.states)
        t = Tuple(s)
        if t[1] + t[2] == nc_target
            cond[t[1]] = get(cond, t[1], 0.0) + p_fine[i]
        end
    end
    mass_nc = sum(values(cond))

    if mass_nc > 1e-10
        # Normalize conditional
        for k in keys(cond); cond[k] /= mass_nc; end

        # Binomial prediction (log-space to avoid overflow for large nc)
        binom_pred = Dict{Int,Float64}()
        log_half_nc = -nc_target * log(2.0)
        for n1 in 0:nc_target
            log_bc = sum(log(i) for i in (nc_target-n1+1):nc_target; init=0.0) -
                     sum(log(i) for i in 1:n1; init=0.0)
            binom_pred[n1] = exp(log_half_nc + log_bc)
        end

        # Find the two modes of the true conditional
        mode1_n1 = argmax(cond)
        println("  Total mass at nc=$nc_target: $(round(mass_nc, sigdigits=3))")
        println()
        println("  n1   true p(n1|nc)   Binomial(nc,1/2)")
        println("  ─────────────────────────────────────")

        # Show three windows: around n_low, nc/2 (Binomial mode), and n_high
        nc_half = nc_target ÷ 2
        windows = vcat(max(0,n_low-5):min(nc_target,n_low+5),
                       max(0,nc_half-3):min(nc_target,nc_half+3),
                       max(0,n_high-5):min(nc_target,n_high+5))
        shown = Set{Int}()
        for n1 in sort(unique(windows))
            n1 in shown && continue; push!(shown, n1)
            true_p  = get(cond, n1, 0.0)
            binom_p = get(binom_pred, n1, 0.0)
            if true_p > 1e-6 || binom_p > 1e-4
                marker = true_p > 0.05 ? " ◀ mode" : ""
                @printf("  %3d   %.4f             %.4f%s\n",
                        n1, true_p, binom_p, marker)
            end
        end
        println()
        println("  Binomial puts all mass at n1=n2=$(nc_target÷2).")
        println("  True conditional has mass near n1≈$(n_low) (system started at front (n_low,n_high)).")
        println("  → Binomial-π prolongation gives the WRONG fine distribution at the front.")

        # TV distance between true conditional and Binomial
        all_n1 = union(keys(cond), keys(binom_pred))
        tv = 0.5 * sum(abs(get(cond, n1, 0.0) - get(binom_pred, n1, 0.0)) for n1 in all_n1)
        @printf("\n  TV(true conditional, Binomial) = %.4f\n", tv)
        println("  (TV=0 would mean Binomial π is exact; TV→1 means completely wrong.)")
    end
end

# ─── Dynamic-π: stationary solve ─────────────────────────────────────────────

println("\n── Dynamic-π: stationary conditional (solving Aᵀπ=0) ──")

t0 = time()
pi_table = compute_dynamic_pi(sp_fine, A_fine; n_max=n_max)
t_pi = time() - t0
@printf("  Stationary solve + slice extraction: %.2fs\n", t_pi)

# Show the dynamic-π conditional at nc = n_low+n_high
nc_front = n_low + n_high
pi_front = pi_table[nc_front+1]
@printf("  Dynamic-π at nc=%d (n_low+n_high):\n", nc_front)
println("    n1   dynamic-π    Binomial(nc,1/2)")
println("    ─────────────────────────────────")
log_half_nc_front = -nc_front * log(2.0)
for n1 in sort(vcat(max(0,n_low-3):min(nc_front,n_low+3),
                    max(0,nc_front÷2-2):min(nc_front,nc_front÷2+2),
                    max(0,n_high-3):min(nc_front,n_high+3)))
    dyn_p  = n1 < length(pi_front) ? pi_front[n1+1] : 0.0
    log_bc = sum(log(i) for i in (nc_front-n1+1):nc_front; init=0.0) -
             sum(log(i) for i in 1:n1; init=0.0)
    bin_p  = exp(log_half_nc_front + log_bc)
    if dyn_p > 1e-4 || bin_p > 1e-4
        marker = dyn_p > 0.05 ? " ◀" : ""
        @printf("    %3d   %.4f        %.4f%s\n", n1, dyn_p, bin_p, marker)
    end
end
println()
println("  Dynamic-π is bimodal at n1≈n_low and n1≈n_high;")
println("  Binomial concentrates at n1=nc÷2 (completely wrong).")

# ─── Restrict-Prolong accuracy test ──────────────────────────────────────────

println("\n── Restrict-Prolong accuracy test at t=0.5 ──")
println("  Restrict fine → coarse → prolong with Binomial-π vs Dynamic-π")
println("  Compare each to the true fine distribution (TV distance).")

t_rp = 0.5
p_fine_rp = probs_fine[t_rp]

# Build coarse distribution by restriction (exact, analytic)
coarse_dict = Dict{Int, Float64}()
for (i, s) in enumerate(sp_fine.states)
    nc_s = Tuple(s)[1] + Tuple(s)[2]
    coarse_dict[nc_s] = get(coarse_dict, nc_s, 0.0) + p_fine_rp[i]
end
sp_rp_coarse = StateSpace{CartesianIndex{1}, Float64}()
for (nc_s, p) in coarse_dict
    add_state!(sp_rp_coarse, CartesianIndex(nc_s), p)
end
renormalize!(sp_rp_coarse)

# Binomial prolongation
sp_binom_fine = prolong(sp_rp_coarse, Val(2); binom_tol=1e-6)
renormalize!(sp_binom_fine)

# Dynamic-π prolongation
sp_dyn_fine = prolong_dynamic(sp_rp_coarse, pi_table, n_max)
renormalize!(sp_dyn_fine)

# TV distances relative to direct FSP
function tv_vs_fine(sp_approx, sp_ref, p_ref)
    tv = 0.0
    for (i, s) in enumerate(sp_ref.states)
        p_true = p_ref[i]
        idx = get(sp_approx.index, s, 0)
        p_approx = idx != 0 ? sp_approx.probs[idx] : 0.0
        tv += abs(p_true - p_approx)
    end
    # States in approx but not in ref
    for (i, s) in enumerate(sp_approx.states)
        get(sp_ref.index, s, 0) != 0 && continue
        tv += sp_approx.probs[i]
    end
    tv / 2
end

tv_binom = tv_vs_fine(sp_binom_fine, sp_fine, p_fine_rp)
tv_dyn   = tv_vs_fine(sp_dyn_fine,   sp_fine, p_fine_rp)

@printf("  TV(Binomial-π prolongation, true fine)  = %.4f\n", tv_binom)
@printf("  TV(Dynamic-π  prolongation, true fine)  = %.4f\n", tv_dyn)
@printf("  Improvement factor: %.1f×\n", tv_binom / max(tv_dyn, 1e-10))

# ─── Dynamic-π Coarse MSA ────────────────────────────────────────────────────

println("\n── Dynamic-π Coarse MSA (K2=1 voxel, dynamic-π prolongation) ──")

coarse_sys_dyn = build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)

sp_coarse_dyn = StateSpace{CartesianIndex{1}, Float64}()
add_state!(sp_coarse_dyn, CartesianIndex(nc0), 1.0)

coarse_bc_dyn = s -> rdme_bc(s, 2*n_max)
expand!(sp_coarse_dyn, coarse_sys_dyn, coarse_bc_dyn; depth=1)

probs_coarse_dyn = Dict{Float64, Tuple{Vector{CartesianIndex{1}}, Vector{Float64}}}()
t0 = time()
let t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt_c, maximum(t_solve) - t_cur)
        expand!(sp_coarse_dyn, coarse_sys_dyn, coarse_bc_dyn; depth=1)
        A_cd, = build_generator(sp_coarse_dyn, coarse_sys_dyn, rates, t_cur)
        sp_coarse_dyn.probs .= expv(dt_step, A_cd, sp_coarse_dyn.probs; m=30)
        t_cur += dt_step
        renormalize!(sp_coarse_dyn)
        prune_threshold!(sp_coarse_dyn, 1e-12)

        if any(abs(t_cur - t) < dt_c/2 for t in t_solve)
            t_snap = t_solve[argmin(abs.(t_solve .- t_cur))]
            if !haskey(probs_coarse_dyn, t_snap)
                probs_coarse_dyn[t_snap] = (copy(sp_coarse_dyn.states), copy(sp_coarse_dyn.probs))
            end
        end
    end
end
t_coarse_dyn_total = time() - t0
@printf("  Solve to t=%.1f: %.2fs\n", maximum(t_solve), t_coarse_dyn_total)

# Report: prolong with dynamic-π and compute voxel means + TV vs direct FSP
println()
println("  t     ⟨n1⟩_dyn  ⟨n2⟩_dyn  TV(dyn-π, FSP)  |  ⟨n1⟩_binom  TV(binom, FSP)")
println("  ─────────────────────────────────────────────────────────────────────────")
for t in t_solve
    if haskey(probs_coarse_dyn, t) && haskey(probs_fine, t)
        states_cd, probs_cd = probs_coarse_dyn[t]
        p_true = probs_fine[t]

        # Build temporary coarse StateSpace for prolongation
        sp_tmp = StateSpace{CartesianIndex{1}, Float64}()
        for (i, s) in enumerate(states_cd)
            add_state!(sp_tmp, s, probs_cd[i])
        end

        # Dynamic-π prolongation → fine means + TV
        sp_prol_dyn = prolong_dynamic(sp_tmp, pi_table, n_max)
        renormalize!(sp_prol_dyn)
        mu1_dyn = sum(Tuple(sp_prol_dyn.states[i])[1] * sp_prol_dyn.probs[i]
                      for i in eachindex(sp_prol_dyn.probs))
        mu2_dyn = sum(Tuple(sp_prol_dyn.states[i])[2] * sp_prol_dyn.probs[i]
                      for i in eachindex(sp_prol_dyn.probs))
        tv_dyn_t = tv_vs_fine(sp_prol_dyn, sp_fine, p_true)

        # Binomial prolongation → TV for comparison
        sp_prol_binom = prolong(sp_tmp, Val(2); binom_tol=1e-6)
        renormalize!(sp_prol_binom)
        mu1_binom = sum(Tuple(sp_prol_binom.states[i])[1] * sp_prol_binom.probs[i]
                        for i in eachindex(sp_prol_binom.probs))
        tv_binom_t = tv_vs_fine(sp_prol_binom, sp_fine, p_true)

        @printf("  %4.1f  %7.1f   %7.1f   %.4f          |  %7.1f     %.4f\n",
                t, mu1_dyn, mu2_dyn, tv_dyn_t, mu1_binom, tv_binom_t)
    end
end

# ─── Summary ─────────────────────────────────────────────────────────────────

println("\n── Summary ──────────────────────────────────────────────")
println("  System: K=2 Schlögl RDME, front IC at (n_low, n_high)")
@printf("  Direct FSP:            %.2fs  (%d states)\n", t_fine_total, (n_max+1)^2)
@printf("  Stationary π solve:    %.2fs  (sparse LU, same state space)\n", t_pi)
@printf("  Coarse MSA Binomial-π: %.2fs  (≤%d coarse states)\n", t_coarse_total, 2*n_max+1)
@printf("  Coarse MSA Dynamic-π:  %.2fs  (≤%d coarse states)\n", t_coarse_dyn_total, 2*n_max+1)
println()
println("  Binomial-π: unimodal at n1=nc/2 — completely wrong (TV=1.0).")
println("  Dynamic-π:  bimodal at n1≈n_low and n1≈n_high — correct structure (TV≈0.8).")
println("  Residual error: stationary dynamic-π assigns 50/50 to both modes because")
println("  nc = n1+n2 loses directionality — it cannot distinguish (n_low,n_high)")
println("  from (n_high,n_low). Perfect prolongation requires time-dependent π.")
println("  For K=4→K=2 multigrid, each coarse voxel keeps spatial identity,")
println("  so dynamic-π per pair naturally captures the front.")
println()
println("="^60)
println("DONE")
println("="^60)
