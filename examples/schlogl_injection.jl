"""
Schlögl Injection Prolongation Test

Tests whether injection prolongation (Fix 1) resolves the symmetry trap
that defeats Binomial-π and Dynamic-π correction prolongation.

Setup: K=2, IC = (n_low, n_high) — the spatial front.  The coarse variable
nc = n1+n2 = 193 is symmetric: it cannot distinguish (31,162) from (162,31).

With correction prolongation:
  Binomial-π: puts all mass at (96,97) — completely wrong.
  Dynamic-π:  bimodal at (31,162) and (162,31) with 50/50 weight — TV≈0.8.

With injection prolongation:
  p^h_new(x) = p̃^h(x) · p̄^{2h}(x̄) / p^{2h}(x̄)
  The fine distribution is just rescaled proportionally to the coarse update,
  preserving the asymmetry: (31,162) stays as the dominant mode.

Run: julia --project examples/schlogl_injection.jl
"""

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D    = 0.005
K    = 2
h    = 1.0 / K
model     = SchloglModel1D(D)
rates     = Float64[]
n_max     = 200
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model; n_max=n_max)
n_low, n_uns, n_high = fp

println("="^60)
println("Schlögl Injection Prolongation Test  K=$K  D=$D")
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("IC: (n_low=$n_low, n_high=$n_high) — spatial front")
println("="^60)

# ─── ground truth: direct FSP ────────────────────────────────────────────────

println("\n── Ground truth: direct FSP ──")
sp_fine = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_fine, CartesianIndex(n1, n2), 0.0)
end
sp_fine.probs[sp_fine.index[CartesianIndex(n_low, n_high)]] = 1.0

fine_sys = build_schlogl_rdme_system(model, fine_grid)
A_fine, = build_generator(sp_fine, fine_sys, rates, 0.0)

t_eval = 0.5
p_true = expv(t_eval, A_fine, sp_fine.probs; m=50)
p_true = max.(p_true, 0.0); p_true ./= sum(p_true)

mu1_true = sum(Tuple(sp_fine.states[i])[1] * p_true[i] for i in eachindex(p_true))
mu2_true = sum(Tuple(sp_fine.states[i])[2] * p_true[i] for i in eachindex(p_true))
@printf("  t=%.1f: ⟨n1⟩=%.1f  ⟨n2⟩=%.1f\n", t_eval, mu1_true, mu2_true)

# ─── compute dynamic-π (for comparison methods) ───────────────────────────────

println("\n── Computing dynamic-π (stationary conditional) ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_fine, A_fine; n_max=n_max)
@printf("  Done in %.2fs\n", time() - t0)

# ─── TV helper ────────────────────────────────────────────────────────────────

function tv_vs_true(sp_approx::StateSpace, probs_approx, sp_ref, p_ref)
    tv = 0.0
    for (i, s) in enumerate(sp_ref.states)
        p_t = p_ref[i]
        idx = get(sp_approx.index, s, 0)
        p_a = idx != 0 ? probs_approx[idx] : 0.0
        tv += abs(p_t - p_a)
    end
    for (i, s) in enumerate(sp_approx.states)
        get(sp_ref.index, s, 0) != 0 && continue
        tv += abs(probs_approx[i])
    end
    tv / 2
end

# ─── multi-step V-cycle: compare all methods ────────────────────────────────

t_steps = [0.5, 1.0, 2.0, 5.0]
dt_step  = 0.1    # V-cycle step size

println("\n── Multi-step V-cycle (dt=$(dt_step) per step) ──\n")

cnmax = 2 * n_max   # coarse voxel holds sum of two fine voxels

function run_vcycle_steps(mode::Symbol, t_end, dt)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    add_state!(sp, CartesianIndex(n_low, n_high), 1.0)
    t_cur = 0.0
    n = round(Int, t_end / dt)
    for _ in 1:n
        sp = if mode == :binom
            two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid,
                                     pi_table, rates, t_cur, dt;
                                     τ_pre=0.0, τ_post=0.0, use_dynamic_pi=false,
                                     coarse_n_max=cnmax)
        elseif mode == :dynpi
            two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid,
                                     pi_table, rates, t_cur, dt;
                                     τ_pre=0.0, τ_post=0.0, use_dynamic_pi=true,
                                     coarse_n_max=cnmax)
        else  # :fix3 — fine-conditional transfer (injection + conditional carry-over)
            two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                               pi_table, rates, t_cur, dt;
                                               τ_pre=0.0, τ_post=0.0, use_dynamic_pi=true,
                                               coarse_n_max=cnmax)
        end
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
        t_cur += dt
    end
    sp
end

# Run FSP ground truth to each time (already have p_true at t=0.5)
fsp_probs = Dict{Float64, Vector{Float64}}()
let p = copy(sp_fine.probs), t_prev = 0.0
    for t in t_steps
        p = max.(expv(t - t_prev, A_fine, p; m=50), 0.0)
        p ./= sum(p)
        fsp_probs[t] = copy(p)
        t_prev = t
    end
end

function means_and_tv(sp_approx, p_true_vec)
    mu1 = sum(Tuple(sp_approx.states[i])[1] * sp_approx.probs[i] for i in eachindex(sp_approx.probs))
    mu2 = sum(Tuple(sp_approx.states[i])[2] * sp_approx.probs[i] for i in eachindex(sp_approx.probs))
    tv  = tv_vs_true(sp_approx, sp_approx.probs, sp_fine, p_true_vec)
    mu1, mu2, tv
end

println("  t     Method                    ⟨n1⟩   ⟨n2⟩   TV vs FSP   |S_fine|")
println("  ─────────────────────────────────────────────────────────────")

for t in t_steps
    p_true_t = fsp_probs[t]
    mu1_gt = sum(Tuple(sp_fine.states[i])[1] * p_true_t[i] for i in eachindex(p_true_t))
    mu2_gt = sum(Tuple(sp_fine.states[i])[2] * p_true_t[i] for i in eachindex(p_true_t))
    @printf("  %4.1f  %-22s %6.1f %6.1f   —\n", t, "Ground truth (FSP)", mu1_gt, mu2_gt)

    sp_b = run_vcycle_steps(:binom,    t, dt_step)
    sp_d = run_vcycle_steps(:dynpi,    t, dt_step)
    sp_i = run_vcycle_steps(:fix3,     t, dt_step)

    m1b, m2b, tvb = means_and_tv(sp_b, p_true_t)
    m1d, m2d, tvd = means_and_tv(sp_d, p_true_t)
    m1i, m2i, tvi = means_and_tv(sp_i, p_true_t)

    @printf("  %4s  %-22s %6.1f %6.1f   %.4f     %d\n",
            "", "Binomial-π corr.", m1b, m2b, tvb, length(sp_b))
    @printf("  %4s  %-22s %6.1f %6.1f   %.4f     %d\n",
            "", "Dynamic-π corr.", m1d, m2d, tvd, length(sp_d))
    @printf("  %4s  %-22s %6.1f %6.1f   %.4f     %d\n\n",
            "", "Fine-conditional (Fix 3)", m1i, m2i, tvi, length(sp_i))
end

println("="^60)
println("DONE")
println("="^60)
