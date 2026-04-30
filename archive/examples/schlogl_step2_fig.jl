"""
Generate figure: Schlögl K=4 asymmetric front — GT vs Binom-π vs Mult. prolong..
Tracks ⟨n₁⟩ and |S| at every time step.
Saves paper/figures/fig7_schlogl_front.{pdf,png}
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

# ─── parameters ───────────────────────────────────────────────────────────────

D     = 0.005
K     = 4
h     = 1.0 / K
n_max = 180

model       = SchloglModel1D(D)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp            = schlogl_fixed_points(model)
n_low, _, n_high = fp

u0     = CartesianIndex(n_low, n_high, n_high, n_high)
dt     = 0.2
t_max  = 2.0
nsteps = round(Int, t_max / dt)
cnmax  = 2 * n_max

# ─── dynamic-π table ──────────────────────────────────────────────────────────

println("Pre-computing dynamic-π (K=2 at h=$(h))...")
pair_grid = VoxelGrid(2, h, 0)
pair_sys  = build_schlogl_rdme_system(model, pair_grid)
sp_pair   = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_pair, CartesianIndex(n1, n2), 0.0)
end
sp_pair.probs[sp_pair.index[CartesianIndex(n_low, n_low)]] = 1.0
A_pair, = build_generator(sp_pair, pair_sys, rates, 0.0)
pi_table = compute_dynamic_pi(sp_pair, A_pair; n_max=n_max)
println("  Done.")

# ─── helper ───────────────────────────────────────────────────────────────────

function voxel_mean1(states, probs)
    mu = 0.0
    for (i, s) in enumerate(states)
        mu += probs[i] * Tuple(s)[1]
    end
    mu
end

# ─── ground truth ─────────────────────────────────────────────────────────────

println("Running GT (K=4 adaptive FSP)...")
fine_sys = build_schlogl_rdme_system(model, fine_grid)
fine_bc  = state -> rdme_bc(state, n_max)

sp_gt = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp_gt, u0, 1.0)

ts_gt   = Float64[]
mu1_gt  = Float64[]
sz_gt   = Int[]
t_gt = @elapsed let t_cur = 0.0
    for step in 1:nsteps
        dt_step = min(dt, t_max - t_cur)
        expand!(sp_gt, fine_sys, fine_bc; depth=1)
        A, = build_generator(sp_gt, fine_sys, rates, t_cur)
        sp_gt.probs .= max.(expv(dt_step, A, sp_gt.probs; m=50), 0.0)
        t_cur += dt_step
        renormalize!(sp_gt)
        prune_threshold!(sp_gt, 1e-14)
        push!(ts_gt,  t_cur)
        push!(mu1_gt, voxel_mean1(sp_gt.states, sp_gt.probs))
        push!(sz_gt,  length(sp_gt))
    end
end
@printf("  Done: %.1fs  final |S|=%d\n", t_gt, length(sp_gt))

# ─── V-cycle runner (function avoids soft-scope issues) ───────────────────────

function run_vcycle_tracked(mode::Symbol)
    sp     = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    ts     = Float64[]
    mu1s   = Float64[]
    szs    = Int[]
    t_cur  = 0.0
    for step in 1:nsteps
        dt_step = min(dt, t_max - t_cur)
        sp = if mode == :binom
            two_level_vcycle_schlogl(sp, model, fine_grid, coarse_grid,
                                     pi_table, rates, t_cur, dt_step;
                                     use_dynamic_pi=false, coarse_n_max=cnmax)
        else
            two_level_vcycle_schlogl_injection(sp, model, fine_grid, coarse_grid,
                                               pi_table, rates, t_cur, dt_step;
                                               coarse_n_max=cnmax)
        end
        t_cur += dt_step
        renormalize!(sp)
        prune_threshold!(sp, 1e-12)
        push!(ts,   t_cur)
        push!(mu1s, voxel_mean1(sp.states, sp.probs))
        push!(szs,  length(sp))
    end
    ts, mu1s, szs, length(sp)
end

println("Running V-cycle Binom-π (may take ~60s)...")
t_binom = @elapsed (ts_b, mu1_b, sz_b, sz_b_final) = run_vcycle_tracked(:binom)
@printf("  Done: %.1fs  final |S|=%d\n", t_binom, sz_b_final)

println("Running V-cycle Mult. prolong....")
t_mult  = @elapsed (ts_f, mu1_f, sz_f, sz_f_final) = run_vcycle_tracked(:multiplicative)
@printf("  Done: %.1fs  final |S|=%d\n", t_mult, sz_f_final)

# ─── plot ─────────────────────────────────────────────────────────────────────

fig = Figure(size=(800, 300), fontsize=13)

# left panel: ⟨n₁⟩ vs t
ax1 = Axis(fig[1, 1],
    xlabel = L"time $t$",
    ylabel = L"\langle n_1 \rangle \;(\text{front voxel})",
    title  = "Front voxel mean",
    yticks = [0, 31, 64, 96, 128, 162])

hlines!(ax1, [n_low],  linestyle=:dot, color=(:gray50, 0.8), linewidth=1.2)
hlines!(ax1, [n_high], linestyle=:dot, color=(:gray50, 0.8), linewidth=1.2)

lines!(ax1, ts_gt,  mu1_gt, color=:black,      linewidth=2.2, label="GT (K=4 FSP)")
lines!(ax1, ts_b,   mu1_b,  color=:firebrick,  linewidth=1.8,
       linestyle=:dash,  label="V-cyc Binom-π")
lines!(ax1, ts_f,   mu1_f,  color=:royalblue,  linewidth=2.2, label="V-cyc Mult. prolong.")

axislegend(ax1, position=:rt, framevisible=true, labelsize=11)

# right panel: |S| vs t (log scale)
ax2 = Axis(fig[1, 2],
    xlabel = L"time $t$",
    ylabel = L"|{\mathcal{S}}| \;(\text{fine states})",
    title  = "State space size",
    yscale = log10)

lines!(ax2, ts_gt, Float64.(sz_gt), color=:black,     linewidth=2.2)
lines!(ax2, ts_b,  Float64.(sz_b),  color=:firebrick,  linewidth=1.8, linestyle=:dash)
lines!(ax2, ts_f,  Float64.(sz_f),  color=:royalblue,  linewidth=2.2)

out_dir = joinpath(@__DIR__, "..", "paper", "figures")
save(joinpath(out_dir, "fig7_schlogl_front.pdf"), fig)
save(joinpath(out_dir, "fig7_schlogl_front.png"), fig, px_per_unit=2)
println("\nSaved fig7_schlogl_front.pdf / .png")
