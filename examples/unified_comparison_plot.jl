# Unified FP-FSP vs Krylov-FSP-SSA comparison figure
# Three failure modes in one 3-panel figure:
#   (a) Oregonator   — state-space explosion
#   (b) Bottleneck   — wavefront truncation
#   (c) Robertson    — stiffness-induced τ-lock
#
# Oregonator and Robertson results loaded from pre-computed .jls files.
# Bottleneck run inline (~2 s).
#
# Pre-requisites (run once each):
#   julia --project examples/oregonator_fpfsp.jl      → results/oregonator_fpfsp.jls
#   julia --project examples/oregonator_krylov.jl     → results/oregonator_krylov.jls
#   julia --project examples/robertson.jl             → results/robertson_fpfsp.jls
#   julia --project examples/robertson_krylov.jl      → results/robertson_krylov.jls
#
# Run: julia --project examples/unified_comparison_plot.jl
# Expected wall time: < 30 s (Bottleneck only)

using DiscStochSim, Catalyst, Serialization
using Expokit: expmv
using Plots

gr()
default(
    fontfamily  = "Computer Modern",
    titlefontsize  = 11,
    guidefontsize  = 10,
    tickfontsize   = 9,
    legendfontsize = 8,
    linewidth = 2,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

const RESDIR  = joinpath(@__DIR__, "results")
const OUTDIR  = joinpath(@__DIR__, "output", "comparison")
mkpath(OUTDIR)

# ─── Load pre-computed Oregonator results ───────────────────────────────────
println("Loading Oregonator results...")
oreg_fp = deserialize(joinpath(RESDIR, "oregonator_fpfsp.jls"))
oreg_k  = deserialize(joinpath(RESDIR, "oregonator_krylov.jls"))

# FP-FSP: per-iteration size_log + t_log (thin to ≤2000 pts)
sz_fp  = oreg_fp["size_log"];  t_fp  = oreg_fp["t_log"]
stride = max(1, length(sz_fp) ÷ 2000)
idx    = 1:stride:length(sz_fp)

# Krylov: per-iteration size_log + time_log
sz_k  = oreg_k["size_log"];  t_k  = oreg_k["time_log"]
println("  FP-FSP: $(length(sz_fp)) iters, max|S|=$(maximum(sz_fp))")
println("  KrylovFSP: $(length(sz_k)) iters, t_reached=$(oreg_k["t_reached"]), max|S|=$(oreg_k["S_max"])")

# ─── Bottleneck: run both methods inline (~2 s) ──────────────────────────────
println("Running Bottleneck...")

rn_bottle = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

k1_b, k2_b = 1e-6, 1e-1
prob_bottle = FSPProblem(rn_bottle, CartesianIndex(1, 0, 0), (0.0, 1e4),
                         [k1_b, k2_b]; bounds = (0, 1500))
E_C_exact = t -> (k2_b / k1_b) * (k1_b * t - 1 + exp(-k1_b * t))

# FP-FSP: fixed δt=10 s step (adaptive dt would give dt~1e6 >> tf at t=0)
function _fpfsp_fixed_step(prob::FSPProblem{E,T};
                            δt = 10.0, depth = 1,
                            quantile = 0.9, flux_tol = 1e-6,
                            save_every = 100) where {E,T}
    model = prob.model; rates = prob.rates; bc = prob.bc
    t0, tf = prob.tspan; ns = length(Tuple(prob.u0))
    sp = StateSpace{E, Float64}(); add_state!(sp, prob.u0, 1.0)
    sol_t = [t0]; snaps = [(copy(sp.states), copy(sp.probs))]; sizes = [length(sp)]
    t = Float64(t0); iter = 0
    while t < tf
        iter += 1
        expand!(sp, model, bc; depth)
        A, _, _, _ = build_generator(sp, model, rates, t)
        dt = min(δt, tf - t)
        sp.probs = expmv(dt, A, sp.probs)
        compress!(sp, model, rates, t, quantile; flux_tolerance = flux_tol)
        renormalize!(sp); t += dt
        if iter % save_every == 0 || t >= tf
            push!(sol_t, t); push!(snaps, (copy(sp.states), copy(sp.probs))); push!(sizes, length(sp))
        end
    end
    FSPSolution{E}(sol_t, snaps, sizes, ns)
end

t0 = time()
sol_b_fp = _fpfsp_fixed_step(prob_bottle)
t_b_fp   = time() - t0
traj_b_fp = mean_trajectory(sol_b_fp)
println("  FP-FSP: $(round(t_b_fp; digits=2))s, ⟨C⟩=$(round(traj_b_fp[end,3]; sigdigits=4))")

using Random; Random.seed!(1)
t0 = time()
sol_b_k, _ = solve(prob_bottle, KrylovFSP(ε=1e-2, τ_init=100.0, r=1, save_interval=1, max_iter=100))
t_b_k    = time() - t0
traj_b_k = mean_trajectory(sol_b_k)
println("  KrylovFSP: $(round(t_b_k; digits=2))s, ⟨C⟩=$(round(traj_b_k[end,3]; sigdigits=4))")

# ─── Load pre-computed Robertson results ────────────────────────────────────
println("Loading Robertson results...")
rober_fp = deserialize(joinpath(RESDIR, "robertson_fpfsp.jls"))
rober_k  = deserialize(joinpath(RESDIR, "robertson_krylov.jls"))

dt_r_fp_full  = rober_fp["dt_log"];  t_r_fp_full  = rober_fp["t_log"]
acc_τ_k_all  = rober_k["τ_log"];   acc_t_k_all = rober_k["time_log"]

# Log-subsample FP-FSP: pick ~2000 points log-spaced in time
# (stride-based would massively oversample the burst region)
let t_full = t_r_fp_full, dt_full = dt_r_fp_full, n_save = 2000
    t_grid  = 10 .^ range(log10(t_full[1]), log10(t_full[end]); length=n_save)
    raw_idx = searchsortedfirst.(Ref(t_full), t_grid)
    log_idx = unique(clamp.(raw_idx, 1, length(t_full)))
    global dt_r_fp = dt_full[log_idx]
    global t_r_fp  = t_full[log_idx]
end

# Filter Krylov to valid (positive) time points for log-scale
let mask = acc_t_k_all .> 0
    global acc_τ_k = acc_τ_k_all[mask]
    global acc_t_k = acc_t_k_all[mask]
end

println("  FP-FSP:   $(length(dt_r_fp)) pts (subsampled), dt=[$(round(minimum(dt_r_fp);sigdigits=2)), $(round(maximum(dt_r_fp);sigdigits=2))]")
println("  KrylovFSP: $(length(acc_τ_k)) accepted steps, t_reached=$(rober_k["t_reached"])")

# ─── Build unified 3-panel figure ───────────────────────────────────────────
println("Building figure...")

# Panel (a): Oregonator — |S| vs simulated time
p1 = plot(xlabel = "Time", ylabel = "|S|", title = "(a) Oregonator")
plot!(p1, t_fp[idx], sz_fp[idx],
      label = "FP-FSP", color = :black, linestyle = :solid, linewidth = 2)
plot!(p1, t_k, sz_k,
      label = "Krylov-FSP-SSA", color = :black, linestyle = :dash, linewidth = 2)

# Panel (b): Bottleneck — ⟨C⟩ vs time
t_ref   = range(0.0, 1e4; length = 400)
t_stars = range(0.0, 1e4; length = 10)
p2 = plot(xlabel = "Time", ylabel = "\u27E8C\u27E9", title = "(b) Bottleneck")
plot!(p2, t_ref, E_C_exact.(t_ref),
      label = "", color = :black, linestyle = :dot, linewidth = 1.2)
scatter!(p2, t_stars, E_C_exact.(t_stars),
         label = "Analytical", color = :black,
         markershape = :star5, markersize = 9,
         markerstrokewidth = 0.8, markerstrokecolor = :black)
plot!(p2, sol_b_fp.t, traj_b_fp[:, 3],
      label = "FP-FSP", color = :black, linestyle = :solid, linewidth = 2)
plot!(p2, sol_b_k.t, traj_b_k[:, 3],
      label = "Krylov-FSP-SSA", color = :black, linestyle = :dash, linewidth = 2)

# Panel (c): Robertson — dt/τ vs simulated time (log scale)
# FP-FSP: explicitly stored per-step dt; Krylov: explicitly stored accepted τ
p3 = plot(xlabel = "Time", ylabel = "Step size", title = "(c) Robertson",
          xscale = :log10, yscale = :log10)
plot!(p3, t_r_fp, dt_r_fp,
      label = "FP-FSP dt", color = :black, linestyle = :solid, linewidth = 1.5)
plot!(p3, acc_t_k, acc_τ_k,
      label = "Krylov-FSP-SSA \u03C4", color = :black, linestyle = :dash, linewidth = 1.5)
# Mark where Krylov stops
vline!(p3, [acc_t_k[end]], color = :gray, linestyle = :dot, linewidth = 1, label = "")

fig = plot(p1, p2, p3; layout = (1, 3), size = (1200, 360), margin = 5Plots.mm)
savefig(fig, joinpath(OUTDIR, "unified_comparison.png"))
savefig(fig, joinpath(OUTDIR, "unified_comparison.pdf"))
println("Saved unified_comparison.{png,pdf}")

# Copy to paper/
paper_dir = joinpath(@__DIR__, "..", "paper")
if isdir(paper_dir)
    cp(joinpath(OUTDIR, "unified_comparison.pdf"),
       joinpath(paper_dir, "unified_comparison.pdf"); force = true)
    println("Copied → paper/unified_comparison.pdf")
end
println("Done.")
