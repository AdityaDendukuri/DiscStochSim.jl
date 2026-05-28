# Robertson comparison plot
# Loads results saved by robertson.jl and robertson_krylov.jl,
# produces robertson_comparison.pdf in paper/.
#
# Run after both data scripts have completed:
#   julia --project examples/robertson.jl        # ~9 min
#   julia --project examples/robertson_krylov.jl # ~1 min
#   julia --project examples/robertson_plot.jl   # fast
using Plots, Serialization

gr()
default(
    fontfamily = "Computer Modern",
    titlefontsize  = 13,
    guidefontsize  = 11,
    tickfontsize   = 9,
    legendfontsize = 9,
    linewidth = 1.5,
    framestyle = :box,
    grid = false,
    dpi  = 300,
)

results_dir = joinpath(@__DIR__, "results")
fpfsp  = deserialize(joinpath(results_dir, "robertson_fpfsp.jls"))
krylov = deserialize(joinpath(results_dir, "robertson_krylov.jls"))

# ── Panel (a): FP-FSP dt vs simulation time t ─────────────────────────────
dt_log = fpfsp["dt_log"]
t_log  = fpfsp["t_log"]   # same length as dt_log

# Downsample for plotting (keep ≤ 3000 points)
stride = max(1, length(dt_log) ÷ 3000)
idx = 1:stride:length(dt_log)

p1 = plot(t_log[idx], dt_log[idx];
    title  = "(a) FP-FSP: time step vs. time",
    xlabel = "Simulation time t",
    ylabel = "\u03B4t",
    yscale = :log10,
    color  = :black,
    linestyle = :solid,
    label  = "FP-FSP \u03B4t",
    legend = :topleft,
)

# ── Panel (b): KrylovFSP τ vs iteration ───────────────────────────────────
τ_log = krylov["τ_log"]
n_τ   = length(τ_log)

p2 = plot(1:n_τ, τ_log;
    title  = "(b) Krylov-FSP-SSA: time step vs. iteration",
    xlabel = "Iteration",
    ylabel = "\u03C4",
    yscale = :log10,
    color  = :black,
    linestyle = :dash,
    label  = "Krylov-FSP-SSA \u03C4",
    legend = :topright,
)

fig = plot(p1, p2; layout=(1,2), size=(900, 350), margin=5Plots.mm)

# Print summary
tf = 1e4
println("FP-FSP:     $(fpfsp["total_iters"]) iters, " *
        "wall=$(round(fpfsp["wall_time"]; digits=1))s, " *
        "dt∈[$(round(minimum(dt_log); sigdigits=2)), $(round(maximum(dt_log); sigdigits=2))], " *
        "⟨A⟩=$(round(fpfsp["mean_A"]; digits=1)), ⟨C⟩=$(round(fpfsp["mean_C"]; digits=1)), " *
        "|S|=$(fpfsp["S_final"])")
println("KrylovFSP:  $(krylov["total_iters"]) iters, " *
        "$(krylov["reject_count"]) rejections, " *
        "wall=$(round(krylov["wall_time"]; digits=2))s, " *
        "t=$(round(krylov["t_reached"]; sigdigits=4))/$(tf) " *
        "($(round(100*krylov["t_reached"]/tf; digits=3))%)")
if n_τ > 0
    println("  τ∈[$(round(minimum(τ_log); sigdigits=2)), $(round(maximum(τ_log); sigdigits=2))]")
end

# Save
out_dir = joinpath(@__DIR__, "output", "comparison")
mkpath(out_dir)
savefig(fig, joinpath(out_dir, "robertson_comparison.png"))
savefig(fig, joinpath(out_dir, "robertson_comparison.pdf"))
println("Saved robertson_comparison.{png,pdf}")

paper_dir = joinpath(@__DIR__, "..", "paper")
if isdir(paper_dir)
    cp(joinpath(out_dir, "robertson_comparison.pdf"),
       joinpath(paper_dir, "robertson_comparison.pdf"); force=true)
    println("Copied → paper/robertson_comparison.pdf")
end
