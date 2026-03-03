# Oregonator comparison plot
# Loads results from oregonator_fpfsp.jls and oregonator_krylov.jls,
# produces oregonator_comparison.pdf in paper/.
#
# Run after both data scripts have completed:
#   julia --project examples/oregonator_fpfsp.jl    # ~40 s
#   julia --project examples/oregonator_krylov.jl   # ~2 min
#   julia --project examples/oregonator_plot.jl     # fast
using Plots, Serialization

gr()
default(
    fontfamily  = "Computer Modern",
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
fpfsp  = deserialize(joinpath(results_dir, "oregonator_fpfsp.jls"))
krylov = deserialize(joinpath(results_dir, "oregonator_krylov.jls"))

# ── Panel (a): |S| vs simulation time ─────────────────────────────────────
p1 = plot(title="(a) State space size",
          xlabel="Simulation time t", ylabel="|S|")

# FP-FSP: per-step size_log with t_log (downsample to ≤3000 pts)
dt_log   = fpfsp["dt_log"]
t_log    = fpfsp["t_log"]
size_log = fpfsp["size_log"]
stride = max(1, length(t_log) ÷ 3000)
idx = 1:stride:length(t_log)
plot!(p1, t_log[idx], size_log[idx];
    label="FP-FSP", color=:black, linestyle=:solid)

# KrylovFSP: size_log vs time_log
if !isempty(krylov["size_log"]) && !isempty(krylov["time_log"])
    plot!(p1, krylov["time_log"], krylov["size_log"];
        label="Krylov-FSP-SSA", color=:black, linestyle=:dash)
end

# ── Panel (b): time step vs simulation time (log scale) ──────────────────
p2 = plot(title="(b) Time step",
          xlabel="Simulation time t", ylabel="dt / \u03C4", yscale=:log10)

plot!(p2, t_log[idx], dt_log[idx];
    label="FP-FSP dt", color=:black, linestyle=:solid)

if !isempty(krylov["τ_log"]) && !isempty(krylov["time_log"])
    plot!(p2, krylov["time_log"], krylov["τ_log"];
        label="Krylov-FSP-SSA \u03C4", color=:black, linestyle=:dash)
end

fig = plot(p1, p2; layout=(1,2), size=(900, 350), margin=5Plots.mm)

# Print summary
tf = 0.11
println("FP-FSP:    $(fpfsp["total_iters"]) iters, " *
        "wall=$(round(fpfsp["wall_time"]; digits=1))s, " *
        "max|S|=$(fpfsp["S_max"]), final|S|=$(fpfsp["S_final"])")
println("           dt∈[$(round(minimum(dt_log[dt_log .> 0]); sigdigits=2)), " *
        "$(round(maximum(dt_log); sigdigits=2))]")
println("KrylovFSP: $(krylov["total_iters"]) iters, " *
        "$(krylov["reject_count"]) rejections, " *
        "wall=$(round(krylov["wall_time"]; digits=1))s, " *
        "t=$(round(krylov["t_reached"]; sigdigits=4))/$tf " *
        "($(round(100*krylov["t_reached"]/tf; digits=2))%), " *
        "max|S|=$(krylov["S_max"])")

out_dir = joinpath(@__DIR__, "output", "comparison")
mkpath(out_dir)
savefig(fig, joinpath(out_dir, "oregonator_comparison.png"))
savefig(fig, joinpath(out_dir, "oregonator_comparison.pdf"))
println("Saved oregonator_comparison.{png,pdf}")

paper_dir = joinpath(@__DIR__, "..", "paper")
if isdir(paper_dir)
    cp(joinpath(out_dir, "oregonator_comparison.pdf"),
       joinpath(paper_dir, "oregonator_comparison.pdf"); force=true)
    println("Copied → paper/oregonator_comparison.pdf")
end
