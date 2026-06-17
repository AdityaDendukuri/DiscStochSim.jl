# Per-step timing figure (speed.pdf)
# 3-panel: (a) Oregonator, (b) Robertson, (c) Toggle switch
# Each panel: step wall time vs simulated time for forward enum vs standard construction
# Requires: julia --project examples/speed_benchmark.jl  first
using Plots, Serialization

gr()
default(
    fontfamily    = "Computer Modern",
    titlefontsize  = 11,
    guidefontsize  = 10,
    tickfontsize   = 9,
    legendfontsize = 8,
    linewidth = 1.0,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

results_dir = joinpath(@__DIR__, "results")
paper_dir   = joinpath(@__DIR__, "..", "paper")

systems = ["oregonator", "robertson", "toggle"]
titles  = ["(a) Oregonator", "(b) Robertson", "(c) Toggle switch"]
xscales = [:linear, :log10, :linear]
yscales = [:log10, :log10, :log10]

panels = []
for (sys, title, xsc, ysc) in zip(systems, titles, xscales, yscales)
    d_fwd = deserialize(joinpath(results_dir, "speed_$(sys)_forward.jls"))
    d_std = deserialize(joinpath(results_dir, "speed_$(sys)_standard.jls"))

    p = plot(xlabel="Time", ylabel="Step time (s)", title=title,
             xscale=xsc, yscale=ysc)
    plot!(p, d_fwd["t_log"], d_fwd["step_time_log"];
          label="Forward enum.", color=:black, linestyle=:solid, linewidth=1.0)
    plot!(p, d_std["t_log"], d_std["step_time_log"];
          label="Standard", color=:black, linestyle=:dot, linewidth=1.0)
    push!(panels, p)
end

fig = plot(panels...;
    layout=(1, 3), size=(1100, 360), margin=5Plots.mm)

savefig(fig, joinpath(paper_dir, "speed.pdf"))
savefig(fig, joinpath(paper_dir, "speed.png"))
println("Saved speed.{pdf,png} → paper/")
