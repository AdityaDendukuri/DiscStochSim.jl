# Robertson simulation results figure (simulation_results.pdf)
# 4-panel: (a) mean trajectories, (b) |S| vs time,
#           (c) dt vs time (log-log), (d) flux vs time
# Requires: julia --project examples/robertson.jl   (~19 min)
using Plots, Serialization

gr()
default(
    fontfamily   = "Computer Modern",
    titlefontsize  = 11,
    guidefontsize  = 10,
    tickfontsize   = 9,
    legendfontsize = 8,
    linewidth = 1.5,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

d = deserialize(joinpath(@__DIR__, "results", "robertson_fpfsp.jls"))

snap_t    = d["snap_t"]
snap_traj = d["snap_traj"]   # rows = snapshots, cols = species
snap_sz   = d["snap_sizes"]
t_sub     = d["t_log"]
dt_sub    = d["dt_log"]
sz_sub    = d["size_log"]
fl_sub    = d["flux_log"]

# Only plot snapshot times > 0 (first entry is t=0)
t_snap_pos = snap_t[snap_t .> 0]
traj_pos   = snap_traj[snap_t .> 0, :]
sz_snap_pos = snap_sz[snap_t .> 0]

# Panel (a): mean trajectories
p1 = plot(xlabel = "Time (s)", ylabel = "Mean population",
          title = "(a) Robertson: mean trajectories",
          xscale = :log10)
plot!(p1, t_snap_pos, traj_pos[:, 1],
      label = "⟨A⟩", color = :black, linestyle = :solid)
plot!(p1, t_snap_pos, traj_pos[:, 2],
      label = "⟨B⟩", color = :black, linestyle = :dash)
plot!(p1, t_snap_pos, traj_pos[:, 3],
      label = "⟨C⟩", color = :black, linestyle = :dot)

# Panel (b): state space size
p2 = plot(xlabel = "Time (s)", ylabel = "|S|",
          title = "(b) State space size",
          xscale = :log10)
plot!(p2, t_snap_pos, sz_snap_pos,
      label = "", color = :black, linestyle = :solid)

# Panel (c): adaptive time step (log-log)
p3 = plot(xlabel = "Time (s)", ylabel = "δt",
          title = "(c) Adaptive time step",
          xscale = :log10, yscale = :log10)
plot!(p3, t_sub, dt_sub,
      label = "", color = :black, linestyle = :solid)

# Panel (d): total flux
p4 = plot(xlabel = "Time (s)", ylabel = "Φ",
          title = "(d) Total flux",
          xscale = :log10, yscale = :log10)
plot!(p4, t_sub, fl_sub,
      label = "", color = :black, linestyle = :solid)

fig = plot(p1, p2, p3, p4; layout = (2, 2), size = (900, 600), margin = 5Plots.mm)

paper_dir = joinpath(@__DIR__, "..", "paper")
savefig(fig, joinpath(paper_dir, "simulation_results.pdf"))
savefig(fig, joinpath(paper_dir, "simulation_results.png"))
println("Saved simulation_results.{pdf,png} → paper/")

println("total_iters=$(d["total_iters"]), wall=$(round(d["wall_time"];digits=0))s")
println("dt: [$(round(minimum(dt_sub);sigdigits=2)), $(round(maximum(dt_sub);sigdigits=2))]")
println("⟨A⟩=$(round(d["mean_A"];digits=0)), ⟨C⟩=$(round(d["mean_C"];digits=0)), |S|_final=$(d["S_final"])")
