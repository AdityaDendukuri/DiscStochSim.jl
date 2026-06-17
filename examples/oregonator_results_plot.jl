# Oregonator simulation results figure (orag.pdf)
# 6-panel: (a) dt, (b) |S|, (c) mean trajectories X/Y/Z, (d-f) phase projections
# Requires: julia --project examples/oregonator_fpfsp.jl  (~40 s)
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

d = deserialize(joinpath(@__DIR__, "results", "oregonator_fpfsp.jls"))

snap_t    = d["snap_t"]
snap_traj = d["snap_traj"]   # rows = snapshots, cols = X/Y/Z
snap_sz   = d["snap_sizes"]
t_sub     = d["t_log"]
dt_sub    = d["dt_log"]
sz_sub    = d["size_log"]

# Drop t=0 snapshot for log-scale plots
pos = snap_t .> 0
t_pos     = snap_t[pos]
traj_pos  = snap_traj[pos, :]
sz_pos    = snap_sz[pos]

mean_X = traj_pos[:, 1]
mean_Y = traj_pos[:, 2]
mean_Z = traj_pos[:, 3]

# Panel (a): adaptive time step
p1 = plot(t_sub, dt_sub;
          xlabel="Time", ylabel="δt", title="(a)",
          label="", color=:black)

# Panel (b): state space size
p2 = plot(t_pos, sz_pos;
          xlabel="Time", ylabel="|S|", title="(b)",
          label="", color=:black)

# Panel (c): mean trajectories
p3 = plot(xlabel="Time", ylabel="Mean population", title="(c)")
plot!(p3, t_pos, mean_X; label="X", color=:black, linestyle=:solid)
plot!(p3, t_pos, mean_Y; label="Y", color=:black, linestyle=:dash)
plot!(p3, t_pos, mean_Z; label="Z", color=:black, linestyle=:dot)

# Panel (d): X-Y phase projection
p4 = plot(mean_X, mean_Y;
          xlabel="X", ylabel="Y", title="(d)",
          label="", color=:black)

# Panel (e): X-Z phase projection
p5 = plot(mean_X, mean_Z;
          xlabel="X", ylabel="Z", title="(e)",
          label="", color=:black)

# Panel (f): Y-Z phase projection
p6 = plot(mean_Y, mean_Z;
          xlabel="Y", ylabel="Z", title="(f)",
          label="", color=:black)

fig = plot(p1, p2, p3, p4, p5, p6;
           layout=(2, 3), size=(1100, 600), margin=5Plots.mm)

paper_dir = joinpath(@__DIR__, "..", "paper")
savefig(fig, joinpath(paper_dir, "orag.pdf"))
savefig(fig, joinpath(paper_dir, "orag.png"))
println("Saved orag.{pdf,png} → paper/")
println("total_iters=$(d["total_iters"]), wall=$(round(d["wall_time"]; digits=1))s")
println("max|S|=$(d["S_max"]), final|S|=$(d["S_final"])")
println("X=$(round(d["mean_X"]; digits=1)), Y=$(round(d["mean_Y"]; digits=1)), Z=$(round(d["mean_Z"]; digits=1))")
