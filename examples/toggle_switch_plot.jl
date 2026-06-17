# Toggle switch results figure (toggle_switch_grid.pdf)
# 3×4 grid: rows = configs, cols = t≈2 / t≈15 / t≈30 / |S| trajectory
# Requires: julia --project examples/toggle_switch.jl  (~10 min total)
using Plots, Serialization

gr()
default(
    fontfamily    = "Computer Modern",
    titlefontsize  = 10,
    guidefontsize  = 9,
    tickfontsize   = 8,
    legendfontsize = 8,
    linewidth = 1.0,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

results_dir = joinpath(@__DIR__, "results")
paper_dir   = joinpath(@__DIR__, "..", "paper")

# Load all 3 configs
data = [deserialize(joinpath(results_dir, "toggle_switch_$(i).jls")) for i in 1:3]

NGRID = data[1]["NGRID"]  # 300

row_labels = [
    "α=0.3, no flux",
    "α=0.3, flux=1e-7",
    "α=0.001, no flux",
]

target_times = data[1]["target_times"]   # [2.0, 15.0, 30.0]
time_labels  = ["t ≈ 2", "t ≈ 15", "t ≈ 30"]

# Build one panel per (config, time snapshot): log₁₀(p) contours + dotted boundary
function make_prob_panel(grid_entry, row_idx, col_idx)
    P_prob   = grid_entry.P_prob    # Float32 300×300
    P_binary = grid_entry.P_binary  # Float32 300×300, 1 where active
    t_actual = grid_entry.t

    # log₁₀ probability, masking zeros
    logP = fill(NaN32, NGRID, NGRID)
    for i in 1:NGRID, j in 1:NGRID
        p = P_prob[i, j]
        p > 0f0 && (logP[i, j] = log10(p))
    end

    finite_vals = filter(isfinite, logP)
    vmin = isempty(finite_vals) ? -10.0 : max(minimum(finite_vals), -10.0)
    vmax = isempty(finite_vals) ?   0.0 : maximum(finite_vals)

    xlabel_str = col_idx == 2 ? "U" : ""
    ylabel_str = row_idx == 2 ? "V" : ""
    title_str  = row_idx == 1 ? time_labels[col_idx] : ""

    p = heatmap(1:NGRID, 1:NGRID, logP';
        clims=(vmin, vmax),
        color=:viridis,
        xlabel=xlabel_str, ylabel=ylabel_str,
        title=title_str,
        colorbar=false,
        aspect_ratio=:equal,
        xlims=(1, NGRID), ylims=(1, NGRID),
    )

    # Dotted truncation boundary: contour of P_binary at 0.5
    contour!(p, 1:NGRID, 1:NGRID, P_binary';
        levels=[0.5],
        color=:white, linestyle=:dot, linewidth=1.0,
        colorbar=false, label="",
    )

    # Row annotation on leftmost column
    if col_idx == 1
        annotate!(p, NGRID*0.02, NGRID*0.95,
            text(row_labels[row_idx], :white, :left, 7))
    end

    p
end

# Build right-column size trajectory panels
function make_size_panel(d, row_idx)
    t_log    = d["t_log"]
    size_log = d["size_log"]

    ylabel_str = row_idx == 2 ? "|S|" : ""
    title_str  = row_idx == 1 ? "|S| trajectory" : ""
    xlabel_str = row_idx == 3 ? "Time" : ""

    plot(t_log, size_log;
        xlabel=xlabel_str, ylabel=ylabel_str, title=title_str,
        label="", color=:black, linewidth=1.0,
    )
end

panels = []
for row in 1:3
    d = data[row]
    grids = d["grids"]  # Vector of (t, P_binary, P_prob) NamedTuples
    for col in 1:3
        push!(panels, make_prob_panel(grids[col], row, col))
    end
    push!(panels, make_size_panel(d, row))
end

fig = plot(panels...;
    layout=(3, 4),
    size=(1100, 750),
    margin=3Plots.mm,
)

savefig(fig, joinpath(paper_dir, "toggle_switch_grid.pdf"))
savefig(fig, joinpath(paper_dir, "toggle_switch_grid.png"))
println("Saved toggle_switch_grid.{pdf,png} → paper/")
