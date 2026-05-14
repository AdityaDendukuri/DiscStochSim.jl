"""
K=8 Schlögl persistent front: state-space sparsity + mean heatmap.

Two-panel figure:
  (a) |S| vs time — grows then plateaus: only the interface is active
  (b) Mean voxel count E[nₖ] heatmap — hi on left, lo on right, narrow interface

Run: julia --project examples/fig_k8_front.jl
"""

using DiscStochSim
using CairoMakie

# ── model ──────────────────────────────────────────────────────────────────────
K = 8;  D = 0.005;  h = 1.0 / K
model = SchloglModel1D(D)
grid  = VoxelGrid(K, h, 0)
n_low, n_unstable, n_high = schlogl_fixed_points(model)
println("Fixed points: lo=$n_low  uns=$n_unstable  hi=$n_high")
println("Diffusion rate d = $(round(D/h^2, digits=2))")

# Persistent front: left half hi, right half lo
u0 = CartesianIndex(n_high, n_high, n_high, n_high, n_low, n_low, n_low, n_low)

alg = RDMEMultigridFSP(
    dt            = 0.5,
    n_levels      = 1,
    n_max         = 250,
    save_every    = 1,
    max_states    = 300_000,
    expand_depth  = 0,
    coarse_r_depth = 1,
    coarse_d_depth = 1,
    prune_tol     = 1e-6,
    weight_tol    = 1e-6,
    τ_pre         = 0.5,
)

println("Running K=$K n_levels=1 persistent front ...")
@time sol = solve_rdme_multigrid(model, grid, u0, (0.0, 8.0), Float64[], alg)
println("Done. $(length(sol.t)) snapshots, peak |S|=$(maximum(sol.state_space_sizes))")

mu = mean_voxel_counts(sol)   # (n_t × K) matrix

# ── figure layout ──────────────────────────────────────────────────────────────
set_theme!(Theme(
    fontsize    = 11,
    Axis        = (spinewidth = 0.7, xgridvisible = false, ygridvisible = false,
                   ticksize = 3, tickwidth = 0.6f0),
    Lines       = (linewidth = 1.6,),
))

fig = Figure(size = (1120, 310))

# ── panel (a): |S| vs time ────────────────────────────────────────────────────
ax_s = Axis(fig[1, 1];
    xlabel = "time  t",
    ylabel = "fine states  |S|",
    title  = "(a)  State-space: grows then plateaus",
    titlesize = 11,
    ytickformat = v -> ["$(round(Int,x/1000))k" for x in v],
)

lines!(ax_s,  sol.t, Float64.(sol.state_space_sizes); color = :steelblue4)
scatter!(ax_s, sol.t, Float64.(sol.state_space_sizes);
         color = :steelblue4, markersize = 5)

# Plateau annotation
plateau = maximum(sol.state_space_sizes)
hlines!(ax_s, [Float64(plateau)];
        color = (:gray50, 0.55), linestyle = :dash, linewidth = 1.2)

# Annotate the two phases
t_turn = sol.t[findfirst(x -> x > plateau * 0.85, sol.state_space_sizes)]
text!(ax_s, t_turn * 0.3, plateau * 0.55;
      text = "front\nestablishes",
      fontsize = 9, color = :gray40, align = (:center, :center))
text!(ax_s, sol.t[end] * 0.45, plateau * 0.86;
      text = "plateau:\nbulk voxels\nadd no new states",
      fontsize = 9, color = :gray40, align = (:left, :top))

ylims!(ax_s, 0, plateau * 1.22)

# ── panel (b): mean count heatmap ────────────────────────────────────────────
ax_h = Axis(fig[1, 2];
    xlabel = "voxel  k",
    ylabel = "time  t",
    title  = "(b)  Mean count E[nₖ]:  K = $K, D = $D",
    titlesize = 11,
    xticks = (1:K, ["v$k" for k in 1:K]),
)

# Build (K × n_t) matrix for heatmap(xs=voxels, ys=times, zs[k,ti])
mu_mat = Matrix(mu')   # (K × n_t)

hm = heatmap!(ax_h, 1:K, sol.t, mu_mat;
    colormap   = Reverse(:RdBu),
    colorrange = (n_low - 8, n_high + 8),
    lowclip    = :dodgerblue4,
    highclip   = :firebrick3,
)

# Dashed line at interface
vlines!(ax_h, [4.5]; color = (:white, 0.7), linewidth = 1.4, linestyle = :dash)
text!(ax_h, 4.55, sol.t[end] * 0.93;
      text = "interface", fontsize = 8, color = (:white, 0.8))

cb = Colorbar(fig[1, 3], hm;
    label  = "E[nₖ]",
    ticks  = ([n_low, n_unstable, n_high],
              ["$n_low (lo)", "$n_unstable (uns)", "$n_high (hi)"]),
    width  = 14, ticksize = 4, labelsize = 10,
)

# ── panel (c): marginal distributions at t=8 for each voxel ──────────────────
ti_final = length(sol.t)
t_final  = sol.t[ti_final]

function voxel_marginal(sol, ti, vox, n_max)
    states, probs = sol.snapshots[ti]
    P = zeros(n_max + 1)
    for (s, p) in zip(states, probs)
        n = Tuple(s)[vox]; n <= n_max && (P[n+1] += p)
    end
    P
end

n_show = 200
ax_c = Axis(fig[1, 4];
    xlabel    = "molecules  n",
    ylabel    = "P(nₖ = n)",
    title     = "(c)  Marginals at t = $(round(t_final, digits=0))",
    titlesize = 11,
)

colors_vox = Makie.resample_cmap(Reverse(:RdBu), K)

for k in 1:K
    P = voxel_marginal(sol, ti_final, k, n_show)
    n_vals = 0:n_show
    lines!(ax_c, n_vals, P; color = colors_vox[k], linewidth = 1.3, label = "v$k")
end

# Reference lines at stable points
vlines!(ax_c, [n_low];      color = (:steelblue, 0.5), linestyle = :dash, linewidth = 1.0)
vlines!(ax_c, [n_high];     color = (:firebrick, 0.5), linestyle = :dash, linewidth = 1.0)
vlines!(ax_c, [n_unstable]; color = (:gray60, 0.4),    linestyle = :dot,  linewidth = 0.8)

text!(ax_c, n_low + 2,    0.055; text = "nₗₒ", fontsize = 8, color = :steelblue4)
text!(ax_c, n_high + 2,   0.055; text = "nₕᵢ", fontsize = 8, color = :firebrick4)

xlims!(ax_c, 0, n_show)
axislegend(ax_c; position = :rt, labelsize = 8, framevisible = false, nbanks = 2)

colgap!(fig.layout, 1, 22)
colgap!(fig.layout, 2, 6)
colgap!(fig.layout, 3, 18)
colsize!(fig.layout, 1, Relative(0.28))
colsize!(fig.layout, 2, Relative(0.30))

# ── save ───────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "paper", "figures", "fig_k8_front.pdf")
save(out, fig)
println("\nSaved → $out")
fig
