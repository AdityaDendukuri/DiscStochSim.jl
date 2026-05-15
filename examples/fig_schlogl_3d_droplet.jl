"""
Schlögl RDME — 3D Critical Droplet on 2×2×2 Grid
==================================================

Two initial conditions on a 3D 2×2×2 voxel lattice (8 voxels, 12 edges):

  SUBCRITICAL  — 4 voxels at n_uns (tipping point), 4 at n_low
                 → droplet COLLAPSES: P(all-low) = 0.72 by t=30

  SUPERCRITICAL — 4 voxels at n_high (stable), 4 at n_low
                 → front PERSISTS: P(m=4 high) = 1.0 at t=30

Visualization: colored 3D voxel cubes, color = E[n] per voxel.
The multigrid FSP (solve_rdme_graph) works on arbitrary VoxelGraph,
enabling 3D spatial RDME computation.

Run: julia --project examples/fig_schlogl_3d_droplet.jl
"""

using DiscStochSim
using CairoMakie, Printf

e(k,d) = parse(typeof(d), get(ENV, k, string(d)))

const D         = e("D",         0.1)
const t_end     = e("T_END",     30.0)
const dt_fsp    = e("DT_FSP",    1.0)
const prune_tol = e("PRUNE_TOL", 1e-5)

model = SchloglModel1D(D)
n_low, n_uns, n_high = schlogl_fixed_points(model)
g3 = grid_3d(2, 2, 2)

@printf("grid_3d(2,2,2): %d voxels, %d edges\n", g3.n_voxels, length(g3.edges))
@printf("n_low=%d  n_uns=%d  n_high=%d  D=%.2f\n", n_low, n_uns, n_high, D)

snap_ts = unique([0.0; filter(<(t_end), [8.0, 18.0]); t_end])

alg = RDMEMultigridFSP(
    dt             = dt_fsp,
    n_levels       = 1,
    n_max          = 220,
    save_every     = 1,
    max_states     = 400_000,
    expand_depth   = 0,
    coarse_r_depth = 1,
    coarse_d_depth = 1,
    prune_tol      = prune_tol,
    weight_tol     = prune_tol,
    τ_pre          = 0.3,
)

# ICs: voxel ordering matches grid_3d(2,2,2)
# idx(i,j,k) = (i-1)*4 + (j-1)*2 + k
# Voxels 1-4: i=1 face;  5-8: i=2 face
u0_sup = CartesianIndex(n_high, n_high, n_high, n_high, n_low, n_low, n_low, n_low)
u0_sub = CartesianIndex(n_uns,  n_uns,  n_uns,  n_uns,  n_low, n_low, n_low, n_low)

println("Running supercritical FSP (4×n_high | 4×n_low)..."); flush(stdout)
@time sol_sup = solve_rdme_graph(model, g3, u0_sup, (0.0, t_end), alg)

println("Running subcritical FSP (4×n_uns | 4×n_low)..."); flush(stdout)
@time sol_sub = solve_rdme_graph(model, g3, u0_sub, (0.0, t_end), alg)

μ_sup = mean_voxel_counts(sol_sup)
μ_sub = mean_voxel_counts(sol_sub)

# Config probabilities over time
function all_cfg(sol, n_uns)
    K = sol.graph.n_voxels
    [[begin
        states, probs = sol.snapshots[ti]
        p = zeros(K+1)
        for (s,pr) in zip(states,probs)
            m = count(n -> n > n_uns, Tuple(s)); p[m+1] += pr
        end
        p
    end for ti in 1:length(sol.t)]...]
end

cfg_sup = all_cfg(sol_sup, n_uns)
cfg_sub = all_cfg(sol_sub, n_uns)

# Voxel positions in 3D (from graph)
pos = g3.positions   # 8×3 matrix: col 1=x, 2=y, 3=z

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11))
fig = Figure(size=(1400, 820))

Label(fig[0, 1:2],
    "Schlögl RDME — 3D Critical Droplet  |  2×2×2 grid  D=$D  " *
    "n_low=$n_low / n_uns=$n_uns / n_high=$n_high";
    fontsize=12, tellwidth=false)

c_lo = Float64(n_low); c_hi = Float64(n_high)
snap_idxs = [argmin(abs.(sol_sup.t .- ts)) for ts in snap_ts]

# ── Rows 1–2: 3D voxel snapshots ─────────────────────────────────────────────
for (row, (sol, mu, ic_label)) in enumerate([
    (sol_sup, μ_sup, "Supercritical IC: 4 voxels at n_high=$n_high (red), 4 at n_low=$n_low (blue)"),
    (sol_sub, μ_sub, "Subcritical IC:   4 voxels at n_uns=$n_uns (tipping point), 4 at n_low=$n_low"),
])
    for (col, ti) in enumerate(snap_idxs)
        ts = sol.t[ti]
        ax = Axis3(fig[row, col];
                   title="t = $(round(ts; digits=0))",
                   titlesize=10,
                   aspect=:equal,
                   azimuth  = 0.9π,
                   elevation = 0.18π,
                   xticksvisible=false, yticksvisible=false, zticksvisible=false,
                   xticklabelsvisible=false, yticklabelsvisible=false, zticklabelsvisible=false,
                   xlabelvisible=false, ylabelvisible=false, zlabelvisible=false)

        means_t = mu[ti, :]
        colors  = [(m - c_lo)/(c_hi - c_lo) for m in means_t]   # 0=blue, 1=red
        clrs    = [RGBAf(c, 0.1, 1.0-c, 0.85) for c in colors]

        # Draw each voxel as a sphere scaled by n
        xs = pos[:,1]; ys = pos[:,2]; zs = pos[:,3]
        meshscatter!(ax, xs, ys, zs;
                     markersize = 0.42,
                     color      = clrs,
                     shading    = FastShading)

        # Annotate E[n] on each voxel
        for k in 1:g3.n_voxels
            text!(ax, pos[k,1], pos[k,2], pos[k,3]+0.52;
                  text   = "$(round(Int, means_t[k]))",
                  fontsize = 8,
                  align  = (:center, :bottom),
                  color  = :black)
        end

        col == 1 && text!(ax, 0.0, 0.0, 2.6;
                          text=ic_label, fontsize=9, color=:gray30)
    end
    # Colorbar for the row
    Colorbar(fig[row, length(snap_ts)+1];
             colormap = Reverse(:RdBu),
             colorrange = (c_lo, c_hi),
             label  = "E[n]",
             width  = 14,
             ticks  = ([c_lo, Float64(n_uns), c_hi],
                        ["n_lo=$n_low", "n_uns=$n_uns", "n_hi=$n_high"]))
end

# ── Row 3: Config probabilities over time ─────────────────────────────────────
K = g3.n_voxels
cfg_clrs = [:navy, :steelblue4, :darkorange, :crimson,
             :forestgreen, :purple, :hotpink, :gray40, :black]

ax_sub = Axis(fig[3, 1:2];
              xlabel="time", ylabel="P(m voxels > n_uns)",
              title="Subcritical: droplet collapses — P(all-low) → 1",
              titlesize=11)
for m in 0:K
    p_m = [cfg_sub[ti][m+1] for ti in eachindex(sol_sub.t)]
    any(>(0.01), p_m) || continue
    lines!(ax_sub, sol_sub.t, p_m; color=cfg_clrs[m+1], linewidth=2.5,
           label="m=$m $(m==0 ? "(all-lo)" : m==K ? "(all-hi)" : "")")
end
axislegend(ax_sub; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_sub, -0.02, 1.05)

ax_sup = Axis(fig[3, 3:4];
              xlabel="time", ylabel="P(m voxels > n_uns)",
              title="Supercritical: front persists — P(m=4) = 1",
              titlesize=11)
for m in 0:K
    p_m = [cfg_sup[ti][m+1] for ti in eachindex(sol_sup.t)]
    any(>(0.01), p_m) || continue
    lines!(ax_sup, sol_sup.t, p_m; color=cfg_clrs[m+1], linewidth=2.5,
           label="m=$m $(m==0 ? "(all-lo)" : m==K ? "(all-hi)" : "")")
end
axislegend(ax_sup; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_sup, -0.02, 1.05)

# ── Row 4: |S| and mean at each voxel ────────────────────────────────────────
ax_ss = Axis(fig[4, 1:2]; xlabel="time", ylabel="|S|",
             title="FSP state space size  (graph V-cycle keeps it tractable in 3D)",
             titlesize=10)
lines!(ax_ss, sol_sub.t, Float64.(sol_sub.state_space_sizes);
       color=:steelblue4, linewidth=2, label="Subcritical")
lines!(ax_ss, sol_sup.t, Float64.(sol_sup.state_space_sizes);
       color=:crimson, linewidth=2, label="Supercritical")
axislegend(ax_ss; position=:rt, labelsize=9, framevisible=false)

ax_mu = Axis(fig[4, 3:4]; xlabel="time", ylabel="E[n]",
             title="Voxel means over time (solid=supercritical, dash=subcritical)",
             titlesize=10)
vclrs = [:crimson, :tomato, :darksalmon, :lightsalmon,
          :steelblue4, :cornflowerblue, :deepskyblue, :lightblue]
for k in 1:K
    lines!(ax_mu, sol_sup.t, μ_sup[:,k]; color=vclrs[k], linewidth=2,
           label="v$k")
    lines!(ax_mu, sol_sub.t, μ_sub[:,k]; color=vclrs[k], linewidth=1.5,
           linestyle=:dash)
end
hlines!(ax_mu, [Float64(n_high), Float64(n_low)]; color=(:black,0.2), linestyle=:dash, linewidth=1)
hlines!(ax_mu, [Float64(n_uns)]; color=(:gray50,0.4), linestyle=:dot, linewidth=1)
axislegend(ax_mu; position=:rc, labelsize=8, framevisible=false, nbanks=2)

rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 5); rowgap!(fig.layout, 3, 8)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_schlogl_3d_droplet.png")
save(out, fig; px_per_unit=2)
println("Saved → $out")
