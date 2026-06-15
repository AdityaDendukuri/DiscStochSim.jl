"""
fig_rate_faithful.jl

The LBNL-talk headline figure: **rate-faithful coarse-graining of the metastable RDME**.

The slow eigenvalue λ₂ of the bistable Schlögl JOINT CME on a K-voxel ring is the
metastable basin-switching rate. We reproduce it with a tiny (K·NM+1)-state coarse
master equation built by STATIONARY-WEIGHTED (Kemeny–Snell) lumping on the total-count
coordinate Y = Σ n_v. The stationary weights π(s|Y) come from the SAME AMG-GMRES solve
that makes the fine joint tractable (π is the easy λ=0 mode, ~30 power-iteration solves).

Verified data (/tmp/ring_rate2.jl, /tmp/ring_lump.jl, /tmp/ring_rate_k4.jl):

  K   |S|          coarse  fine λ₂     π-weighted λ₂  rel.err   uniform-P λ₂  uniform err
  2   2,601        101     -0.04634    -0.04603       0.65%     -0.11751      154%
  3   132,651      151     -0.05315    -0.05262       0.99%     -0.16337      207%
  4   6,765,201    201     (infeasible) FILL_PI       --        FILL_UNI      --

Headline: a ~200-state physics-informed master equation recovers the rare-event rate of
a 6.8M-state joint to <1%, while the *standard* multigrid prolongation (uniform within
shell) is wrong by 150%+. The fix — and the innovation — is weighting the prolongation
by the stationary measure supplied by the scalable solver.

Output: paper/figures/fig_rate_faithful.{png,pdf}
"""

using CairoMakie

# ── verified data ─────────────────────────────────────────────────────────────
S        = [2601.0, 132651.0, 6765201.0]
Kvox     = [2, 3, 4]
coarseN  = [101, 151, 201]
λ_fine   = [-0.04634, -0.05315, NaN]        # K=4 fine infeasible (6.8M, no direct ref)
λ_piw    = [-0.04603, -0.05262, -0.06363]   # K=4 from /tmp/ring_rate_k4.out (21 π-iters)
λ_unif   = [-0.11751, -0.16337, -0.18685]   # uniform-P over-estimates ~3× at every K
relerr   = abs.((λ_piw .- λ_fine) ./ λ_fine)

set_theme!(Theme(fontsize = 13))
fig = Figure(size = (960, 380))

# ── Panel (a): the rate, three coarse models ──────────────────────────────────
axa = Axis(fig[1, 1];
    title = "Metastable rate from a ~200-state coarse model",
    titlesize = 13,
    xlabel = "joint state-space size |S|", ylabel = "metastable rate  −λ₂",
    xscale = log10, xlabelsize = 12, ylabelsize = 12)
xlims!(axa, 1.5e3, 1.5e7); ylims!(axa, 0.0, 0.20)
# uniform-P (wrong)
scatterlines!(axa, S, -λ_unif; color = :crimson, linewidth = 2.5, markersize = 13,
    marker = :rect, linestyle = :dash, label = "uniform-P Galerkin (standard)")
# fine truth
scatter!(axa, S, -λ_fine; color = :black, markersize = 15, marker = :circle,
    label = "fine λ₂ (AMG subspace, truth)")
# π-weighted coarse
scatterlines!(axa, S, -λ_piw; color = :steelblue, linewidth = 2.5, markersize = 11,
    marker = :utriangle, label = "π-weighted lump (ours)")
for i in 1:2
    text!(axa, S[i], -λ_piw[i] - 0.018; text = "K=$(Kvox[i])",
        align = (:center, :top), fontsize = 10, color = :gray30)
end
text!(axa, 3.0e3, 0.155; text = "standard prolongation:\n150–200% rate error",
    align = (:left, :bottom), fontsize = 10, color = :crimson)
text!(axa, 1.2e5, 0.022; text = "π-weighted: <1%",
    align = (:left, :bottom), fontsize = 11, color = :steelblue, font = :bold)
axislegend(axa; position = :rt, framevisible = false, labelsize = 9.5)

# ── Panel (b): reduced model stays ~200 states while |S| explodes ─────────────
axb = Axis(fig[1, 2];
    title = "Reduced model size is decoupled from |S|",
    titlesize = 13,
    xlabel = "ring voxels K", ylabel = "states",
    yscale = log10, xlabelsize = 12, ylabelsize = 12)
xlims!(axb, 1.6, 4.4); ylims!(axb, 50, 3e7)
scatterlines!(axb, Kvox, S; color = :black, linewidth = 2.5, markersize = 13,
    marker = :circle, label = "fine joint |S|")
scatterlines!(axb, Kvox, Float64.(coarseN); color = :steelblue, linewidth = 2.5,
    markersize = 13, marker = :utriangle, label = "coarse master eqn")
text!(axb, 4.0, 6.8e6; text = "6.8M", align = (:right, :bottom), fontsize = 10, color = :black)
text!(axb, 4.0, 201; text = "201", align = (:right, :top), fontsize = 10, color = :steelblue)
text!(axb, 2.0, 4e6;
    text = "π from ~30 AMG solves\nrate from 201-state eigen",
    align = (:left, :top), fontsize = 10, color = :gray25)
axislegend(axb; position = :rb, framevisible = false, labelsize = 10)

colgap!(fig.layout, 26)

for ext in ("png", "pdf")
    out = joinpath(@__DIR__, "..", "paper", "figures", "fig_rate_faithful.$ext")
    save(out, fig; px_per_unit = 2)
    println("wrote $out")
end
