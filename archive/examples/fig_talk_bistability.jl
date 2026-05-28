"""
fig_talk_bistability.jl

Talk figure: single-voxel Schlögl bistability.
Two panels:
  (a) Rate equation f(n) = net production — three zero crossings show bistability
  (b) Stationary distribution P_stat(n) — bimodal, peaks at n_lo and n_hi

Fast: no time integration needed.
"""

using DiscStochSim
using CairoMakie
using LinearAlgebra
using SparseArrays
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)
save_fig(name, fig) = begin
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=3)
    println("  Saved: $name.{pdf,png}")
end

# ── colors ─────────────────────────────────────────────────────────────────────
const C_LO   = RGBf(0.208, 0.475, 0.694)   # blue  — low state
const C_HI   = RGBf(0.835, 0.369, 0.000)   # orange — high state
const C_UN   = RGBf(0.600, 0.600, 0.600)   # gray  — unstable
const C_RATE = RGBf(0.122, 0.122, 0.122)

set_theme!(Theme(
    fontsize = 12,
    Axis = (spinewidth=0.8, xgridvisible=false, ygridvisible=false,
            xticklabelsize=10, yticklabelsize=10, xlabelsize=11, ylabelsize=11,
            titlesize=11),
))

model = SchloglModel1D(0.005)
fp    = schlogl_fixed_points(model; n_max=250)
n_lo, n_un, n_hi = fp
@printf("Fixed points: n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

n_max = 220

# ── (a) rate equation ──────────────────────────────────────────────────────────
ns = 0:0.5:n_max
fs = [schlogl_rate_eq(model, n) for n in ns]

# ── (b) stationary distribution via null vector of A ──────────────────────────
println("Computing stationary distribution...")
K_stat = 1
h_stat = 1.0
grid1  = VoxelGrid(K_stat, h_stat, 0)
sys1   = build_schlogl_rdme_system(model, grid1)

sp1 = StateSpace{CartesianIndex{1}, Float64}()
for n in 0:n_max
    add_state!(sp1, CartesianIndex(n), 0.0)
end
sp1.probs[sp1.index[CartesianIndex(n_hi)]] = 1.0

A1, = build_generator(sp1, sys1, Float64[], 0.0)
# Long-time evolution to approximate stationary distribution
using ExponentialUtilities
pstat = copy(sp1.probs)
for _ in 1:20
    global pstat
    pstat = expv(5.0, A1, pstat; m=50)
    pstat ./= sum(pstat)
end
ns_stat = [Tuple(sp1.states[i])[1] for i in eachindex(sp1.states)]
ps_stat = pstat
println("  done  Σp=$(round(sum(ps_stat),digits=5))")

# ── figure ─────────────────────────────────────────────────────────────────────
fig = Figure(size=(700, 300))

# Panel (a): rate equation
ax1 = Axis(fig[1,1];
    xlabel = "molecule count  n",
    ylabel = "f(n)  (net production rate)",
    title  = "(a)  Schlögl rate equation")
hlines!(ax1, [0.0]; color=(:black, 0.3), linewidth=0.8, linestyle=:dash)
lines!(ax1, collect(ns), fs; color=C_RATE, linewidth=2.0)
# Mark fixed points
scatter!(ax1, [n_lo]; color=C_LO,   markersize=10, marker=:circle, label="n_lo (stable)")
scatter!(ax1, [n_un]; color=C_UN,   markersize=8,  marker=:utriangle, label="n_un (unstable)")
scatter!(ax1, [n_hi]; color=C_HI,   markersize=10, marker=:circle, label="n_hi (stable)")
axislegend(ax1; position=:rt, framevisible=false, labelsize=9, rowgap=1)
ylims!(ax1, -8, 8)

# Panel (b): stationary distribution
ax2 = Axis(fig[1,2];
    xlabel = "molecule count  n",
    ylabel = "P_stat(n)",
    title  = "(b)  Stationary distribution")
barplot!(ax2, ns_stat, ps_stat; color=(:steelblue, 0.6), strokewidth=0, gap=0)
vlines!(ax2, [n_lo, n_un, n_hi]; color=[C_LO, C_UN, C_HI],
        linewidth=1.5, linestyle=:dash)
text!(ax2, n_lo-2, maximum(ps_stat)*0.85; text="n_lo", color=C_LO,
      fontsize=10, align=(:right,:center))
text!(ax2, n_hi+2, maximum(ps_stat)*0.85; text="n_hi", color=C_HI,
      fontsize=10, align=(:left,:center))
xlims!(ax2, 0, n_max)

save_fig("fig_talk_bistability", fig)
println("Done.")
