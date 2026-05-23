"""
2D A+B→∅ Annihilation — Mixed-Resolution Visualization
========================================================

A produced in top-left corner block, B in bottom-right corner block.
They diffuse through a 20×20 grid and annihilate wherever they meet.

Color scheme:  RED = E[nA],  BLUE = E[nB],  PURPLE = both present (reaction zone)

Two layers:
  Mean-field ODE  : exact means for the full 20×20 grid (cheap)
  Mixed-resolution: identifies the "active strip" where the stochastic FSP
                    would be applied to capture fluctuations at the front.
                    Outside the strip: mean-field is sufficient.

Scaling argument (shown in final panel):
  Full 2D FSP on K×K grid: CartesianIndex{2K²} — intractable beyond K≈3
  Mixed-resolution (active strip width W):
    1D along diagonal, W~√(D/k_rxn) voxels wide:
    Ka ≈ W × K_y active voxels.  Even for K=64:
    Ka ≈ 4 × 64 = 256 — still intractable for joint FSP.
  Column-product factoring: treat each active column independently,
    each a 1D FSP with Ka≈2-4. |S| per column ≈ (n_max+1)^4 ≈ thousands.
    This is the natural 2D extension of the 1D mixed-resolution approach.

Run: julia --project examples/mixed_annihilation_2d.jl
"""

using DiscStochSim
using CairoMakie
using Printf
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# 2D mean-field ODE for A+B→∅
# ─────────────────────────────────────────────────────────────────────────────

function run_ode_2d(K_x, K_y, D_A, D_B, k_A, k_B, k_rxn,
                   src_A, src_B,    # BitMatrix K_x×K_y
                   T_END, DT)
    dA = D_A; dB = D_B   # dx=1 voxel unit; D is in voxels²/time

    μA = zeros(K_x, K_y)
    μB = zeros(K_x, K_y)
    # Rough initial condition: source voxels near SS
    for i in 1:K_x, j in 1:K_y
        src_A[i,j] && (μA[i,j] = 0.5 * k_A / k_rxn)
        src_B[i,j] && (μB[i,j] = 0.5 * k_B / k_rxn)
    end

    max_rate = 4*max(dA,dB) + k_rxn * max(maximum(μA), maximum(μB), 1.0)
    n_sub    = max(1, ceil(Int, max_rate * DT * 2))
    dt_sub   = DT / n_sub

    n_steps = round(Int, T_END / DT)
    snaps = Tuple{Float64, Matrix{Float64}, Matrix{Float64}}[]
    push!(snaps, (0.0, copy(μA), copy(μB)))

    ΔA = zeros(K_x, K_y); ΔB = zeros(K_x, K_y)
    for step in 1:n_steps
        for _ in 1:n_sub
            fill!(ΔA, 0.0); fill!(ΔB, 0.0)
            for i in 1:K_x, j in 1:K_y
                rxn = k_rxn * μA[i,j] * μB[i,j]
                ΔA[i,j] = (src_A[i,j] ? k_A : 0.0) - rxn
                ΔB[i,j] = (src_B[i,j] ? k_B : 0.0) - rxn
                # 4-connected diffusion (no-flux BC)
                i > 1   && (ΔA[i,j] += dA*(μA[i-1,j]-μA[i,j]); ΔB[i,j] += dB*(μB[i-1,j]-μB[i,j]))
                i < K_x && (ΔA[i,j] += dA*(μA[i+1,j]-μA[i,j]); ΔB[i,j] += dB*(μB[i+1,j]-μB[i,j]))
                j > 1   && (ΔA[i,j] += dA*(μA[i,j-1]-μA[i,j]); ΔB[i,j] += dB*(μB[i,j-1]-μB[i,j]))
                j < K_y && (ΔA[i,j] += dA*(μA[i,j+1]-μA[i,j]); ΔB[i,j] += dB*(μB[i,j+1]-μB[i,j]))
            end
            μA .= max.(0.0, μA .+ ΔA .* dt_sub)
            μB .= max.(0.0, μB .+ ΔB .* dt_sub)
        end
        push!(snaps, (step * DT, copy(μA), copy(μB)))
    end
    snaps
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

const K_X   = 20
const K_Y   = 20
# D in VOXEL units (dx=1).  Must be O(K²·D_physical) to get macroscopic transport.
# For D_physical=0.01 on [0,1] with K=20: D_vox = 0.01*20² = 4.
const D_A   = 5.0      # voxels²/time
const D_B   = 5.0
const K_A   = 3.0      # ∅→A rate per source voxel
const K_B   = 3.0
const K_RXN = 0.20     # A+B→∅  (weak reaction → wide, visible front)
const T_END = 20.0     # long enough for species to meet at center
const DT    = 0.2

# Source regions: top-left 4×4 block for A; bottom-right 4×4 for B
const SRC_WIDTH = 4
src_A = falses(K_X, K_Y);  src_A[1:SRC_WIDTH,         1:SRC_WIDTH]           .= true
src_B = falses(K_X, K_Y);  src_B[K_X-SRC_WIDTH+1:K_X, K_Y-SRC_WIDTH+1:K_Y]  .= true

# Active-strip threshold: voxel is "active" (needs FSP) when both species present
const ACTIVE_THRESH = 0.05

println("Running 2D ODE  ($(K_X)×$(K_Y) grid)...")
@time snaps = run_ode_2d(K_X, K_Y, D_A, D_B, K_A, K_B, K_RXN,
                          src_A, src_B, T_END, DT)

println("Done ($(length(snaps)) snapshots).")

# ─────────────────────────────────────────────────────────────────────────────
# Build RGB image: red channel = μA, blue channel = μB
# ─────────────────────────────────────────────────────────────────────────────

# Find global max across all snapshots for consistent colorscale
max_A = max(1e-9, maximum(maximum(s[2]) for s in snaps))
max_B = max(1e-9, maximum(maximum(s[3]) for s in snaps))
# Use same scale for A and B so equal amounts → symmetric purple
scale = max(max_A, max_B)

function ab_image(μA, μB, scale; gamma=0.45)
    K_x, K_y = size(μA)
    # gamma < 1 lifts dim reaction-zone pixels so the purple front is visible
    # alongside bright red/blue source regions
    [RGBAf(clamp(μA[i, K_y+1-j] / scale, 0f0, 1f0)^Float32(gamma),
           0f0,
           clamp(μB[i, K_y+1-j] / scale, 0f0, 1f0)^Float32(gamma),
           1f0)
     for i in 1:K_x, j in 1:K_y]
end

function active_mask(μA, μB, thresh)
    [μA[i,j] > thresh && μB[i,j] > thresh for i in 1:size(μA,1), j in 1:size(μA,2)]
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure: 4 time snapshots + active-strip panel + scaling panel
# ─────────────────────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

snap_ts  = [0.0, 4.0, 10.0, T_END]
snap_idx = [argmin(abs(s[1] - tt) for s in snaps) for tt in snap_ts]

n_steps = round(Int, T_END / DT)
times   = [s[1] for s in snaps]
total_A = [sum(s[2]) for s in snaps]
total_B = [sum(s[3]) for s in snaps]
# Reaction rate density at center voxel over time
center_rxn = [K_RXN * snaps[i][2][K_X÷2, K_Y÷2] * snaps[i][3][K_X÷2, K_Y÷2]
              for i in eachindex(snaps)]

# Active strip voxel count over time
n_active_t = [sum(active_mask(s[2], s[3], ACTIVE_THRESH)) for s in snaps]

fig = Figure(size=(1600, 1050), backgroundcolor=:black)

# ── Row 1: four time snapshots ────────────────────────────────────────────────
snap_axes = Axis[]
for (ci, si) in enumerate(snap_idx)
    t_s  = snap_ts[ci]
    μA_s = snaps[si][2]; μB_s = snaps[si][3]
    img  = ab_image(μA_s, μB_s, scale)
    mask = active_mask(μA_s, μB_s, ACTIVE_THRESH)

    ax = Axis(fig[1, ci];
              title      = "t = $(round(t_s, digits=1))",
              titlecolor = :white, titlesize = 14,
              aspect     = DataAspect(),
              backgroundcolor = :black,
              xgridvisible=false, ygridvisible=false,
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false,
              leftspinevisible=false, rightspinevisible=false,
              topspinevisible=false, bottomspinevisible=false)

    image!(ax, img)

    # Active-strip overlay: white contour where BOTH species > threshold
    # Draw it as a heatmap overlay with alpha
    mask_mat = Float64.(mask)
    heatmap!(ax, 1:K_X, 1:K_Y, mask_mat;
             colormap = [RGBAf(0,0,0,0), RGBAf(1,1,1,0.25)],
             colorrange=(0,1))
    push!(snap_axes, ax)
end

# ── Custom colorbar ───────────────────────────────────────────────────────────
cb_ax = Axis(fig[1, 5];
             title      = "color key",
             titlecolor = :white, titlesize = 12,
             aspect     = DataAspect(),
             backgroundcolor = :black,
             xgridvisible=false, ygridvisible=false,
             xticksvisible=false, yticksvisible=false,
             xticklabelsvisible=false, yticklabelsvisible=false,
             leftspinevisible=false, rightspinevisible=false,
             topspinevisible=false, bottomspinevisible=false)
# Build a gradient: top=pure A (red), bottom=pure B (blue), center=equal (purple)
N_cb = 100
cb_img = [RGBAf(Float32(1 - r/N_cb), 0f0, Float32(r/N_cb), 1f0)
          for r in 0:N_cb-1, _ in 1:15]
image!(cb_ax, cb_img)
text!(cb_ax, 7, 2;    text="B (blue)",    color=:white, fontsize=10, align=(:center,:bottom))
text!(cb_ax, 7, 50;   text="both\n(purple)",color=:white, fontsize=10, align=(:center,:center))
text!(cb_ax, 7, 98;   text="A (red)",    color=:white, fontsize=10, align=(:center,:top))
# White border
Ns=N_cb; lw=2
lines!(cb_ax,[1,15,15,1,1],[1,1,Ns,Ns,1];color=:white,linewidth=lw)

# ── Row 2: reaction dynamics ──────────────────────────────────────────────────
ax_tot = Axis(fig[2, 1:2];
              title      = "Total molecules over time",
              xlabel     = "time", ylabel = "Σ E[nA] or Σ E[nB]",
              titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
              xticklabelcolor = :white, yticklabelcolor = :white,
              backgroundcolor = :black,
              xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))
lines!(ax_tot, times, total_A; color=:red,   linewidth=2.5, label="Σ E[nA]")
lines!(ax_tot, times, total_B; color=:blue,  linewidth=2.5, label="Σ E[nB]")
for tt in snap_ts
    vlines!(ax_tot, [tt]; color=RGBAf(1,1,1,0.4), linewidth=1, linestyle=:dash)
end
axislegend(ax_tot; labelcolor=:white, backgroundcolor=:transparent, framecolor=:transparent, labelsize=11)

ax_active = Axis(fig[2, 3:4];
                 title      = "Active strip size (voxels with both A,B > $(ACTIVE_THRESH))",
                 xlabel     = "time", ylabel = "active voxels",
                 titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
                 xticklabelcolor = :white, yticklabelcolor = :white,
                 backgroundcolor = :black,
                 xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))
lines!(ax_active, times, Float64.(n_active_t); color=:white, linewidth=2.5)
for tt in snap_ts
    vlines!(ax_active, [tt]; color=RGBAf(1,1,1,0.4), linewidth=1, linestyle=:dash)
end

# ── Row 2 far right: scaling comparison ──────────────────────────────────────
ax_scale = Axis(fig[2, 5];
                title      = "State-space scaling",
                xlabel     = "grid size K (K×K)",
                ylabel     = "log₁₀ |S|",
                titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
                xticklabelcolor = :white, yticklabelcolor = :white,
                backgroundcolor = :black,
                xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))

Ks      = [2, 3, 4, 5, 6, 8, 10, 12]
n_max   = 8
K_strip = 4   # active strip width (voxels)
Ka_col  = 3   # active voxels per column (1D FSP)

log_full  = [min(2*K^2 * log10(n_max+1), 30.0) for K in Ks]   # full 2D FSP
log_strip = [min((K_strip * K) * log10(n_max+1), 30.0) for K in Ks]  # joint strip
log_col   = [log10((n_max+1)^(2*Ka_col)) + log10(K_strip * K) for K in Ks]  # column-product

lines!(ax_scale, Float64.(Ks), log_full;
       color=:red,   linewidth=2.5, label="Full 2D FSP (intractable >K≈3)")
lines!(ax_scale, Float64.(Ks), log_strip;
       color=:orange,linewidth=2,   linestyle=:dash, label="Joint strip FSP")
lines!(ax_scale, Float64.(Ks), log_col;
       color=:lime,  linewidth=2.5, label="Column-product (1D FSPs, tractable!)")
hlines!(ax_scale, [6.0]; color=RGBAf(1,1,1,0.5), linewidth=1,
        linestyle=:dot, label="~10⁶ (memory limit)")
axislegend(ax_scale; labelcolor=:white, backgroundcolor=:transparent,
           framecolor=:transparent, labelsize=9, position=:lt)
ylims!(ax_scale, 0, 32)

# ── Labels ────────────────────────────────────────────────────────────────────
Label(fig[0, 1:5],
      "2D A+B→∅  |  A source: top-left corner,  B source: bottom-right corner  |  " *
      "D=$(D_A), k_rxn=$(K_RXN)  |  White overlay = active FSP strip (both species present)",
      color=:white, fontsize=13)

Label(fig[3, 1:5],
      "Red = A only   ·   Blue = B only   ·   Purple/magenta = reaction zone (A∩B)   ·   " *
      "Black = empty   ·   White outline = active FSP strip (needs stochastic tracking)",
      color=RGBAf(0.8,0.8,0.8,1), fontsize=11)

png_path = joinpath(outdir, "mixed_annihilation_2d.png")
save(png_path, fig; px_per_unit=3)
println("Saved: $png_path")
println("Done.")
