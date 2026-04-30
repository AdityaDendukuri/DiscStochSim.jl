"""
fig_prob_diffusion.jl

Shows how the probability distribution evolves on the K=2 Schlögl state space.
State = (n_1, n_2); axes are voxel molecule counts.

Two initial conditions:
  (a) Uniform pair:  IC = (n_lo, n_lo) — both voxels in the same stable state
  (b) Front pair:    IC = (n_lo, n_hi) — voxels in opposite stable states

Visualizes:
  - Diffusion:  moves probability ALONG the constraint lines n_1+n_2 = const
                (molecules swap between voxels, total conserved)
  - Reactions:  move probability ACROSS constraint lines
                (molecules created/destroyed, total changes)

The key contrast:
  - Uniform pair:  diffusion and reactions cooperate; conditional → Binomial
  - Front pair:    bistable reactions fight diffusion; probability stays at corners

Run from repo root:
  julia --project examples/fig_prob_diffusion.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

function save_fig(name, fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=2)
    println("  Saved: $name.{pdf,png}")
end

set_theme!(theme_light(); fontsize=11,
    Axis=(spinewidth=0.8, xtickwidth=0.7, ytickwidth=0.7,
          xgridvisible=false, ygridvisible=false))

# ── Model: K=2, D=0.005 ───────────────────────────────────────────────────────
D     = 0.005
K     = 2
h     = 1.0 / K        # h = 0.5, so d = D/h² = 0.02
model = SchloglModel1D(D)
grid  = VoxelGrid(K, h, 0)
sys   = build_schlogl_rdme_system(model, grid)

fp = schlogl_fixed_points(model; n_max=250)
n_lo, n_uns, n_hi = fp
@printf("Fixed points: n_lo=%d  n_uns=%d  n_hi=%d\n", n_lo, n_uns, n_hi)
@printf("d = D/h² = %.4f\n", D/h^2)

# ── Build full (n_1, n_2) state space up to n_max ────────────────────────────
# We build the full dense grid once so the generator A is fixed.
n_max = 175   # covers both fixed points (31 and 162) with margin
println("\nBuilding full state space ($(n_max+1)² = $((n_max+1)^2) states)...")

sp_base = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_base, CartesianIndex(n1, n2), 0.0)
end

A, = build_generator(sp_base, sys, Float64[], 0.0)
println("  Generator built: size=$(size(A,1))×$(size(A,2))  nnz=$(nnz(A))")

# ── Evolve from a given IC ────────────────────────────────────────────────────
function evolve(u0::CartesianIndex{2}, times::Vector{Float64})
    p = zeros(length(sp_base.states))
    p[sp_base.index[u0]] = 1.0
    snaps = Dict{Float64, Matrix{Float64}}()
    # t=0 snapshot
    mat = zeros(n_max+1, n_max+1)
    for (i,s) in enumerate(sp_base.states)
        t = Tuple(s)
        mat[t[1]+1, t[2]+1] = p[i]
    end
    snaps[0.0] = mat

    t_prev = 0.0
    for t in sort(times[times .> 0])
        p .= max.(expv(t - t_prev, A, p; m=40), 0.0)
        p ./= sum(p)
        mat = zeros(n_max+1, n_max+1)
        for (i,s) in enumerate(sp_base.states)
            tv = Tuple(s)
            mat[tv[1]+1, tv[2]+1] = p[i]
        end
        snaps[t] = mat
        t_prev = t
        @printf("  t=%.1f  Σp=%.6f  peak at (%d,%d)\n", t, sum(p),
                Tuple(sp_base.states[argmax(p)])[1],
                Tuple(sp_base.states[argmax(p)])[2])
    end
    snaps
end

T_PLOT = [0.0, 1.0, 5.0, 20.0]

println("\n── Uniform IC: ($n_lo, $n_lo) ──")
snaps_uni  = evolve(CartesianIndex(n_lo, n_lo),  T_PLOT)

println("\n── Front IC: ($n_lo, $n_hi) ──")
snaps_front = evolve(CartesianIndex(n_lo, n_hi), T_PLOT)

# ── Plotting ──────────────────────────────────────────────────────────────────
# Zoom window for display: show [0, n_zoom] × [0, n_zoom]
n_zoom = 170
idx = 1:n_zoom+1   # indices into the n_max+1 array

# For the front case, annotate the two corners and the Binomial peak
n_bar = n_lo + n_hi   # = 193; Binomial peak at n_bar/2

println("\n━━━ Plotting ━━━")

fig = Figure(size = (920, 520))

titles_uni  = ["(a)  t = $t" for t in T_PLOT]
titles_front = ["(b)  t = $t" for t in T_PLOT]

# Shared color range: use log scale, clip very small values
function prob_to_color(mat, idx)
    m = mat[idx, idx]
    # log-scale: map 0→min, max→max
    ε = 1e-8
    log.(m .+ ε)
end

clims_uni   = let vals = vcat([prob_to_color(snaps_uni[t],   idx) for t in T_PLOT]...)
                  (minimum(vals[isfinite.(vals)]), maximum(vals))
              end
clims_front = let vals = vcat([prob_to_color(snaps_front[t], idx) for t in T_PLOT]...)
                  (minimum(vals[isfinite.(vals)]), maximum(vals))
              end

for (row, snaps, clims, row_label) in [
        (1, snaps_uni,   clims_uni,   "uniform IC  ($n_lo, $n_lo)"),
        (2, snaps_front, clims_front, "front IC  ($n_lo, $n_hi)"),
    ]

    for (col, t) in enumerate(T_PLOT)
        ax = Axis(fig[row, col];
            xlabel = col == 1 ? "n₁" : "",
            ylabel = row  == 1 ? "n₂" : "",
            title  = "t = $t",
            aspect = 1,
            xticks = [0, 50, 100, 150],
            yticks = [0, 50, 100, 150],
        )

        mat_log = prob_to_color(snaps[t], idx)
        heatmap!(ax, 0:n_zoom, 0:n_zoom, mat_log;
                 colorrange = clims, colormap = :inferno)

        # Constraint lines: n_1 + n_2 = const (iso-total lines, slope -1)
        for nc in [n_lo, n_bar, n_hi + n_lo, 2*n_hi]
            nc > 0 && nc <= 2*n_zoom || continue
            x0, x1 = max(0, nc - n_zoom), min(nc, n_zoom)
            y0, y1 = nc - x0, nc - x1
            lines!(ax, [Float64(x0), Float64(x1)], [Float64(y0), Float64(y1)];
                   color = (:white, 0.35), linewidth = 0.8, linestyle = :dash)
        end

        # Mark fixed points and Binomial peak for front case
        if row == 2 && col == 1
            scatter!(ax, [Float64(n_lo)], [Float64(n_hi)];
                     color = :cyan, markersize = 8, marker = :circle)
            scatter!(ax, [Float64(n_hi)], [Float64(n_lo)];
                     color = :cyan, markersize = 8, marker = :circle)
            scatter!(ax, [Float64(n_bar÷2)], [Float64(n_bar - n_bar÷2)];
                     color = :yellow, markersize = 6, marker = :diamond)
        end

        # Row label on first column only
        if col == 1
            text!(ax, 4, n_zoom - 10; text = row_label,
                  fontsize = 8, color = :white)
        end
    end
end

# Column headers
for (col, t) in enumerate(T_PLOT)
    Label(fig[0, col], "t = $t"; fontsize = 11, tellwidth = false)
end

# Row labels
Label(fig[1, 0], "uniform\nIC (n_lo, n_lo)";
      fontsize = 10, tellheight = false, rotation = π/2)
Label(fig[2, 0], "front\nIC (n_lo, n_hi)";
      fontsize = 10, tellheight = false, rotation = π/2)

# Annotation: diffusion vs reaction directions (on first panel of each row)
text!(contents(fig[1,1])[1], n_lo + 15.0, n_lo - 12.0;
      text = "← diffusion\n    (along diagonal)",
      fontsize = 7, color = (:white, 0.8))
text!(contents(fig[1,1])[1], n_lo - 8.0, n_lo + 20.0;
      text = "reactions\n(across diagonals) ↗",
      fontsize = 7, color = (:white, 0.8))

save_fig("fig_prob_diffusion", fig)
println("\nDone.")
