"""
2D Schlögl Wavefront — Mean-Field ODE + Active-Strip Visualization
===================================================================

Bistable reaction-diffusion on a K×K grid:
  ∅  →(c3)           X       (production)
  X  →(c4·n)         ∅       (death)
  2X →(c1·n(n-1)/2)  3X      (autocatalysis)
  3X →(c2·n(n-1)(n-2)/6) 2X  (suppression)
  + diffusion D between adjacent voxels (no-flux BC)

Two stable states per voxel: n_lo ≈ 31 and n_hi ≈ 162.

IC: left half at n_hi, right half at n_lo → sharp vertical front.
∫f(n)dn > 0 → n_hi state is thermodynamically favored → front moves right.

Active strip (where FSP tracking is needed):
  n_lo < μ < n_hi  (voxels in the bistable transition zone)
  Outside: unimodal distribution → mean-field sufficient.

Run: julia --project examples/schlogl_2d.jl
"""

using DiscStochSim
using CairoMakie
using LinearAlgebra
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────

const C1  = 0.028
const C2  = 3.2e-4
const C3  = 19.5
const C4  = 1.0
const D   = 5.0

model = SchloglModel1D(D, C1, C2, C3, C4)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Fixed points:  n_lo=%d  n_unstable=%d  n_hi=%d\n", n_lo, n_un, n_hi)

# Per-voxel deterministic rate (used in the mean-field ODE)
f_schlogl(μ) = C1*μ*(μ-1)/2 - C2*μ*(μ-1)*(μ-2)/6 + C3 - C4*μ

# ─────────────────────────────────────────────────────────────────────────────
# 2D mean-field ODE solver
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_schlogl_ode_2d(K_x, K_y, D, T_END, DT, μ0) -> Vector{(t, μ)}

Explicit-Euler mean-field ODE on a K_x × K_y grid with no-flux BC.
dμ[i,j]/dt = f(μ[i,j]) + D · Δ²μ[i,j]
where Δ² is the 4-connected discrete Laplacian (dx = 1 voxel unit).

Sub-steps are chosen automatically for explicit-Euler stability.
"""
function run_schlogl_ode_2d(K_x::Int, K_y::Int, D_val::Float64,
                             T_END::Float64, DT::Float64,
                             μ0::Matrix{Float64};
                             n_snap::Int = 4)
    # Stability: dt_sub < 1 / (4D + |f'|_max)  — use conservative estimate
    d          = D_val   # dx = 1
    max_react  = abs(f_schlogl(0.0)) + 2.0   # rough bound on |f'| near n_lo, n_hi
    dt_stable  = 1.0 / (4*d + max_react + 1.0)
    n_sub      = max(1, ceil(Int, DT / dt_stable))
    dt_sub     = DT / n_sub

    μ  = copy(μ0)
    Δ  = zeros(K_x, K_y)

    n_steps    = round(Int, T_END / DT)
    snap_every = max(1, n_steps ÷ (n_snap - 1))
    snaps      = Tuple{Float64, Matrix{Float64}}[(0.0, copy(μ))]

    for step in 1:n_steps
        for _ in 1:n_sub
            fill!(Δ, 0.0)
            for i in 1:K_x, j in 1:K_y
                Δ[i,j] = f_schlogl(μ[i,j])
                i > 1   && (Δ[i,j] += d*(μ[i-1,j] - μ[i,j]))
                i < K_x && (Δ[i,j] += d*(μ[i+1,j] - μ[i,j]))
                j > 1   && (Δ[i,j] += d*(μ[i,j-1] - μ[i,j]))
                j < K_y && (Δ[i,j] += d*(μ[i,j+1] - μ[i,j]))
            end
            μ .+= Δ .* dt_sub
        end
        if step % snap_every == 0 || step == n_steps
            push!(snaps, (step * DT, copy(μ)))
        end
    end

    # Trim to exactly n_snap snapshots
    snaps[1:min(n_snap, length(snaps))]
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters & IC
# ─────────────────────────────────────────────────────────────────────────────

const K    = 40          # square grid K × K
const T_END = 24.0
const DT    = 0.1
const N_SNAP = 4

# Flat front: bottom half n_hi, top half n_lo → front sweeps upward
# (n_hi state is thermodynamically favored: ∫f(n)dn > 0, so it expands)
μ0 = Float64[j <= K÷2 ? n_hi : n_lo for i in 1:K, j in 1:K]

println("Running 2D Schlögl ODE on $(K)×$(K) grid ...")
@time snaps = run_schlogl_ode_2d(K, K, D, T_END, DT, μ0; n_snap=N_SNAP)
println("Done ($(length(snaps)) snapshots).")

# ─────────────────────────────────────────────────────────────────────────────
# Colormap: blue (n_lo) → white (n_unstable) → red (n_hi)
# ─────────────────────────────────────────────────────────────────────────────

# Clamp display range slightly beyond stable states
const vmin = Float32(n_lo - 15)
const vmax = Float32(n_hi + 15)
const cmap = :RdBu   # blue=low, red=high; we'll reverse so n_hi=red

# ─────────────────────────────────────────────────────────────────────────────
# Compute single-voxel stationary distribution (shows bistability directly)
# ─────────────────────────────────────────────────────────────────────────────

function schlogl_stationary(c1, c2, c3, c4; n_max=250, tol=1e-14)
    # Build tridiagonal generator for a single voxel
    n_states = n_max + 1
    Q = zeros(n_states, n_states)
    for n in 0:n_max
        prod_rate = c3 + c1 * max(0, n*(n-1)) / 2
        deg_rate  = c4 * n + c2 * max(0, n*(n-1)*(n-2)) / 6
        if n < n_max
            Q[n+2, n+1] += prod_rate   # n → n+1
            Q[n+1, n+1] -= prod_rate
        end
        if n > 0
            Q[n,   n+1] += deg_rate    # n → n-1
            Q[n+1, n+1] -= deg_rate
        end
    end
    # Null vector = stationary distribution
    vals, vecs = eigen(Q)
    idx  = argmin(abs.(real.(vals)))
    π    = real.(vecs[:, idx])
    π   .= abs.(π)
    π  ./= sum(π)
    π[π .< tol] .= 0.0
    π ./= sum(π)
    0:n_max, π
end

println("Computing single-voxel stationary distribution ...")
ns_stat, π_stat = schlogl_stationary(C1, C2, C3, C4; n_max=250)

# ─────────────────────────────────────────────────────────────────────────────
# Figure: 4 time snapshots + bistability panel
# ─────────────────────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

fig = Figure(size = (1600, 460), backgroundcolor = :black)

axes = Axis[]
for (ci, (t_s, μ_s)) in enumerate(snaps)
    ax = Axis(fig[1, ci];
              title           = @sprintf("t = %.0f", t_s),
              titlecolor      = :white,
              titlesize       = 14,
              aspect          = DataAspect(),
              backgroundcolor = :black,
              xgridvisible    = false,
              ygridvisible    = false,
              xticksvisible   = false,
              yticksvisible   = false,
              xticklabelsvisible = false,
              yticklabelsvisible = false,
              leftspinevisible   = false,
              rightspinevisible  = false,
              topspinevisible    = false,
              bottomspinevisible = false)

    hm = heatmap!(ax, 1:K, 1:K, μ_s;
                  colormap   = Reverse(:RdBu),
                  colorrange = (vmin, vmax))

    # Active-strip overlay: where n_lo + margin < μ < n_hi - margin
    # These voxels need FSP tracking; the rest can use mean-field.
    margin  = 20
    active  = [Float64(n_lo + margin < μ_s[i,j] < n_hi - margin)
               for i in 1:K, j in 1:K]
    heatmap!(ax, 1:K, 1:K, active;
             colormap   = [RGBAf(0,0,0,0), RGBAf(1,1,1,0.35)],
             colorrange = (0,1))

    push!(axes, ax)
end

# Shared colorbar
Colorbar(fig[1, N_SNAP+1];
         colormap   = Reverse(:RdBu),
         limits     = (vmin, vmax),
         label      = "E[nₖ]",
         labelcolor = :white,
         ticklabelcolor = :white,
         tickcolor  = :white,
         height     = Relative(0.85),
         tellheight = false)

# ── Bistability panel: single-voxel stationary distribution P(n) ─────────────
ax_bi = Axis(fig[1, N_SNAP+2];
             title      = "P(n) — single voxel at SS",
             xlabel     = "molecule count n",
             ylabel     = "probability",
             titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
             xticklabelcolor = :white, yticklabelcolor = :white,
             backgroundcolor = :black,
             xgridcolor = RGBAf(1,1,1,0.08), ygridcolor = RGBAf(1,1,1,0.08),
             titlesize = 11, xlabelsize = 10, ylabelsize = 10)

# Fill under peaks with species colors
band!(ax_bi, Float64.(ns_stat), zeros(length(ns_stat)), Float64.(π_stat);
      color = (:gray, 0.4))
lines!(ax_bi, Float64.(ns_stat), Float64.(π_stat); color = :white, linewidth = 1.5)

# Mark the stable states and unstable point
vlines!(ax_bi, [Float64(n_lo)]; color = :royalblue, linewidth = 2,
        linestyle = :dash, label = "n_lo=$(n_lo)")
vlines!(ax_bi, [Float64(n_un)]; color = :yellow,    linewidth = 1.5,
        linestyle = :dash, label = "n_un=$(n_un)")
vlines!(ax_bi, [Float64(n_hi)]; color = :crimson,   linewidth = 2,
        linestyle = :dash, label = "n_hi=$(n_hi)")

# Shade the two modes
peak_lo = ns_stat[argmax(π_stat[1:n_un])]
peak_hi = ns_stat[n_un+1 + argmax(π_stat[n_un+2:end])]
band!(ax_bi, Float64.(ns_stat[1:n_un+1]), zeros(n_un+1),
      Float64.(π_stat[1:n_un+1]); color = (:royalblue, 0.25))
band!(ax_bi, Float64.(ns_stat[n_un+2:end]), zeros(length(ns_stat)-n_un-1),
      Float64.(π_stat[n_un+2:end]); color = (:crimson, 0.25))

axislegend(ax_bi; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 9, position = :rt)
xlims!(ax_bi, 0, 220)
text!(ax_bi, n_lo - 5,  maximum(π_stat)*0.7;  text="bimodal!",
      color=:white, fontsize=10, align=(:right,:center))

Label(fig[0, 1:N_SNAP+2],
      "2D Schlögl wavefront  |  n_lo=$(n_lo)  n_unstable=$(n_un)  n_hi=$(n_hi)  D=$(D)  " *
      "Blue = n_lo state ($(n_lo) molecules) · Red = n_hi state ($(n_hi) molecules) · White overlay = FSP active strip",
      color    = :white,
      fontsize = 12)

Label(fig[2, 1:N_SNAP+2],
      "IC: bottom half n_hi=$(n_hi), top half n_lo=$(n_lo).  " *
      "n_hi state thermodynamically favored (∫f(n)dn > 0): front sweeps upward.  " *
      "Front width ≈ √(D/|f′(n_unstable)|) ≈ $(round(sqrt(D/0.22), digits=1)) voxels  " *
      "→ narrow strip of K_act ≈ K × width voxels tracked by FSP; rest mean-field.",
      color    = :white,
      fontsize = 10)

png_path = joinpath(outdir, "schlogl_2d.png")
save(png_path, fig; px_per_unit = 3)
println("Saved: $png_path")

# ─────────────────────────────────────────────────────────────────────────────
# Supplementary: 1D profile through the center row at each snapshot
# ─────────────────────────────────────────────────────────────────────────────

fig2 = Figure(size = (900, 380), backgroundcolor = :black)

ax_prof = Axis(fig2[1,1];
               title      = "Mean-field profile E[nₖ] through center row",
               xlabel     = "voxel j", ylabel     = "E[nⱼ]",
               titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
               xticklabelcolor = :white, yticklabelcolor = :white,
               backgroundcolor = :black,
               xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))

snap_colors = [:white, :cyan, :orange, :red]
for (ci, (t_s, μ_s)) in enumerate(snaps)
    center_row = μ_s[K÷2, :]
    lines!(ax_prof, 1:K, center_row;
           color = snap_colors[ci], linewidth = 2,
           label = "t=$(round(t_s, digits=0))")
end
hlines!(ax_prof, [Float64(n_lo)]; color=:blue,  linewidth=1, linestyle=:dash, label="n_lo")
hlines!(ax_prof, [Float64(n_un)]; color=:yellow, linewidth=1, linestyle=:dash, label="n_unstable")
hlines!(ax_prof, [Float64(n_hi)]; color=:red,   linewidth=1, linestyle=:dash, label="n_hi")
# Shade the active strip
band!(ax_prof, 1:K,
      fill(Float64(n_lo + 20), K), fill(Float64(n_hi - 20), K);
      color = RGBAf(1,1,1,0.10))

axislegend(ax_prof; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 11)
ylims!(ax_prof, n_lo - 20, n_hi + 20)

ax_size = Axis(fig2[1,2];
               title      = "Active-strip voxel count vs time",
               xlabel     = "time", ylabel     = "active voxels K_act",
               titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
               xticklabelcolor = :white, yticklabelcolor = :white,
               backgroundcolor = :black,
               xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))

# Track K_act over all time steps (re-run lightweight)
n_steps_full = round(Int, T_END / DT)
μ_track = copy(μ0)
times_all = Float64[]
kact_all  = Int[]
Δ_track   = zeros(K, K)
d_val = D
dt_sub_track = 0.1 / max(1, ceil(Int, 0.1 / (1/(4*D+3))))

for step in 0:n_steps_full
    active_n = sum(n_lo + 20 < μ_track[i,j] < n_hi - 20
                   for i in 1:K, j in 1:K)
    push!(times_all, step * DT)
    push!(kact_all, active_n)
    step == n_steps_full && break
    # Advance one DT
    n_sub_t = max(1, ceil(Int, DT / (1/(4*d_val+4))))
    dt_st   = DT / n_sub_t
    for _ in 1:n_sub_t
        fill!(Δ_track, 0.0)
        for i in 1:K, j in 1:K
            Δ_track[i,j] = f_schlogl(μ_track[i,j])
            i>1   && (Δ_track[i,j] += d_val*(μ_track[i-1,j]-μ_track[i,j]))
            i<K   && (Δ_track[i,j] += d_val*(μ_track[i+1,j]-μ_track[i,j]))
            j>1   && (Δ_track[i,j] += d_val*(μ_track[i,j-1]-μ_track[i,j]))
            j<K   && (Δ_track[i,j] += d_val*(μ_track[i,j+1]-μ_track[i,j]))
        end
        μ_track .+= Δ_track .* dt_st
    end
end

lines!(ax_size, times_all, Float64.(kact_all);
       color = :white, linewidth = 2)
hlines!(ax_size, [2.0*K]; color = :cyan, linewidth = 1, linestyle = :dash,
        label = "2K (1 col = 2 active)")
axislegend(ax_size; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 11)
ylims!(ax_size, 0, K*K*0.6)

Label(fig2[0, 1:2],
      "2D Schlögl ($(K)×$(K))  |  n_lo=$(n_lo)  n_hi=$(n_hi)  D=$(D)  " *
      "Shaded band = active FSP strip  |  K_act peak ≈ $(K) rows × front_width voxels",
      color = :white, fontsize = 12)

png_path2 = joinpath(outdir, "schlogl_2d_profile.png")
save(png_path2, fig2; px_per_unit = 3)
println("Saved: $png_path2")
println("Done.")
