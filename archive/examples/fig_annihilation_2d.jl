"""
2D A+B→∅ Annihilation — Spatial Diffusion and Reaction Demo
=============================================================

Two species start localized in opposite corners of a 4×4 spatial grid:
  A: concentrated in top-left  voxel (1,1)
  B: concentrated in bottom-right voxel (4,4)

They diffuse toward each other; when they meet, A+B→∅ fires.
The FSP is solved exactly at the 2×2 coarse level (4 coarse voxels),
then means are broadcast back to the 4×4 fine grid for visualization.

Output: animated GIF + final PNG showing spatial heatmaps of E[nA], E[nB],
        and a reaction progress trace E[nA+nB](t).

Run: julia --project examples/fig_annihilation_2d.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

# ─── Parameters ───────────────────────────────────────────────────────────────

const KF_x  = 4        # fine grid rows
const KF_y  = 4        # fine grid cols
const KC_x  = 2        # coarse grid rows  (2×2 coarsening)
const KC_y  = 2        # coarse grid cols
const KF    = KF_x * KF_y   # 16 fine voxels
const KC    = KC_x * KC_y   #  4 coarse voxels
const NDIM  = 2 * KC         #  8 state dimensions at coarse level

const D_A   = 1.0      # diffusion coefficient for A
const D_B   = 1.0      # diffusion coefficient for B
const K_R   = 3.0      # A + B → ∅  bimolecular rate
const N_MOL = 5        # molecules per species at t=0
const N_MAX = N_MOL    # hard cap per species per coarse voxel (for BC)

const T_END    = 3.0
const DT       = 0.05
const N_STEPS  = round(Int, T_END / DT)
const SAVE_EV  = 1     # save every step for animation

# Fine grid: dx = 1/KF_x = 0.25
const DX_FINE   = 1.0 / KF_x
# Coarse grid: each coarse voxel covers a 2×2 patch of fine voxels → dx_c = 2*dx_f
const DX_COARSE = DX_FINE * 2

# Diffusion rates at coarse level: D / dx_c²
const DA_C = D_A / DX_COARSE^2
const DB_C = D_B / DX_COARSE^2

println("=" ^ 60)
println("A+B→∅  Annihilation  $(KF_x)×$(KF_y) fine  →  $(KC_x)×$(KC_y) coarse")
println("=" ^ 60)
@printf("  D_A=D_B=%.1f,  k_r=%.1f,  N_mol=%d per species\n", D_A, K_R, N_MOL)
@printf("  Coarse diffusion rate: %.2f\n", DA_C)
@printf("  Fast-diffusion ratio:  D/dx²/k_r = %.1f\n\n", DA_C / K_R)

# ─── Build coarse-level RDME system ───────────────────────────────────────────
#
# State:  CartesianIndex{NDIM}  = (nA₁, nB₁, nA₂, nB₂, nA₃, nB₃, nA₄, nB₄)
# Coarse voxel layout (row-major, 2×2):
#   k=1 (1,1)  k=2 (1,2)
#   k=3 (2,1)  k=4 (2,2)

_unit(k::Int, ::Val{N}) where {N} =
    CartesianIndex(ntuple(j -> j == k ? 1 : 0, Val(N)))

function build_annihilation_system_2d()
    stoichs      = CartesianIndex{NDIM}[]
    propensities = Function[]

    for k in 1:KC
        let k=k
            iA = 2k-1; iB = 2k
            eA = _unit(iA, Val(NDIM)); eB = _unit(iB, Val(NDIM))

            # A + B → ∅
            push!(stoichs, -eA - eB)
            push!(propensities, (x, rv, t) -> K_R * max(0, Tuple(x)[iA]) * max(0, Tuple(x)[iB]))
        end
    end

    # 4-connected diffusion on 2×2 coarse grid
    for i in 1:KC_x, j in 1:KC_y
        let i=i, j=j
            k  = (i-1)*KC_y + j
            iAk = 2k-1; iBk = 2k
            eAk = _unit(iAk, Val(NDIM)); eBk = _unit(iBk, Val(NDIM))

            # Right neighbor
            if j < KC_y
                k2 = k + 1; iAk2=2k2-1; iBk2=2k2
                eAk2=_unit(iAk2,Val(NDIM)); eBk2=_unit(iBk2,Val(NDIM))
                push!(stoichs, -eAk+eAk2); push!(propensities, (x,rv,t)->DA_C*max(0,Tuple(x)[iAk]))
                push!(stoichs,  eAk-eAk2); push!(propensities, (x,rv,t)->DA_C*max(0,Tuple(x)[iAk2]))
                push!(stoichs, -eBk+eBk2); push!(propensities, (x,rv,t)->DB_C*max(0,Tuple(x)[iBk]))
                push!(stoichs,  eBk-eBk2); push!(propensities, (x,rv,t)->DB_C*max(0,Tuple(x)[iBk2]))
            end
            # Down neighbor
            if i < KC_x
                k2 = k + KC_y; iAk2=2k2-1; iBk2=2k2
                eAk2=_unit(iAk2,Val(NDIM)); eBk2=_unit(iBk2,Val(NDIM))
                push!(stoichs, -eAk+eAk2); push!(propensities, (x,rv,t)->DA_C*max(0,Tuple(x)[iAk]))
                push!(stoichs,  eAk-eAk2); push!(propensities, (x,rv,t)->DA_C*max(0,Tuple(x)[iAk2]))
                push!(stoichs, -eBk+eBk2); push!(propensities, (x,rv,t)->DB_C*max(0,Tuple(x)[iBk]))
                push!(stoichs,  eBk-eBk2); push!(propensities, (x,rv,t)->DB_C*max(0,Tuple(x)[iBk2]))
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{NDIM}}(stoichs, propensities)
end

# Boundary condition: all counts in [0, N_MAX]
annihilation_bc(s::CartesianIndex{NDIM}) = all(c -> 0 ≤ c ≤ N_MAX, Tuple(s))

# ─── Initial condition ────────────────────────────────────────────────────────
# A molecules (N_MOL) in coarse voxel 1 (top-left = index k=1)
# B molecules (N_MOL) in coarse voxel 4 (bottom-right = index k=4)
# State layout: (nA1, nB1, nA2, nB2, nA3, nB3, nA4, nB4)
#                A in k=1 →  nA1=N_MOL, all others 0
#                B in k=4 →  nB4=N_MOL, all others 0

u0 = CartesianIndex(N_MOL, 0, 0, 0, 0, 0, 0, N_MOL)
sp = StateSpace{CartesianIndex{NDIM}, Float64}()
add_state!(sp, u0, 1.0)

sys = build_annihilation_system_2d()
rates = Float64[]

# ─── Extract coarse-level means ───────────────────────────────────────────────

function coarse_means(sp::StateSpace{CartesianIndex{NDIM}, Float64})
    μA = zeros(KC_x, KC_y)
    μB = zeros(KC_x, KC_y)
    for (idx, state) in enumerate(sp.states)
        t = Tuple(state); p = sp.probs[idx]
        p > 0.0 || continue
        for k in 1:KC
            i = (k-1) ÷ KC_y + 1
            j = (k-1) % KC_y + 1
            μA[i,j] += p * t[2k-1]
            μB[i,j] += p * t[2k]
        end
    end
    μA, μB
end

function fine_means(μA_c::Matrix, μB_c::Matrix)
    # Each coarse voxel (i,j) covers fine voxels (2i-1:2i, 2j-1:2j)
    # Under fast diffusion: uniform within patch → mean per fine voxel = coarse mean / 4
    patch_size = (KF_x ÷ KC_x) * (KF_y ÷ KC_y)   # = 4
    μA = zeros(KF_x, KF_y)
    μB = zeros(KF_x, KF_y)
    for ci in 1:KC_x, cj in 1:KC_y
        for fi in (2ci-1):2ci, fj in (2cj-1):2cj
            μA[fi,fj] = μA_c[ci,cj] / patch_size
            μB[fi,fj] = μB_c[ci,cj] / patch_size
        end
    end
    μA, μB
end

total_molecules(sp) = begin
    tot = 0.0
    for (idx, state) in enumerate(sp.states)
        t = Tuple(state); p = sp.probs[idx]
        p > 0.0 || continue
        tot += p * sum(t)
    end
    tot
end

# ─── Run FSP ──────────────────────────────────────────────────────────────────

println("Running coarse FSP ($N_STEPS steps × dt=$DT)...")

snapshots = Tuple{Float64, Matrix{Float64}, Matrix{Float64}}[]
times_total = Float64[]
totals      = Float64[]

push!(snapshots, (0.0, fine_means(coarse_means(sp)...)...))
push!(times_total, 0.0); push!(totals, total_molecules(sp))

t_wall = @elapsed for step in 1:N_STEPS
    expand!(sp, sys, annihilation_bc; depth=1)
    A_mat, = build_generator(sp, sys, rates, (step-1)*DT)
    sp.probs .= expv(Float64(DT), A_mat, sp.probs; m=30)
    renormalize!(sp)
    prune_threshold!(sp, 1e-12)

    t_now = step * DT
    μA_c, μB_c = coarse_means(sp)
    μA_f, μB_f = fine_means(μA_c, μB_c)
    push!(snapshots, (t_now, μA_f, μB_f))
    push!(times_total, t_now); push!(totals, total_molecules(sp))

    step % 10 == 0 && @printf("  t=%.2f  |S|=%5d  E[nA+nB]=%.2f\n",
                               t_now, length(sp), totals[end])
end
@printf("\nWall time: %.2fs  Final |S|=%d\n\n", t_wall, length(sp))

# ─── Visualization ────────────────────────────────────────────────────────────

max_A = max(maximum(maximum(s[2]) for s in snapshots), 1e-9)
max_B = max(maximum(maximum(s[3]) for s in snapshots), 1e-9)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

# --- Animated GIF ---
gif_path = joinpath(outdir, "annihilation_2d.gif")
println("Rendering animation → $gif_path")

fig = Figure(size=(820, 340), backgroundcolor=:black)

ax_A   = Axis(fig[1,1]; title="E[nA]",   titlecolor=:white, titlesize=14,
              aspect=DataAspect(), backgroundcolor=:black,
              xgridvisible=false, ygridvisible=false,
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
ax_B   = Axis(fig[1,2]; title="E[nB]",   titlecolor=:white, titlesize=14,
              aspect=DataAspect(), backgroundcolor=:black,
              xgridvisible=false, ygridvisible=false,
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
ax_tot = Axis(fig[1,3]; title="Reaction progress",
              xlabel="time", ylabel="E[nA+nB]",
              titlecolor=:white, xlabelcolor=:white, ylabelcolor=:white,
              xticklabelcolor=:white, yticklabelcolor=:white,
              titlesize=13, xlabelsize=11, ylabelsize=11,
              backgroundcolor=:black,
              xgridcolor=RGBAf(1,1,1,0.1), ygridcolor=RGBAf(1,1,1,0.1))

colsize!(fig.layout, 1, Relative(0.3))
colsize!(fig.layout, 2, Relative(0.3))
colsize!(fig.layout, 3, Relative(0.4))

# Observables
t_obs  = Observable(0.0)
μA_obs = Observable(snapshots[1][2])
μB_obs = Observable(snapshots[1][3])
ttrace = Observable(Float64[0.0])
ntrace = Observable(Float64[totals[1]])

heatmap!(ax_A, μA_obs; colormap=:reds,   colorrange=(0, max_A/4))
heatmap!(ax_B, μB_obs; colormap=:blues,  colorrange=(0, max_B/4))
lines!(ax_tot, ttrace, ntrace; color=:white, linewidth=2)
xlims!(ax_tot, 0, T_END)
ylims!(ax_tot, 0, 2*N_MOL + 0.5)

title_label = Label(fig[0, 1:3],
    @lift("A+B→∅ Annihilation   t = $(round($t_obs, digits=2))"),
    color=:white, fontsize=15)

record(fig, gif_path, 1:length(snapshots); framerate=20) do frame
    t_now, μA_f, μB_f = snapshots[frame]
    t_obs[]  = t_now
    μA_obs[] = μA_f
    μB_obs[] = μB_f
    ttrace[] = times_total[1:frame]
    ntrace[] = totals[1:frame]
end
println("Saved: $gif_path")

# --- Static 4-panel figure at key times ---
key_times_target = [0.0, 0.5, 1.5, 3.0]
key_frames = [argmin(abs(s[1] - tt) for s in snapshots) for tt in key_times_target]

fig2 = Figure(size=(900, 420), backgroundcolor=:black)

for (col, frame) in enumerate(key_frames)
    t_now, μA_f, μB_f = snapshots[frame]

    # Overlay: red=A, blue=B, purple=both
    img = [RGBAf(clamp(μA_f[i,j] / (max_A/4), 0, 1),
                 0f0,
                 clamp(μB_f[i,j] / (max_B/4), 0, 1),
                 1f0)
           for j in 1:KF_y, i in 1:KF_x]

    ax = Axis(fig2[1, col];
              title = @sprintf("t = %.1f", t_now),
              titlecolor=:white, titlesize=13,
              aspect=DataAspect(), backgroundcolor=:black,
              xgridvisible=false, ygridvisible=false,
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false,
              leftspinevisible=false, rightspinevisible=false,
              topspinevisible=false, bottomspinevisible=false)
    image!(ax, img)
end

# Reaction progress in last column
ax_p = Axis(fig2[1, 5];
            title="E[nA+nB](t)", titlecolor=:white, titlesize=13,
            xlabel="time", ylabel="total molecules",
            xlabelcolor=:white, ylabelcolor=:white,
            xticklabelcolor=:white, yticklabelcolor=:white,
            backgroundcolor=:black,
            xgridcolor=RGBAf(1,1,1,0.1), ygridcolor=RGBAf(1,1,1,0.1))
lines!(ax_p, times_total, totals; color=:white, linewidth=2)
for tt in key_times_target
    vlines!(ax_p, [tt]; color=:yellow, linewidth=1, linestyle=:dash)
end
xlims!(ax_p, 0, T_END); ylims!(ax_p, 0, 2*N_MOL + 0.5)

Label(fig2[2, 1:5],
      "Red = E[nA]   Blue = E[nB]   Purple = overlap   " *
      "4×4 fine grid → 2×2 coarse FSP (CartesianIndex{8})\n" *
      "A starts top-left, B starts bottom-right, A+B→∅ on contact",
      color=:white, fontsize=10)

png_path = joinpath(outdir, "annihilation_2d_panels.png")
save(png_path, fig2; px_per_unit=3)
println("Saved: $png_path")
println("Done.")
