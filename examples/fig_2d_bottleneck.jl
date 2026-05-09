"""
2D RDME Multigrid FSP — Bottleneck Visualization

Two-species A→B reaction-diffusion on an 8×8 spatial grid.
A is produced along the left edge (j=1), diffuses fast, slowly converts to B.
B diffuses slowly → B stays concentrated near the source.

Key result: Direct FSP on 8×8 (CartesianIndex{128}) is intractable.
Multigrid coarsening 8×8 → 4×4 → 2×2 reduces to CartesianIndex{8} (~100k states).

Color scheme: Red channel ∝ E[nA_k], Blue channel ∝ E[nB_k].

Run: julia --project examples/fig_2d_bottleneck.jl
"""

using DiscStochSim
using CairoMakie
using Printf

# ─── Parameters ───────────────────────────────────────────────────────────────

const K_x   = 8       # fine grid: 8 rows
const K_y   = 8       # fine grid: 8 columns
const D_A   = 1.0     # fast A diffusion — ensures within-2×2-patch equilibrium
const D_B   = 0.05    # slow B diffusion — B stays near source
const k_prod = 0.05   # ∅ → A in source voxels
const k_AB  = 0.1     # A → B (bottleneck: 10× slower than intra-patch diffusion)
const k_deg = 1.0     # B → ∅
const T_END = 60.0
const N_MAX = 8       # per species per fine voxel (conservative)

# Source: left column (j=1, all rows i)
const SOURCE = [(i-1)*K_y + 1 for i in 1:K_x]

# Analytical steady-state means at source voxels (rough, ignoring diffusion)
const μ_A_src = k_prod / k_AB   # = 0.5 per source fine voxel
const μ_B_src = k_prod / k_deg  # = 0.05 per source fine voxel

model     = BottleneckModel1D(D_A, D_B, k_prod, k_AB, k_deg)
fine_grid = Grid2D(K_x, K_y, 1.0/K_x, 0)
rates     = Float64[]

# CoarseningMaps: 8×8 → 4×4 → 2×2
cmap1 = CoarseningMap2D(K_x, K_y)             # level 0→1: 64→16 voxels, patch_size=4
cmap2 = CoarseningMap2D(K_x÷2, K_y÷2)        # level 1→2: 16→4 voxels,  patch_size=4

grids = build_hierarchy2d(fine_grid, 2)       # [fine, 4×4, 2×2]
cmaps = [cmap1, cmap2]

# Stiffness diagnostics
d_intra = D_A / fine_grid.dx^2
d_inter = D_A / (2*fine_grid.dx)^2
cfl     = 1.0 / (2 * d_intra * N_MAX)

println("=" ^ 65)
println("2D Bottleneck RDME  $(K_x)×$(K_y) fine  →  $(K_x÷4)×$(K_y÷4) coarsest")
println("=" ^ 65)
@printf("  D_A = %.2f,  D_B = %.3f,  k_AB = %.2f,  k_deg = %.2f\n", D_A, D_B, k_AB, k_deg)
@printf("  d_intra/k_AB = %.0f×  (multigrid condition)\n", d_intra / k_AB)
@printf("  Explicit CFL ≈ %.2e;  V-cycle dt = 1.0  (%.0f× larger)\n\n", cfl, 1.0/cfl)
println("  State space:  fine CartesianIndex{$(2*K_x*K_y)} → coarsest CartesianIndex{$(2*(K_x÷4)*(K_y÷4))}")
println("  Direct FSP on fine grid: intractable")
println("  V-cycle solves at 2×2 coarsest (~100k states)")
println()

# ─── IC: all zeros (start from empty) ────────────────────────────────────────

u0 = CartesianIndex(ntuple(_ -> 0, Val(2*K_x*K_y)))

# ─── Run multigrid ────────────────────────────────────────────────────────────

dt      = 1.0
n_steps = round(Int, T_END / dt)
save_at = [0, 5, 15, 30, 60]   # save snapshots at these timesteps

println("Running multigrid V-cycle ($n_steps steps × dt=$dt)...")

sp = StateSpace{CartesianIndex{2*K_x*K_y}, Float64}()
add_state!(sp, u0, 1.0)

snapshots = Dict{Int, StateSpace{CartesianIndex{2*K_x*K_y}, Float64}}()
0 in save_at && (snapshots[0] = deepcopy(sp))

t_wall = @elapsed for step in 1:n_steps
    global sp
    t_cur = (step-1) * dt
    sp = multi_level_vcycle_2d_2s(sp, model, grids, cmaps, rates, t_cur, dt;
                                   τ_pre          = 0.0,
                                   coarse_n_max   = 2 * N_MAX,
                                   coarse_r_depth = 2,
                                   coarse_d_depth = 1,
                                   expand_coarse  = true,
                                   binom_tol      = 1e-8,
                                   weight_tol     = 1e-12,
                                   prune_tol      = 1e-11,
                                   source_voxels  = SOURCE)
    renormalize!(sp)
    prune_threshold!(sp, 1e-11)
    t_now = step * dt
    if round(Int, t_now) in save_at
        snapshots[round(Int, t_now)] = deepcopy(sp)
        @printf("  t=%4.0f  |S|=%6d\n", t_now, length(sp))
    end
end
@printf("\n  Total wall time: %.1fs\n\n", t_wall)

# ─── Extract mean counts per fine voxel ───────────────────────────────────────

function mean_counts_2d(sp::StateSpace{CartesianIndex{N}, Float64},
                         K_x::Int, K_y::Int) where {N}
    μA = zeros(K_x, K_y)
    μB = zeros(K_x, K_y)
    for (idx, state) in enumerate(sp.states)
        t = Tuple(state); p = sp.probs[idx]
        p > 0.0 || continue
        for i in 1:K_x, j in 1:K_y
            k = (i-1)*K_y + j
            μA[i,j] += p * t[2k-1]
            μB[i,j] += p * t[2k]
        end
    end
    μA, μB
end

# ─── Visualization ────────────────────────────────────────────────────────────

times_to_show = sort(collect(keys(snapshots)))
n_panels = length(times_to_show)

fig = Figure(size = (160 * n_panels, 220), backgroundcolor = :black)

# Determine global color scale across all snapshots
all_μA = [mean_counts_2d(snapshots[t], K_x, K_y)[1] for t in times_to_show]
all_μB = [mean_counts_2d(snapshots[t], K_x, K_y)[2] for t in times_to_show]
max_A  = max(maximum(maximum(m) for m in all_μA), 1e-6)
max_B  = max(maximum(maximum(m) for m in all_μB), 1e-6)

for (col, t) in enumerate(times_to_show)
    μA, μB = all_μA[col], all_μB[col]

    # Build RGB image: R=A, G=0, B=B (black background)
    img = [RGBAf(clamp(μA[i,j]/max_A, 0, 1),
                 0f0,
                 clamp(μB[i,j]/max_B, 0, 1),
                 1f0)
           for j in 1:K_y, i in 1:K_x]   # CairoMakie image: [x=col, y=row]

    ax = Axis(fig[1, col];
              title        = "t = $t",
              titlecolor   = :white,
              titlesize    = 13,
              backgroundcolor = :black,
              aspect       = DataAspect(),
              xgridvisible = false, ygridvisible = false,
              xticksvisible = false, yticksvisible = false,
              xticklabelsvisible = false, yticklabelsvisible = false,
              leftspinevisible = false, rightspinevisible = false,
              topspinevisible = false,  bottomspinevisible = false)

    image!(ax, img)
end

# Color legend
Label(fig[2, 1:n_panels],
      "  Red → E[nA]  (A species, max=$(round(max_A, digits=2)))    " *
      "Blue → E[nB]  (B species, max=$(round(max_B, digits=2)))    " *
      "Purple → both present",
      color = :white, fontsize = 11)

rowgap!(fig.layout, 4)

outpath = joinpath(@__DIR__, "output", "fig_2d_bottleneck.png")
mkpath(dirname(outpath))
save(outpath, fig; px_per_unit = 4)
println("Saved: $outpath")

display(fig)
println("Done.")
