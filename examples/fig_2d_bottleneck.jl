"""
2D RDME Coarse-only FSP — Bottleneck Visualization

Two-species A→B reaction-diffusion on an 8×8 spatial grid.
A is produced along the left edge (j=1), diffuses fast, slowly converts to B.
B diffuses slowly → B stays concentrated near the source.

Strategy: solve FSP entirely at the 2×2 coarsest level (CartesianIndex{8}, ~100k states).
Fine-grid means recovered at save times: E[nA_fine[i,j]] = E[nA_coarse[J]] / patch_size.
This is exact under fast-diffusion (d_intra/k_AB = 640 >> 1).

Color scheme: Red channel ∝ E[nA_k], Blue channel ∝ E[nB_k].

Run: julia --project examples/fig_2d_bottleneck.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

# ─── Parameters ───────────────────────────────────────────────────────────────

const K_x    = 8       # fine grid: 8 rows
const K_y    = 8       # fine grid: 8 columns
const D_A    = 1.0     # fast A diffusion — ensures within-2×2-patch equilibrium
const D_B    = 0.05    # slow B diffusion — B stays near source
const k_prod = 0.05    # ∅ → A in source voxels
const k_AB   = 0.1     # A → B (bottleneck)
const k_deg  = 1.0     # B → ∅
const T_END  = 60.0
const N_MAX  = 8       # per species per fine voxel

# Coarsest grid: 2×2 = 4 voxels
const KC_x = K_x ÷ 4
const KC_y = K_y ÷ 4

# CoarseningMaps: 8×8 → 4×4 → 2×2
const cmap1 = CoarseningMap2D(K_x, K_y)          # 64→16 voxels, patch_size=4
const cmap2 = CoarseningMap2D(K_x÷2, K_y÷2)      # 16→4  voxels, patch_size=4

# Composed fine→coarsest map: fine voxel k → coarsest voxel J
const fine_to_coarsest = [cmap2.fine_to_coarse[cmap1.fine_to_coarse[k]] for k in 1:K_x*K_y]

# Source voxels: left column (j=1, all rows), mapped to coarsest level
const SOURCE_FINE     = [(i-1)*K_y + 1 for i in 1:K_x]
const SOURCE_COARSEST = unique(fine_to_coarsest[k] for k in SOURCE_FINE)

model       = BottleneckModel1D(D_A, D_B, k_prod, k_AB, k_deg)
fine_grid   = Grid2D(K_x, K_y, 1.0/K_x, 0)
coarse_grid = Grid2D(KC_x, KC_y, fine_grid.dx * 4, 2)   # 2×2 grid, dx=4×fine

# Scale k_prod by patch_size^2 = 16 (two levels of 4× coarsening)
model_coarse = coarsen_model(coarsen_model(model, 4.0), 4.0)

rates = Float64[]

# Stiffness diagnostics
d_intra = D_A / fine_grid.dx^2
cfl     = 1.0 / (2 * d_intra * N_MAX)

println("=" ^ 65)
println("2D Bottleneck RDME  $(K_x)×$(K_y) fine  →  $(KC_x)×$(KC_y) coarsest")
println("=" ^ 65)
@printf("  D_A = %.2f,  D_B = %.3f,  k_AB = %.2f,  k_deg = %.2f\n", D_A, D_B, k_AB, k_deg)
@printf("  d_intra/k_AB = %.0f×  (multigrid fast-diffusion condition)\n", d_intra / k_AB)
@printf("  Explicit CFL ≈ %.2e;  coarse-only dt = 1.0  (%.0f× larger)\n\n", cfl, 1.0/cfl)
println("  Fine state space:    CartesianIndex{$(2*K_x*K_y)} — intractable")
println("  Coarsest level:      CartesianIndex{$(2*KC_x*KC_y)} — ~100k states")
println("  Source coarse voxels: $SOURCE_COARSEST")
println()

# ─── IC: empty system at coarsest level ──────────────────────────────────────

const N_COARSE = 2 * KC_x * KC_y   # = 8

u0 = CartesianIndex(ntuple(_ -> 0, Val(N_COARSE)))
sp = StateSpace{CartesianIndex{N_COARSE}, Float64}()
add_state!(sp, u0, 1.0)

# ─── Adaptive FSP at coarsest level ──────────────────────────────────────────

dt      = 1.0
n_steps = round(Int, T_END / dt)
save_at = Set([0, 5, 15, 30, 60])

snapshots = Dict{Int, StateSpace{CartesianIndex{N_COARSE}, Float64}}()
0 in save_at && (snapshots[0] = deepcopy(sp))

sys_coarse = build_bottleneck_system_2d(model_coarse, coarse_grid;
                                         source_voxels = SOURCE_COARSEST)
bc_coarse  = s -> all(c -> 0 ≤ c ≤ 2*N_MAX, Tuple(s))

println("Running coarse-only FSP ($n_steps steps × dt=$dt)...")

t_wall = @elapsed for step in 1:n_steps
    t_cur = (step-1) * dt

    expand!(sp, sys_coarse, bc_coarse; depth=1)
    A, = build_generator(sp, sys_coarse, rates, t_cur)
    sp.probs .= expv(Float64(dt), A, sp.probs; m=30)
    renormalize!(sp)
    prune_threshold!(sp, 1e-11)

    t_now = step * dt
    if round(Int, t_now) in save_at
        snapshots[round(Int, t_now)] = deepcopy(sp)
        @printf("  t=%4.0f  |S|=%6d\n", t_now, length(sp))
    end
end
@printf("\n  Total wall time: %.1fs\n\n", t_wall)

# ─── Extract fine-grid means from coarse distribution ────────────────────────
# Under fast diffusion, molecules are uniformly distributed within each coarse
# patch, so E[nA_fine[i,j]] = E[nA_coarse[J]] / patch_size  (patch_size = 16).

function mean_counts_fine(sp::StateSpace{CartesianIndex{N_COARSE}, Float64})
    μA_c = zeros(KC_x * KC_y)
    μB_c = zeros(KC_x * KC_y)
    for (idx, state) in enumerate(sp.states)
        t = Tuple(state); p = sp.probs[idx]
        p > 0.0 || continue
        for J in 1:KC_x*KC_y
            μA_c[J] += p * t[2J-1]
            μB_c[J] += p * t[2J]
        end
    end
    # Broadcast to fine grid
    μA = zeros(K_x, K_y)
    μB = zeros(K_x, K_y)
    for i in 1:K_x, j in 1:K_y
        k = (i-1)*K_y + j
        J = fine_to_coarsest[k]
        μA[i,j] = μA_c[J] / 16.0
        μB[i,j] = μB_c[J] / 16.0
    end
    μA, μB
end

# ─── Visualization ────────────────────────────────────────────────────────────

times_to_show = sort(collect(keys(snapshots)))
n_panels      = length(times_to_show)

all_μA = [mean_counts_fine(snapshots[t])[1] for t in times_to_show]
all_μB = [mean_counts_fine(snapshots[t])[2] for t in times_to_show]
max_A  = max(maximum(maximum(m) for m in all_μA), 1e-6)
max_B  = max(maximum(maximum(m) for m in all_μB), 1e-6)

fig = Figure(size = (160 * n_panels, 220), backgroundcolor = :black)

for (col, t) in enumerate(times_to_show)
    μA, μB = all_μA[col], all_μB[col]

    img = [RGBAf(clamp(μA[i,j]/max_A, 0, 1),
                 0f0,
                 clamp(μB[i,j]/max_B, 0, 1),
                 1f0)
           for j in 1:K_y, i in 1:K_x]

    ax = Axis(fig[1, col];
              title            = "t = $t",
              titlecolor       = :white,
              titlesize        = 13,
              backgroundcolor  = :black,
              aspect           = DataAspect(),
              xgridvisible     = false, ygridvisible     = false,
              xticksvisible    = false, yticksvisible    = false,
              xticklabelsvisible = false, yticklabelsvisible = false,
              leftspinevisible  = false, rightspinevisible = false,
              topspinevisible   = false, bottomspinevisible = false)

    image!(ax, img)
end

Label(fig[2, 1:n_panels],
      "  Red → E[nA]  (A species, max=$(round(max_A, digits=3)))    " *
      "Blue → E[nB]  (B species, max=$(round(max_B, digits=3)))    " *
      "Purple → both present\n" *
      "  Coarse-only FSP: 2×2 grid (CartesianIndex{8}) → fine means by uniform-within-patch approximation",
      color = :white, fontsize = 10)

rowgap!(fig.layout, 4)

outpath = joinpath(@__DIR__, "output", "fig_2d_bottleneck.png")
mkpath(dirname(outpath))
save(outpath, fig; px_per_unit = 4)
println("Saved: $outpath")

display(fig)
println("Done.")
