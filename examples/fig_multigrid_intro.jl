"""
fig_multigrid_intro.jl

Beautiful intro figure for high school research mentorship.
Shows the birth-death RDME wavefront and multigrid level structure.

Layout: 2 rows × 3 columns  (dark background, plasma colourmap)
  Row 1: plasma ⟨n_k⟩ snapshots at t = t1, t2, t3
  Row 2: multigrid level-structure maps at the same three times
          dark = empty | grey = equil | red = wavefront edge |
          orange = L1 pairs | teal = L2 quads | green = L3 blocks

Output: paper/figures/fig_multigrid_intro.png
"""

using DiscStochSim, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const k_b = 0.5;  const k_d = 0.1;  const D = 1.0
const K   = 60;   const N0  = 3;    const dt = 0.4
const N_LEV = 3
const n_max = 14
v_kpp = 2*sqrt(D*(k_b - k_d))
@printf("BD RDME %d×%d  v_KPP=%.3f  MG L=%d\n", K, K, v_kpp, N_LEV)

# snapshot times: wave at ~20%, ~45%, ~70% of grid radius
const T_SNAPS = [5.0, 11.0, 18.0]

# ── Colour palette (light professional) ────────────────────────────────────────
const BG         = RGBf(0.97, 0.97, 0.99)   # figure background (near-white)
const COL_EMPTY  = RGBf(0.91, 0.91, 0.94)   # not-yet-reached — light grey
const COL_EQUIL  = RGBf(0.74, 0.76, 0.84)   # settled interior — slate
const COL_EDGE   = RGBf(0.84, 0.12, 0.18)   # wavefront edge   — crimson
const COL_L1     = RGBf(0.94, 0.58, 0.08)   # L1 pairs         — amber
const COL_L2     = RGBf(0.20, 0.58, 0.86)   # L2 quads         — sky blue
const COL_L3     = RGBf(0.18, 0.70, 0.62)   # L3 blocks        — teal
const COL_UNPAIR = RGBf(0.72, 0.72, 0.40)   # orphan interior  — olive

# Wavefront heatmap: empty → BG, occupied → magma (deep purple → orange → yellow)
const MAGMA = CairoMakie.to_colormap(:magma)
magma_c(t) = let c = MAGMA[clamp(floor(Int, t*(length(MAGMA)-1))+1, 1, length(MAGMA))];
             RGBf(c.r, c.g, c.b) end

@inline function gnb(idx::CartesianIndex{2}, Kx, Ky)
    i, j = Tuple(idx)
    ns = CartesianIndex{2}[]
    i > 1  && push!(ns, CartesianIndex(i-1, j))
    i < Kx && push!(ns, CartesianIndex(i+1, j))
    j > 1  && push!(ns, CartesianIndex(i, j-1))
    j < Ky && push!(ns, CartesianIndex(i, j+1))
    ns
end

# ── Multigrid hierarchy (from bd_spatial_fsp.jl logic) ─────────────────────────
function mg_structure(s, n_lev, row_pass)
    active_keys = collect(keys(s.dists))
    active_set  = Set(active_keys)
    equil_set   = Set(keys(s.equil_dists))
    interior = CartesianIndex{2}[]
    edge     = CartesianIndex{2}[]
    for idx in active_keys
        if all(nb -> nb in active_set || nb in equil_set, gnb(idx, s.K_x, s.K_y))
            push!(interior, idx)
        else
            push!(edge, idx)
        end
    end
    Δ_par  = row_pass ? CartesianIndex(0,1) : CartesianIndex(1,0)
    Δ_perp = row_pass ? CartesianIndex(1,0) : CartesianIndex(0,1)
    Δdir(ℓ) = isodd(ℓ) ? Δ_par : Δ_perp
    pair_delta(ℓ) = CartesianIndex(Tuple(Δdir(ℓ+1)) .* 2^((ℓ-1)÷2))
    level_reps   = [interior]
    level_voxels = [[[v] for v in interior]]
    level_pairs  = Vector{Tuple{Int,Int}}[]
    for ℓ in 1:n_lev
        prev = level_reps[end]
        isempty(prev) && break
        ridx = Dict(r => i for (i, r) in enumerate(prev))
        δ    = pair_delta(ℓ)
        used = Set{Int}()
        new_reps = CartesianIndex{2}[]
        new_vox  = Vector{Vector{CartesianIndex{2}}}()
        pairs_ℓ  = Tuple{Int,Int}[]
        for (i, rep) in enumerate(prev)
            i in used && continue
            j = get(ridx, rep+δ, 0)
            (j == 0 || j in used) && continue
            push!(pairs_ℓ, (i, j)); push!(used, i); push!(used, j)
            push!(new_reps, min(rep, rep+δ))
            push!(new_vox, vcat(level_voxels[end][i], level_voxels[end][j]))
        end
        isempty(pairs_ℓ) && break
        push!(level_pairs,  pairs_ℓ)
        push!(level_reps,   new_reps)
        push!(level_voxels, new_vox)
    end
    edge, interior, level_voxels, level_pairs
end

function build_vox_depth(level_voxels, level_pairs)
    n_built = length(level_pairs)
    vox_depth = Dict{CartesianIndex{2}, Int}()
    for ℓ in 1:n_built
        covered = Set{Int}()
        ℓ < n_built && for (i2, j2) in level_pairs[ℓ+1]; push!(covered, i2); push!(covered, j2); end
        for g in 1:length(level_voxels[ℓ+1])
            g in covered && continue
            for v in level_voxels[ℓ+1][g]
                vox_depth[v] = ℓ
            end
        end
    end
    vox_depth
end

# ── Build one plasma image ──────────────────────────────────────────────────────
mu_peak = Ref(1.0)
function make_plasma(s)
    mg  = mean_grid(s)
    cur = maximum(mg; init=1e-3)
    mu_peak[] = max(mu_peak[], cur * 1.02)
    sc  = max(mu_peak[], 1e-3)
    # empty voxels (n≈0) → BG colour; occupied → magma 0.08..1.0
    [let v = clamp(sqrt(mg[i,j] / sc), 0.0, 1.0)
         v < 0.015 ? BG : magma_c(clamp(0.08 + 0.92*v, 0.0, 1.0))
     end for j in 1:K, i in 1:K]
end

# ── Build one level-structure image ────────────────────────────────────────────
function make_struct(s, row_pass)
    img = fill(COL_EMPTY, K, K)
    for idx in keys(s.equil_dists)
        i, j = Tuple(idx); img[j, i] = COL_EQUIL
    end
    edge, interior, level_voxels, level_pairs = mg_structure(s, N_LEV, row_pass)
    vox_depth = build_vox_depth(level_voxels, level_pairs)
    for idx in edge
        i, j = Tuple(idx); img[j, i] = COL_EDGE
    end
    for idx in interior
        i, j = Tuple(idx)
        lv = get(vox_depth, idx, 0)
        img[j, i] = lv == 1 ? COL_L1 :
                    lv == 2 ? COL_L2 :
                    lv == 3 ? COL_L3 : COL_UNPAIR
    end
    img
end

# ── Run simulation ──────────────────────────────────────────────────────────────
model = BranchingDeathRDME(k_b, k_d, D, 1.0)
state = SpatialFSP(model, K, K;
    n_max, ε_expand=0.02, ε_equil=0.15,
    use_multigrid=true, multigrid_levels=N_LEV, use_block_vcycle=false)

ci, cj = K÷2+1, K÷2+1
for di in -2:2, dj in -2:2
    di^2+dj^2 <= 4 || continue
    1 <= ci+di <= K && 1 <= cj+dj <= K || continue
    set_ic!(state, CartesianIndex(ci+di, cj+dj), N0)
end

snap_plasma = Dict{Float64, Matrix{RGBf}}()
snap_struct = Dict{Float64, Matrix{RGBf}}()
remaining   = Set(T_SNAPS)

# capture t=0 only if in snaps
0.0 in T_SNAPS && (snap_plasma[0.0] = make_plasma(state);
                   snap_struct[0.0] = make_struct(state, true))

@printf("Simulating to t=%.1f …\n", maximum(T_SNAPS)); flush(stdout)
t0_wall = time()
row_pass = true
for step_i in 1:round(Int, maximum(T_SNAPS)/dt)
    global row_pass
    step!(state, dt)
    t = round(step_i * dt, digits=6)
    for ts in collect(remaining)
        if t >= ts - 1e-9
            snap_plasma[ts] = make_plasma(state)
            snap_struct[ts] = make_struct(state, row_pass)
            delete!(remaining, ts)
            @printf("  t=%.1f  active=%d  equil=%d\n",
                    ts, n_active(state), n_equilibrated(state))
            flush(stdout)
        end
    end
    isempty(remaining) && break
    row_pass = !row_pass
end
@printf("Done in %.1fs\n", time()-t0_wall)

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(backgroundcolor=:white, textcolor=:gray20,
    Axis=(backgroundcolor=:white,
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false,
          leftspinevisible=false, rightspinevisible=false,
          topspinevisible=false, bottomspinevisible=false)))

fig = Figure(size=(1100, 780))
hide_ax(ax) = (hidedecorations!(ax); hidespines!(ax);
               ax.backgroundcolor[] = RGBAf(BG.r, BG.g, BG.b, 1f0))

t_labels = [@sprintf("t = %.0f", t) for t in T_SNAPS]

# ── Overall title ──────────────────────────────────────────────────────────────
Label(fig[1, 1:3];
      text="Adaptive Multigrid FSP: tracking a stochastic birth–death wave",
      fontsize=16, color=RGBf(0.15,0.15,0.20), font=:bold,
      halign=:center, tellwidth=false)

# ── Row A: wavefront heatmaps ─────────────────────────────────────────────────
Label(fig[2, 0]; text="mean particle count  ⟨n⟩",
      fontsize=11, color=RGBf(0.35,0.35,0.45),
      rotation=π/2, tellheight=false, tellwidth=false)
for (c, ts) in enumerate(T_SNAPS)
    ax = Axis(fig[2, c]; title=t_labels[c],
              titlecolor=RGBf(0.25,0.25,0.35), titlesize=13, yreversed=true)
    hide_ax(ax)
    image!(ax, 0.5..K+0.5, 0.5..K+0.5, snap_plasma[ts])
end

# ── Row B: multigrid structure ─────────────────────────────────────────────────
Label(fig[3, 0]; text="multigrid level structure",
      fontsize=11, color=RGBf(0.35,0.35,0.45),
      rotation=π/2, tellheight=false, tellwidth=false)
for (c, ts) in enumerate(T_SNAPS)
    ax = Axis(fig[3, c]; yreversed=true)
    hide_ax(ax)
    image!(ax, 0.5..K+0.5, 0.5..K+0.5, snap_struct[ts])
end

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_items = [
    (COL_EMPTY, "not yet reached"),
    (COL_EQUIL, "settled  (Poisson approx.)"),
    (COL_L3,    "L3 blocks  (8 voxels)"),
    (COL_L2,    "L2 blocks  (4 voxels)"),
    (COL_L1,    "L1 pairs  (2 voxels)"),
    (COL_EDGE,  "wavefront  (per-voxel exact)"),
]
legend_els  = [PolyElement(color=c, strokecolor=:gray70, strokewidth=0.5)
               for (c,_) in legend_items]
legend_labs = [l for (_,l) in legend_items]
Legend(fig[4, 1:3], legend_els, legend_labs;
       orientation=:horizontal, framevisible=false,
       labelsize=11, labelcolor=RGBf(0.2,0.2,0.3),
       patchsize=(16, 16), colgap=18, rowgap=2,
       tellwidth=false, halign=:center)

# ── Sizes ──────────────────────────────────────────────────────────────────────
rowsize!(fig.layout, 2, Fixed(300))
rowsize!(fig.layout, 3, Fixed(300))
colsize!(fig.layout, 0, Fixed(22))
rowgap!(fig.layout, 6); colgap!(fig.layout, 6)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out = joinpath(outdir, "fig_multigrid_intro.png")
save(out, fig; px_per_unit=3)
@printf("Saved → %s\n", out)
