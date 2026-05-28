"""
anim_nucleus_cell.jl — VISUALIZATION TOOL (not an FSP simulation)
Demonstrates what the adaptive FSP grid structure would look like for a
nucleus/cytoplasm geometry.  Uses mean-field ODE to drive the phase classification
(valid here because mRNA reactions are linear → mean-field is exact for means).

anim_nucleus_cell.jl — Level 1: mRNA gradient, circular cell + nucleus

Active voxels: Cartesian grid inside circular cell (R_cell=25, ~1976 voxels).
Nucleus = inner circle (R_nuc=10).  Nuclear membrane blocks diffusion except at pore.
Gene at nucleus center → mRNA → exits through pore → diffuses + degrades.

MULTIGRID CRITERION (inverse): merge where max(μ_u, μ_v) < FINE_THRESHOLD.
  t=0 : whole domain empty → 1–2 giant super-voxels  (Kc tiny)
  t>0 : wavefront spreads from pore → voxels split off → Kc grows
  SS  : near-pore / nucleus fine,  far cytoplasm still merged

Panels:
  Left  — concentration heatmap  (:magma)
  Middle — log₂(super-voxel size) spatial map  (white=fine, blue=coarse)
  Right — Kc × Kc block matrix image

Run: julia --project examples/anim_nucleus_cell.jl
"""

using CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const N_BB     = 52
const CX       = 26.5;  const CY = 26.5
const R_CELL   = 25.0
const R_NUC    = 10.0
const PORE_ANG  = 0.0          # east
const PORE_HALF = π / 6        # 30° half-width
const D        = 1.5
const k_b      = 30.0
const k_d      = 0.05          # λ = √(D/k_d) ≈ 5.5 voxels
const T_tot    = 120.0
const N_FRM    = 120
const IMGSIZE   = 250
const ABS_FLOOR = 0.05         # μ < this  →  "empty"   phase (merge)
const EQ_TOL    = 0.10         # |μ−μ_SS|/μ_SS < this  →  "equil" phase (merge)
# "front" phase (everything else) → singletons

# ── Geometry ───────────────────────────────────────────────────────────────────
in_cell(i, j)    = (i - CX)^2 + (j - CY)^2 < R_CELL^2
in_nucleus(i, j) = (i - CX)^2 + (j - CY)^2 < R_NUC^2

const PORE_ROW = CX + R_NUC * sin(PORE_ANG)
const PORE_COL = CY + R_NUC * cos(PORE_ANG)

function is_pore_edge(i1, j1, i2, j2)
    mi = (i1 + i2) / 2 - CX;  mj = (j1 + j2) / 2 - CY
    r  = sqrt(mi^2 + mj^2);   r < 1e-6 && return false
    cosθ = (mi * sin(PORE_ANG) + mj * cos(PORE_ANG)) / r
    cosθ > cos(PORE_HALF)
end

# ── Active voxels ordered by distance from pore ────────────────────────────────
const _IJ_ALL  = [(i, j) for i in 1:N_BB for j in 1:N_BB if in_cell(i, j)]
const _PDIST   = [sqrt((ij[1]-PORE_ROW)^2 + (ij[2]-PORE_COL)^2) for ij in _IJ_ALL]
const _ORDER   = sortperm(_PDIST)
const VOX_IJ   = _IJ_ALL[_ORDER]
const VOX_PDIST = _PDIST[_ORDER]          # pore distance of voxel k
const K        = length(VOX_IJ)
const IJ_TO_K  = Dict(ij => k for (k, ij) in enumerate(VOX_IJ))
const VOX_NUC  = [in_nucleus(ij...) for ij in VOX_IJ]

const GENE_K = let
    best_k = 0;  best_d = Inf
    for k in 1:K
        !VOX_NUC[k] && continue
        di = VOX_IJ[k][1] - CX;  dj = VOX_IJ[k][2] - CY
        d  = di^2 + dj^2
        if d < best_d;  best_d = d;  best_k = k;  end
    end
    best_k
end

# ── Edges ──────────────────────────────────────────────────────────────────────
const EDGE_U = Int[];  const EDGE_V = Int[];  const EDGE_D = Float64[]

for k in 1:K
    i, j = VOX_IJ[k]
    for (i2, j2) in ((i, j+1), (i+1, j))
        !haskey(IJ_TO_K, (i2, j2)) && continue
        l = IJ_TO_K[(i2, j2)]
        crosses = (VOX_NUC[k] != VOX_NUC[l])
        d_rate  = crosses ? (is_pore_edge(i, j, i2, j2) ? D : 0.0) : D
        push!(EDGE_U, k);  push!(EDGE_V, l);  push!(EDGE_D, d_rate)
    end
end

const ACT_MASK = findall(>(0), EDGE_D)
const ACT_U = EDGE_U[ACT_MASK];  const ACT_V = EDGE_V[ACT_MASK];  const ACT_DR = EDGE_D[ACT_MASK]

# ── Neighbour lists ────────────────────────────────────────────────────────────
const NBRS = [Tuple{Int,Float64}[] for _ in 1:K]
for e in eachindex(EDGE_U)
    EDGE_D[e] == 0.0 && continue
    u, v = EDGE_U[e], EDGE_V[e]
    push!(NBRS[u], (v, EDGE_D[e]));  push!(NBRS[v], (u, EDGE_D[e]))
end

const BIRTH_V = [k == GENE_K ? k_b : 0.0 for k in 1:K]
const DEATH_V = [VOX_NUC[k]  ? 0.0 : k_d  for k in 1:K]

# ── Mean-field ODE ─────────────────────────────────────────────────────────────
function mf_step!(mu, dmu, dt)
    for k in 1:K
        flux = 0.0
        for (l, d) in NBRS[k];  flux += d * (mu[l] - mu[k]);  end
        dmu[k] = BIRTH_V[k] - DEATH_V[k] * mu[k] + flux
    end
    mu .+= dt .* dmu
end

function advance_mf!(mu, Δt; dt_inner = 0.05)
    dmu = similar(mu);  t = 0.0
    while t < Δt
        s = min(dt_inner, Δt - t);  mf_step!(mu, dmu, s);  t += s
    end
end

# ── Phase-based flood-fill coarsening ─────────────────────────────────────────
# Each voxel is classified into one of three phases:
#   :empty  μ < ABS_FLOOR                        → not yet reached, merge
#   :equil  |μ − μ_SS| / μ_SS < EQ_TOL           → at steady state, merge
#   :front  everything else                       → wavefront, keep as singleton
#
# BFS expands connected same-phase regions.  Front voxels are always singletons.
# Result: dark-blue blobs (empty + equil), white ring (front singletons).
function flood_coarsen(mu, mu_ss)
    # Phase assignment (1=empty, 2=equil, 3=front)
    ph = Vector{Int8}(undef, K)
    for k in 1:K
        if mu[k] < ABS_FLOOR
            ph[k] = 1
        elseif mu_ss[k] > ABS_FLOOR && abs(mu[k] - mu_ss[k]) / mu_ss[k] < EQ_TOL
            ph[k] = 2
        else
            ph[k] = 3
        end
    end

    labels = zeros(Int, K)
    comp   = 0
    for k in 1:K
        labels[k] != 0 && continue
        comp += 1;  labels[k] = comp
        ph[k] == 3 && continue   # front → singleton, no BFS expansion
        queue = [k];  head = 1
        while head <= length(queue)
            v = queue[head];  head += 1
            for (u, _) in NBRS[v]
                labels[u] != 0 && continue
                ph[u] != ph[v]  && continue   # different phase → boundary
                ph[u] == 3      && continue   # never pull front into a blob
                labels[u] = comp;  push!(queue, u)
            end
        end
    end
    Kc = comp

    sv_sizes = zeros(Int, Kc)
    for k in 1:K;  sv_sizes[labels[k]] += 1;  end

    # Order by mean pore distance (near-pore first → ordered block matrix)
    sv_dist = zeros(Kc)
    for k in 1:K;  sv_dist[labels[k]] += VOX_PDIST[k];  end
    sv_dist ./= sv_sizes

    perm  = sortperm(sv_dist)
    remap = zeros(Int, Kc)
    for (new_id, old_id) in enumerate(perm);  remap[old_id] = new_id;  end

    asgn_f  = [remap[labels[k]] for k in 1:K]
    sv_sz_f = sv_sizes[perm]

    fin_u = Int[];  fin_v = Int[]
    seen  = Set{Tuple{Int,Int}}()
    for e in eachindex(ACT_U)
        gu, gv = asgn_f[ACT_U[e]], asgn_f[ACT_V[e]]
        gu == gv && continue
        key = gu < gv ? (gu, gv) : (gv, gu)
        key ∈ seen && continue
        push!(seen, key);  push!(fin_u, key[1]);  push!(fin_v, key[2])
    end

    asgn_f, Kc, sv_sz_f, fin_u, fin_v
end

# ── Grid builders (N_BB × N_BB, NaN outside cell) ─────────────────────────────
# heatmap!(ax, 1:N_BB, 1:N_BB, grid) needs grid[col, row]
function mu_grid(mu)
    g = fill(NaN32, N_BB, N_BB)
    for k in 1:K;  i, j = VOX_IJ[k];  g[j, i] = Float32(mu[k]);  end
    g
end

function size_grid(asgn, sv_sizes)
    g = fill(NaN32, N_BB, N_BB)
    for k in 1:K
        i, j = VOX_IJ[k]
        g[j, i] = Float32(log2(max(1, sv_sizes[asgn[k]])))
    end
    g
end

# ── Block matrix renderer ──────────────────────────────────────────────────────
const COL_DIAG  = RGBAf(0.27, 0.51, 0.71, 1.0)
const COL_COUPL = RGBAf(0.89, 0.45, 0.13, 1.0)
const COL_ZERO  = RGBAf(0.94, 0.94, 0.94, 1.0)
const COL_GRID  = RGBAf(0.45, 0.45, 0.45, 1.0)

function render_block(Kc, fin_u, fin_v)
    img = fill(COL_ZERO, IMGSIZE, IMGSIZE)
    cp  = IMGSIZE / Kc
    for i in 1:Kc
        r0 = floor(Int,(i-1)*cp)+1;  r1 = min(IMGSIZE, floor(Int,i*cp))
        for r in r0:r1, c in r0:r1;  img[r, c] = COL_DIAG;  end
    end
    for (u, v) in zip(fin_u, fin_v)
        r0u=floor(Int,(u-1)*cp)+1; r1u=min(IMGSIZE,floor(Int,u*cp))
        r0v=floor(Int,(v-1)*cp)+1; r1v=min(IMGSIZE,floor(Int,v*cp))
        for r in r0u:r1u, c in r0v:r1v;  img[r,c]=COL_COUPL;  end
        for r in r0v:r1v, c in r0u:r1u;  img[r,c]=COL_COUPL;  end
    end
    if cp >= 2
        for k in 0:Kc
            px = min(IMGSIZE, floor(Int,k*cp)+1)
            for q in 1:IMGSIZE;  img[px,q]=COL_GRID;  img[q,px]=COL_GRID;  end
        end
    end
    img
end

# ── Circle overlays: x=col, y=row in plot coords ──────────────────────────────
function circle_pts(R, n=500; gap_ang=nothing, gap_half=0.0)
    pts = Point2f[]
    for s in 0:n-1
        θ = 2π*s/n
        if !isnothing(gap_ang)
            Δ = abs(θ - gap_ang);  Δ = min(Δ, 2π-Δ)
            Δ < gap_half && continue
        end
        push!(pts, Point2f(CY + R*cos(θ), CX + R*sin(θ)))
    end
    pts
end
const NUC_PTS  = circle_pts(R_NUC;  gap_ang=PORE_ANG, gap_half=PORE_HALF)
const CELL_PTS = circle_pts(R_CELL)

# ── Pre-compute mean-field frames ──────────────────────────────────────────────
println("K=$K active voxels  (nucleus=$(count(VOX_NUC)), cytoplasm=$(K-count(VOX_NUC)))")
println("λ = $(round(sqrt(D/k_d),digits=2)) voxels   SS total mRNA ≈ $(round(k_b/k_d,digits=0))")

t_frames  = collect(range(0.0, T_tot; length=N_FRM))
mu_frames = [zeros(K) for _ in 1:N_FRM]   # IC = all zero → Kc=1 at frame 1

println("Pre-computing $N_FRM mean-field frames …")
for f in 2:N_FRM
    mu_frames[f] .= mu_frames[f-1]
    advance_mf!(mu_frames[f], t_frames[f] - t_frames[f-1])
    f % 30 == 0 && @printf("  frame %3d  t=%6.1f  max⟨n⟩=%5.2f\n",
        f, t_frames[f], maximum(mu_frames[f]))
end

println("SS max⟨n⟩=$(round(maximum(mu_frames[end]),digits=2))  ABS_FLOOR=$ABS_FLOOR  EQ_TOL=$EQ_TOL")

# ── Pre-compute coarsenings ────────────────────────────────────────────────────
println("Computing coarsenings …")
asgn_frames  = Vector{Vector{Int}}(undef, N_FRM)
kc_frames    = zeros(Int, N_FRM)
szgrid_frames = Vector{Matrix{Float32}}(undef, N_FRM)
mat_frames   = Vector{Matrix{RGBAf}}(undef, N_FRM)

const MAX_SZ_LOG2 = 11.0   # 2^11 = 2048 ≥ K; cap for colormap

for f in 1:N_FRM
    asgn, Kc, sv_sizes, fu, fv = flood_coarsen(mu_frames[f], mu_frames[end])
    asgn_frames[f]   = asgn
    kc_frames[f]     = Kc
    szgrid_frames[f] = size_grid(asgn, sv_sizes)
    mat_frames[f]    = render_block(Kc, fu, fv)
    f % 30 == 0 && @printf("  frame %3d  Kc=%d\n", f, Kc)
end
println("Kc trajectory: $(kc_frames[1:12:end])")

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1450, 560))

MU_MAX = maximum(mu_frames[end]) * 1.05

mu_obs  = Observable(mu_grid(mu_frames[1]))
sz_obs  = Observable(szgrid_frames[1])
mat_obs = Observable(mat_frames[1])

axis_kw = (aspect=DataAspect(),
    xticks=(5:10:N_BB, string.(5:10:N_BB)),
    yticks=(5:10:N_BB, string.(5:10:N_BB)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=7, yticklabelsize=7)

# Panel A — concentration
ax_mu = Axis(fig[1,1]; title="mRNA concentration", axis_kw...)
heatmap!(ax_mu, 1:N_BB, 1:N_BB, mu_obs;
    colormap=:magma, colorrange=(0, MU_MAX), nan_color=(:black,0.0))
scatter!(ax_mu, NUC_PTS;  color=:white, markersize=1.5, marker=:circle)
scatter!(ax_mu, CELL_PTS; color=:white, markersize=1.0, marker=:circle)
gi, gj = VOX_IJ[GENE_K]
scatter!(ax_mu, [gj], [gi]; marker=:star5, markersize=13,
    color=:yellow, strokewidth=0.4, strokecolor=:white)
scatter!(ax_mu, [CY+R_NUC*cos(PORE_ANG)+1.0], [CX+R_NUC*sin(PORE_ANG)];
    marker=:utriangle, markersize=11, color=:cyan, strokewidth=0)
Colorbar(fig[1,1][1,2]; colormap=:magma, colorrange=(0,MU_MAX),
    label="⟨mRNA⟩", width=10, labelsize=8, ticklabelsize=6)

# Panel B — coarsening size map
# Colormap: white (log2=0, singleton=fine) → dark blue (log2 large, coarse)
COARSE_CMAP = cgrad([:white, :royalblue4])
ax_sz = Axis(fig[1,2]; title="Multigrid resolution  (white = fine,  blue = coarse)", axis_kw...)
heatmap!(ax_sz, 1:N_BB, 1:N_BB, sz_obs;
    colormap=COARSE_CMAP, colorrange=(0, MAX_SZ_LOG2), nan_color=(:black,0.0))
scatter!(ax_sz, NUC_PTS;  color=(:white,0.5), markersize=1.5, marker=:circle)
scatter!(ax_sz, CELL_PTS; color=(:white,0.5), markersize=1.0, marker=:circle)
Colorbar(fig[1,2][1,2]; colormap=COARSE_CMAP, colorrange=(0,MAX_SZ_LOG2),
    label="log₂(super-voxel size)", width=10, labelsize=8, ticklabelsize=6,
    ticks=(0:2:10, ["1","4","16","64","256","1024"]))

# Panel C — block matrix
ax_mx = Axis(fig[1,3]; title="RDME block matrix  (K_c super-voxels)",
    aspect=DataAspect(),
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    yreversed=true, xaxisposition=:top, titlesize=11)
image!(ax_mx, 0.5..N_BB+0.5, 0.5..N_BB+0.5, mat_obs)

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))

t_label  = Label(fig[0,1:2], "t = 0.0"; fontsize=13, halign=:left)
kc_label = Label(fig[0,3], "K_c = $K"; fontsize=12, halign=:center)

Legend(fig[2,1:3],
    [PolyElement(color=RGBf(0.27,0.51,0.71)),
     PolyElement(color=RGBf(0.89,0.45,0.13)),
     PolyElement(color=RGBf(0.94,0.94,0.94), strokecolor=:gray60, strokewidth=1),
     MarkerElement(marker=:star5,     color=:yellow, markersize=10),
     MarkerElement(marker=:utriangle, color=:cyan,   markersize=10)],
    ["Same super-voxel (diagonal)",
     "Coupled super-voxels (off-diagonal)",
     "No coupling",
     "Gene (mRNA source)",
     "Nuclear pore"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(16,11))

rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,4)
colgap!(fig.layout,1,8); colgap!(fig.layout,2,8)

# ── Record ─────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "paper", "figures", "anim_nucleus_cell.gif")
println("\nRendering → $out")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]        = mu_grid(mu_frames[f])
    sz_obs[]        = szgrid_frames[f]
    mat_obs[]       = mat_frames[f]
    t_label.text[]  = "t = $(round(t_frames[f], digits=1))"
    kc_label.text[] = "K_c = $(kc_frames[f]) super-voxels"
end
println("Saved → $out")
println("Kc: $(kc_frames[1]) → $(maximum(kc_frames))")
