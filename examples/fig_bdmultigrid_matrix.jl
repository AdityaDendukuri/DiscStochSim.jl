"""
fig_bdmultigrid_matrix.jl — 2D birth-death RDME: adaptive mesh + block matrix structure.

Shows how flux-adaptive voxel pairing removes inter-block diffusion coupling
from the RDME generator matrix, reducing state-space dimensionality.

Run: julia --project examples/fig_bdmultigrid_matrix.jl
"""

using CairoMakie

# ── Parameters ────────────────────────────────────────────────────────────────
const K_x, K_y = 4, 4; const K = K_x * K_y
const D   = 1.0; const k_b = 0.5; const k_d = 0.5   # n_ss = 1.0
const N_ic   = 50    # IC molecules at voxel (1,1)
const T_snap = 0.7   # snapshot time
const dt_ref = 0.25  # reference τ for flux-weight computation

# ── Indexing ──────────────────────────────────────────────────────────────────
vid(i,j) = (i-1)*K_y + j
vij(k)   = ((k-1)÷K_y + 1, mod1(k, K_y))   # → (row i, col j)

# ── Edges ─────────────────────────────────────────────────────────────────────
const fine_edges = let es = Tuple{Int,Int}[]
    for i in 1:K_x, j in 1:K_y
        j < K_y && push!(es, (vid(i,j), vid(i,j+1)))
        i < K_x && push!(es, (vid(i,j), vid(i+1,j)))
    end
    es
end

# ── Mean-field ODE (exact for linear reactions) ───────────────────────────────
function mean_field(T; dt=0.005)
    mu = zeros(K_x, K_y); mu[1,1] = float(N_ic)
    t = 0.0
    while t < T
        s = min(dt, T - t); Δ = similar(mu)
        for i in 1:K_x, j in 1:K_y
            Δ[i,j] = k_b - k_d * mu[i,j]
            i > 1   && (Δ[i,j] += D*(mu[i-1,j] - mu[i,j]))
            i < K_x && (Δ[i,j] += D*(mu[i+1,j] - mu[i,j]))
            j > 1   && (Δ[i,j] += D*(mu[i,j-1] - mu[i,j]))
            j < K_y && (Δ[i,j] += D*(mu[i,j+1] - mu[i,j]))
        end
        mu .+= s .* Δ; t += s
    end
    mu
end

# ── Greedy flux matching ──────────────────────────────────────────────────────
function flux_match(mu_vec, edges)
    w    = [D * min(mu_vec[i], mu_vec[j]) * dt_ref for (i,j) in edges]
    used = fill(false, length(mu_vec))
    pairs = Tuple{Int,Int}[]
    for e in sortperm(w; rev=true)
        i, j = edges[e]
        used[i] || used[j] || (push!(pairs, (i,j)); used[i] = used[j] = true)
    end
    pairs, findall(.!used)
end

# ── Coarsen graph ─────────────────────────────────────────────────────────────
function coarsen_graph(K, edges, pairs, singles)
    groups = vcat([[a,b] for (a,b) in pairs], [[s] for s in singles])
    Kc = length(groups)
    nm = zeros(Int, K)
    for (g,ks) in enumerate(groups), k in ks; nm[k] = g; end
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in edges
        gi, gj = nm[i], nm[j]
        gi ≠ gj && push!(ced, gi < gj ? (gi,gj) : (gj,gi))
    end
    groups, nm, sort(collect(ced))
end

# ── Block adjacency (0=zero, 1=diagonal/reaction, 2=diffusion coupling) ───────
function block_adj(Kv, edges)
    B = zeros(Int, Kv, Kv)
    for k in 1:Kv; B[k,k] = 1; end
    for (i,j) in edges; B[i,j] = 2; B[j,i] = 2; end
    B
end

# ── Compute ───────────────────────────────────────────────────────────────────
println("Computing mean field at t=$T_snap …")
mu2d   = mean_field(T_snap)
muvec  = [mu2d[vij(k)...] for k in 1:K]

pairs, singles = flux_match(muvec, fine_edges)
groups, nmap, cedges = coarsen_graph(K, fine_edges, pairs, singles)
Kc = length(groups)

println("K=$K fine → Kc=$Kc coarse  ($(length(pairs)) pairs, $(length(singles)) singletons)")
for (pi,(a,b)) in enumerate(pairs)
    ia,ja = vij(a); ib,jb = vij(b)
    println("  Pair C$pi: ($ia,$ja)↔($ib,$jb)  μ=$(round(muvec[a],digits=2)), $(round(muvec[b],digits=2))")
end

Bf = block_adj(K,  fine_edges)
Bc = block_adj(Kc, cedges)

ew_raw  = [D * min(muvec[i], muvec[j]) for (i,j) in fine_edges]
ew_norm = ew_raw ./ max(maximum(ew_raw), 1e-9)

println("Fine matrix nonzeros: $(count(Bf .> 0)) / $(K^2)  ",
        "($(count(Bf.==2)) coupling blocks)")
println("Coarse matrix nonzeros: $(count(Bc .> 0)) / $(Kc^2)  ",
        "($(count(Bc.==2)) coupling blocks)")

# ── Pair colors (one per coarse super-voxel) ──────────────────────────────────
PCOLORS = [:royalblue3, :firebrick2, :forestgreen, :darkorchid3,
           :chocolate3, :teal, :goldenrod3, :deeppink3]
pcolor(pi) = PCOLORS[mod1(pi, length(PCOLORS))]

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1020, 700))

# Colormap for block matrix: 0→white, 1→steelblue, 2→darkorange
CMAP_BLK  = cgrad([:white, :steelblue, :darkorange])
CRANGE    = (0, 2)
CMAP_HEAT = :YlOrRd_9
HEAT_MAX  = max(N_ic * 0.65, maximum(muvec))

# ─── helper: draw grid lines ──────────────────────────────────────────────────
function grid_lines!(ax, nx, ny; lw=0.5, alpha=0.25)
    for k in 0:nx; hlines!(ax, [k+0.5]; color=(:black,alpha), linewidth=lw); end
    for k in 0:ny; vlines!(ax, [k+0.5]; color=(:black,alpha), linewidth=lw); end
end

# ─── Panel A: fine spatial grid ──────────────────────────────────────────────
ax_A = Axis(fig[1,1];
    title="Fine grid: $(K_x)×$(K_y) voxels  [⟨n⟩ at t = $T_snap]",
    aspect=DataAspect(), xlabel="column j", ylabel="row i",
    xticks=(1:K_y, string.(1:K_y)), yticks=(1:K_x, string.(1:K_x)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=9, yticklabelsize=9)

hm = heatmap!(ax_A, 1:K_y, 1:K_x, mu2d';
    colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX))

# Edge lines: opacity+width ∝ flux weight
for (e,(vk,vl)) in enumerate(fine_edges)
    ik,jk = vij(vk); il,jl = vij(vl)
    linesegments!(ax_A, [Point2f(jk,ik), Point2f(jl,il)];
        color=(:dodgerblue, 0.15 + 0.75*ew_norm[e]),
        linewidth=1.2 + 5.0*ew_norm[e])
end

grid_lines!(ax_A, K_x, K_y)

for k in 1:K
    ik,jk = vij(k)
    text!(ax_A, Float64(jk), Float64(ik);
        text="$(round(muvec[k],digits=1))",
        align=(:center,:center), fontsize=7.5, color=:black)
end

Colorbar(fig[1,1][1,2], hm; label="⟨n⟩", width=9, labelsize=9, ticklabelsize=7)

# ─── Panel B: fine block matrix ──────────────────────────────────────────────
ax_B = Axis(fig[1,2];
    title="Fine RDME matrix  ($(K)×$(K) blocks)",
    aspect=DataAspect(), xlabel="voxel j", ylabel="voxel i",
    xticks=(1:K, string.(1:K)), yticks=(1:K, string.(1:K)),
    xticklabelsize=6, yticklabelsize=6,
    yreversed=true, xaxisposition=:top, titlesize=11)

heatmap!(ax_B, 1:K, 1:K, Bf'; colormap=CMAP_BLK, colorrange=CRANGE)
for k in 0:K
    hlines!(ax_B, [k+0.5]; color=(:gray60,0.35), linewidth=0.3)
    vlines!(ax_B, [k+0.5]; color=(:gray60,0.35), linewidth=0.3)
end

# Colour-coded row/col bands showing which voxels will be merged
for k in 1:K
    pi = nmap[k]; clr = pcolor(pi)
    # top band (x-axis side): thin colored strip above row k
    poly!(ax_B, Point2f[(k-0.5,0.5),(k+0.5,0.5),(k+0.5,0.0),(k-0.5,0.0)];
          color=(clr, 0.55), strokewidth=0)
    # left band (y-axis side): thin strip left of column k
    poly!(ax_B, Point2f[(0.0,k-0.5),(0.5,k-0.5),(0.5,k+0.5),(0.0,k+0.5)];
          color=(clr, 0.55), strokewidth=0)
end

# ─── Panel C: coarse spatial grid ────────────────────────────────────────────
ax_C = Axis(fig[2,1];
    title="Coarse: $Kc super-voxels  (flux-adaptive pairing)",
    aspect=DataAspect(), xlabel="column j", ylabel="row i",
    xticks=(1:K_y, string.(1:K_y)), yticks=(1:K_x, string.(1:K_x)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=9, yticklabelsize=9)

heatmap!(ax_C, 1:K_y, 1:K_x, mu2d';
    colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX))
grid_lines!(ax_C, K_x, K_y)

# Pair borders (solid) + coarse label
for (pi,(vk,vl)) in enumerate(pairs)
    ik,jk = vij(vk); il,jl = vij(vl); clr = pcolor(pi)
    x1, x2 = min(jk,jl)-0.46, max(jk,jl)+0.46
    y1, y2 = min(ik,il)-0.46, max(ik,il)+0.46
    lines!(ax_C, [x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1]; color=(clr,0.95), linewidth=3.5)
    text!(ax_C, (x1+x2)/2, (y1+y2)/2;
        text="C$pi", align=(:center,:center), fontsize=9, color=clr, font=:bold)
end

# Singleton borders (dashed)
for (si,vk) in enumerate(singles)
    ik,jk = vij(vk); pi = length(pairs)+si; clr = pcolor(pi)
    lines!(ax_C, [jk-0.46,jk+0.46,jk+0.46,jk-0.46,jk-0.46],
                 [ik-0.46,ik-0.46,ik+0.46,ik+0.46,ik-0.46];
           color=(clr,0.95), linewidth=3.0, linestyle=:dash)
    text!(ax_C, Float64(jk), Float64(ik);
        text="C$pi", align=(:center,:center), fontsize=9, color=clr, font=:bold)
end

# ─── Panel D: coarse block matrix ────────────────────────────────────────────
ctick_labels = ["C$k" for k in 1:Kc]
ax_D = Axis(fig[2,2];
    title="Coarse RDME matrix  ($Kc×$Kc blocks)",
    aspect=DataAspect(), xlabel="super-voxel j", ylabel="super-voxel i",
    xticks=(1:Kc, ctick_labels), yticks=(1:Kc, ctick_labels),
    xticklabelsize=8, yticklabelsize=8,
    yreversed=true, xaxisposition=:top, titlesize=11)

heatmap!(ax_D, 1:Kc, 1:Kc, Bc'; colormap=CMAP_BLK, colorrange=CRANGE)
for k in 0:Kc
    hlines!(ax_D, [k+0.5]; color=(:gray60,0.5), linewidth=0.5)
    vlines!(ax_D, [k+0.5]; color=(:gray60,0.5), linewidth=0.5)
end

# Colour-coded diagonal cells matching spatial panel
for k in 1:Kc
    clr = pcolor(k)
    poly!(ax_D, Point2f[(k-0.5,0.5),(k+0.5,0.5),(k+0.5,0.0),(k-0.5,0.0)];
          color=(clr,0.55), strokewidth=0)
    poly!(ax_D, Point2f[(0.0,k-0.5),(0.5,k-0.5),(0.5,k+0.5),(0.0,k+0.5)];
          color=(clr,0.55), strokewidth=0)
end

# ─── Legend ───────────────────────────────────────────────────────────────────
n_coup_fine   = count(Bf .== 2) ÷ 2
n_coup_coarse = count(Bc .== 2) ÷ 2
Legend(fig[3,1:2],
    [PolyElement(color=:steelblue),
     PolyElement(color=:darkorange),
     PolyElement(color=:white, strokecolor=:gray60, strokewidth=1)],
    ["Aₖ — per-voxel reactions (diagonal)",
     "Dₖₖ' — diffusion coupling (off-diagonal)  fine: $n_coup_fine → coarse: $n_coup_coarse edges",
     "Zero (no coupling)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(20,12))

# ─── Title ────────────────────────────────────────────────────────────────────
Label(fig[0,1:2],
    "Birth-death RDME (D=$D, λ=$k_b, μ=$k_d, n̄=$( k_b/k_d )) — " *
    "flux-adaptive coarsening: $K voxels → $Kc super-voxels",
    fontsize=12, tellwidth=false)

rowgap!(fig.layout, 1, 10)
rowgap!(fig.layout, 2, 6)
colgap!(fig.layout, 1, 14)

# ─── Save ─────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "paper", "figures", "fig_bdmultigrid_matrix.pdf")
save(out, fig)
println("Saved → $out")
