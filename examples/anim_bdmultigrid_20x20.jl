"""
anim_bdmultigrid_20x20.jl — Animated 20×20 diffusion RDME: wavefront + matrix collapse.

Point source at center (11,11) diffuses across a 20×20 grid (pure diffusion).
Left panel : mean-field heatmap + white super-voxel boundary lines.
Right panel: K_c×K_c coarse block matrix rendered as a fixed 300×300 image.
             Each cell GROWS as K_c shrinks → literally "collapsing into one".

Run: julia --project examples/anim_bdmultigrid_20x20.jl
"""

using CairoMakie

# ── Parameters ────────────────────────────────────────────────────────────────
const K_x, K_y = 20, 20;  const K = K_x * K_y
const D   = 1.2          # diffusion coefficient
const k_b = 0.0;  const k_d = 0.0   # pure diffusion — wavefront from IC only
const N_ic   = 300        # IC molecules at center
const T_tot  = 36.0
const N_FRM  = 120
const dt_ref = 0.5
const mu_eq  = float(N_ic) / K    # ≈ 0.75 per voxel at equilibrium
const MIN_MU = 0.62 * mu_eq       # pair only when ≥ 62% of equilibrium
const FLUX_MIN = D * MIN_MU * dt_ref

const CENTER = (K_x÷2 + 1, K_y÷2 + 1)   # (11,11)
const IMGSIZE = 300                        # matrix image pixels

# ── Indexing + edges ──────────────────────────────────────────────────────────
vid(i,j) = (i-1)*K_y + j
vij(k)   = ((k-1)÷K_y + 1, mod1(k, K_y))

const EDGES = let es = Tuple{Int,Int}[]
    for i in 1:K_x, j in 1:K_y
        j < K_y && push!(es, (vid(i,j), vid(i,j+1)))
        i < K_x && push!(es, (vid(i,j), vid(i+1,j)))
    end
    es
end

# ── Mean-field ODE step ───────────────────────────────────────────────────────
function mf_step!(mu, dt)
    Δ = similar(mu)
    for i in 1:K_x, j in 1:K_y
        Δ[i,j] = k_b - k_d * mu[i,j]
        i > 1   && (Δ[i,j] += D*(mu[i-1,j]-mu[i,j]))
        i < K_x && (Δ[i,j] += D*(mu[i+1,j]-mu[i,j]))
        j > 1   && (Δ[i,j] += D*(mu[i,j-1]-mu[i,j]))
        j < K_y && (Δ[i,j] += D*(mu[i,j+1]-mu[i,j]))
    end
    mu .+= dt .* Δ
end

function advance_mf(mu0, Δt; dt_inner=0.01)
    mu = copy(mu0); t = 0.0
    while t < Δt; s = min(dt_inner, Δt-t); mf_step!(mu, s); t += s; end
    mu
end

# ── Greedy flux matching with threshold (min-means criterion) ─────────────────
function flux_match(mu_min_vec, edges)
    w    = [D * min(mu_min_vec[i], mu_min_vec[j]) * dt_ref for (i,j) in edges]
    used = fill(false, length(mu_min_vec)); pairs = Tuple{Int,Int}[]
    for e in sortperm(w; rev=true)
        w[e] < FLUX_MIN && break
        i,j = edges[e]
        used[i]||used[j]||(push!(pairs,(i,j)); used[i]=used[j]=true)
    end
    pairs, findall(.!used)
end

function coarsen_one(Kv, edges, pairs, singles)
    groups = vcat([[a,b] for (a,b) in pairs], [[s] for s in singles])
    nm = zeros(Int, Kv)
    for (g,ks) in enumerate(groups), k in ks; nm[k]=g; end
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in edges
        gi,gj=nm[i],nm[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    groups, nm, sort(collect(ced)), length(groups)
end

# ── Recursive coarsening using MINIMUM fine mean (prevents cascade) ───────────
function recursive_coarsen(mu_vec)
    cur_min   = copy(mu_vec)
    cur_edges = copy(EDGES)
    cur_K     = K
    c2o       = [[k] for k in 1:K]
    assign    = collect(1:K)

    while cur_K > 1
        pairs, singles = flux_match(cur_min, cur_edges)
        isempty(pairs) && break
        groups, nm, new_edges, new_K = coarsen_one(cur_K, cur_edges, pairs, singles)
        new_c2o  = [vcat([c2o[g] for g in gs]...) for gs in groups]
        new_min  = [minimum(mu_vec[k] for k in new_c2o[j]) for j in 1:new_K]
        for (j,gs) in enumerate(groups), g in gs
            for k in c2o[g]; assign[k]=j; end
        end
        c2o=new_c2o; cur_K=new_K; cur_edges=new_edges; cur_min=new_min
    end

    uid = unique(assign); mp = Dict(v=>i for (i,v) in enumerate(uid))
    asgn = [mp[assign[k]] for k in 1:K]; Kc = length(uid)
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in EDGES
        gi,gj=asgn[i],asgn[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    asgn, Kc, sort(collect(ced))
end

# ── Reorder coarse super-voxels by centroid distance from center ──────────────
# Puts inner (merged) super-voxels first so connections stay near the diagonal.
function reorder_by_centroid(asgn, Kc)
    si = zeros(Kc); sj = zeros(Kc); cnt = zeros(Int, Kc)
    for k in 1:K
        g = asgn[k]; i,j = vij(k)
        si[g] += i; sj[g] += j; cnt[g] += 1
    end
    ci, cj = float(CENTER[1]), float(CENTER[2])
    dists = [sqrt((si[g]/cnt[g]-ci)^2 + (sj[g]/cnt[g]-cj)^2) for g in 1:Kc]
    order = sortperm(dists)               # inner super-voxels get low indices
    remap = invperm(order)
    [remap[asgn[k]] for k in 1:K]
end

# ── Build K_c×K_c block adjacency ─────────────────────────────────────────────
function block_adj(Kc, cedges)
    B = zeros(Int8, Kc, Kc)
    for k in 1:Kc; B[k,k] = 1; end
    for (i,j) in cedges; B[i,j] = 2; B[j,i] = 2; end
    B
end

# ── Render K_c×K_c matrix to fixed IMGSIZE×IMGSIZE RGBA image ────────────────
# Each cell is IMGSIZE/K_c pixels wide — cells GROW as K_c shrinks.
const COL_DIAG   = RGBAf(0.27, 0.51, 0.71, 1.0)   # steelblue  — same super-voxel
const COL_COUPL  = RGBAf(0.89, 0.45, 0.13, 1.0)   # darkorange — active coupling
const COL_ZERO   = RGBAf(0.96, 0.96, 0.96, 1.0)   # near-white — no coupling
const COL_GRID   = RGBAf(0.55, 0.55, 0.55, 1.0)   # grid line color

function render_matrix(Bc, Kc)
    img = fill(COL_ZERO, IMGSIZE, IMGSIZE)
    cp  = IMGSIZE / Kc
    for i in 1:Kc, j in 1:Kc
        clr = Bc[i,j] == 1 ? COL_DIAG : Bc[i,j] == 2 ? COL_COUPL : COL_ZERO
        r_lo = floor(Int,(i-1)*cp)+1; r_hi = min(IMGSIZE, floor(Int,i*cp))
        c_lo = floor(Int,(j-1)*cp)+1; c_hi = min(IMGSIZE, floor(Int,j*cp))
        for r in r_lo:r_hi, c in c_lo:c_hi; img[r,c] = clr; end
    end
    if cp >= 3
        for k in 0:Kc
            px = min(IMGSIZE, floor(Int, k*cp) + 1)
            for q in 1:IMGSIZE; img[px,q] = COL_GRID; img[q,px] = COL_GRID; end
        end
    end
    img
end

# ── Boundary segments for spatial panel ───────────────────────────────────────
function boundary_segs(asgn)
    segs = Point2f[]
    for (vk,vl) in EDGES
        asgn[vk]==asgn[vl] && continue
        ik,jk=vij(vk); il,jl=vij(vl)
        if ik==il
            x=(jk+jl)/2.0; push!(segs, Point2f(x,ik-0.5), Point2f(x,ik+0.5))
        else
            y=(ik+il)/2.0; push!(segs, Point2f(jk-0.5,y), Point2f(jk+0.5,y))
        end
    end
    segs
end

# ── Pre-compute frames ────────────────────────────────────────────────────────
println("Pre-computing $N_FRM mean-field frames on $(K_x)×$(K_y) grid …")
t_frames  = collect(range(0.0, T_tot; length=N_FRM))
mu_frames = Vector{Matrix{Float64}}(undef, N_FRM)
mu_frames[1] = zeros(K_x, K_y)
mu_frames[1][CENTER...] = float(N_ic)

for f in 2:N_FRM
    mu_frames[f] = advance_mf(mu_frames[f-1], t_frames[f]-t_frames[f-1])
end

println("Computing recursive coarsenings …")
asgn_frames = Vector{Vector{Int}}(undef, N_FRM)
kc_frames   = zeros(Int, N_FRM)
Bc_frames   = Vector{Matrix{Int}}(undef, N_FRM)
bseg_frames = Vector{Vector{Point2f}}(undef, N_FRM)
img_frames  = Vector{Matrix{RGBAf}}(undef, N_FRM)

for f in 1:N_FRM
    muvec = [mu_frames[f][vij(k)...] for k in 1:K]
    asgn_raw, Kc, _ = recursive_coarsen(muvec)
    asgn = reorder_by_centroid(asgn_raw, Kc)  # inner super-voxels get low indices
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in EDGES
        gi,gj=asgn[i],asgn[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    asgn_frames[f] = asgn; kc_frames[f] = Kc
    Bc_frames[f]   = block_adj(Kc, sort(collect(ced)))
    bseg_frames[f] = boundary_segs(asgn)
    img_frames[f]  = render_matrix(Bc_frames[f], Kc)
    f % 20 == 0 && println("  frame $f / $N_FRM  t=$(round(t_frames[f],digits=1))  Kc=$Kc")
end
println("Kc trajectory: $(kc_frames[1:12:end])")

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1050, 530))

CMAP_HEAT = :inferno
HEAT_MAX  = mu_eq * 2.2

# Observables
mu_obs   = Observable(mu_frames[1]')
bseg_obs = Observable(bseg_frames[1])
img_obs  = Observable(img_frames[1])

# ── Panel A: spatial ──────────────────────────────────────────────────────────
ax_sp = Axis(fig[1,1];
    title="$(K_x)×$(K_y) diffusion RDME  (D=$D,  μ_eq=$(round(mu_eq,digits=2)))",
    aspect=DataAspect(), xlabel="column j", ylabel="row i",
    xticks=(1:2:K_y, string.(1:2:K_y)),
    yticks=(1:2:K_x, string.(1:2:K_x)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=8, yticklabelsize=8)

hm = heatmap!(ax_sp, 1:K_y, 1:K_x, mu_obs;
    colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX),
    highclip=:white, lowclip=:black)
linesegments!(ax_sp, bseg_obs; color=:white, linewidth=1.8)
scatter!(ax_sp, [CENTER[2]], [CENTER[1]]; marker=:cross, markersize=10,
    color=:cyan, strokewidth=0)
for i in 0:K_x; hlines!(ax_sp,[i+0.5]; color=(:white,0.06), linewidth=0.3); end
for j in 0:K_y; vlines!(ax_sp,[j+0.5]; color=(:white,0.06), linewidth=0.3); end

# Colorbar in its own narrow column so it doesn't steal width from panel A
Colorbar(fig[1,3]; colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX),
    highclip=:white, lowclip=:black,
    label="⟨n⟩", width=12, labelsize=9, ticklabelsize=7)

# ── Panel B: matrix image ─────────────────────────────────────────────────────
ax_mx = Axis(fig[1,2];
    title="RDME block matrix",
    aspect=DataAspect(),
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    yreversed=true, xaxisposition=:top, titlesize=11)

image!(ax_mx, 0.5..K_y+0.5, 0.5..K_x+0.5, img_obs)

# Make both data columns equal width; colorbar column is fixed narrow
colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))

# ── Labels ────────────────────────────────────────────────────────────────────
t_label  = Label(fig[0,1], "t = 0.00"; fontsize=13, halign=:left)
kc_label = Label(fig[0,2],
    "K_c = $K super-voxels"; fontsize=12, halign=:center)

Legend(fig[2,1:3],
    [PolyElement(color=:steelblue),
     PolyElement(color=:darkorange),
     PolyElement(color=RGBf(0.96,0.96,0.96), strokecolor=:gray60, strokewidth=1)],
    ["Same super-voxel (merged — diagonal block)",
     "Adjacent super-voxels (diffusion coupling)",
     "No coupling  (white lines = super-voxel boundaries)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(18,12))

rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 4)
colgap!(fig.layout, 1, 12); colgap!(fig.layout, 2, 4)

# ── Record ────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "paper", "figures", "anim_bdmultigrid_20x20.gif")
println("Rendering animation …")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]   = mu_frames[f]'
    bseg_obs[] = bseg_frames[f]
    img_obs[]  = img_frames[f]
    t_label.text[]  = "t = $(round(t_frames[f], digits=1))"
    kc_label.text[] = "K_c = $(kc_frames[f]) super-voxels"
end
println("Saved → $out")
println("Kc: $(maximum(kc_frames)) → $(minimum(kc_frames))")
