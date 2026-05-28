"""
anim_schlogl_transition.jl — Schlögl phase transition on a 20×20 RDME grid.

IC: entire grid at low stable state (n_lo≈31).
    Centre 9×9 patch set to high stable state (n_hi≈162) — supercritical droplet.
The high-state front propagates outward until the whole grid transitions.

Coarsening threshold MIN_MU > n_lo so ONLY high-state voxels merge.
Block matrix: starts as one blue block (centre droplet) surrounded by diagonal
singletons; blue block grows with the front and finally collapses to 1×1.

Run: julia --project examples/anim_schlogl_transition.jl
"""

using CairoMakie

# ── Schlögl parameters (default SchloglModel1D(1.0)) ─────────────────────────
const c1 = 0.028; const c2 = 0.00032; const c3 = 19.5; const c4 = 1.0
const D  = 1.0
const n_lo  = 31.0;  const n_uns = 72.0;  const n_hi = 162.0

# ── Grid ──────────────────────────────────────────────────────────────────────
const K_x, K_y = 20, 20;  const K = K_x * K_y
const CENTER    = (K_x÷2 + 1, K_y÷2 + 1)   # (11,11)
const PATCH_R   = 6   # half-width of the initial high-state square patch (13×13)

# ── Animation ─────────────────────────────────────────────────────────────────
const T_tot  = 25.0
const N_FRM  = 100
const dt_ref = 0.4

# Coarsening: only voxels well above n_lo merge (singletons stay as singletons)
const MIN_MU   = 50.0            # > n_lo=31, << n_hi=162
const FLUX_MIN = D * MIN_MU * dt_ref
const IMGSIZE  = 300

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

# ── Schlögl mean-field ODE ────────────────────────────────────────────────────
@inline schlogl_f(n) = c3 - c4*n + c1*n*(n-1)/2 - c2*n*(n-1)*(n-2)/6

function mf_step!(mu, dt)
    Δ = similar(mu)
    for i in 1:K_x, j in 1:K_y
        n = mu[i,j]
        Δ[i,j] = schlogl_f(n)
        i > 1   && (Δ[i,j] += D*(mu[i-1,j]-mu[i,j]))
        i < K_x && (Δ[i,j] += D*(mu[i+1,j]-mu[i,j]))
        j > 1   && (Δ[i,j] += D*(mu[i,j-1]-mu[i,j]))
        j < K_y && (Δ[i,j] += D*(mu[i,j+1]-mu[i,j]))
    end
    mu .+= dt .* Δ
end

function advance_mf(mu0, Δt; dt_inner=0.005)
    mu = copy(mu0); t = 0.0
    while t < Δt; s = min(dt_inner, Δt-t); mf_step!(mu, s); t += s; end
    mu
end

# ── Greedy flux matching with min-means recursive coarsening ──────────────────
function flux_match(mu_min, edges)
    w    = [D * min(mu_min[i], mu_min[j]) * dt_ref for (i,j) in edges]
    used = fill(false, length(mu_min)); pairs = Tuple{Int,Int}[]
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

function recursive_coarsen(mu_vec)
    cur_min = copy(mu_vec); cur_edges = copy(EDGES); cur_K = K
    c2o = [[k] for k in 1:K]; assign = collect(1:K)
    while cur_K > 1
        pairs, singles = flux_match(cur_min, cur_edges)
        isempty(pairs) && break
        groups, nm, new_edges, new_K = coarsen_one(cur_K, cur_edges, pairs, singles)
        new_c2o = [vcat([c2o[g] for g in gs]...) for gs in groups]
        new_min = [minimum(mu_vec[k] for k in new_c2o[j]) for j in 1:new_K]
        for (j,gs) in enumerate(groups), g in gs
            for k in c2o[g]; assign[k]=j; end
        end
        c2o=new_c2o; cur_K=new_K; cur_edges=new_edges; cur_min=new_min
    end
    uid=unique(assign); mp=Dict(v=>i for (i,v) in enumerate(uid))
    asgn=[mp[assign[k]] for k in 1:K]; Kc=length(uid)
    ced=Set{Tuple{Int,Int}}()
    for (i,j) in EDGES
        gi,gj=asgn[i],asgn[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    asgn, Kc, sort(collect(ced))
end

# ── Reorder super-voxels by centroid distance from centre ─────────────────────
function reorder_by_centroid(asgn, Kc)
    si=zeros(Kc); sj=zeros(Kc); cnt=zeros(Int,Kc)
    for k in 1:K
        g=asgn[k]; i,j=vij(k); si[g]+=i; sj[g]+=j; cnt[g]+=1
    end
    ci,cj = float(CENTER[1]),float(CENTER[2])
    dists=[sqrt((si[g]/cnt[g]-ci)^2+(sj[g]/cnt[g]-cj)^2) for g in 1:Kc]
    remap=invperm(sortperm(dists))
    [remap[asgn[k]] for k in 1:K]
end

# ── Block adjacency + rendering ───────────────────────────────────────────────
function block_adj(Kc, cedges)
    B=zeros(Int,Kc,Kc); for k in 1:Kc; B[k,k]=1; end
    for (i,j) in cedges; B[i,j]=2; B[j,i]=2; end; B
end

const COL_DIAG  = RGBAf(0.27,0.51,0.71,1.0)
const COL_COUPL = RGBAf(0.89,0.45,0.13,1.0)
const COL_ZERO  = RGBAf(0.96,0.96,0.96,1.0)
const COL_GRID  = RGBAf(0.55,0.55,0.55,1.0)

function render_matrix(Bc, Kc)
    img=fill(COL_ZERO,IMGSIZE,IMGSIZE); cp=IMGSIZE/Kc
    for i in 1:Kc, j in 1:Kc
        clr=Bc[i,j]==1 ? COL_DIAG : Bc[i,j]==2 ? COL_COUPL : COL_ZERO
        r_lo=floor(Int,(i-1)*cp)+1; r_hi=min(IMGSIZE,floor(Int,i*cp))
        c_lo=floor(Int,(j-1)*cp)+1; c_hi=min(IMGSIZE,floor(Int,j*cp))
        for r in r_lo:r_hi, c in c_lo:c_hi; img[r,c]=clr; end
    end
    if cp>=3
        for k in 0:Kc
            px=min(IMGSIZE,floor(Int,k*cp)+1)
            for q in 1:IMGSIZE; img[px,q]=COL_GRID; img[q,px]=COL_GRID; end
        end
    end
    img
end

function boundary_segs(asgn)
    segs=Point2f[]
    for (vk,vl) in EDGES
        asgn[vk]==asgn[vl] && continue
        ik,jk=vij(vk); il,jl=vij(vl)
        if ik==il; x=(jk+jl)/2.0; push!(segs,Point2f(x,ik-0.5),Point2f(x,ik+0.5))
        else;      y=(ik+il)/2.0; push!(segs,Point2f(jk-0.5,y),Point2f(jk+0.5,y)); end
    end
    segs
end

# ── Pre-compute frames ────────────────────────────────────────────────────────
println("Setting up IC and computing $N_FRM Schlögl mean-field frames …")
t_frames = collect(range(0.0, T_tot; length=N_FRM))

mu0 = fill(n_lo, K_x, K_y)   # start at low stable state
ci,cj = CENTER
for i in max(1,ci-PATCH_R):min(K_x,ci+PATCH_R),
    j in max(1,cj-PATCH_R):min(K_y,cj+PATCH_R)
    mu0[i,j] = n_hi           # supercritical droplet at centre
end

mu_frames = Vector{Matrix{Float64}}(undef, N_FRM)
mu_frames[1] = mu0
for f in 2:N_FRM
    mu_frames[f] = advance_mf(mu_frames[f-1], t_frames[f]-t_frames[f-1])
end

println("Computing coarsenings …")
asgn_frames = Vector{Vector{Int}}(undef, N_FRM)
kc_frames   = zeros(Int, N_FRM)
Bc_frames   = Vector{Matrix{Int}}(undef, N_FRM)
bseg_frames = Vector{Vector{Point2f}}(undef, N_FRM)
img_frames  = Vector{Matrix{RGBAf}}(undef, N_FRM)

for f in 1:N_FRM
    muvec = [mu_frames[f][vij(k)...] for k in 1:K]
    asgn_raw, Kc, _ = recursive_coarsen(muvec)
    asgn = reorder_by_centroid(asgn_raw, Kc)
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in EDGES
        gi,gj=asgn[i],asgn[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    asgn_frames[f]=asgn; kc_frames[f]=Kc
    Bc_frames[f]  = block_adj(Kc, sort(collect(ced)))
    bseg_frames[f] = boundary_segs(asgn)
    img_frames[f]  = render_matrix(Bc_frames[f], Kc)
    f%20==0 && println("  frame $f  t=$(round(t_frames[f],digits=1))  Kc=$Kc")
end
println("Kc: $(kc_frames[1]) → $(minimum(kc_frames))")

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1050, 530))

CMAP_HEAT = :inferno
HEAT_MAX  = 175.0   # just above n_hi=162

mu_obs   = Observable(mu_frames[1]')
bseg_obs = Observable(bseg_frames[1])
img_obs  = Observable(img_frames[1])

# Panel A: spatial
ax_sp = Axis(fig[1,1];
    title="Schlögl 20×20  (D=$D,  n_lo=$(Int(n_lo)), n_hi=$(Int(n_hi)))",
    aspect=DataAspect(), xlabel="column j", ylabel="row i",
    xticks=(1:2:K_y, string.(1:2:K_y)), yticks=(1:2:K_x, string.(1:2:K_x)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=8, yticklabelsize=8)

hm = heatmap!(ax_sp, 1:K_y, 1:K_x, mu_obs;
    colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX), lowclip=:black)
linesegments!(ax_sp, bseg_obs; color=(:white,0.85), linewidth=1.8)
for i in 0:K_x; hlines!(ax_sp,[i+0.5]; color=(:white,0.05), linewidth=0.3); end
for j in 0:K_y; vlines!(ax_sp,[j+0.5]; color=(:white,0.05), linewidth=0.3); end

Colorbar(fig[1,3]; colormap=CMAP_HEAT, colorrange=(0,HEAT_MAX), lowclip=:black,
    ticks=([0,n_lo,n_uns,n_hi], ["0","n_lo=$(Int(n_lo))","n_uns=$(Int(n_uns))","n_hi=$(Int(n_hi))"]),
    label="⟨n⟩", width=12, labelsize=9, ticklabelsize=7)

# Panel B: matrix
ax_mx = Axis(fig[1,2];
    title="RDME block matrix",
    aspect=DataAspect(),
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    yreversed=true, xaxisposition=:top, titlesize=11)

image!(ax_mx, 0.5..K_y+0.5, 0.5..K_x+0.5, img_obs)

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))

t_label  = Label(fig[0,1], "t = 0.00"; fontsize=13, halign=:left)
kc_label = Label(fig[0,2], "K_c = $(kc_frames[1]) super-voxels"; fontsize=12, halign=:center)

Legend(fig[2,1:3],
    [PolyElement(color=:steelblue), PolyElement(color=:darkorange),
     PolyElement(color=RGBf(0.96,0.96,0.96), strokecolor=:gray60, strokewidth=1)],
    ["Same super-voxel (high-state merged)",
     "Adjacent super-voxels (wavefront coupling)",
     "Singleton (low-state or interface — not yet merged)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(18,12))

rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,4)
colgap!(fig.layout,1,12); colgap!(fig.layout,2,4)

# ── Record ────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__,"..","paper","figures","anim_schlogl_transition.gif")
println("Rendering …")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]   = mu_frames[f]'
    bseg_obs[] = bseg_frames[f]
    img_obs[]  = img_frames[f]
    t_label.text[]  = "t = $(round(t_frames[f],digits=1))"
    kc_label.text[] = "K_c = $(kc_frames[f]) super-voxels"
end
println("Saved → $out")
