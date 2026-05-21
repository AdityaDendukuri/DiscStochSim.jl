"""
anim_bdmultigrid.jl — Animated 2D diffusion RDME: adaptive coarsening collapse.

Point source at (1,1) diffuses across a 6×6 grid (pure diffusion, no reactions).
Flux-adaptive RECURSIVE coarsening merges voxels once both neighbours have
enough molecules. The block matrix fills with blue and collapses to 1×1 as
the wavefront sweeps the domain and the distribution equilibrates.

Left panel : 2D mean-field heatmap + white boundary lines between super-voxels.
Right panel: K×K block-coarsening matrix — blue = same super-voxel, orange = coupled.

Run: julia --project examples/anim_bdmultigrid.jl
"""

using CairoMakie

# ── Parameters ────────────────────────────────────────────────────────────────
const K_x, K_y = 6, 6;  const K = K_x * K_y
const D   = 0.8          # diffusion coefficient (dx=1) — slow enough to show spreading
const k_b = 0.0;  const k_d = 0.0   # pure diffusion — no reactions
const N_ic  = 200        # IC molecules at (1,1)
const T_tot = 30.0       # long enough for wavefront to equilibrate far corner
const N_FRM = 120
const dt_ref  = 0.5      # reference τ for flux weights
# Pair only when BOTH fine voxels have ≥ MIN_MU molecules (75% of equilibrium μ)
const mu_eq   = float(N_ic) / K  # ≈ 5.56 molecules per voxel at equilibrium
const MIN_MU  = 0.70 * mu_eq     # ≈ 3.9 — high bar so collapse tracks wavefront
const FLUX_MIN = D * MIN_MU * dt_ref   # absolute flux threshold

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

function advance_mf(mu0, Δt; dt_inner=0.005)
    mu = copy(mu0); t = 0.0
    while t < Δt; s = min(dt_inner, Δt-t); mf_step!(mu, s); t += s; end
    mu
end

# ── Greedy flux matching (threshold) ─────────────────────────────────────────
function flux_match(mu_vec, edges; min_w=FLUX_MIN)
    w    = [D * min(mu_vec[i], mu_vec[j]) * dt_ref for (i,j) in edges]
    used = fill(false, length(mu_vec)); pairs = Tuple{Int,Int}[]
    for e in sortperm(w; rev=true)
        w[e] < min_w && break
        i,j = edges[e]
        used[i]||used[j]||(push!(pairs,(i,j)); used[i]=used[j]=true)
    end
    pairs, findall(.!used)
end

function coarsen_one(Kv, edges, pairs, singles)
    groups = vcat([[a,b] for (a,b) in pairs], [[s] for s in singles])
    nm     = zeros(Int, Kv)
    for (g,ks) in enumerate(groups), k in ks; nm[k] = g; end
    ced = Set{Tuple{Int,Int}}()
    for (i,j) in edges
        gi,gj = nm[i],nm[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    groups, nm, sort(collect(ced)), length(groups)
end

# ── Recursive coarsening using MINIMUM fine mean (prevents cascade) ───────────
#
# At every recursion level the pairing criterion is based on the minimum
# original-voxel mean inside each super-voxel, NOT the average.  This stops
# the cascade: a super-voxel that contains even one "cold" fine voxel will not
# merge further until that voxel warms up.  The collapse therefore tracks the
# wavefront accurately.
function recursive_coarsen(mu_vec)
    cur_min   = copy(mu_vec)   # minimum fine mean within each current super-voxel
    cur_edges = copy(EDGES)
    cur_K     = K
    c2o       = [[k] for k in 1:K]   # current super-voxel → original voxels
    assign    = collect(1:K)

    while cur_K > 1
        # flux weight = D * min(min_mean_i, min_mean_j) * dt_ref
        w    = [D * min(cur_min[i], cur_min[j]) * dt_ref for (i,j) in cur_edges]
        used = fill(false, cur_K); pairs = Tuple{Int,Int}[]
        for e in sortperm(w; rev=true)
            w[e] < FLUX_MIN && break
            i,j = cur_edges[e]
            used[i]||used[j]||(push!(pairs,(i,j)); used[i]=used[j]=true)
        end
        isempty(pairs) && break

        singles = findall(.!used)
        groups, nm, new_edges, new_K = coarsen_one(cur_K, cur_edges, pairs, singles)

        new_c2o   = [vcat([c2o[g] for g in gs]...) for gs in groups]
        # new minimum = min fine mean over all original voxels in the super-voxel
        new_min   = [minimum(mu_vec[k] for k in new_c2o[j]) for j in 1:new_K]

        for (j,gs) in enumerate(groups), g in gs
            for k in c2o[g]; assign[k] = j; end
        end

        c2o = new_c2o; cur_K = new_K; cur_edges = new_edges; cur_min = new_min
    end

    uid  = unique(assign); mp = Dict(v=>i for (i,v) in enumerate(uid))
    asgn = [mp[assign[k]] for k in 1:K]
    Kc   = length(uid)

    ced = Set{Tuple{Int,Int}}()
    for (i,j) in EDGES
        gi,gj = asgn[i],asgn[j]; gi≠gj && push!(ced, gi<gj ? (gi,gj) : (gj,gi))
    end
    asgn, Kc, sort(collect(ced))
end

# ── Build K×K "coarsened block" matrix ────────────────────────────────────────
# 0 = no coupling, 1 = same super-voxel (merged/diagonal), 2 = adjacent coupling
function build_Bc(asgn, coarse_edges)
    ced_set = Set(coarse_edges)
    B = zeros(Int, K, K)
    for k1 in 1:K, k2 in 1:K
        g1,g2 = asgn[k1], asgn[k2]
        if g1 == g2; B[k1,k2] = 1
        elseif (min(g1,g2),max(g1,g2)) ∈ ced_set; B[k1,k2] = 2
        end
    end
    B
end

# ── Boundary lines between super-voxels (for spatial panel) ───────────────────
function boundary_segs(asgn)
    segs = Point2f[]
    for (vk,vl) in EDGES
        asgn[vk] == asgn[vl] && continue
        ik,jk = vij(vk); il,jl = vij(vl)
        if ik == il   # horizontal edge → vertical boundary line
            x = (jk+jl)/2.0
            push!(segs, Point2f(x,ik-0.5), Point2f(x,ik+0.5))
        else          # vertical edge → horizontal boundary line
            y = (ik+il)/2.0
            push!(segs, Point2f(jk-0.5,y), Point2f(jk+0.5,y))
        end
    end
    segs
end

# ── Pre-compute frames ────────────────────────────────────────────────────────
println("Pre-computing $N_FRM mean-field frames …")
t_frames  = collect(range(0.0, T_tot; length=N_FRM))
mu_frames = Vector{Matrix{Float64}}(undef, N_FRM)
mu_frames[1] = zeros(K_x, K_y); mu_frames[1][1,1] = float(N_ic)
for f in 2:N_FRM
    mu_frames[f] = advance_mf(mu_frames[f-1], t_frames[f]-t_frames[f-1])
end

println("Computing recursive coarsenings …")
asgn_frames = Vector{Vector{Int}}(undef, N_FRM)
kc_frames   = zeros(Int, N_FRM)
Bc_frames   = Vector{Matrix{Int}}(undef, N_FRM)
bseg_frames = Vector{Vector{Point2f}}(undef, N_FRM)

for f in 1:N_FRM
    muvec = [mu_frames[f][vij(k)...] for k in 1:K]
    asgn, Kc, ced = recursive_coarsen(muvec)
    asgn_frames[f] = asgn; kc_frames[f] = Kc
    Bc_frames[f]   = build_Bc(asgn, ced)
    bseg_frames[f] = boundary_segs(asgn)
    f % 20 == 0 && println("  frame $f  t=$(round(t_frames[f],digits=1))  Kc=$Kc")
end
println("Kc: $(maximum(kc_frames)) → $(minimum(kc_frames))")

# ── Figure layout ─────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(980, 520))

CMAP_BLK  = cgrad([:white, :steelblue, :darkorange])
CMAP_HEAT = :YlOrRd_9
HEAT_MAX  = mu_eq * 2.5   # 2.5× equilibrium; IC saturates via highclip

# Observables
mu_obs   = Observable(mu_frames[1]')          # K_y × K_x
bseg_obs = Observable(bseg_frames[1])
Bc_obs   = Observable(Bc_frames[1]')          # K × K

# ── Panel A: spatial ──────────────────────────────────────────────────────────
ax_sp = Axis(fig[1,1];
    title="Fine: $(K_x)×$(K_y) voxels  (D=$D, μ_eq=$(round(mu_eq,digits=1)))",
    aspect=DataAspect(), xlabel="column j", ylabel="row i",
    xticks=(1:K_y, string.(1:K_y)), yticks=(1:K_x, string.(1:K_x)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=9, yticklabelsize=9)

hm_sp = heatmap!(ax_sp, 1:K_y, 1:K_x, mu_obs;
    colormap=CMAP_HEAT, colorrange=(0, HEAT_MAX),
    highclip=:darkred, lowclip=:lightyellow)
linesegments!(ax_sp, bseg_obs; color=:white, linewidth=2.8)
for i in 0:K_x; hlines!(ax_sp,[i+0.5]; color=(:black,0.12), linewidth=0.4); end
for j in 0:K_y; vlines!(ax_sp,[j+0.5]; color=(:black,0.12), linewidth=0.4); end
Colorbar(fig[1,1][1,2], hm_sp; label="⟨n⟩", width=9, labelsize=9, ticklabelsize=7)

# ── Panel B: block matrix ─────────────────────────────────────────────────────
ax_mx = Axis(fig[1,2];
    aspect=DataAspect(), xlabel="voxel j", ylabel="voxel i",
    xticks=(1:2:K, string.(1:2:K)), yticks=(1:2:K, string.(1:2:K)),
    xticklabelsize=6, yticklabelsize=6, yreversed=true, xaxisposition=:top)

heatmap!(ax_mx, 1:K, 1:K, Bc_obs; colormap=CMAP_BLK, colorrange=(0,2))
for k in 0:K
    hlines!(ax_mx,[k+0.5]; color=(:gray55,0.28), linewidth=0.22)
    vlines!(ax_mx,[k+0.5]; color=(:gray55,0.28), linewidth=0.22)
end

# ── Labels ────────────────────────────────────────────────────────────────────
t_label  = Label(fig[0,1], "t = 0.00"; fontsize=13, halign=:left)
kc_label = Label(fig[0,2], "K_c = $K super-voxels"; fontsize=13, halign=:center)

Legend(fig[2,1:2],
    [PolyElement(color=:steelblue),
     PolyElement(color=:darkorange),
     PolyElement(color=:white, strokecolor=:gray60, strokewidth=1)],
    ["Same super-voxel (merged — no coupling needed)",
     "Adjacent super-voxels (active diffusion coupling)",
     "No coupling"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(18,12))

rowgap!(fig.layout, 1, 6); rowgap!(fig.layout, 2, 4)
colgap!(fig.layout, 1, 12)

# ── Record animation ──────────────────────────────────────────────────────────
out_gif = joinpath(@__DIR__, "..", "paper", "figures", "anim_bdmultigrid.gif")
out_mp4 = joinpath(@__DIR__, "..", "paper", "figures", "anim_bdmultigrid.mp4")

println("Rendering animation …")
@time record(fig, out_gif, 1:N_FRM; framerate=15) do f
    mu_obs[]   = mu_frames[f]'
    bseg_obs[] = bseg_frames[f]
    Bc_obs[]   = Bc_frames[f]'
    t_label.text[]  = "t = $(round(t_frames[f], digits=2))"
    kc_label.text[] = "K_c = $(kc_frames[f]) super-voxels"
end

println("Saved GIF → $out_gif")
println("Kc trajectory: ", kc_frames[1:10:end])
