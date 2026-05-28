using DiscStochSim, CairoMakie, Printf

# ── Model ─────────────────────────────────────────────────────────────────────
const k_b  = 0.5; const k_d = 0.1; const D = 1.0
const K    = 200; const N0 = 3
const dt   = 0.4; const T_END = 200.0; const N_FRM = 150
const n_max = 14
const IMGSIZE = 300   # matrix thumbnail resolution

model  = BranchingDeathRDME(k_b, k_d, D, 1.0)
v_wave = 2*sqrt(D*(k_b-k_d))
@printf("Model: A→2A (k_b=%.1f) A→∅ (k_d=%.1f) D=%.1f  v_wave≈%.2f\n",
        k_b, k_d, D, v_wave)

state = SpatialFSP(model, K, K;
                   n_max, ε_expand=0.02, ε_equil=0.15,
                   use_block_vcycle=true)

# IC: disc of radius 4 voxels — visible initial bright spot
ci, cj = K÷2+1, K÷2+1
for di in -4:4, dj in -4:4
    di^2+dj^2 <= 16 || continue
    1<=ci+di<=K && 1<=cj+dj<=K || continue
    set_ic!(state, CartesianIndex(ci+di, cj+dj), N0)
end

# ── Colour helpers ─────────────────────────────────────────────────────────────
const _PLASMA  = CairoMakie.to_colormap(:plasma)
_plasma(t)     = let c=_PLASMA[clamp(floor(Int,t*(length(_PLASMA)-1))+1,1,length(_PLASMA))];
                 RGBf(c.r,c.g,c.b) end
const _mu_max  = Ref(1.0)

const C_EMPTY  = RGBf(0.04, 0.02, 0.06)
const C_ACT    = RGBf(0.88, 0.48, 0.22)   # orange — active wave front
const C_EQ     = RGBf(0.27, 0.45, 0.72)   # blue  — equil interior
const COL_DIAG = RGBAf(0.27,0.51,0.71,1.0)   # per-voxel CME blocks
const COL_COUP = RGBAf(0.89,0.45,0.13,1.0)   # diffusion coupling blocks
const COL_ZERO = RGBAf(0.96,0.96,0.96,1.0)   # empty entries

# ── Frame builder ──────────────────────────────────────────────────────────────
"""
Returns (merged_img, matrix_img, n_active).

merged_img (Change 1): heatmap with orange tint for active ring, blue tint for
equil interior.  Single panel; grid structure visible through colour tinting.

matrix_img (Change 2): block-level RDME sparsity — K_active × K_active pixels:
  diagonal (per-voxel CME) = blue, adjacent off-diagonal (diffusion) = orange.
"""
function make_frames(s)
    mg = mean_grid(s)
    sg = status_grid(s)

    # Auto-scale colormap to current peak
    cur_max = maximum(mg; init=1e-3)
    _mu_max[] = max(_mu_max[], cur_max * 1.05)
    mu_scale = _mu_max[]

    # ── Heatmap: pure plasma colormap (grid overlay added separately as lines) ─
    merged = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        v = clamp(sqrt(mg[i,j] / max(mu_scale,1e-3)), 0.0, 1.0)
        merged[j,i] = _plasma(v)
    end

    # ── Change 2: multiscale RDME matrix ──────────────────────────────────────
    # Ordering: equil voxels first (big merged block), then active ring (fine).
    # Pixel size proportional to voxel count:
    #   eq_pix  = IMGSIZE × N_eq/(N_eq+N_act)  → one large equil block
    #   act_pix = IMGSIZE × N_act/(N_eq+N_act) → individual active blocks
    act  = [(i,j) for i in 1:K, j in 1:K if sg[i,j]==1]
    eq   = [(i,j) for i in 1:K, j in 1:K if sg[i,j]==2]
    n_a  = length(act); n_e = length(eq); total = n_a + n_e
    mat  = fill(COL_ZERO, IMGSIZE, IMGSIZE)

    total == 0 && return merged, mat, sg, n_a

    # Pixel regions
    eq_pix  = round(Int, IMGSIZE * n_e / total)   # top-left equil block
    act_pix = IMGSIZE - eq_pix                     # bottom-right active block

    # ── Equil block (top-left) ─────────────────────────────────────────────
    # Individual equil voxels have their own CME but are shown as a merged unit.
    # Internal structure: near-diagonal (each voxel coupled to 4 neighbours).
    if eq_pix > 0 && n_e > 0
        eidx = Dict(v=>k for (k,v) in enumerate(eq))
        eq_adj = Set{Tuple{Int,Int}}()
        for (k,(i,j)) in enumerate(eq)
            for (ni,nj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
                l = get(eidx,(ni,nj),0)
                l>0 && push!(eq_adj, minmax(k,l))
            end
        end
        for r in 1:eq_pix, c in 1:eq_pix
            ei = clamp(floor(Int,(r-0.5)/eq_pix*n_e)+1, 1, n_e)
            ej = clamp(floor(Int,(c-0.5)/eq_pix*n_e)+1, 1, n_e)
            if ei==ej
                mat[r,c] = COL_DIAG
            elseif (minmax(ei,ej)) ∈ eq_adj
                mat[r,c] = COL_COUP
            end
        end
        # Draw separator line between equil and active blocks
        for p in 1:IMGSIZE
            mat[min(eq_pix+1,IMGSIZE),p] = RGBAf(0.3,0.3,0.3,1)
            mat[p,min(eq_pix+1,IMGSIZE)] = RGBAf(0.3,0.3,0.3,1)
        end
    end

    # ── Active block (bottom-right) ────────────────────────────────────────
    if act_pix > 0 && n_a > 0
        aidx = Dict(v=>k for (k,v) in enumerate(act))
        act_adj = Set{Tuple{Int,Int}}()
        for (k,(i,j)) in enumerate(act)
            for (ni,nj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
                l = get(aidx,(ni,nj),0)
                l>0 && push!(act_adj, minmax(k,l))
            end
        end
        off = eq_pix   # pixel offset into the active sub-region
        for r in 1:act_pix, c in 1:act_pix
            bi = clamp(floor(Int,(r-0.5)/act_pix*n_a)+1, 1, n_a)
            bj = clamp(floor(Int,(c-0.5)/act_pix*n_a)+1, 1, n_a)
            pr = off+r; pc = off+c
            if bi==bj
                mat[pr,pc] = COL_DIAG
            elseif (minmax(bi,bj)) ∈ act_adj
                mat[pr,pc] = COL_COUP
            end
        end
    end

    # ── Interface coupling: equil boundary ↔ adjacent active ──────────────
    if eq_pix > 0 && act_pix > 0 && n_e > 0 && n_a > 0
        eidx_m = Dict(v=>k for (k,v) in enumerate(eq))
        aidx_m = Dict(v=>k for (k,v) in enumerate(act))
        for (qi,(ei,ej)) in enumerate(eq)
            for (ni,nj) in ((ei-1,ej),(ei+1,ej),(ei,ej-1),(ei,ej+1))
                ai = get(aidx_m,(ni,nj),0); ai==0 && continue
                # Equil row qi → active col ai in off-diagonal rectangles
                r_eq  = clamp(round(Int,(qi-0.5)*eq_pix/n_e)+1, 1, eq_pix)
                c_act = eq_pix + clamp(round(Int,(ai-0.5)*act_pix/n_a)+1, 1, act_pix)
                for dr in -1:1, dc in -1:1
                    rr = clamp(r_eq+dr, 1, eq_pix)
                    cc = clamp(c_act+dc, eq_pix+1, IMGSIZE)
                    mat[rr,cc] = COL_COUP   # equil→active coupling
                    mat[cc,rr] = COL_COUP   # active→equil (symmetric)
                end
            end
        end
    end

    merged, mat, sg, n_a
end

# ── Pre-compute frames + capture snapshot for paper figure ─────────────────────
"""
    compute_overlay_segs(sg) → (boundary_segs, fine_segs, coarse_segs)

Three sets of line segments for the quad-tree overlay:
  boundary_segs: edges at active↔empty or active↔equil — marks the wave front ring
  fine_segs:     edges between adjacent active voxels — shows fine 1×1 grid in ring
  coarse_segs:   8×8 block edges fully inside equil interior — shows coarse resolution
"""
function compute_overlay_segs(sg)
    bound_x = Float32[]; bound_y = Float32[]   # wave front boundary
    fine_x  = Float32[]; fine_y  = Float32[]   # fine grid within ring
    coarse_x= Float32[]; coarse_y= Float32[]   # coarse equil grid

    # Draw a horizontal edge between row i and i+1 at column j
    function hadd!(xs, ys, i, j)
        push!(xs, j-0.5f0, j+0.5f0, NaN32)
        push!(ys, Float32(i)+0.5f0, Float32(i)+0.5f0, NaN32)
    end
    # Draw a vertical edge between col j and j+1 at row i
    function vadd!(xs, ys, i, j)
        push!(xs, Float32(j)+0.5f0, Float32(j)+0.5f0, NaN32)
        push!(ys, i-0.5f0, i+0.5f0, NaN32)
    end

    for i in 1:K, j in 1:K
        s = sg[i,j]
        # Horizontal edge below voxel (i,j) — border with (i+1,j)
        if i < K
            s2 = sg[i+1,j]
            if s==1 && s2!=1       # active → boundary
                hadd!(bound_x, bound_y, i, j)
            elseif s==1 && s2==1   # active-active → fine grid
                hadd!(fine_x, fine_y, i, j)
            end
        end
        # Vertical edge right of (i,j) — border with (i,j+1)
        if j < K
            s2 = sg[i,j+1]
            if s==1 && s2!=1
                vadd!(bound_x, bound_y, i, j)
            elseif s==1 && s2==1
                vadd!(fine_x, fine_y, i, j)
            end
        end
        # Coarse equil grid: draw 8×8 block edges within equil interior
        if s==2
            if i%8==0 && i<K && sg[i+1,j]==2
                hadd!(coarse_x, coarse_y, i, j)
            end
            if j%8==0 && j<K && sg[i,j+1]==2
                vadd!(coarse_x, coarse_y, i, j)
            end
        end
    end
    # Also add outer boundary of active ring (active voxels at the grid edge)
    for j in 1:K
        sg[1,j]==1   && hadd!(bound_x, bound_y, 0, j)
        sg[K,j]==1   && hadd!(bound_x, bound_y, K, j)
    end
    for i in 1:K
        sg[i,1]==1   && vadd!(bound_x, bound_y, i, 0)
        sg[i,K]==1   && vadd!(bound_x, bound_y, i, K)
    end

    (Point2f.(bound_x, bound_y),
     Point2f.(fine_x,  fine_y),
     Point2f.(coarse_x,coarse_y))
end

t_frames   = collect(range(0.0, T_END; length=N_FRM))
mu_frames  = Vector{Matrix{RGBf}}(undef, N_FRM)
mat_frames = Vector{Matrix{RGBAf}}(undef, N_FRM)
sg_frames  = Vector{Matrix{Int}}(undef, N_FRM)   # status grids for overlay
nact_frames = zeros(Int, N_FRM)

# State snapshot for the paper figure (captured at ~40% through the simulation)
const PAPER_FRAME = round(Int, N_FRM * 0.40)
paper_state_dists   = Dict{CartesianIndex{2},Vector{Float64}}()
paper_state_equil   = Dict{CartesianIndex{2},Vector{Float64}}()
paper_mg = zeros(K, K); paper_sg = zeros(Int, K, K); paper_t = 0.0

mu_frames[1], mat_frames[1], sg_frames[1], nact_frames[1] = make_frames(state)
println("Pre-computing $N_FRM frames…"); flush(stdout)

let fi = 2
    for step in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step*dt
        if fi <= N_FRM && t >= t_frames[fi]-1e-9
            mu_frames[fi], mat_frames[fi], sg_frames[fi], nact_frames[fi] = make_frames(state)
            # Snapshot for paper figure
            if fi == PAPER_FRAME
                global paper_t = t
                for (idx,p) in state.dists;       paper_state_dists[idx] = copy(p); end
                for (idx,p) in state.equil_dists; paper_state_equil[idx] = copy(p); end
                paper_mg .= mean_grid(state)
                paper_sg .= status_grid(state)
            end
            fi%10==0 && @printf("  frame %3d  t=%5.1f  act=%4d  eq=%5d\n",
                                fi, t, n_active(state), n_equilibrated(state))
            fi += 1
        end
        fi > N_FRM && break
    end
end
@printf("Peak active: %d\n\n", maximum(nact_frames))

# ── Animation figure (2 panels) ────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8,xgridvisible=false,ygridvisible=false)))
fig = Figure(size=(1100, 560))

mu_obs  = Observable(mu_frames[1])
mat_obs = Observable(mat_frames[1])

ax_kw = (aspect=DataAspect(), yreversed=true, xaxisposition=:top, titlesize=10,
         xticksvisible=false, yticksvisible=false,
         xticklabelsvisible=false, yticklabelsvisible=false)

# Left: merged heatmap + grid overlay
ax_mu = Axis(fig[1,1];
    title="FSP ⟨n_k⟩  ($(K)×$(K), A→2A)  |  quad-tree: fine=active ring, coarse=equil, none=empty",
    ax_kw...)
image!(ax_mu, 0.5..K+0.5, 0.5..K+0.5, mu_obs)

# Quad-tree line overlay: three levels of resolution
init_segs = compute_overlay_segs(sg_frames[1])
bound_obs  = Observable(init_segs[1])   # wave front boundary (thick)
fine_obs   = Observable(init_segs[2])   # fine grid in active ring (thin)
coarse_obs = Observable(init_segs[3])   # coarse equil grid (subtle)

linesegments!(ax_mu, coarse_obs; color=(:white,0.20), linewidth=0.4)
linesegments!(ax_mu, fine_obs;   color=(:white,0.45), linewidth=0.6)
linesegments!(ax_mu, bound_obs;  color=(:white,0.90), linewidth=1.4)

# Right: RDME block-matrix sparsity
ax_mat = Axis(fig[1,2];
    title="RDME matrix structure  (top-left=equil merged block, bottom-right=active ring)",
    ax_kw...)
image!(ax_mat, 0.5..IMGSIZE+0.5, 0.5..IMGSIZE+0.5, mat_obs)

for col in 1:2; colsize!(fig.layout, col, Aspect(1,1.0)); end
t_label  = Label(fig[0,1], "t = 0.0"; fontsize=13, halign=:left)
na_label = Label(fig[0,2], "";         fontsize=10, halign=:left)
Legend(fig[2,1:2],
    [LineElement(color=(:white,0.9),linewidth=2),
     LineElement(color=(:white,0.5),linewidth=1),
     LineElement(color=(:white,0.2),linewidth=0.5),
     PolyElement(color=RGBf(COL_DIAG.r,COL_DIAG.g,COL_DIAG.b)),
     PolyElement(color=RGBf(COL_COUP.r,COL_COUP.g,COL_COUP.b))],
    ["Wave front boundary (thick)","Active ring fine grid","Equil coarse grid (8×8)",
     "CME block (diagonal)","Diffusion coupling (off-diagonal)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(14,10))
rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,4); colgap!(fig.layout,1,8)

outdir = joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)
out    = joinpath(outdir,"anim_branching_200x200.gif")
println("Rendering animation → $out"); flush(stdout)
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]  = mu_frames[f]
    mat_obs[] = mat_frames[f]
    # Update quad-tree line overlay
    segs = compute_overlay_segs(sg_frames[f])
    bound_obs[]  = segs[1]
    fine_obs[]   = segs[2]
    coarse_obs[] = segs[3]
    t_label.text[]  = "t = $(round(t_frames[f],digits=1))"
    na_label.text[] = "active: $(nact_frames[f])"
end
println("Saved animation → $out\n")

# ── Change 3: paper-style snapshot figure ─────────────────────────────────────
println("Generating paper figure at t=$(round(paper_t,digits=1))…")

# Choose three representative voxels from the snapshot
# 1. Equil interior (behind the wave front, distribution peaked at ~7)
# 2. Active front (wave front, transitioning distribution)
# 3. Newly active (just entered ring, low-n distribution)
function pick_voxels(sg, dists, equil_dists)
    act_voxels  = [(i,j) for i in 1:K, j in 1:K if sg[i,j]==1]
    equil_voxels = [(i,j) for i in 1:K, j in 1:K if sg[i,j]==2]
    isempty(act_voxels) || isempty(equil_voxels) && return nothing

    # Equil interior: pick one with highest mean (well-saturated)
    eq_means = [(voxel_mean(equil_dists[CartesianIndex(i,j)]), (i,j))
                for (i,j) in equil_voxels if haskey(equil_dists, CartesianIndex(i,j))]
    isempty(eq_means) && return nothing
    sort!(eq_means; rev=true)
    v_equil = eq_means[max(1,length(eq_means)÷2)][2]   # median equil

    # Active front: pick one with mid-range mean (in transition)
    act_means = [(voxel_mean(dists[CartesianIndex(i,j)]), (i,j))
                 for (i,j) in act_voxels if haskey(dists, CartesianIndex(i,j))]
    isempty(act_means) && return nothing
    sort!(act_means)
    # "Mid-ring": voxel that has been active a while — building up toward equil
    v_front = act_means[max(1, round(Int,length(act_means)*0.75))][2]

    # "Leading edge": freshest/sparsest — shows the stochastic sparse front
    v_new = act_means[max(1, round(Int,length(act_means)*0.15))][2]

    v_equil, v_front, v_new
end

vox_pick = pick_voxels(paper_sg, paper_state_dists, paper_state_equil)

if vox_pick !== nothing
    v_eq, v_fr, v_nw = vox_pick

    get_dist(v) = begin
        ci = CartesianIndex(v...)
        if haskey(paper_state_dists, ci);  return paper_state_dists[ci]
        elseif haskey(paper_state_equil, ci); return paper_state_equil[ci]
        else; p=zeros(n_max+1); p[1]=1.0; return p
        end
    end

    dist_eq = get_dist(v_eq)
    dist_fr = get_dist(v_fr)
    dist_nw = get_dist(v_nw)

    # ── Build paper figure ─────────────────────────────────────────────────────
    set_theme!(Theme(fontsize=10,
        Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
              topspinevisible=false, rightspinevisible=false)))

    pfig = Figure(size=(800, 700))

    # Main heatmap (reuse snapshot)
    ax_main = Axis(pfig[1,1];
        aspect=DataAspect(), yreversed=true,
        xaxisposition=:top, titlesize=9,
        title="A→2A RDME, t=$(round(paper_t,digits=0)),  $(K)×$(K) grid",
        xticksvisible=false, yticksvisible=false,
        xticklabelsvisible=false, yticklabelsvisible=false)

    # Recompute merged image at the snapshot time using stored mg/sg
    snap_img = zeros(RGBf, K, K)
    snap_max = maximum(paper_mg; init=1e-3)
    for j in 1:K, i in 1:K
        v2 = clamp(sqrt(paper_mg[i,j]/max(snap_max,1e-3)), 0.0, 1.0)
        bc = _plasma(v2)
        if paper_sg[i,j]==1
            α=0.40
            snap_img[j,i]=RGBf((1-α)*bc.r+α*C_ACT.r,(1-α)*bc.g+α*C_ACT.g,(1-α)*bc.b+α*C_ACT.b)
        elseif paper_sg[i,j]==2
            α=0.18
            snap_img[j,i]=RGBf((1-α)*bc.r+α*C_EQ.r,(1-α)*bc.g+α*C_EQ.g,(1-α)*bc.b+α*C_EQ.b)
        else
            snap_img[j,i]=bc
        end
    end
    image!(ax_main, 0.5..K+0.5, 0.5..K+0.5, snap_img)

    # Mark the three chosen voxels on the heatmap
    vox_colors = [:gold, :tomato, :cyan]
    vox_labels = ["A: equil", "B: mid-ring", "C: leading edge"]
    for (v, col, lbl) in zip([v_eq, v_fr, v_nw], vox_colors, vox_labels)
        i, j = v
        # Small square marker at voxel center
        poly!(ax_main, Rect(j-0.5, i-0.5, 1.0, 1.0);
              color=(:transparent,0), strokecolor=col, strokewidth=2.5)
        text!(ax_main, j+1.5, i; text=lbl, fontsize=7, color=col, align=(:left,:center))
    end

    # ── Three inset panels: full P(n) distributions ────────────────────────────
    ns = 0:n_max
    inset_data = [(dist_eq, vox_colors[1], "Voxel A: equilibrated (QSD)"),
                  (dist_fr, vox_colors[2], "Voxel B: mid-ring (building up)"),
                  (dist_nw, vox_colors[3], "Voxel C: leading edge (sparse)")]

    # Row below the main figure
    for (col, (dist, col_c, lbl)) in enumerate(inset_data)
        ax_ins = Axis(pfig[2,col];
            title=lbl, titlesize=9, titlecolor=col_c,
            xlabel="n (molecules)", ylabel="P(n)",
            xlabelsize=8, ylabelsize=8,
            xticklabelsize=7, yticklabelsize=7,
            spinewidth=0.6, ygridvisible=true, ygridcolor=(:gray, 0.3))
        barplot!(ax_ins, ns, dist; color=(col_c, 0.75), strokewidth=0.5, gap=0)
        ylims!(ax_ins, 0, max(maximum(dist)*1.15, 0.01))
        xlims!(ax_ins, -0.5, n_max+0.5)
        # Expected mean annotation
        μ = sum(n*p for (n,p) in zip(ns,dist))
        text!(ax_ins, n_max*0.65, maximum(dist)*0.9;
              text="⟨n⟩=$(round(μ,digits=1))", fontsize=7, color=:gray40)
    end

    rowsize!(pfig.layout, 1, Aspect(1, 0.9))
    rowgap!(pfig.layout, 1, 12)
    colgap!(pfig.layout, 1, 8)
    colgap!(pfig.layout, 2, 8)

    # Colorbar for main heatmap
    Colorbar(pfig[1,2];
             colormap=:plasma, limits=(0, snap_max),
             label="√-scaled ⟨n_k⟩", width=12,
             labelsize=8, ticklabelsize=7)

    pout = joinpath(outdir, "fig_wave_snapshot.pdf")
    save(pout, pfig)
    @printf("Saved paper figure → %s\n", pout)

    # Also as PNG for quick viewing
    save(replace(pout, ".pdf"=>".png"), pfig; px_per_unit=3)
    println("Saved PNG → ", replace(pout, ".pdf"=>".png"))
else
    println("Could not pick voxels for paper figure (snapshot not yet populated)")
end
