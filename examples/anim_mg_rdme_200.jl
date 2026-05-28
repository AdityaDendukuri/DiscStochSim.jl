"""
anim_mg_rdme_200.jl — BranchingDeath RDME 200×200 + MG block matrix

Left  — plasma heatmap of ⟨n_k⟩ with white wavefront boundary.
         Cyan rectangle = fixed 20-voxel center-row strip used for right panel.
Right — RDME K_c×K_c block matrix for the strip.
         Diagonal blocks = CME for each MG group (coloured by level).
         Off-diagonal bands = diffusion coupling between adjacent groups.
         Blocks GROW as the equil region expands and K_c shrinks.
"""

using DiscStochSim, CairoMakie, Printf

const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4; const T_END=80.0
const N_FRM=80; const n_max=14; const N_LEV=3
const ci=K÷2+1; const cj=K÷2+1
const IMGSIZE=300

# Strip: 20 voxels in center row, starting 35 voxels ahead of IC center
# (wave KPP ≈ 1.265 → arrives at strip t ≈ 28)
const STRIP_J0  = K÷2+1+35
const STRIP_J1  = STRIP_J0+19
const STRIP_VXS = [CartesianIndex(ci,j) for j in STRIP_J0:STRIP_J1]
const N_STRIP   = length(STRIP_VXS)

# ── Colours ────────────────────────────────────────────────────────────────────
const PLASMA = CairoMakie.to_colormap(:plasma)
plasma_col(t) = let c=PLASMA[clamp(floor(Int,t*(length(PLASMA)-1))+1,1,length(PLASMA))];
                RGBAf(c.r,c.g,c.b,1f0) end
const mu_peak = Ref(1.0)

# Matrix diagonal: coloured by MG level
const MDIAG_EQUIL  = RGBAf(0.28,0.28,0.32,1.0)
const MDIAG_EDGE   = RGBAf(0.95,0.18,0.25,1.0)
const MDIAG_L1     = RGBAf(0.98,0.58,0.10,1.0)
const MDIAG_L2     = RGBAf(0.10,0.82,0.78,1.0)
const MDIAG_L3     = RGBAf(0.20,0.92,0.28,1.0)
const MDIAG_UNPAIR = RGBAf(0.75,0.72,0.08,1.0)
const MAT_COUPL    = RGBAf(0.89,0.45,0.13,1.0)   # darkorange — diffusion
const MAT_GRID     = RGBAf(0.30,0.30,0.30,1.0)
const MAT_BG       = RGBAf(0.04,0.02,0.06,1.0)

mlev_col(lv) = lv==0 ? MDIAG_EQUIL : lv==-1 ? MDIAG_EDGE :
               lv==1  ? MDIAG_L1   : lv==2  ? MDIAG_L2   :
               lv==3  ? MDIAG_L3   : MDIAG_UNPAIR

# ── Inline neighbors ──────────────────────────────────────────────────────────
@inline function gnb(idx::CartesianIndex{2}, Kx, Ky)
    i,j=Tuple(idx); ns=CartesianIndex{2}[]
    i>1  && push!(ns,CartesianIndex(i-1,j)); i<Kx && push!(ns,CartesianIndex(i+1,j))
    j>1  && push!(ns,CartesianIndex(i,j-1)); j<Ky && push!(ns,CartesianIndex(i,j+1))
    ns
end

# ── MG hierarchy builder ──────────────────────────────────────────────────────
function mg_structure(s, n_lev, row_pass)
    aks=collect(keys(s.dists)); as=Set(aks); eq=Set(keys(s.equil_dists))
    interior=CartesianIndex{2}[]; edge=CartesianIndex{2}[]
    for idx in aks
        all(nb->nb in as||nb in eq, gnb(idx,s.K_x,s.K_y)) ?
            push!(interior,idx) : push!(edge,idx)
    end
    Δp=row_pass ? CartesianIndex(0,1) : CartesianIndex(1,0)
    Δq=row_pass ? CartesianIndex(1,0) : CartesianIndex(0,1)
    Δd(ℓ)=isodd(ℓ) ? Δp : Δq
    pd(ℓ)=CartesianIndex(Tuple(Δd(ℓ+1)).*(2^((ℓ-1)÷2)))
    lr=[interior]; lv=[[[v] for v in interior]]; lp=Vector{Tuple{Int,Int}}[]
    for ℓ in 1:n_lev
        prev=lr[end]; isempty(prev) && break
        ridx=Dict(r=>i for (i,r) in enumerate(prev)); δ=pd(ℓ)
        used=Set{Int}(); prs=Tuple{Int,Int}[]; nr=CartesianIndex{2}[]
        nv=Vector{Vector{CartesianIndex{2}}}()
        for (i,rep) in enumerate(prev)
            i in used && continue; j=get(ridx,rep+δ,0)
            (j==0||j in used) && continue
            push!(prs,(i,j)); push!(used,i); push!(used,j)
            push!(nr,min(rep,rep+δ)); push!(nv,vcat(lv[end][i],lv[end][j]))
        end
        isempty(prs) && break
        push!(lp,prs); push!(lr,nr); push!(lv,nv)
    end
    edge, lv, lp
end

function build_vox_maps(lv, lp)
    nb=length(lp)
    depth=Dict{CartesianIndex{2},Int}(); grp=Dict{CartesianIndex{2},Int}(); gc=Ref(0)
    for ℓ in 1:nb
        cov=Set{Int}()
        ℓ<nb && (for (a,b) in lp[ℓ+1]; push!(cov,a); push!(cov,b); end)
        for g in 1:length(lv[ℓ+1])
            g in cov && continue; gc[]+=1; gid=gc[]
            for v in lv[ℓ+1][g]; depth[v]=ℓ; grp[v]=gid; end
        end
    end
    depth, grp
end

# ── Strip → K_c groups + adjacency ────────────────────────────────────────────
function get_strip_info(s, edge_set, depth, grp)
    as=Set(keys(s.dists)); eq=Set(keys(s.equil_dists))
    vkey=Dict{CartesianIndex{2},Any}()
    for v in STRIP_VXS
        v in eq && (vkey[v]=:equil; continue)
        v in edge_set && (vkey[v]=(:e,v); continue)
        v in as && (vkey[v]=haskey(grp,v) ? (:i,grp[v]) : (:u,v))
    end
    ko=Dict{Any,Int}(); gs=Vector{Vector{CartesianIndex{2}}}(); lvs=Int[]
    for v in STRIP_VXS
        haskey(vkey,v)||continue; k=vkey[v]
        if !haskey(ko,k)
            lv=k==:equil ? 0 : (k isa Tuple&&k[1]==:e ? -1 : get(depth,v,0))
            push!(gs,CartesianIndex{2}[]); push!(lvs,lv); ko[k]=length(gs)
        end
        push!(gs[ko[k]],v)
    end
    K_c=length(gs)
    adj=Set{Tuple{Int,Int}}()
    for k in 1:N_STRIP-1
        v1,v2=STRIP_VXS[k],STRIP_VXS[k+1]
        haskey(vkey,v1)&&haskey(vkey,v2)||continue
        g1,g2=ko[vkey[v1]],ko[vkey[v2]]
        g1≠g2 && push!(adj,(min(g1,g2),max(g1,g2)))
    end
    K_c, lvs, sort(collect(adj))
end

# ── Render K_c×K_c → IMGSIZE×IMGSIZE ─────────────────────────────────────────
function render_matrix(K_c, lvs, adj)
    img=fill(MAT_BG, IMGSIZE, IMGSIZE)
    K_c==0 && return img
    cp=IMGSIZE/K_c
    for i in 1:K_c
        r1=floor(Int,(i-1)*cp)+1; r2=min(IMGSIZE,floor(Int,i*cp))
        col=mlev_col(lvs[i])
        for r in r1:r2, c in r1:r2; img[r,c]=col; end
    end
    for (i,j) in adj
        r1=floor(Int,(i-1)*cp)+1; r2=min(IMGSIZE,floor(Int,i*cp))
        c1=floor(Int,(j-1)*cp)+1; c2=min(IMGSIZE,floor(Int,j*cp))
        for r in r1:r2, c in c1:c2; img[r,c]=MAT_COUPL; img[c,r]=MAT_COUPL; end
    end
    if cp>=2
        for k in 0:K_c
            px=min(IMGSIZE,floor(Int,k*cp)+1)
            for q in 1:IMGSIZE; img[px,q]=MAT_GRID; img[q,px]=MAT_GRID; end
        end
    end
    img
end

# ── FrameData ─────────────────────────────────────────────────────────────────
struct FrameData
    heat    :: Matrix{RGBAf}
    bnd     :: Vector{Point2f}
    mat_img :: Matrix{RGBAf}
    K_c     :: Int
end

function make_frame(s)
    # Heatmap
    mg=mean_grid(s); cur=maximum(mg;init=1e-3)
    mu_peak[]=max(mu_peak[],cur*1.02); sc=max(mu_peak[],1e-3)
    Kx,Ky=s.K_x,s.K_y
    heat=[plasma_col(clamp(sqrt(mg[i,j]/sc),0f0,1f0)) for j in 1:Ky, i in 1:Kx]

    # Wavefront boundary
    sg=status_grid(s); xs=Float32[]; ys=Float32[]
    for i in 1:Kx, j in 1:Ky
        sg[i,j]==0 && continue
        i>1   && sg[i-1,j]==0 && (push!(xs,j-.5f0,j+.5f0,NaN32); push!(ys,i-.5f0,i-.5f0,NaN32))
        i<Kx  && sg[i+1,j]==0 && (push!(xs,j-.5f0,j+.5f0,NaN32); push!(ys,i+.5f0,i+.5f0,NaN32))
        j>1   && sg[i,j-1]==0 && (push!(xs,j-.5f0,j-.5f0,NaN32); push!(ys,i-.5f0,i+.5f0,NaN32))
        j<Ky  && sg[i,j+1]==0 && (push!(xs,j+.5f0,j+.5f0,NaN32); push!(ys,i-.5f0,i+.5f0,NaN32))
    end

    # MG block matrix for strip
    edge, lv, lp = mg_structure(s, N_LEV, true)
    edge_set = Set(edge)
    depth, grp = build_vox_maps(lv, lp)
    K_c, lvs, adj = get_strip_info(s, edge_set, depth, grp)
    mat = render_matrix(K_c, lvs, adj)

    FrameData(heat, Point2f.(xs,ys), mat, K_c)
end

# ── State ─────────────────────────────────────────────────────────────────────
model = BranchingDeathRDME(k_b,k_d,D,1.0)
state = SpatialFSP(model,K,K;
    n_max, ε_expand=999.0, ε_equil=0.15,
    organic_activation=true, ε_organic=1e-6,
    use_multigrid=true, multigrid_levels=N_LEV, use_block_vcycle=false)
for di in -4:4, dj in -4:4
    di^2+dj^2<=16||continue; 1<=ci+di<=K&&1<=cj+dj<=K||continue
    set_ic!(state,CartesianIndex(ci+di,cj+dj),N0)
end

v_kpp=2*sqrt(D*(k_b-k_d))
@printf("Strip: row %d  cols %d‥%d  v_KPP≈%.3f  wave arrives t≈%.1f\n\n",
        ci, STRIP_J0, STRIP_J1, v_kpp, (STRIP_J0-cj)/v_kpp)

# ── Pre-compute ───────────────────────────────────────────────────────────────
t_frames=collect(range(0.0,T_END;length=N_FRM))
frames=Vector{FrameData}(undef,N_FRM)
frames[1]=make_frame(state)

println("Pre-computing $N_FRM frames…"); flush(stdout)
t0=time()
let fi=2
    for si in 1:round(Int,T_END/dt)
        step!(state,dt); t=si*dt
        if fi<=N_FRM && t>=t_frames[fi]-1e-9
            frames[fi]=make_frame(state)
            if fi%10==0
                @printf("  frame %3d/%d  t=%5.1f  Kc=%2d  active=%5d  equil=%5d  %.0fs\n",
                        fi,N_FRM,t,frames[fi].K_c,length(state.dists),length(state.equil_dists),time()-t0)
                flush(stdout)
            end
            fi+=1
        end
        fi>N_FRM && break
    end
end
@printf("Done in %.1fs\n\n",time()-t0)

# ── Figure ────────────────────────────────────────────────────────────────────
fig=Figure(size=(1060,560), backgroundcolor=:black)

ax_sp=Axis(fig[1,1]; aspect=DataAspect(), yreversed=true, backgroundcolor=:black,
    xticksvisible=false, yticksvisible=false, xticklabelsvisible=false,
    yticklabelsvisible=false, leftspinevisible=false, rightspinevisible=false,
    topspinevisible=false, bottomspinevisible=false)
ax_mx=Axis(fig[1,2]; aspect=DataAspect(), yreversed=true,
    backgroundcolor=RGBf(0.04,0.02,0.06),
    xticksvisible=false, yticksvisible=false, xticklabelsvisible=false,
    yticklabelsvisible=false, leftspinevisible=false, rightspinevisible=false,
    topspinevisible=false, bottomspinevisible=false)

heat_obs = Observable(frames[1].heat)
bnd_obs  = Observable(frames[1].bnd)
mat_obs  = Observable(frames[1].mat_img)

image!(ax_sp, 0.5..K+0.5, 0.5..K+0.5, heat_obs)
linesegments!(ax_sp, bnd_obs; color=(:white,0.75), linewidth=0.8)
linesegments!(ax_sp, Point2f[
    (STRIP_J0-.5f0,ci-1.5f0),(STRIP_J1+.5f0,ci-1.5f0),
    (STRIP_J1+.5f0,ci-1.5f0),(STRIP_J1+.5f0,ci+1.5f0),
    (STRIP_J1+.5f0,ci+1.5f0),(STRIP_J0-.5f0,ci+1.5f0),
    (STRIP_J0-.5f0,ci+1.5f0),(STRIP_J0-.5f0,ci-1.5f0)];
    color=:cyan, linewidth=1.8)

image!(ax_mx, 0.5..IMGSIZE+0.5, 0.5..IMGSIZE+0.5, mat_obs)

Label(fig[0,1], "⟨n_k⟩  heatmap  (cyan = center-row strip)";
      fontsize=10, color=:gray90, halign=:center)
Label(fig[0,2],
    "RDME block matrix  (diag=CME · orange=diffusion · blocks grow as equil expands)";
      fontsize=9, color=:gray70, halign=:center)

t_lbl  = Label(fig[2,1], "t = 0.0"; fontsize=11, color=:gray70, halign=:left)
kc_lbl = Label(fig[2,2], @sprintf("K_c = %d", frames[1].K_c);
               fontsize=11, color=:gray70, halign=:center)

Legend(fig[3,1:2],
    [PolyElement(color=RGBf(0.28,0.28,0.32)), PolyElement(color=RGBf(0.95,0.18,0.25)),
     PolyElement(color=RGBf(0.98,0.58,0.10)), PolyElement(color=RGBf(0.10,0.82,0.78)),
     PolyElement(color=RGBf(0.20,0.92,0.28)), PolyElement(color=RGBf(0.89,0.45,0.13))],
    ["equil","edge (per-voxel)","L=1 pair","L=2 quad","L=3 group","diffusion coupling"],
    orientation=:horizontal, framevisible=false,
    labelsize=8, patchsize=(14,10), labelcolor=:gray80)

rowgap!(fig.layout,1,2); rowgap!(fig.layout,2,4); rowgap!(fig.layout,3,4)
colgap!(fig.layout,1,10)
colsize!(fig.layout,1,Aspect(1,1.0)); colsize!(fig.layout,2,Aspect(1,1.0))

outdir=joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)
out=joinpath(outdir,"anim_mg_rdme_200.gif")
println("Rendering → $out"); flush(stdout)

@time record(fig,out,1:N_FRM;framerate=15) do f
    fr=frames[f]
    heat_obs[]=fr.heat; bnd_obs[]=fr.bnd; mat_obs[]=fr.mat_img
    t_lbl.text[]  = @sprintf("t = %4.1f", t_frames[f])
    kc_lbl.text[] = @sprintf("K_c = %d  groups", fr.K_c)
end
println("Saved → $out")
