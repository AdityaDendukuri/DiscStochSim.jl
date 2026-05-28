"""
anim_nucleus_fsp_100.jl  —  A+B→2A catalytic invasion, graph-mesh FSP

Model (species-agnostic solver, spatial heterogeneity only in parameters):
  ∅ → A   at rate k_bA_source  (nonzero only at source patch)
  A → ∅   at rate k_dA × n_A
  A+B→2A  at rate k_cat × n_A × n_B  (B always consumed; A created if room)
  Diffusion A: D_A,   Diffusion B: D_B

IC: B pre-loaded inside nuclear boundary (geometry only, not a solver branch).
    A starts as single point at source patch in cytoplasm.

FSP three-phase adaptive state space (uniform criterion, no species checks):
  EMPTY  → accumulated flux > ε_expand  →  ACTIVE
  ACTIVE → TV(p_new,p_old) < ε_equil   →  EQUIL  (stored μA, μB)
  EQUIL  → re-activated when A flux arrives at a B-containing equil'd voxel

Nuclear voxels equil immediately at rest (TV≈0, uniform B).
They re-activate when the A invasion front arrives (TV spikes).
After A converts all B: equil again at new SS.

Nuclear membrane = absent voxels.  Pore = gap in absent ring.
"""

using DiscStochSim, ExponentialUtilities, SparseArrays, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const K        = 150
const CX       = K÷2+1;  const CY = K÷2+1
const R_CELL   = 67.0
const R_NUC    = 22.0
const MEMBRANE = 1.5
const PORE_ANG  = 0.0
const PORE_HALF = π/5

const k_bA_source = 5.0   # birth rate at source patch
const k_dA        = 0.001  # λ = √(D_A/k_dA) = √500 ≈ 22 vox
const k_cat       = 0.02   # A+B→2A (nuclear wave visible over ~500s)
const D_A         = 0.5
const D_B         = 0.0
const N_B0        = 3

const PATCH_ROW   = CX
const PATCH_COL   = CY - round(Int, R_CELL*0.55)  # west cytoplasm, far from pore
const PATCH_R     = 3

const n_max   = 5
const NV      = n_max+1;  const NS = NV*NV;  const NS4 = NV^4

const ε_expand   = 1e-5   # near-zero: all cytoplasm activates → smooth gradient
const ε_equil    = 0.01
const ε_react    = 1e-4
const τ_pre      = 0.05
const VCYCLE_MAX = 0
const dt         = 2.0
const T_END      = 1200.0
const N_FRM      = 100
const KRYLOV_M   = 20

# ── State encoding ─────────────────────────────────────────────────────────────
lin(nA,nB)   = nA*NV + nB + 1
lj(a,b,c,d) = a*NV^3 + b*NV^2 + c*NV + d + 1

function voxel_means(p)
    μA=μB=0.0
    for nA in 0:n_max, nB in 0:n_max
        pr=p[lin(nA,nB)]; μA+=nA*pr; μB+=nB*pr
    end
    μA,μB
end

# ── Graph (membrane = absent voxels) ──────────────────────────────────────────
function voxel_exists(i,j)
    r=sqrt((i-CX)^2+(j-CY)^2); r>=R_CELL&&return false
    if abs(r-R_NUC)<=MEMBRANE
        mi=i-CX; mj=j-CY
        θ=atan(Float64(mi),Float64(mj)); Δ=abs(θ-PORE_ANG); Δ=min(Δ,2π-Δ)
        return Δ<PORE_HALF
    end
    true
end

const VOXELS_SET = Set(CartesianIndex(i,j)
    for i in 1:K, j in 1:K if voxel_exists(i,j))
const NBRS = Dict(v => CartesianIndex{2}[
    nb for (di,dj) in ((-1,0),(1,0),(0,-1),(0,1))
    for nb in [CartesianIndex(v[1]+di,v[2]+dj)]
    if nb ∈ VOXELS_SET] for v in VOXELS_SET)

# ── Per-voxel birth rate (only at source patch) ────────────────────────────────
bA_at(i,j) = (i-PATCH_ROW)^2+(j-PATCH_COL)^2 <= PATCH_R^2 ? k_bA_source : 0.0
in_nuc(i,j)  = sqrt((i-CX)^2+(j-CY)^2) < R_NUC-MEMBRANE   # IC only

@printf("Voxels=%d  source→pore≈%d  λ≈%.1f  NS=%d\n",
        length(VOXELS_SET),
        round(Int,sqrt((PATCH_ROW-CX)^2+(PATCH_COL-CY-R_NUC)^2)),
        sqrt(D_A/k_dA), NS)

# ── Per-voxel 2-species generator ─────────────────────────────────────────────
function build_gen(bA, μA_in, μB_in, dA_diff, dB_diff)
    Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(r,c,v)=(push!(Is,r);push!(Js,c);push!(Vs,v))
    for nA in 0:n_max, nB in 0:n_max
        j=lin(nA,nB); diag=0.0
        arr=bA+μA_in
        nA<n_max&&arr>0&&(add!(lin(nA+1,nB),j,arr);diag-=arr)
        nA>0&&(r=k_dA*nA;add!(lin(nA-1,nB),j,r);diag-=r)
        if nA>0&&nB>0
            r=k_cat*nA*nB; diag-=r
            nA<n_max ? add!(lin(nA+1,nB-1),j,r) : add!(lin(nA,nB-1),j,r)
        end
        nB<n_max&&μB_in>0&&(add!(lin(nA,nB+1),j,μB_in);diag-=μB_in)
        nA>0&&(r=dA_diff*nA;add!(lin(nA-1,nB),j,r);diag-=r)
        nB>0&&(r=dB_diff*nB;add!(lin(nA,nB-1),j,r);diag-=r)
        add!(j,j,diag)
    end
    sparse(Is,Js,Vs,NS,NS)
end

# ── Joint generator for V-cycle ────────────────────────────────────────────────
function build_joint_gen(v1,v2,μA_all,μB_all)
    ext1=[nb for nb in NBRS[v1] if nb≠v2]
    ext2=[nb for nb in NBRS[v2] if nb≠v1]
    μAin1=D_A*sum(get(μA_all,nb,0.0) for nb in ext1;init=0.0)
    μBin1=D_B*sum(get(μB_all,nb,0.0) for nb in ext1;init=0.0)
    μAin2=D_A*sum(get(μA_all,nb,0.0) for nb in ext2;init=0.0)
    μBin2=D_B*sum(get(μB_all,nb,0.0) for nb in ext2;init=0.0)
    bv1=bA_at(v1[1],v1[2]); bv2=bA_at(v2[1],v2[2])
    dAd1=D_A*length(ext1); dBd1=D_B*length(ext1)
    dAd2=D_A*length(ext2); dBd2=D_B*length(ext2)
    Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(r,c,v)=(push!(Is,r);push!(Js,c);push!(Vs,v))
    for a1 in 0:n_max,b1 in 0:n_max,a2 in 0:n_max,b2 in 0:n_max
        j=lj(a1,b1,a2,b2); diag=0.0
        r1=bv1+μAin1; a1<n_max&&r1>0&&(add!(lj(a1+1,b1,a2,b2),j,r1);diag-=r1)
        a1>0&&(r=k_dA*a1;add!(lj(a1-1,b1,a2,b2),j,r);diag-=r)
        if a1>0&&b1>0; r=k_cat*a1*b1; diag-=r
            a1<n_max ? add!(lj(a1+1,b1-1,a2,b2),j,r) : add!(lj(a1,b1-1,a2,b2),j,r); end
        r2=bv2+μAin2; a2<n_max&&r2>0&&(add!(lj(a1,b1,a2+1,b2),j,r2);diag-=r2)
        a2>0&&(r=k_dA*a2;add!(lj(a1,b1,a2-1,b2),j,r);diag-=r)
        if a2>0&&b2>0; r=k_cat*a2*b2; diag-=r
            a2<n_max ? add!(lj(a1,b1,a2+1,b2-1),j,r) : add!(lj(a1,b1,a2,b2-1),j,r); end
        b1<n_max&&μBin1>0&&(add!(lj(a1,b1+1,a2,b2),j,μBin1);diag-=μBin1)
        b2<n_max&&μBin2>0&&(add!(lj(a1,b1,a2,b2+1),j,μBin2);diag-=μBin2)
        a1>0&&(r=dAd1*a1;add!(lj(a1-1,b1,a2,b2),j,r);diag-=r)
        b1>0&&(r=dBd1*b1;add!(lj(a1,b1-1,a2,b2),j,r);diag-=r)
        a2>0&&(r=dAd2*a2;add!(lj(a1,b1,a2-1,b2),j,r);diag-=r)
        b2>0&&(r=dBd2*b2;add!(lj(a1,b1,a2,b2-1),j,r);diag-=r)
        a1>0&&a2<n_max&&(r=D_A*a1;add!(lj(a1-1,b1,a2+1,b2),j,r);diag-=r)
        a2>0&&a1<n_max&&(r=D_A*a2;add!(lj(a1+1,b1,a2-1,b2),j,r);diag-=r)
        b1>0&&b2<n_max&&(r=D_B*b1;add!(lj(a1,b1-1,a2,b2+1),j,r);diag-=r)
        b2>0&&b1<n_max&&(r=D_B*b2;add!(lj(a1,b1+1,a2,b2-1),j,r);diag-=r)
        add!(j,j,diag)
    end
    sparse(Is,Js,Vs,NS4,NS4)
end

# ── Adaptive FSP state ─────────────────────────────────────────────────────────
const dists = Dict{CartesianIndex{2},Vector{Float64}}()
const eq_μA = Dict{CartesianIndex{2},Float64}()
const eq_μB = Dict{CartesianIndex{2},Float64}()
const flux  = Dict{CartesianIndex{2},Float64}()

# IC: nuclear voxels start active (B loaded, equil quickly);
#     cytoplasm starts EMPTY (activates via A flux from source)
for v in VOXELS_SET
    i,j=Tuple(v)
    in_nuc(i,j) || continue   # cytoplasm: empty (activated by A flux)
    p=zeros(NS); p[lin(0,N_B0)]=1.0; dists[v]=p
end
# Source: 1 A molecule at one patch voxel
let vlist=collect(VOXELS_SET)
    idx=findfirst(v->(v[1]-PATCH_ROW)^2+(v[2]-PATCH_COL)^2<=PATCH_R^2, vlist)
    if !isnothing(idx)
        sv=vlist[idx]; p=zeros(NS); p[lin(1,0)]=1.0; dists[sv]=p
    end
end

# ── V-cycle helpers ────────────────────────────────────────────────────────────
function greedy_pairs(ks)
    used=Set{CartesianIndex{2}}(); pairs=Tuple{CartesianIndex{2},CartesianIndex{2}}[]
    for v in sort(ks;by=Tuple); v∈used&&continue; paired=false
        for nb in NBRS[v]; (nb∈keys(dists))&&(nb∉used)||continue
            push!(pairs,(v,nb)); push!(used,v); push!(used,nb); paired=true; break
        end
    end; pairs
end

function presmooth!(μA_all,μB_all)
    for (v1,v2) in greedy_pairs(collect(keys(dists)))
        p1=dists[v1]; p2=dists[v2]
        pj=zeros(NS4)
        for a in 0:n_max,b in 0:n_max,c in 0:n_max,d in 0:n_max
            pj[lj(a,b,c,d)]=p1[lin(a,b)]*p2[lin(c,d)]
        end
        Qj=build_joint_gen(v1,v2,μA_all,μB_all)
        pjn=expv(τ_pre,Qj,pj;m=min(KRYLOV_M,NS4))
        pjn.=max.(0,pjn); pjn./=sum(pjn)
        p1n=zeros(NS); p2n=zeros(NS)
        for a in 0:n_max,b in 0:n_max,c in 0:n_max,d in 0:n_max
            pr=pjn[lj(a,b,c,d)]; p1n[lin(a,b)]+=pr; p2n[lin(c,d)]+=pr
        end
        dists[v1]=p1n; dists[v2]=p2n
    end
end

# ── FSP step ───────────────────────────────────────────────────────────────────
function fsp_step!(dt_s)
    μA_all=Dict{CartesianIndex{2},Float64}()
    μB_all=Dict{CartesianIndex{2},Float64}()
    for (v,p) in dists; a,b=voxel_means(p); μA_all[v]=a; μB_all[v]=b; end
    for v in keys(eq_μA); μA_all[v]=eq_μA[v]; μB_all[v]=eq_μB[v]; end

    length(dists)<=VCYCLE_MAX && presmooth!(μA_all,μB_all)

    # Per-voxel solve — equil on TV only, no species-specific thresholds
    to_eq=CartesianIndex{2}[]
    for (v,p_old) in dists
        i,j=Tuple(v); nbs=NBRS[v]; n_nb=length(nbs)
        μAin=D_A*sum(get(μA_all,nb,0.0) for nb in nbs;init=0.0)
        μBin=D_B*sum(get(μB_all,nb,0.0) for nb in nbs;init=0.0)
        Q=build_gen(bA_at(i,j),μAin,μBin,D_A*n_nb,D_B*n_nb)
        p_new=expv(dt_s,Q,p_old;m=min(KRYLOV_M,NS))
        p_new.=max.(0,p_new); p_new./=sum(p_new)
        if sum(abs.(p_new.-p_old))/2 < ε_equil
            a,b=voxel_means(p_new); push!(to_eq,v); eq_μA[v]=a; eq_μB[v]=b
        else
            dists[v]=p_new
        end
    end
    for v in to_eq; delete!(dists,v); end

    # Flux propagation: empty → activate; B-equil'd → re-activate when A arrives
    for (v,μa) in μA_all
        μa<1e-9&&continue
        for nb in NBRS[v]
            haskey(dists,nb)&&continue
            if haskey(eq_μA,nb)
                # Re-activate B-containing equil'd voxel when A arrives
                eq_μB[nb]>0.1&&D_A*μa*dt_s>ε_react || continue
                # Initialize from Poisson(eq_μB) so fractional B depletion accumulates
                μb0=clamp(eq_μB[nb],0.0,Float64(n_max))
                p=zeros(NS)
                for nB in 0:n_max; p[lin(0,nB)]=exp(-μb0)*μb0^nB/factorial(nB); end
                p./=sum(p)
                dists[nb]=p; delete!(eq_μA,nb); delete!(eq_μB,nb)
            else
                flux[nb]=get(flux,nb,0.0)+D_A*μa*dt_s
            end
        end
    end
    for (v,f) in flux
        haskey(dists,v)||haskey(eq_μA,v)&&continue; f>ε_expand||continue
        p=zeros(NS); p[lin(0,0)]=1.0; dists[v]=p; delete!(flux,v)
    end
end

# ── Rendering ──────────────────────────────────────────────────────────────────
const C_WALL=RGBf(0.15,0.10,0.05); const C_OUT=RGBf(0,0,0)
const C_EMPTY_CELL=RGBf(0.0,0.0,0.0)   # black: blends seamlessly with zero-concentration
const C_ACT=RGBf(0.88,0.48,0.22); const C_EQ=RGBf(0.27,0.45,0.72)
const COL_DIAG=RGBAf(0.27,0.51,0.71,1.0); const COL_COUP=RGBAf(0.89,0.45,0.13,1.0)
const COL_ZERO=RGBAf(0.94,0.94,0.94,1.0); const COL_SEP=RGBAf(0.20,0.20,0.20,1.0)
const _MAGMA=CairoMakie.to_colormap(:magma)
_magma(t)=let c=_MAGMA[clamp(floor(Int,Float64(t)*(length(_MAGMA)-1))+1,1,length(_MAGMA))]; RGBf(c.r,c.g,c.b) end
const _mu_maxA=Ref(0.1); const _mu_maxB=Ref(Float64(N_B0))

function make_img(μA_all,μB_all)
    img=fill(C_OUT,K,K)
    for j in 1:K,i in 1:K
        sqrt((i-CX)^2+(j-CY)^2)>=R_CELL&&continue
        v=CartesianIndex(i,j)
        if v∉VOXELS_SET; img[j,i]=C_WALL; continue; end
        μa=get(μA_all,v,0.0); μb=get(μB_all,v,0.0)
        if μa<1e-6&&μb<1e-6; img[j,i]=C_EMPTY_CELL; continue; end
        a=Float32(clamp(sqrt(μa/_mu_maxA[]),0,1))  # √ gamma: reveals gradient
        b=Float32(clamp(μb/_mu_maxB[],0,1))
        img[j,i]=RGBf(a,0f0,b)
    end
    img
end

function make_ph()
    img=fill(C_OUT,K,K)
    for j in 1:K,i in 1:K
        sqrt((i-CX)^2+(j-CY)^2)>=R_CELL&&continue
        v=CartesianIndex(i,j)
        img[j,i]=v∉VOXELS_SET ? C_WALL :
                 haskey(dists,v) ? C_ACT :
                 haskey(eq_μA,v) ? C_EQ : C_EMPTY_CELL
    end
    img
end

function render_mat()
    act=sort!(collect(keys(dists));by=Tuple)
    n_a=length(act); has_eq=!isempty(eq_μA)
    total=n_a*NS+(has_eq ? 1 : 0); total==0&&return fill(COL_ZERO,250,250)
    aidx=Dict(v=>i for (i,v) in enumerate(act))
    stov=Vector{Int}(undef,total)
    for (vi,_) in enumerate(act); for s in (vi-1)*NS+1:vi*NS; stov[s]=vi; end; end
    has_eq&&(stov[total]=0)
    adj=Set{NTuple{2,Int}}(); eq_b=Set{Int}()
    for (vi,v) in enumerate(act), nb in NBRS[v]
        nidx=get(aidx,nb,nothing)
        !isnothing(nidx)&&push!(adj, vi<nidx ? (vi,nidx) : (nidx,vi))
        haskey(eq_μA,nb)&&push!(eq_b,vi)
    end
    img=fill(COL_ZERO,250,250)
    for r in 1:250
        sr=clamp(floor(Int,(r-.5)/250*total)+1,1,total); vr=stov[sr]
        for c in 1:250
            sc=clamp(floor(Int,(c-.5)/250*total)+1,1,total); vc=stov[sc]
            if vr==vc; img[r,c]=COL_DIAG
            elseif vr==0||vc==0; av=vr==0 ? vc : vr; av>0&&av∈eq_b&&(img[r,c]=COL_COUP)
            else (min(vr,vc),max(vr,vc))∈adj&&(img[r,c]=COL_COUP)
            end
        end
    end
    n_a>0&&has_eq&&(sep=clamp(round(Int,n_a*NS/total*250),1,250);
        for q in 1:250; img[sep,q]=COL_SEP; img[q,sep]=COL_SEP; end)
    img
end

# ── Pre-compute frames ─────────────────────────────────────────────────────────
t_frames=collect(range(0.0,T_END;length=N_FRM))
ab_frames=Vector{Matrix{RGBf}}(undef,N_FRM)
ph_frames=Vector{Matrix{RGBf}}(undef,N_FRM)
mat_frames=Vector{Matrix{RGBAf}}(undef,N_FRM)
nact_frames=zeros(Int,N_FRM)

function snap!(f)
    μA_all=Dict{CartesianIndex{2},Float64}(); μB_all=Dict{CartesianIndex{2},Float64}()
    for (v,p) in dists; a,b=voxel_means(p); μA_all[v]=a; μB_all[v]=b; end
    for v in keys(eq_μA); μA_all[v]=eq_μA[v]; μB_all[v]=eq_μB[v]; end
    _mu_maxA[]=isempty(μA_all) ? 0.1 : max(maximum(values(μA_all))*1.1, 0.1)
    ab_frames[f]=make_img(μA_all,μB_all)
    ph_frames[f]=make_ph(); mat_frames[f]=render_mat(); nact_frames[f]=length(dists)
end

snap!(1)
println("Pre-computing $N_FRM frames…  Voxels=$(length(VOXELS_SET))  NS=$NS")
let fi=2
    for step in 1:round(Int,T_END/dt)
        fsp_step!(dt); t=step*dt
        if fi<=N_FRM&&t>=t_frames[fi]-1e-9
            snap!(fi)
            fi%10==0&&@printf("  frame %3d  t=%6.1f  active=%4d  equil=%5d\n",
                              fi,t,nact_frames[fi],length(eq_μA))
            fi+=1
        end
        fi>N_FRM&&break
    end
end
println("Peak active: $(maximum(nact_frames))")

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,Axis=(spinewidth=0.8,xgridvisible=false,ygridvisible=false)))
fig=Figure(size=(1200,530))
ab_obs=Observable(ab_frames[1]); ph_obs=Observable(ph_frames[1]); mat_obs=Observable(mat_frames[1])
ax_kw=(aspect=DataAspect(),yreversed=true,xaxisposition=:top,titlesize=10,
       xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false)
ax_ab=Axis(fig[1,1];title="A (red) diffuses → finds B (blue) → A+B→2A catalytic takeover",ax_kw...)
image!(ax_ab,0.5..K+0.5,0.5..K+0.5,ab_obs)
ax_ph=Axis(fig[1,2];title="Adaptive FSP state (orange=active wavefront, blue=equil)",ax_kw...)
image!(ax_ph,0.5..K+0.5,0.5..K+0.5,ph_obs)
ax_mx=Axis(fig[1,3];title="RDME block matrix (active=NS=$NS states, equil=1)",
    aspect=DataAspect(),yreversed=true,xaxisposition=:top,titlesize=10,
    xticksvisible=false,yticksvisible=false,xticklabelsvisible=false,yticklabelsvisible=false)
image!(ax_mx,0.5..K+0.5,0.5..K+0.5,mat_obs)
for col in 1:3; colsize!(fig.layout,col,Aspect(1,1.0)); end
t_label=Label(fig[0,1],"t=0";fontsize=13,halign=:left)
kc_label=Label(fig[0,2:3],"";fontsize=11,halign=:center)
Legend(fig[2,1:3],
  [PolyElement(color=c) for c in [C_WALL,C_EMPTY_CELL,C_ACT,C_EQ,
   RGBf(0.27,0.51,0.71),RGBf(0.89,0.45,0.13)]],
  ["Membrane","Empty cell","Active (CME+V-cycle)","Equil (1 merged state)",
   "Diagonal","Coupling"],
  orientation=:horizontal,framevisible=false,labelsize=9,patchsize=(13,10))
rowgap!(fig.layout,1,5);rowgap!(fig.layout,2,4);colgap!(fig.layout,1,8);colgap!(fig.layout,2,8)

out=joinpath(@__DIR__,"..","paper","figures","anim_nucleus_fsp_100.gif")
println("Rendering → $out")
@time record(fig,out,1:N_FRM;framerate=15) do f
    ab_obs[]=ab_frames[f]; ph_obs[]=ph_frames[f]; mat_obs[]=mat_frames[f]
    t_label.text[]="t=$(round(t_frames[f],digits=0))"
    na=nact_frames[f]
    kc_label.text[]="active: $na ($(na*NS) FSP states) | equil: $(length(eq_μA))"
end
println("Saved | Peak active: $(maximum(nact_frames))")
