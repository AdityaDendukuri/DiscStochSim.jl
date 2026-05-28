"""
anim_annihilation_500.jl — A+B→∅ on 500×500, correct adaptive per-voxel FSP

Implements the SAME three-phase framework as SpatialFSP (single-species):

  EMPTY  — no flux accumulated; zero cost
  ACTIVE — per-voxel 2-species CME via expv (unconditionally stable)
           V-cycle: adjacent active pairs pre-smoothed with JOINT 4-species CME
           capturing intra-pair diffusion correlations before the full solve
  EQUIL  — TV(p_new,p_old) < ε_equil; collapsed to stored (μ_A,μ_B); 1 merged state

Classification is purely flux/TV driven — no species-specific hardcoding.

Birth model sustains visible blobs. Warmup pre-equilibrates source interior
so only the wavefront and interface remain active during the animation.
"""

using ExponentialUtilities, SparseArrays, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const K       = 500
const k_b     = 0.3      # birth rate in source patches
const k_d     = 0.1      # death rate  →  μ_SS = 3,  λ = √(D/k_d) = 10 voxels
const k_ann   = 0.5      # A + B → ∅
const D       = 10.0
const PATCH_R = 5        # patch half-width: (2R+1)² = 121 voxels
const CY      = K÷2      # patch row centre (250)
const CXA     = K÷2 - 15 # patch col centre for A (235)
const CXB     = K÷2 + 15 # patch col centre for B (265)
# patches 30 voxels centre-to-centre, 20 edge-to-edge
# λ=10: blobs extend ~30 voxels from each patch and OVERLAP in the middle

const n_max   = 7        # per-species truncation  (Poisson(3)@n>7: <0.2% truncated)
const NV      = n_max+1  # 8 states per species
const NS      = NV*NV    # 64 states per voxel
const NS4     = NV^4     # 4096 states for joint pair (V-cycle pre-smooth)

const ε_expand = 0.05    # activate empty voxel when accumulated flux > this
const ε_equil  = 0.005   # collapse active voxel when TV < this
const τ_pre    = 0.05    # V-cycle pre-smooth time (equilibrates fast intra-pair diffusion)
const WARMUP   = 20      # warmup steps (not recorded) to equil source interior
const dt_warm  = 3.0     # large dt for warmup (fast equilibration, expv stable)
const dt_fsp   = 1.5     # FSP step for animation
const T_END    = 60.0
const N_FRM    = 100
const KRYLOV_M = 20
const IMGSIZE  = 250

# ── State encoding ─────────────────────────────────────────────────────────────
lin(nA, nB) = nA*NV + nB + 1                               # 1..NS
linj(a,b,c,d) = a*NV^3 + b*NV^2 + c*NV + d + 1            # 1..NS4

function voxel_means(p::AbstractVector)
    μA = μB = 0.0
    for nA in 0:n_max, nB in 0:n_max
        pr = p[lin(nA,nB)];  μA += nA*pr;  μB += nB*pr
    end
    μA, μB
end

function poisson2d(μA, μB)
    p = zeros(NS)
    for nA in 0:n_max, nB in 0:n_max
        p[lin(nA,nB)] = exp(-μA)*μA^nA/factorial(nA) *
                        exp(-μB)*μB^nB/factorial(nB)
    end
    p ./= sum(p)
end

# 8-connected neighbours with isotropic rates:
#   orthogonal (4): D/2,  diagonal (4): D/4
#   D_eff_x = 4*(D/2)*(1)² / 2 + 4*(D/4)*(1/√2)² / 2 ... no — simpler:
#   E[Δx²]/2Δt = D/2*1 + D/2*1 + D/4*(1/2) + D/4*(1/2) = D  ✓
const D_ORTHO = D / 2
const D_DIAG  = D / 4

const _NB8_OFFSETS = ((-1,0,D_ORTHO),(1,0,D_ORTHO),(0,-1,D_ORTHO),(0,1,D_ORTHO),
                      (-1,-1,D_DIAG),(-1,1,D_DIAG),(1,-1,D_DIAG),(1,1,D_DIAG))

# Returns list of ((ni,nj), rate) for all valid 8-connected neighbours
function valid_nb8(i, j)
    nbs = Tuple{Tuple{Int,Int},Float64}[]
    for (di,dj,r) in _NB8_OFFSETS
        ni,nj=i+di,j+dj
        1≤ni≤K&&1≤nj≤K&&push!(nbs,((ni,nj),r))
    end
    nbs
end

const BIRTH_A = [(i-CY)^2+(j-CXA)^2 ≤ PATCH_R^2 ? k_b : 0.0
                 for i in 1:K, j in 1:K]
const BIRTH_B = [(i-CY)^2+(j-CXB)^2 ≤ PATCH_R^2 ? k_b : 0.0
                 for i in 1:K, j in 1:K]

# ── Per-voxel 2-species CME generator ─────────────────────────────────────────
function build_gen(bA, bB, μA_in, μB_in, dA_tot, dB_tot)
    Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(r,c,v)=(push!(Is,r);push!(Js,c);push!(Vs,v))
    for nA in 0:n_max, nB in 0:n_max
        j=lin(nA,nB); diag=0.0
        rAin=bA+μA_in; rBin=bB+μB_in
        nA<n_max&&rAin>0&&(add!(lin(nA+1,nB),j,rAin);diag-=rAin)
        nB<n_max&&rBin>0&&(add!(lin(nA,nB+1),j,rBin);diag-=rBin)
        nA>0&&(r=dA_tot*nA;   add!(lin(nA-1,nB),j,r);   diag-=r)
        nB>0&&(r=dB_tot*nB;   add!(lin(nA,nB-1),j,r);   diag-=r)
        nA>0&&nB>0&&(r=k_ann*nA*nB; add!(lin(nA-1,nB-1),j,r); diag-=r)
        add!(j,j,diag)
    end
    sparse(Is,Js,Vs,NS,NS)
end

# ── Joint 4-species CME generator for V-cycle pair (v1,v2) ────────────────────
# State: (nA1,nB1,nA2,nB2). Captures intra-pair diffusion + all reactions.
# External neighbours contribute via mean-field influx.
function build_joint_gen(v1, v2, ext_nbs1, ext_nbs2, μA, μB, pair_rate)
    # ext_nbsX = list of ((ni,nj), rate) for non-pair neighbours of vX
    μAin1 = sum(r*get(μA,nb,0.0) for (nb,r) in ext_nbs1; init=0.0)
    μBin1 = sum(r*get(μB,nb,0.0) for (nb,r) in ext_nbs1; init=0.0)
    μAin2 = sum(r*get(μA,nb,0.0) for (nb,r) in ext_nbs2; init=0.0)
    μBin2 = sum(r*get(μB,nb,0.0) for (nb,r) in ext_nbs2; init=0.0)
    bA1=BIRTH_A[v1...]; bB1=BIRTH_B[v1...]
    bA2=BIRTH_A[v2...]; bB2=BIRTH_B[v2...]
    # external outflow: death + diffusion to non-pair neighbours (8-connected rates)
    dextA1=k_d+sum(r for (_,r) in ext_nbs1; init=0.0); dextB1=dextA1
    dextA2=k_d+sum(r for (_,r) in ext_nbs2; init=0.0); dextB2=dextA2

    Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(r,c,v)=(push!(Is,r);push!(Js,c);push!(Vs,v))

    for nA1 in 0:n_max,nB1 in 0:n_max,nA2 in 0:n_max,nB2 in 0:n_max
        j=linj(nA1,nB1,nA2,nB2); diag=0.0

        # External birth+influx at v1
        rA1=bA1+μAin1; rB1=bB1+μBin1
        nA1<n_max&&rA1>0&&(add!(linj(nA1+1,nB1,nA2,nB2),j,rA1);diag-=rA1)
        nB1<n_max&&rB1>0&&(add!(linj(nA1,nB1+1,nA2,nB2),j,rB1);diag-=rB1)
        # External outflow at v1
        nA1>0&&(r=dextA1*nA1;add!(linj(nA1-1,nB1,nA2,nB2),j,r);diag-=r)
        nB1>0&&(r=dextB1*nB1;add!(linj(nA1,nB1-1,nA2,nB2),j,r);diag-=r)
        # Annihilation at v1
        nA1>0&&nB1>0&&(r=k_ann*nA1*nB1;add!(linj(nA1-1,nB1-1,nA2,nB2),j,r);diag-=r)

        # External birth+influx at v2
        rA2=bA2+μAin2; rB2=bB2+μBin2
        nA2<n_max&&rA2>0&&(add!(linj(nA1,nB1,nA2+1,nB2),j,rA2);diag-=rA2)
        nB2<n_max&&rB2>0&&(add!(linj(nA1,nB1,nA2,nB2+1),j,rB2);diag-=rB2)
        # External outflow at v2
        nA2>0&&(r=dextA2*nA2;add!(linj(nA1,nB1,nA2-1,nB2),j,r);diag-=r)
        nB2>0&&(r=dextB2*nB2;add!(linj(nA1,nB1,nA2,nB2-1),j,r);diag-=r)
        # Annihilation at v2
        nA2>0&&nB2>0&&(r=k_ann*nA2*nB2;add!(linj(nA1,nB1,nA2-1,nB2-1),j,r);diag-=r)

        # Intra-pair diffusion: A hops v1↔v2  (rate = pair_rate, ortho or diag)
        nA1>0&&nA2<n_max&&(r=pair_rate*nA1;add!(linj(nA1-1,nB1,nA2+1,nB2),j,r);diag-=r)
        nA2>0&&nA1<n_max&&(r=pair_rate*nA2;add!(linj(nA1+1,nB1,nA2-1,nB2),j,r);diag-=r)
        # Intra-pair diffusion: B hops v1↔v2
        nB1>0&&nB2<n_max&&(r=pair_rate*nB1;add!(linj(nA1,nB1-1,nA2,nB2+1),j,r);diag-=r)
        nB2>0&&nB1<n_max&&(r=pair_rate*nB2;add!(linj(nA1,nB1,nA2,nB2-1+1),j,r);diag-=r)

        add!(j,j,diag)
    end
    sparse(Is,Js,Vs,NS4,NS4)
end

# ── Adaptive FSP state ─────────────────────────────────────────────────────────
const dists = Dict{Tuple{Int,Int},Vector{Float64}}()
const eq_μA = Dict{Tuple{Int,Int},Float64}()
const eq_μB = Dict{Tuple{Int,Int},Float64}()
const fluxA = Dict{Tuple{Int,Int},Float64}()
const fluxB = Dict{Tuple{Int,Int},Float64}()

# Initialise source patch voxels as active (distribution at (0,0))
for j in 1:K, i in 1:K
    (BIRTH_A[i,j]>0 || BIRTH_B[i,j]>0) || continue
    p=zeros(NS); p[lin(0,0)]=1.0; dists[(i,j)]=p
end
@printf("Init: %d source voxels active.  NS=%d  NS4=%d\n",
        length(dists), NS, NS4)

# ── V-cycle helpers ────────────────────────────────────────────────────────────
# Greedy matching: pair adjacent active voxels (sorted order for reproducibility)
function greedy_pairs(active_keys)
    used   = Set{Tuple{Int,Int}}()
    pairs  = Tuple{Tuple{Int,Int},Tuple{Int,Int}}[]
    singles= Tuple{Int,Int}[]
    for v in sort(active_keys)
        v ∈ used && continue
        paired = false
        for (nb,_) in valid_nb8(v...)
            (nb ∈ keys(dists)) && (nb ∉ used) || continue
            push!(pairs,(v,nb)); push!(used,v); push!(used,nb)
            paired=true; break
        end
        paired || push!(singles,v)
    end
    pairs, singles
end

# Pre-smooth: evolve each active pair jointly for τ, then extract marginals.
# Captures intra-pair diffusion correlations (the "fast" modes for V-cycle).
function presmooth!(μA, μB, τ)
    pairs, _ = greedy_pairs(collect(keys(dists)))
    for (v1,v2) in pairs
        p1=dists[v1]; p2=dists[v2]
        # Joint IC = product of marginals (valid when voxels recently uncorrelated)
        pj = zeros(NS4)
        for nA1 in 0:n_max,nB1 in 0:n_max,nA2 in 0:n_max,nB2 in 0:n_max
            pj[linj(nA1,nB1,nA2,nB2)] = p1[lin(nA1,nB1)] * p2[lin(nA2,nB2)]
        end
        ext1=[(nb,r) for (nb,r) in valid_nb8(v1...) if nb≠v2]
        ext2=[(nb,r) for (nb,r) in valid_nb8(v2...) if nb≠v1]
        # Intra-pair rate: ortho if same row or col, else diagonal
        pr = (v1[1]==v2[1]||v1[2]==v2[2]) ? D_ORTHO : D_DIAG
        Qj = build_joint_gen(v1,v2,ext1,ext2,μA,μB,pr)
        pj_new = expv(τ, Qj, pj; m=min(KRYLOV_M,NS4))
        pj_new .= max.(0.0,pj_new); pj_new ./= sum(pj_new)
        # Extract marginals
        p1n=zeros(NS); p2n=zeros(NS)
        for nA1 in 0:n_max,nB1 in 0:n_max,nA2 in 0:n_max,nB2 in 0:n_max
            pr=pj_new[linj(nA1,nB1,nA2,nB2)]
            p1n[lin(nA1,nB1)]+=pr; p2n[lin(nA2,nB2)]+=pr
        end
        dists[v1]=p1n; dists[v2]=p2n
    end
end

# ── FSP step with V-cycle ──────────────────────────────────────────────────────
function fsp_step!(dt; vcycle=true)
    # 1. Current means of all non-empty voxels
    μA=Dict{Tuple{Int,Int},Float64}(); μB=Dict{Tuple{Int,Int},Float64}()
    for (ij,p) in dists; a,b=voxel_means(p); μA[ij]=a; μB[ij]=b; end
    for ij in keys(eq_μA); μA[ij]=eq_μA[ij]; μB[ij]=eq_μB[ij]; end

    # 2. V-cycle pre-smooth: joint CME for τ_pre on paired active voxels
    vcycle && presmooth!(μA, μB, τ_pre)

    # 3. Per-voxel full solve + TV-based equil detection
    to_equil = Tuple{Int,Int}[]
    for (ij, p_old) in dists
        i,j=ij; nbs8=valid_nb8(i,j)
        μAin=sum(r*get(μA,nb,0.0) for (nb,r) in nbs8; init=0.0)
        μBin=sum(r*get(μB,nb,0.0) for (nb,r) in nbs8; init=0.0)
        dtot=k_d+sum(r for (_,r) in nbs8; init=0.0)
        Q=build_gen(BIRTH_A[i,j],BIRTH_B[i,j],μAin,μBin,dtot,dtot)
        p_new=expv(dt,Q,p_old;m=min(KRYLOV_M,NS))
        p_new.=max.(0.0,p_new); p_new./=sum(p_new)
        tv=sum(abs.(p_new.-p_old))/2
        if tv<ε_equil
            push!(to_equil,ij)
            a,b=voxel_means(p_new); eq_μA[ij]=a; eq_μB[ij]=b
        else
            dists[ij]=p_new
        end
    end
    for ij in to_equil; delete!(dists,ij); end

    # 4. Propagate diffusion flux into empty voxels (from active+equil)
    for (ij,μa) in μA
        μa<1e-9 && continue
        for (nb,r) in valid_nb8(ij...)
            (haskey(dists,nb)||haskey(eq_μA,nb)) && continue
            fluxA[nb]=get(fluxA,nb,0.0)+r*μa*dt
        end
    end
    for (ij,μb) in μB
        μb<1e-9 && continue
        for (nb,r) in valid_nb8(ij...)
            (haskey(dists,nb)||haskey(eq_μA,nb)) && continue
            fluxB[nb]=get(fluxB,nb,0.0)+r*μb*dt
        end
    end

    # 5. Activate empty voxels where accumulated flux > ε_expand
    for ij in union(keys(fluxA),keys(fluxB))
        haskey(dists,ij)||haskey(eq_μA,ij) && continue
        fa=get(fluxA,ij,0.0); fb=get(fluxB,ij,0.0)
        (fa>ε_expand||fb>ε_expand) || continue
        dists[ij]=poisson2d(min(fa,0.5),min(fb,0.5))
        delete!(fluxA,ij); delete!(fluxB,ij)
    end
end

# ── Warmup: equilibrate source interior without recording ─────────────────────
print("Warmup ($WARMUP × dt=$dt_warm = $(WARMUP*dt_warm)s simulation)…  ")
for _ in 1:WARMUP; fsp_step!(dt_warm; vcycle=false); end  # no V-cycle during warmup
@printf("active=%d  equil=%d\n", length(dists), length(eq_μA))

# ── Visualisation helpers ──────────────────────────────────────────────────────
function mean_grids()
    gA=zeros(Float32,K,K); gB=zeros(Float32,K,K)
    for ((i,j),p) in dists; a,b=voxel_means(p); gA[i,j]=a; gB[i,j]=b; end
    for (ij,μa) in eq_μA; gA[ij...]=μa; gB[ij...]=eq_μB[ij]; end
    gA, gB
end

function phase_grid()
    # 0=empty,1=activeA,2=activeB,3=interface(both),4=equilA,5=equilB,6=equilBoth
    g=zeros(Int8,K,K)
    for ((i,j),p) in dists
        a,b=voxel_means(p)
        g[i,j]=b<0.05*max(a,1e-3) ? 1 : a<0.05*max(b,1e-3) ? 2 : 3
    end
    for (ij,μa) in eq_μA
        μb=eq_μB[ij]
        g[ij...]=μb<0.05*max(μa,1e-3) ? 4 : μa<0.05*max(μb,1e-3) ? 5 : 6
    end
    g
end

# RGB: R=A, G=0, B=B; brightness from combined, clamped at MU_CLAMP
const MU_CLAMP = Float32(0.3)
function combined_rgb(gA, gB)
    img=zeros(RGBf,size(gA))
    for i in eachindex(gA)
        a=min(1f0,gA[i]/MU_CLAMP); b=min(1f0,gB[i]/MU_CLAMP)
        img[i]=RGBf(a,0f0,b)
    end
    img
end

const COL_DIAG =RGBAf(0.27,0.51,0.71,1.0)
const COL_COUPL=RGBAf(0.89,0.45,0.13,1.0)
const COL_ZERO =RGBAf(0.94,0.94,0.94,1.0)
const COL_SEP  =RGBAf(0.20,0.20,0.20,1.0)

function render_matrix(ph)
    iface=Tuple{Int,Int}[]; actA=Tuple{Int,Int}[]; actB=Tuple{Int,Int}[]
    has_eqA=false; has_eqB=false
    for j in 1:K,i in 1:K
        ph[i,j]==3&&push!(iface,(i,j))
        ph[i,j]==1&&push!(actA,(i,j))
        ph[i,j]==2&&push!(actB,(i,j))
        ph[i,j] in (4,6)&&(has_eqA=true)
        ph[i,j] in (5,6)&&(has_eqB=true)
    end
    all_act=vcat(iface,actA,actB); n_act=length(all_act)
    n_eq=(has_eqA ? 1 : 0)+(has_eqB ? 1 : 0)
    total=n_act*NS+n_eq; total==0&&return fill(COL_ZERO,IMGSIZE,IMGSIZE)
    act_idx=Dict(v=>i for (i,v) in enumerate(all_act))
    stov=Vector{Int}(undef,total)
    for (vi,_) in enumerate(all_act); for s in (vi-1)*NS+1:vi*NS; stov[s]=vi; end; end
    off=n_act*NS
    has_eqA && (stov[off+1]=-1)
    has_eqB && (stov[off+(has_eqA ? 2 : 1)]=-2)
    adj=Set{NTuple{2,Int}}(); eq_border=Set{Int}()
    for (vi,(i,j)) in enumerate(all_act)
        for (nb,_) in valid_nb8(i,j)
            nidx=get(act_idx,nb,nothing)
            !isnothing(nidx)&&push!(adj, vi<nidx ? (vi,nidx) : (nidx,vi))
            ph[nb...]>=4&&push!(eq_border,vi)
        end
    end
    img=fill(COL_ZERO,IMGSIZE,IMGSIZE)
    for r in 1:IMGSIZE
        sr=clamp(floor(Int,(r-.5)/IMGSIZE*total)+1,1,total); vr=stov[sr]
        for c in 1:IMGSIZE
            sc=clamp(floor(Int,(c-.5)/IMGSIZE*total)+1,1,total); vc=stov[sc]
            if vr==vc; img[r,c]=COL_DIAG
            elseif vr<0||vc<0; av=vr<0 ? vc : vr; av>0&&av∈eq_border&&(img[r,c]=COL_COUPL)
            else (min(vr,vc),max(vr,vc))∈adj&&(img[r,c]=COL_COUPL)
            end
        end
    end
    n_act>0&&n_eq>0&&(sep=clamp(round(Int,n_act*NS/total*IMGSIZE),1,IMGSIZE);
        for q in 1:IMGSIZE; img[sep,q]=COL_SEP; img[q,sep]=COL_SEP; end)
    img
end

# ── Pre-compute frames ─────────────────────────────────────────────────────────
t_frames   =collect(range(0.0,T_END;length=N_FRM))
rgb_frames =Vector{Matrix{RGBf}}(undef,N_FRM)
ph_frames  =Vector{Matrix{Int8}}(undef,N_FRM)
mat_frames =Vector{Matrix{RGBAf}}(undef,N_FRM)
nact_frames=zeros(Int,N_FRM)

function snap!(f)
    gA,gB=mean_grids(); ph=phase_grid()
    rgb_frames[f]=combined_rgb(gA,gB); ph_frames[f]=ph
    mat_frames[f]=render_matrix(ph); nact_frames[f]=length(dists)
end

snap!(1)
println("Pre-computing $N_FRM frames (with V-cycle)…")
for f in 2:N_FRM
    fsp_step!(dt_fsp; vcycle=true)
    snap!(f)
    f%10==0&&@printf("  frame %3d  t=%5.1f  active=%4d  equil=%6d\n",
                     f,t_frames[f],nact_frames[f],length(eq_μA))
end
println("Peak active: $(maximum(nact_frames))  Final equil: $(length(eq_μA))")

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,Axis=(spinewidth=0.8,xgridvisible=false,ygridvisible=false)))
fig=Figure(size=(1200,530))

rgb_obs=Observable(rgb_frames[1]'); ph_obs=Observable(Float32.(ph_frames[1]'))
mat_obs=Observable(mat_frames[1])

ax_kw=(aspect=DataAspect(),yreversed=true,xaxisposition=:top,titlesize=11,
       xticksvisible=false,yticksvisible=false,
       xticklabelsvisible=false,yticklabelsvisible=false)

ax_rgb=Axis(fig[1,1];
    title="A (red) + B (blue) — FSP, exact for A+B→∅, V-cycle within active ring",
    ax_kw...)
image!(ax_rgb,0.5..K+0.5,0.5..K+0.5,rgb_obs)

PCMAP=cgrad([colorant"#111122",colorant"#c0392b",colorant"#2980b9",colorant"#8e44ad",
             colorant"#d08070",colorant"#70a0d0",colorant"#906090"])
ax_ph=Axis(fig[1,2];title="Adaptive zones (flux/TV-driven, species-agnostic)",ax_kw...)
heatmap!(ax_ph,1:K,1:K,ph_obs;colormap=PCMAP,colorrange=(0,6))

ax_mx=Axis(fig[1,3];title="RDME block matrix (active=NS=$NS, equil=1 merged)",
           aspect=DataAspect(),yreversed=true,xaxisposition=:top,titlesize=11,
           xticksvisible=false,yticksvisible=false,
           xticklabelsvisible=false,yticklabelsvisible=false)
image!(ax_mx,0.5..K+0.5,0.5..K+0.5,mat_obs)

for col in 1:3; colsize!(fig.layout,col,Aspect(1,1.0)); end
t_label =Label(fig[0,1],"t = 0.0";fontsize=13,halign=:left)
kc_label=Label(fig[0,2:3],"";fontsize=11,halign=:center)

Legend(fig[2,1:3],
  [PolyElement(color=c) for c in
   [colorant"#111122",colorant"#c0392b",colorant"#2980b9",colorant"#8e44ad",
    colorant"#d08070",colorant"#70a0d0",RGBf(0.27,0.51,0.71),RGBf(0.89,0.45,0.13)]],
  ["Empty","Active A (wavefront)","Active B (wavefront)","Interface A+B (FSP)",
   "Equil A (collapsed)","Equil B (collapsed)","Diagonal block","Coupling"],
  orientation=:horizontal,framevisible=false,labelsize=8,patchsize=(12,10))

rowgap!(fig.layout,1,5);rowgap!(fig.layout,2,4)
colgap!(fig.layout,1,8);colgap!(fig.layout,2,8)

out=joinpath(@__DIR__,"..","paper","figures","anim_annihilation_$(K)x$(K).gif")
println("Rendering → $out")
@time record(fig,out,1:N_FRM;framerate=15) do f
    rgb_obs[]=rgb_frames[f]'; ph_obs[]=Float32.(ph_frames[f]')
    mat_obs[]=mat_frames[f]
    t_label.text[]="t = $(round(t_frames[f],digits=1))"
    na=nact_frames[f]
    kc_label.text[]="active: $na ($(na*NS) states) | equil: $(length(eq_μA))"
end
println("Saved → $out | Peak active: $(maximum(nact_frames))")
