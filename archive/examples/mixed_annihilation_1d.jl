"""
Mixed-Resolution Multigrid — A+B→∅ Annihilation Front
======================================================

Model:
  ∅  →(k_A)   A   [source: left voxels, k in src_A]
  ∅  →(k_B)   B   [source: right voxels, k in src_B]
  A+B →(k_rxn) ∅  [bimolecular annihilation, all voxels]
  + diffusion DA, DB

Steady state: A decays from left, B decays from right, they annihilate
near the center.  The reaction FRONT is a narrow zone (width ≈ sqrt(DA/k_rxn))
where BOTH species coexist and fluctuations matter most.

Mixed-resolution strategy:
  Active FSP:  the reaction front (both μA and μB appreciable)
               → full joint distribution captures A-B correlations and fluctuations
  Mean-field:  A-only zone (left), B-only zone (right)
               → mean-field closure k_rxn·μA·μB is exact when one species ≈ 0

2D visualization: extrude the 1D profile (uniform in y) into a 2D heatmap,
showing the reaction front and active window.

Scaling argument:
  For a K_x × K_y grid with the front perpendicular to x:
  Active zone ≈ 2-4 voxels in x × K_y voxels in y.
  If y-dynamics are homogeneous (sources span full y), each 1D x-slice is
  independent → scale to K_y fine voxels by multiplying k_A, k_B by K_y.
  The mixed FSP size stays O((n_max+1)^(2·Ka)) regardless of K_y!

Run: julia --project examples/mixed_annihilation_1d.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# ODE reference for A+B→∅
# ─────────────────────────────────────────────────────────────────────────────

function run_ode_annihilation(K, model::AnnihilationModel1D, dx,
                              T_END, DT; source_A, source_B)
    dA=model.D_A/dx^2; dB=model.D_B/dx^2
    k_A=model.k_A; k_B=model.k_B; k_rxn=model.k_rxn

    μA=zeros(K); μB=zeros(K)
    # Initialize: uniform across source regions
    for k in source_A; μA[k]=k_A/k_rxn; end   # rough IC at SS for A alone
    for k in source_B; μB[k]=k_B/k_rxn; end

    max_d=max(dA,dB)
    max_rate=2*max_d + k_A + k_B + k_rxn*max(maximum(μA),maximum(μB),1.0)
    n_sub=max(1,ceil(Int,max_rate*DT*2)); dt_sub=DT/n_sub

    n_steps=round(Int,T_END/DT)
    times=Float64[0.0]; μA_hist=[copy(μA)]; μB_hist=[copy(μB)]
    ΔA=zeros(K); ΔB=zeros(K)
    for _ in 1:n_steps
        for _ in 1:n_sub
            fill!(ΔA,0.0); fill!(ΔB,0.0)
            for k in 1:K
                rxn=k_rxn*μA[k]*μB[k]
                ΔA[k]=(k in source_A ? k_A : 0.0)-rxn
                ΔB[k]=(k in source_B ? k_B : 0.0)-rxn
                k>1&&(ΔA[k]+=dA*(μA[k-1]-μA[k]);ΔB[k]+=dB*(μB[k-1]-μB[k]))
                k<K&&(ΔA[k]+=dA*(μA[k+1]-μA[k]);ΔB[k]+=dB*(μB[k+1]-μB[k]))
            end
            μA.=max.(0.0,μA.+ΔA.*dt_sub); μB.=max.(0.0,μB.+ΔB.*dt_sub)
        end
        push!(times,length(times)*DT)
        push!(μA_hist,copy(μA)); push!(μB_hist,copy(μB))
    end
    times, μA_hist, μB_hist
end

# ─────────────────────────────────────────────────────────────────────────────
# Mixed-resolution runner for annihilation
# ─────────────────────────────────────────────────────────────────────────────

function run_mixed_annihilation(K, model::AnnihilationModel1D, dx,
                                 n_max, T_END, DT;
                                 source_A, source_B,
                                 k_lo_init, k_hi_init,
                                 promote_thresh_A=0.1,
                                 promote_thresh_B=0.1,
                                 demote_thresh_A=0.02,
                                 demote_thresh_B=0.02,
                                 min_k_act=2, max_k_act=4,
                                 n_vcycle_levels=1,
                                 prune_tol=1e-12, krylov_m=30)

    μA=zeros(K); μB=zeros(K)
    for k in source_A; μA[k]=model.k_A/model.k_rxn*0.5; end
    for k in source_B; μB[k]=model.k_B/model.k_rxn*0.5; end

    Ka=k_hi_init-k_lo_init+1
    ic=zeros(Int,2Ka)
    # IC: roughly steady-state for the active voxels
    for (ki,k) in enumerate(k_lo_init:k_hi_init)
        ic[2ki-1]=round(Int,μA[k]); ic[2ki]=round(Int,μB[k])
    end
    sp=StateSpace{CartesianIndex{2Ka},Float64}()
    add_state!(sp,CartesianIndex(Tuple(ic)),1.0)
    avs=ActiveVoxelSet1D(k_lo_init,k_hi_init,K); rates=Float64[]

    n_steps=round(Int,T_END/DT)
    times=Float64[0.0]; μA_hist=[copy(μA)]; μB_hist=[copy(μB)]
    sp_sizes=Int[length(sp)]; ka_hist=Int[n_active(avs)]
    alo_hist=Int[avs.k_lo]; ahi_hist=Int[avs.k_hi]
    # Track reaction rate per voxel: k_rxn * E[nA*nB] ≈ k_rxn * μA * μB (mean-field)
    rxn_hist=[k_rxn_profile(μA,μB,model.k_rxn)]

    for step in 1:n_steps
        t=Float64((step-1)*DT); Ka=n_active(avs); Na=2Ka
        lbc=left_inactive(avs); rbc=right_inactive(avs)
        ga=VoxelGrid(Ka,dx,0)
        sys=build_active_annihilation_1d(model,ga,lbc,rbc,μA,μB,Val(Na);
                source_A,source_B,k_lo_full=avs.k_lo)
        expand!(sp,sys,active_bc(n_max);depth=1)

        if n_vcycle_levels>0 && iseven(Ka) && Ka>=2
            nL=min(n_vcycle_levels,trailing_zeros(Ka))
            hier=build_hierarchy(ga,nL)
            # Reuse the bottleneck V-cycle structure (same state layout)
            # For annihilation we use direct expv — V-cycle smoother assumes
            # reactions are slow relative to intra-pair diffusion, which may not hold.
            A,=build_generator(sp,sys,rates,t)
            sp.probs .= expv(DT,A,sp.probs;m=krylov_m)
        else
            A,=build_generator(sp,sys,rates,t)
            sp.probs .= expv(DT,A,sp.probs;m=krylov_m)
        end

        renormalize!(sp); prune_threshold!(sp,prune_tol)
        update_inactive_means_annihilation_1d!(μA,μB,sp,avs,model,dx,DT;
                source_A,source_B)

        sp,avs=adapt_annihilation_1d(sp,avs,μA,μB,n_max;
                    promote_thresh_A,promote_thresh_B,
                    demote_thresh_A,demote_thresh_B,
                    min_k_act,max_k_act)

        push!(times,step*DT); push!(μA_hist,copy(μA)); push!(μB_hist,copy(μB))
        push!(sp_sizes,length(sp)); push!(ka_hist,n_active(avs))
        push!(alo_hist,avs.k_lo); push!(ahi_hist,avs.k_hi)
        push!(rxn_hist,k_rxn_profile(μA,μB,model.k_rxn))
    end

    times,μA_hist,μB_hist,sp_sizes,ka_hist,alo_hist,ahi_hist,rxn_hist
end

k_rxn_profile(μA,μB,k_rxn) = k_rxn .* μA .* μB

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

const K    = 24          # 1D voxels along x
const DX   = 1.0 / K
const K_Y  = 16          # number of y-voxels (for 2D visualization only)
const DA   = 0.005
const DB   = 0.005
const K_A  = 0.5         # ∅→A rate per source voxel
const K_B  = 0.5
const K_RXN = 0.3        # A+B→∅ rate
const N_MAX = 6
const T_END = 6.0
const DT    = 0.2
const MAX_KA = 4

# Source: A from leftmost 2 voxels, B from rightmost 2 voxels
const SRC_A = Set(1:2)
const SRC_B = Set(K-1:K)

# Active region init: center voxels where front forms
const K_LO_INIT = K÷2 - 1
const K_HI_INIT = K÷2 + 2

model = AnnihilationModel1D(DA, DB, K_A, K_B, K_RXN)

println("="^62)
@printf("A+B→∅ front  K=%d  DA=%.4f  d=%.2f  K_rxn=%.2f\n",K,DA,DA/DX^2,K_RXN)
println("="^62)

# ODE reference
print("ODE ............. "); t_ode,μA_ode,μB_ode=@time run_ode_annihilation(
    K,model,DX,T_END,DT;source_A=SRC_A,source_B=SRC_B)

# Mixed multigrid
print("Mixed+MG ........ "); t_mg,μA_mg,μB_mg,sp_sz,ka,alo,ahi,rxn_hist=@time run_mixed_annihilation(
    K,model,DX,N_MAX,T_END,DT;
    source_A=SRC_A, source_B=SRC_B,
    k_lo_init=K_LO_INIT, k_hi_init=K_HI_INIT,
    promote_thresh_A=0.08, promote_thresh_B=0.08,
    demote_thresh_A=0.02,  demote_thresh_B=0.02,
    min_k_act=2, max_k_act=MAX_KA)

@printf("Ka: mean=%.1f  max=%d  |S|_peak=%d\n",mean(ka),maximum(ka),maximum(sp_sz))
# Compare mean reaction rates at steady state (last quarter of sim)
ss_idx=round(Int,3*T_END/(4*DT)):round(Int,T_END/DT)+1
@printf("SS ΣA_mg=%.3f  ΣA_ode=%.3f  (should match for linear-approx inactive)\n",
        mean(sum.(μA_mg[ss_idx])),mean(sum.(μA_ode[ss_idx])))

# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

outdir=joinpath(@__DIR__,"output"); mkpath(outdir)
n_steps=round(Int,T_END/DT)
snap_ts=[0.0, T_END*0.3, T_END*0.6, T_END]
snap_idxs=[round(Int,t/DT)+1 for t in snap_ts]
sclrs=[:black,:royalblue,:crimson,:orange]

fig=Figure(size=(1600,1100))

# ─── 2D heatmaps of μA and μB (extrude 1D profile into K_Y rows) ─────────────
# Each column of the 2D grid = one x-voxel; y is homogeneous.
# Sources: A at x=1,2 (left); B at x=K-1,K (right)
nT=length(t_ode)
matA=zeros(nT,K); matB=zeros(nT,K)
for ti in 1:nT
    matA[ti,:].=μA_ode[ti]; matB[ti,:].=μB_ode[ti]
end

ax_hA=Axis(fig[1,1:2];title="E[nA(x,t)] — ODE mean-field (species A, from left)",
           xlabel="voxel x",ylabel="time")
ax_hB=Axis(fig[1,3:4];title="E[nB(x,t)] — ODE mean-field (species B, from right)",
           xlabel="voxel x",ylabel="time")
hm_A=heatmap!(ax_hA,1:K,t_ode,matA';colormap=:plasma,  colorrange=(0,maximum(matA)+1e-6))
hm_B=heatmap!(ax_hB,1:K,t_ode,matB';colormap=:viridis, colorrange=(0,maximum(matB)+1e-6))
Colorbar(fig[1,2];colormap=:plasma,  limits=(0,maximum(matA)),label="E[nA]",vertical=true,tellheight=false)
Colorbar(fig[1,4];colormap=:viridis, limits=(0,maximum(matB)),label="E[nB]",vertical=true,tellheight=false)
# Overlay active window
for ax in [ax_hA,ax_hB]
    band!(ax,t_mg,Float64.(alo),Float64.(ahi);color=(:white,0.25))
    lines!(ax,t_mg,Float64.(alo);color=:white,linewidth=1.5)
    lines!(ax,t_mg,Float64.(ahi);color=:white,linewidth=1.5)
end

# ─── Spatial profiles at snapshots ──────────────────────────────────────────
ax_pA=Axis(fig[2,1:2];title="E[nA_k] profiles (ODE)",xlabel="voxel k",ylabel="E[nA]")
ax_pB=Axis(fig[2,3:4];title="E[nB_k] profiles (ODE)",xlabel="voxel k",ylabel="E[nB]")
for (si,sn) in enumerate(snap_idxs)
    t_s=snap_ts[si]
    lines!(ax_pA,1:K,μA_ode[sn];color=sclrs[si],linewidth=2,label="t=$(round(t_s,digits=1))")
    lines!(ax_pB,1:K,μB_ode[sn];color=sclrs[si],linewidth=2)
end
# Mark source voxels
vspan!(ax_pA,1,2;color=(:blue,0.1)); vspan!(ax_pA,K-1,K;color=(:green,0.1))
vspan!(ax_pB,1,2;color=(:blue,0.1)); vspan!(ax_pB,K-1,K;color=(:green,0.1))
axislegend(ax_pA;labelsize=10)

# ─── 2D snapshot at SS (extrude to K_Y rows) ──────────────────────────────
ss_snap=snap_idxs[end]
ax_2dA=Axis(fig[3,1:2];title="2D snapshot E[nA] at t=$(round(snap_ts[end],digits=1)) (K_y=$K_Y rows)",
            xlabel="voxel x",ylabel="voxel y")
ax_2dB=Axis(fig[3,3:4];title="2D snapshot E[nB] at t=$(round(snap_ts[end],digits=1))",
            xlabel="voxel x",ylabel="voxel y")
heatmap!(ax_2dA,1:K,1:K_Y,repeat(μA_ode[ss_snap]',K_Y,1)';colormap=:plasma)
heatmap!(ax_2dB,1:K,1:K_Y,repeat(μB_ode[ss_snap]',K_Y,1)';colormap=:viridis)

# ─── K_act, |S|, reaction rate profile ───────────────────────────────────────
ax_ka=Axis(fig[4,1];title="K_act(t) — active voxels",xlabel="time",ylabel="K_act")
ax_ss=Axis(fig[4,2];title="|S_active|(t)",          xlabel="time",ylabel="|S|")
ax_rxn=Axis(fig[4,3:4];title="Reaction rate k_rxn·μA·μB at SS (ODE)",
            xlabel="voxel k",ylabel="rate")

lines!(ax_ka,t_mg,Float64.(ka);color=:royalblue,linewidth=2)
hlines!(ax_ka,[Float64(MAX_KA)];color=:gray,linestyle=:dash,linewidth=1,label="cap")
axislegend(ax_ka;labelsize=10)
lines!(ax_ss,t_mg,Float64.(sp_sz);color=:royalblue,linewidth=2)
lines!(ax_rxn,1:K,rxn_hist[end];color=:crimson,linewidth=2.5,label="SS rate")
vspan!(ax_rxn,minimum(alo),maximum(ahi);color=(:royalblue,0.15),label="active window")
axislegend(ax_rxn;labelsize=10)

png_path=joinpath(outdir,"mixed_annihilation_1d.png")
save(png_path,fig;px_per_unit=3)
println("\nSaved: $png_path\nDone.")
