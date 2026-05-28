"""
Mixed-Resolution Multigrid — 1D Bottleneck, Slow Diffusion
===========================================================

Model:  A →(k_AB) B →(k_deg) ∅,   k_prod=0 (transient bolus of N_MOL A at voxel 1)
        Slow diffusion D_A=D_B=0.001 so the wavefront is narrow (≈2 voxels wide).

At any time the active FSP covers only the ~2-4 voxels at the wavefront.
Behind it: molecules have degraded → empty → demoted (empty criterion).
Ahead of it: molecules haven't arrived → μ≈0 → demoted (empty criterion).

Run: julia --project examples/mixed_multigrid_1d.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# ODE reference (exact means for linear reactions)
# ─────────────────────────────────────────────────────────────────────────────

function run_ode(K, model, dx, N_MOL, T_END, DT;
                 source_voxels=Set{Int}())
    dA=model.D_A/dx^2; dB=model.D_B/dx^2
    k_prod=model.k_prod; k_AB=model.k_AB; k_deg=model.k_deg
    μA=zeros(K); μB=zeros(K); μA[1]=Float64(N_MOL)

    max_rate=2*max(dA,dB)+k_AB+k_deg
    n_sub=max(1,ceil(Int,max_rate*DT*2)); dt_sub=DT/n_sub

    n_steps=round(Int,T_END/DT)
    times=Float64[0.0]; μA_hist=[copy(μA)]; μB_hist=[copy(μB)]
    ΔA=zeros(K); ΔB=zeros(K)
    for _ in 1:n_steps
        for _ in 1:n_sub
            fill!(ΔA,0.0); fill!(ΔB,0.0)
            for k in 1:K
                ΔA[k]=(k in source_voxels ? k_prod : 0.0)-k_AB*μA[k]
                ΔB[k]=k_AB*μA[k]-k_deg*μB[k]
                k>1&&(ΔA[k]+=dA*(μA[k-1]-μA[k]);ΔB[k]+=dB*(μB[k-1]-μB[k]))
                k<K&&(ΔA[k]+=dA*(μA[k+1]-μA[k]);ΔB[k]+=dB*(μB[k+1]-μB[k]))
            end
            μA.=max.(0.0,μA.+ΔA.*dt_sub); μB.=max.(0.0,μB.+ΔB.*dt_sub)
        end
        push!(times,length(times)*DT); push!(μA_hist,copy(μA)); push!(μB_hist,copy(μB))
    end
    times, μA_hist, μB_hist
end

# ─────────────────────────────────────────────────────────────────────────────
# Mixed-resolution multigrid runner
# ─────────────────────────────────────────────────────────────────────────────

function run_mixed_mg(K, model, dx, n_max, T_END, DT;
                      N_MOL=3, k_lo_init=1, k_hi_init=2,
                      source_voxels=Set{Int}(),
                      n_vcycle_levels=1,
                      promote_thresh=0.05,
                      demote_empty_tol=0.01,
                      min_k_act=2, max_k_act=4,
                      prune_tol=1e-12, krylov_m=30)

    μA=zeros(K); μB=zeros(K); μA[k_lo_init]=Float64(N_MOL)
    Ka=k_hi_init-k_lo_init+1; ic=zeros(Int,2Ka); ic[1]=N_MOL
    sp=StateSpace{CartesianIndex{2Ka},Float64}()
    add_state!(sp,CartesianIndex(Tuple(ic)),1.0)
    avs=ActiveVoxelSet1D(k_lo_init,k_hi_init,K); rates=Float64[]

    n_steps=round(Int,T_END/DT)
    times=Float64[0.0]; μA_hist=[copy(μA)]; μB_hist=[copy(μB)]
    sp_sizes=Int[length(sp)]; ka_hist=Int[n_active(avs)]
    alo_hist=Int[avs.k_lo]; ahi_hist=Int[avs.k_hi]

    for step in 1:n_steps
        t=Float64((step-1)*DT); Ka=n_active(avs); Na=2Ka
        lbc=left_inactive(avs); rbc=right_inactive(avs)
        ga=VoxelGrid(Ka,dx,0)
        sys=build_active_bottleneck_1d(model,ga,lbc,rbc,μA,μB,Val(Na);
                source_voxels_full=source_voxels,k_lo_full=avs.k_lo)
        expand!(sp,sys,active_bc(n_max);depth=1)
        if n_vcycle_levels>0 && iseven(Ka) && Ka>=2
            nL=min(n_vcycle_levels,trailing_zeros(Ka))
            hier=build_hierarchy(ga,nL)
            sp=mixed_vcycle_1d_2s(sp,model,hier,lbc,rbc,μA,μB,rates,t,DT;krylov_m,prune_tol)
        else
            A,=build_generator(sp,sys,rates,t)
            sp.probs .= expv(DT,A,sp.probs;m=krylov_m)
        end
        renormalize!(sp); prune_threshold!(sp,prune_tol)
        update_inactive_means_1d!(μA,μB,sp,avs,model,dx,DT;source_voxels_full=source_voxels)
        sp,avs=adapt_active_1d(sp,avs,μA,μB,n_max;
                    μA_ss=Inf,μB_ss=Inf,
                    promote_thresh,demote_empty_tol,min_k_act,max_k_act)
        push!(times,step*DT); push!(μA_hist,copy(μA)); push!(μB_hist,copy(μB))
        push!(sp_sizes,length(sp)); push!(ka_hist,n_active(avs))
        push!(alo_hist,avs.k_lo); push!(ahi_hist,avs.k_hi)
    end

    times, μA_hist, μB_hist, sp_sizes, ka_hist, alo_hist, ahi_hist
end

# ─────────────────────────────────────────────────────────────────────────────
# Parameters: slow diffusion D=0.001 → d·DT ≈ 0.1 for K=16
# ─────────────────────────────────────────────────────────────────────────────

const DA    = 0.001
const DB    = 0.001
const KPROD = 0.0     # transient bolus
const KAB   = 0.5     # A→B (slow: long-lived wavefront)
const KDEG  = 1.0     # B→∅
const N_MOL = 5
const N_MAX = 5
const T_END = 8.0
const DT    = 0.2
const MAX_KA = 4

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

println("="^62)
println("Slow diffusion 1D bottleneck  (k_prod=0, D=0.001)")
println("="^62)

results = Dict{Int,NamedTuple}()
for K in [16, 32]
    dx    = 1.0/K
    model = BottleneckModel1D(DA,DB,KPROD,KAB,KDEG)
    d_val = DA/dx^2
    @printf("K=%d  dx=%.4f  d=%.3f  d·DT=%.3f\n",K,dx,d_val,d_val*DT)

    t_ode,μA_ode,μB_ode = run_ode(K,model,dx,N_MOL,T_END,DT)

    print("  Mixed+MG ... ")
    t_mg,μA_mg,μB_mg,sp_sz,ka,alo,ahi = @time run_mixed_mg(
        K,model,dx,N_MAX,T_END,DT; N_MOL,
        k_lo_init=1,k_hi_init=2,n_vcycle_levels=1,
        promote_thresh=0.05,demote_empty_tol=0.01,
        min_k_act=2,max_k_act=MAX_KA)

    nA_ode=[sum(v) for v in μA_ode]; nA_mg=[sum(v) for v in μA_mg]
    err=maximum(abs.(nA_mg.-nA_ode))
    @printf("  Ka: mean=%.1f  max=%d  |S|_peak=%d  max|ΔΣμA|=%.4f\n",
            mean(ka),maximum(ka),maximum(sp_sz),err)
    results[K]=(t_ode=t_ode,μA_ode=μA_ode,μB_ode=μB_ode,
                nA_ode=nA_ode,t_mg=t_mg,μA_mg=μA_mg,μB_mg=μB_mg,
                nA_mg=nA_mg,sp_sz=sp_sz,ka=ka,alo=alo,ahi=ahi)
end

# ─────────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────────

outdir=joinpath(@__DIR__,"output"); mkpath(outdir)
n_steps=round(Int,T_END/DT)
snap_ts=[0.0, T_END*0.25, T_END*0.5, T_END]
snap_idxs=[round(Int,t/DT)+1 for t in snap_ts]
sclrs=[:black,:royalblue,:crimson,:orange]

# For K=32: 2D heatmap (extrude 1D profile into rows)
r=results[32]; K32=32

fig=Figure(size=(1600,950))

# ─── Row 1: Spatial heatmaps (A and B profiles over time, extruded to 2D) ───
ax_A=Axis(fig[1,1:2];title="K=32: E[nA_k] heatmap (space × time)",
          xlabel="voxel k",ylabel="time")
ax_B=Axis(fig[1,3:4];title="K=32: E[nB_k] heatmap (space × time)",
          xlabel="voxel k",ylabel="time")

# Build matrices: rows=time steps, cols=voxels
nT=length(r.t_ode)
matA=zeros(nT,K32); matB=zeros(nT,K32)
for ti in 1:nT
    matA[ti,:].=r.μA_ode[ti]; matB[ti,:].=r.μB_ode[ti]
end
heatmap!(ax_A,1:K32,r.t_ode,matA';colormap=:plasma)
heatmap!(ax_B,1:K32,r.t_ode,matB';colormap=:viridis)
# Overlay active window
band!(ax_A,r.t_mg,Float64.(r.alo),Float64.(r.ahi);color=(:white,0.3))
band!(ax_B,r.t_mg,Float64.(r.alo),Float64.(r.ahi);color=(:white,0.3))
lines!(ax_A,r.t_mg,Float64.(r.alo);color=:white,linewidth=1.5)
lines!(ax_A,r.t_mg,Float64.(r.ahi);color=:white,linewidth=1.5)
lines!(ax_B,r.t_mg,Float64.(r.alo);color=:white,linewidth=1.5)
lines!(ax_B,r.t_mg,Float64.(r.ahi);color=:white,linewidth=1.5)

# ─── Row 2: Spatial profiles at snapshots ────────────────────────────────────
ax_profA=Axis(fig[2,1:2];title="E[nA_k] profiles at 4 times (K=32)",
              xlabel="voxel k",ylabel="E[nA_k]")
ax_profB=Axis(fig[2,3:4];title="E[nB_k] profiles at 4 times (K=32)",
              xlabel="voxel k",ylabel="E[nB_k]")
for (si,sn) in enumerate(snap_idxs)
    t_s=snap_ts[si]
    lines!(ax_profA,1:K32,r.μA_ode[sn];color=sclrs[si],linewidth=2,
           label="t=$(round(t_s,digits=1))")
    lines!(ax_profB,1:K32,r.μB_ode[sn];color=sclrs[si],linewidth=2)
end
axislegend(ax_profA;labelsize=10)

# ─── Row 3: |S|, K_act, Σ E[nA] for K=16 and K=32 ───────────────────────────
ax_s=Axis(fig[3,1];title="|S_active|(t)",xlabel="time",ylabel="|S|")
ax_k=Axis(fig[3,2];title="K_act(t)",    xlabel="time",ylabel="K_act")
ax_n=Axis(fig[3,3];title="Σ E[nA]: ODE vs mixed",xlabel="time",ylabel="Σ E[nA_k]")
ax_e=Axis(fig[3,4];title="|ΔΣ E[nA]|",xlabel="time",ylabel="abs error")

clrs=Dict(16=>:crimson,32=>:royalblue)
for K in [16,32]
    r2=results[K]
    lines!(ax_s,r2.t_mg,Float64.(r2.sp_sz);color=clrs[K],linewidth=2,label="K=$K")
    lines!(ax_k,r2.t_mg,Float64.(r2.ka);   color=clrs[K],linewidth=2,label="K=$K")
    lines!(ax_n,r2.t_ode,r2.nA_ode;color=:gray,linewidth=3)
    lines!(ax_n,r2.t_mg, r2.nA_mg; color=clrs[K],linewidth=2,linestyle=:dash,label="K=$K mixed")
    err=[abs(a-b) for (a,b) in zip(r2.nA_mg,r2.nA_ode)]
    lines!(ax_e,r2.t_mg,err.+1e-16;color=clrs[K],linewidth=2,label="K=$K")
end
hlines!(ax_k,[Float64(MAX_KA)];color=:gray,linestyle=:dash,linewidth=1,label="cap")
for ax in [ax_s,ax_k,ax_n,ax_e]; axislegend(ax;labelsize=10); end

png_path=joinpath(outdir,"mixed_multigrid_1d.png")
save(png_path,fig;px_per_unit=3)
println("\nSaved: $png_path\nDone.")
