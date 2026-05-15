"""
2D Schlögl RDME — Spontaneous Domain Formation from Uniform Unstable IC
========================================================================

IC: ALL voxels at n_uns≈72 (the unstable fixed point).
Stochastic fluctuations spontaneously break symmetry: some voxels commit
to the HIGH stable state (~162), others to the LOW stable state (~31).
Diffusion correlates neighbours — adjacent voxels tend to co-commit,
forming spatial DOMAINS separated by a stochastic front.

Three rows:
  1. Mean-field ODE     — all voxels stay at n_uns forever (symmetric IC
                          → deterministic system cannot break symmetry)
  2. Per-voxel FSP mean — means drift but no spatial coordination;
                          all voxels evolve identically (mean-field coupling)
  3. Per-voxel FSP Std  — HIGH everywhere (all voxels bimodal), no spatial pattern

Key result (bottom row):
  P(n) at probe voxels: starts as δ_{n_uns}, develops bimodal peaks at
  n_low and n_high over time — stochastic commitment to stable states.

Run: julia --project examples/fig_schlogl_2d_front.jl
     K=6 D=0.1 T_END=30 julia --project examples/fig_schlogl_2d_front.jl
"""

using DiscStochSim
using ExponentialUtilities
using SparseArrays
using CairoMakie, Printf

# ── Parameters ────────────────────────────────────────────────────────────────
e(k,d) = parse(typeof(d), get(ENV, k, string(d)))

const K      = e("K",      8)
const D      = e("D",      0.3)     # fast enough for front to dissolve in MF
const h      = 1.0
const t_end  = e("T_END",  8.0)
const dt_mf  = e("DT_MF",  0.005)
const dt_fsp = e("DT_FSP", 0.2)
const n_max  = e("N_MAX",  220)
const krylov_m = 30

model = SchloglModel1D(D)
n_low, n_uns, n_high = schlogl_fixed_points(model)
d = D / h^2
center = (K÷2+1, K÷2+1)

@printf("Schlögl 2D: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
@printf("Grid %d×%d  D=%.3f  d=%.3f  IC=left-half n_high / right-half n_low  t_end=%.1f\n",
        K, K, D, d, t_end)

snap_ts = unique([0.0; filter(<(t_end), [2.0, 4.0, 6.0]); t_end])

# Probe: deep in high region, at interface, deep in low region
j_iface = K÷2      # last column of high half
probe_voxels = [(K÷2+1, 1), (K÷2+1, j_iface), (K÷2+1, K)]
probe_labels  = ["deep high (col 1)", "interface (col $j_iface)", "deep low (col $K)"]

# ── Neighbours ────────────────────────────────────────────────────────────────
function neighbors_2d(i, j, K)
    nb = Tuple{Int,Int}[]
    i > 1 && push!(nb, (i-1,j)); i < K && push!(nb, (i+1,j))
    j > 1 && push!(nb, (i,j-1)); j < K && push!(nb, (i,j+1))
    nb
end

# ── Schlögl generator ─────────────────────────────────────────────────────────
function build_schlogl_gen(model, d, n_nb, μ_in, n_max)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    N = n_max+1; Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(i,j,v) = (push!(Is,i); push!(Js,j); push!(Vs,v))
    for n in 0:n_max
        j=n+1; diag=0.0
        r_up = c1*n*(n-1)/2 + c3
        n < n_max && r_up > 0 && (add!(j+1,j,r_up); diag -= r_up)
        if n > 0
            r_dn = c2*n*(n-1)*(n-2)/6 + c4*n
            r_dn > 0 && (add!(j-1,j,r_dn); diag -= r_dn)
        end
        n < n_max && μ_in > 0 && (add!(j+1,j,μ_in); diag -= μ_in)
        n > 0 && (r_out = d*n_nb*n; add!(j-1,j,r_out); diag -= r_out)
        add!(j,j,diag)
    end
    sparse(Is,Js,Vs,N,N)
end

pvoxel_mean(p) = sum((i-1)*p[i] for i in eachindex(p))
pvoxel_std(p)  = sqrt(max(0.0, sum((i-1)^2*p[i] for i in eachindex(p)) - pvoxel_mean(p)^2))

# ── IC: left half n_high, right half n_low ───────────────────────────────────
μ_mf = zeros(K, K)
for i in 1:K, j in 1:K
    μ_mf[i,j] = j <= K÷2 ? Float64(n_high) : Float64(n_low)
end

pvoxel = Dict{Tuple{Int,Int}, Vector{Float64}}()
for i in 1:K, j in 1:K
    p = zeros(n_max+1)
    p[j <= K÷2 ? n_high+1 : n_low+1] = 1.0
    pvoxel[(i,j)] = p
end

# ── 1. Mean-field ODE ─────────────────────────────────────────────────────────
function schlogl_mf_rhs!(dμ, μ, K, d, model)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    for i in 1:K, j in 1:K
        n=μ[i,j]
        rxn = c1*n*(n-1)/2 - c2*n*(n-1)*(n-2)/6 + c3 - c4*n
        lap = 0.0
        i>1 && (lap += d*(μ[i-1,j]-n)); i<K && (lap += d*(μ[i+1,j]-n))
        j>1 && (lap += d*(μ[i,j-1]-n)); j<K && (lap += d*(μ[i,j+1]-n))
        dμ[i,j] = rxn + lap
    end
end

dμ = zeros(K,K)
snap_steps_mf = Set(round(Int, ts/dt_mf) for ts in snap_ts)
snap_μ_mf = Dict{Int, Matrix{Float64}}(0 => copy(μ_mf))
for step in 1:round(Int, t_end/dt_mf)
    schlogl_mf_rhs!(dμ, μ_mf, K, d, model)
    μ_mf .+= dt_mf .* dμ
    step in snap_steps_mf && (snap_μ_mf[step] = copy(μ_mf))
end
println("Mean-field done.")

# ── 2. Per-voxel FSP ──────────────────────────────────────────────────────────
function pv_grids(pvoxel, K)
    mg=zeros(K,K); sg=zeros(K,K)
    for i in 1:K, j in 1:K
        p=pvoxel[(i,j)]; mg[i,j]=pvoxel_mean(p); sg[i,j]=pvoxel_std(p)
    end
    mg, sg
end

snap_steps_fsp = Set(round(Int, ts/dt_fsp) for ts in snap_ts)
mg0,sg0 = pv_grids(pvoxel,K)
snap_pv_mean = Dict{Int,Matrix{Float64}}(0=>mg0)
snap_pv_std  = Dict{Int,Matrix{Float64}}(0=>sg0)
probe_snaps  = Dict{Int,Vector{Vector{Float64}}}(
    0 => [copy(pvoxel[pv]) for pv in probe_voxels])

println("Running per-voxel FSP on $(K)×$(K) grid...")
@time for step in 1:round(Int, t_end/dt_fsp)
    global pvoxel
    μs = Dict((i,j)=>pvoxel_mean(pvoxel[(i,j)]) for i in 1:K, j in 1:K)
    pvoxel = Dict((i,j) => begin
        nbs = neighbors_2d(i,j,K); n_nb=length(nbs)
        μ_in = d*sum(μs[nb] for nb in nbs)
        Q = build_schlogl_gen(model, d, n_nb, μ_in, n_max)
        p = expv(Float64(dt_fsp), Q, pvoxel[(i,j)]; m=krylov_m)
        p .= max.(0.0,p); p ./= sum(p); p
    end for i in 1:K, j in 1:K)

    if step in snap_steps_fsp
        mg,sg = pv_grids(pvoxel,K)
        snap_pv_mean[step]=mg; snap_pv_std[step]=sg
        probe_snaps[step] = [copy(pvoxel[pv]) for pv in probe_voxels]
    end
    step % 5 == 0 && @printf("  t=%.1f  E[hi]=%.1f  E[iface]=%.1f  E[lo]=%.1f\n",
        step*dt_fsp, pvoxel_mean(pvoxel[probe_voxels[1]]),
        pvoxel_mean(pvoxel[probe_voxels[2]]), pvoxel_mean(pvoxel[probe_voxels[3]]))
end
println("Done.")

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

n_snaps = length(snap_ts)
fig = Figure(size=(220*n_snaps+200, 1000))

Label(fig[0,1:n_snaps+1],
    "2D Schlögl Bistable Front  |  $(K)×$(K)  D=$D  " *
    "n_low=$n_low / n_uns=$n_uns / n_high=$n_high  |  IC: left half=n_high, right half=n_low";
    fontsize=11, tellwidth=false)

c_lo=Float64(n_low); c_hi=Float64(n_high)
sorted_snap_keys = sort(collect(keys(snap_pv_mean)))
sorted_mf_keys   = sort(collect(keys(snap_μ_mf)))

# ── Row 1: Mean-field E[n] ───────────────────────────────────────────────────
for (col,(step,ts)) in enumerate(zip(sorted_mf_keys, snap_ts))
    ax = Axis(fig[1,col]; aspect=DataAspect(), title="t=$(round(ts,digits=1))",
              titlesize=10, xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
    hm = heatmap!(ax,1:K,1:K,snap_μ_mf[step]';colormap=:inferno,colorrange=(0,c_hi))
    col==1 && text!(ax,0.02,0.98;text="Mean-field",space=:relative,
                    align=(:left,:top),color=:white,fontsize=9,
                    strokecolor=:black,strokewidth=2)
    col==n_snaps && Colorbar(fig[1,col+1],hm;width=12,label="E[n]",
        ticks=([0,Float64(n_uns),c_hi],["0","n_uns","n_hi"]))
end

# ── Row 2: FSP E[n] ──────────────────────────────────────────────────────────
for (col,(step,ts)) in enumerate(zip(sorted_snap_keys, snap_ts))
    ax = Axis(fig[2,col]; aspect=DataAspect(),
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
    hm = heatmap!(ax,1:K,1:K,snap_pv_mean[step]';colormap=:inferno,colorrange=(0,c_hi))
    for (pv,clr) in zip(probe_voxels,[:white,:yellow,:cyan])
        scatter!(ax,[Float32(pv[2])],[Float32(pv[1])];
                 color=clr,marker=:cross,markersize=12,strokewidth=1.5)
    end
    col==1 && text!(ax,0.02,0.98;text="FSP E[n]",space=:relative,
                    align=(:left,:top),color=:white,fontsize=9,
                    strokecolor=:black,strokewidth=2)
    col==n_snaps && Colorbar(fig[2,col+1],hm;width=12,label="E[n]",
        ticks=([0,Float64(n_uns),c_hi],["0","n_uns","n_hi"]))
end

# ── Row 3: FSP Std[n] ────────────────────────────────────────────────────────
std_max = sqrt(c_hi)*1.5
for (col,(step,ts)) in enumerate(zip(sorted_snap_keys, snap_ts))
    ax = Axis(fig[3,col]; aspect=DataAspect(),
              xticksvisible=false, yticksvisible=false,
              xticklabelsvisible=false, yticklabelsvisible=false)
    hm = heatmap!(ax,1:K,1:K,snap_pv_std[step]';colormap=:magma,colorrange=(0,std_max))
    col==1 && text!(ax,0.02,0.98;text="FSP Std[n]\n(bistable ring)",space=:relative,
                    align=(:left,:top),color=:white,fontsize=9,
                    strokecolor=:black,strokewidth=2)
    col==n_snaps && Colorbar(fig[3,col+1],hm;width=12,label="Std[n]",
        ticks=([0,sqrt(c_hi),std_max],["0","√n_hi","$(round(std_max,digits=1))"]))
end

# ── Row 4: P(n) at probe voxels over time ────────────────────────────────────
n_show    = min(n_max, 210)
snap_clrs = Makie.resample_cmap(:plasma, n_snaps+1)[1:end-1]

for (pi, (pv,lbl)) in enumerate(zip(probe_voxels, probe_labels))
    ax = Axis(fig[4,pi]; xlabel="n", ylabel="P(n)",
              title="P(n,t) at $(lbl) voxel $(pv)", titlesize=10)
    for (ci, step) in enumerate(sorted_snap_keys)
        p  = probe_snaps[step][pi]
        ts = snap_ts[ci]
        lines!(ax, 0:n_show, p[1:n_show+1]; color=snap_clrs[ci], linewidth=1.8,
               label="t=$(round(ts,digits=1))")
    end
    vlines!(ax,[Float64(n_low),Float64(n_high)];color=(:black,0.2),linestyle=:dash,linewidth=1)
    vlines!(ax,[Float64(n_uns)];color=(:gray50,0.4),linestyle=:dot,linewidth=1)
    axislegend(ax;position=:rt,labelsize=8,framevisible=false)
    xlims!(ax,0,n_show)
end

# Mean-field vs FSP mean over time for each probe voxel
ax_mu = Axis(fig[4,4:n_snaps]; xlabel="time", ylabel="E[n]",
             title="E[n] over time at probe voxels: FSP (solid) vs mean-field (dashed)",
             titlesize=10)
probe_clrs = [:steelblue4, :crimson, :darkorange]
t_mf_hist = collect(0.0:dt_mf:t_end)
for (pi, (pv,lbl,clr)) in enumerate(zip(probe_voxels, probe_labels, probe_clrs))
    # MF: same for all voxels (symmetric IC)
    mf_mu = [snap_μ_mf[sort(collect(keys(snap_μ_mf)))[1]][pv...]  # constant
             for _ in t_mf_hist]
    # Actually build MF time series from snapshots
    mf_steps = sort(collect(keys(snap_μ_mf)))
    mf_ts = [s * dt_mf for s in mf_steps]
    mf_vs = [snap_μ_mf[s][pv...] for s in mf_steps]
    lines!(ax_mu, mf_ts, mf_vs; color=(clr,0.5), linewidth=1.5, linestyle=:dash)
    # FSP
    fsp_ts = [s * dt_fsp for s in sorted_snap_keys]
    fsp_vs = [snap_pv_mean[s][pv...] for s in sorted_snap_keys]
    lines!(ax_mu, fsp_ts, fsp_vs; color=clr, linewidth=2.5, label=lbl)
end
hlines!(ax_mu,[Float64(n_high)];color=(:black,0.2),linestyle=:dash,linewidth=1)
hlines!(ax_mu,[Float64(n_uns)]; color=(:gray50,0.4),linestyle=:dot,linewidth=1)
hlines!(ax_mu,[Float64(n_low)]; color=(:black,0.2),linestyle=:dash,linewidth=1)
text!(ax_mu,1.0,Float64(n_high)+3;text="n_high=$n_high",fontsize=9)
text!(ax_mu,1.0,Float64(n_uns)+3; text="n_uns=$n_uns",  fontsize=9,color=:gray50)
text!(ax_mu,1.0,Float64(n_low)+3; text="n_low=$n_low",  fontsize=9)
axislegend(ax_mu;position=:rc,labelsize=9,framevisible=false)

rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,5); rowgap!(fig.layout,3,10)

outdir = joinpath(@__DIR__,"..","paper","figures")
mkpath(outdir)
out = joinpath(outdir,"fig_schlogl_2d_front.pdf")
save(out,fig)
println("Saved → $out")
