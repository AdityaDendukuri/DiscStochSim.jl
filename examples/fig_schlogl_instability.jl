"""
Schlögl RDME — Spontaneous Symmetry Breaking from Unstable IC
=============================================================

All K voxels start at the unstable fixed point n_uns≈72.
Watch the system spontaneously bifurcate toward the two stable states.

Three methods compared:

  1. Mean-field ODE
     - All voxels start at n_uns, symmetric IC → stays at n_uns forever.
     - The unstable fixed point is also a mean-field fixed point under
       symmetric coupling → mean-field CANNOT break symmetry.

  2. Per-voxel FSP (mean-field coupling)
     - Each voxel's distribution P(n,t) evolves independently under Schlögl.
     - Correctly develops a bimodal P(n): each voxel bifurcates stochastically.
     - FAILURE: no spatial coordination — coupling uses bimodal mean ≈ n_uns,
       so neighbours can't tell if a voxel went high or low.

  3. Multigrid FSP (V-cycle)
     - Joint P(n₁,…,nₖ,t) captures WHICH state each voxel chose.
     - Spatial correlations: adjacent voxels tend to co-flip due to diffusion.
     - Gives probability of every spatial configuration (HHHH, HHLL, LLLL, …).

ENV parameters:
  K          number of voxels (default 4)
  D          diffusion coefficient (default 0.1)
  T_END      simulation end time (default 15.0)
  DT_FSP     FSP time step (default 0.5)
  N_MAX      CME truncation (default 220)
  PRUNE_TOL  multigrid pruning threshold (default 5e-4)
  TAU_PRE    pre-smoothing time fraction (default 0.3)

Run: julia --project examples/fig_schlogl_instability.jl
"""

using DiscStochSim
using ExponentialUtilities
using SparseArrays
using CairoMakie, Printf

# ── Parameters ────────────────────────────────────────────────────────────────
e(name, default) = parse(typeof(default), get(ENV, name, string(default)))

const K         = e("K",         4)
const D         = e("D",         0.1)
const h         = 1.0
const t_end     = e("T_END",     200.0)
const dt_mf     = e("DT_MF",     0.05)
const dt_fsp    = e("DT_FSP",    5.0)
const n_max     = e("N_MAX",     220)
const prune_tol = e("PRUNE_TOL", 1e-5)
const τ_pre     = e("TAU_PRE",   0.3)
const krylov_m  = 30

model     = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)
n_low, n_uns, n_high = schlogl_fixed_points(model)
d = D / h^2

@printf("Schlögl: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
@printf("D=%.3f  d=%.3f  K=%d  t_end=%.1f  dt_fsp=%.2f  prune=%.0e\n",
        D, d, K, t_end, dt_fsp, prune_tol)
println("IC: all voxels at n_uns=$n_uns (unstable fixed point)")

snap_ts    = unique([0.0; filter(<(t_end), [20.0, 60.0, 120.0]); t_end])
snap_steps = Set(round(Int, ts/dt_fsp) for ts in snap_ts)

# ── 1. Mean-field ODE ─────────────────────────────────────────────────────────
function schlogl_mf_rhs!(dμ, μ, model, K, d)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    for k in 1:K
        n   = μ[k]
        rxn = c1*n*(n-1)/2 - c2*n*(n-1)*(n-2)/6 + c3 - c4*n
        lap = (k>1 ? d*(μ[k-1]-n) : 0.0) + (k<K ? d*(μ[k+1]-n) : 0.0)
        dμ[k] = rxn + lap
    end
end

μ_mf      = fill(Float64(n_uns), K)
dμ        = zeros(K)
t_mf_hist = Float64[0.0]
μ_mf_hist = [copy(μ_mf)]
for _ in 1:round(Int, t_end/dt_mf)
    schlogl_mf_rhs!(dμ, μ_mf, model, K, d)
    μ_mf .+= dt_mf .* dμ
    push!(t_mf_hist, t_mf_hist[end] + dt_mf)
    push!(μ_mf_hist, copy(μ_mf))
end
println("Mean-field done.  Final μ = $(round.(μ_mf_hist[end], digits=1))")

# ── 2. Per-voxel FSP ──────────────────────────────────────────────────────────
function build_schlogl_voxel_gen(model, d, n_nb, μ_in, n_max)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    N = n_max + 1
    Is=Int[]; Js=Int[]; Vs=Float64[]
    add!(i,j,v) = (push!(Is,i); push!(Js,j); push!(Vs,v))
    for n in 0:n_max
        j = n+1; diag = 0.0
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

pvoxel_means(p) = sum((i-1)*p[i] for i in eachindex(p))

# IC: delta at n_uns for all voxels
pvoxel_dists = [begin p=zeros(n_max+1); p[n_uns+1]=1.0; p end for _ in 1:K]

t_pv_hist  = Float64[0.0]
μ_pv_hist  = [[pvoxel_means(p) for p in pvoxel_dists]]
snap_pvdists = Dict{Int, Vector{Vector{Float64}}}(0 => deepcopy(pvoxel_dists))

println("Running per-voxel FSP...")
for step in 1:round(Int, t_end/dt_fsp)
    global pvoxel_dists
    μs = [pvoxel_means(p) for p in pvoxel_dists]
    pvoxel_dists = [begin
        n_nb = (k>1)+(k<K)
        μ_in = d*(k>1 ? μs[k-1] : 0.0) + d*(k<K ? μs[k+1] : 0.0)
        Q    = build_schlogl_voxel_gen(model, d, n_nb, μ_in, n_max)
        p    = expv(Float64(dt_fsp), Q, pvoxel_dists[k]; m=krylov_m)
        p .= max.(0.0, p); p ./= sum(p); p
    end for k in 1:K]
    push!(t_pv_hist, step*dt_fsp)
    push!(μ_pv_hist, [pvoxel_means(p) for p in pvoxel_dists])
    step in snap_steps && (snap_pvdists[step] = deepcopy(pvoxel_dists))
end
println("Per-voxel FSP done.  Final μ = $(round.(μ_pv_hist[end], digits=1))")

# ── 3. Multigrid FSP ──────────────────────────────────────────────────────────
u0  = CartesianIndex(ntuple(_ -> n_uns, K))
alg = RDMEMultigridFSP(
    dt             = dt_fsp,
    n_levels       = 1,
    n_max          = n_max,
    save_every     = 1,
    max_states     = 400_000,
    expand_depth   = 0,
    coarse_r_depth = 1,
    coarse_d_depth = 1,
    prune_tol      = prune_tol,
    weight_tol     = prune_tol,
    τ_pre          = τ_pre,
)
println("Running multigrid FSP...")
@time sol_mg = solve_rdme_multigrid(model, fine_grid, u0, (0.0, t_end), Float64[], alg)
println("Multigrid done. $(length(sol_mg.t)) snapshots")
println("  |S|: ", sol_mg.state_space_sizes)

μ_mg_hist = mean_voxel_counts(sol_mg)
t_mg_hist = sol_mg.t

# ── Configuration probabilities from multigrid at final time ──────────────────
"""
Compute P(exactly m voxels in high state) at a given snapshot.
High state: n > n_uns.  Low state: n ≤ n_uns.
"""
function config_probs(sol, ti, n_uns, K)
    states, probs = sol.snapshots[ti]
    p_m = zeros(K+1)   # p_m[m+1] = P(exactly m voxels in high state)
    for (s, pr) in zip(states, probs)
        m = count(n -> n > n_uns, Tuple(s))
        p_m[m+1] += pr
    end
    p_m
end

mg_snap_idxs = [argmin(abs.(t_mg_hist .- ts)) for ts in snap_ts]
config_snap  = [config_probs(sol_mg, ti, n_uns, K) for ti in mg_snap_idxs]

# ── 4. SSA ground truth ───────────────────────────────────────────────────────
"""
Gillespie SSA for K-voxel 1D Schlögl RDME.
Returns time series of molecule counts sampled at save_times.
"""
function run_ssa(model, K, d, n0::Vector{Int}, t_end, save_times; seed=0)
    rng    = seed > 0 ? Random.MersenneTwister(seed) : Random.default_rng()
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    n      = copy(n0)
    t      = 0.0
    si     = 1
    snaps  = [copy(n)]
    while t < t_end
        # Propensities
        props = Float64[]
        for k in 1:K
            nk = n[k]
            push!(props, c3)                              # ∅→X
            push!(props, c4*nk)                           # X→∅
            push!(props, c1*nk*(nk-1)/2)                  # 2X→3X
            push!(props, c2*nk*(nk-1)*(nk-2)/6)           # 3X→2X
            k > 1 && push!(props, d*nk)                   # X_k→X_{k-1}
            k < K && push!(props, d*nk)                   # X_k→X_{k+1}
        end
        a0 = sum(props)
        a0 == 0.0 && break
        t += -log(rand(rng)) / a0
        r = rand(rng) * a0
        # Select reaction
        cum = 0.0; fired = 0
        for i in eachindex(props)
            cum += props[i]; fired = i
            cum >= r && break
        end
        # Apply reaction — decode index
        idx = 0
        for k in 1:K
            base = idx
            idx += 4 + (k>1 ? 1 : 0) + (k<K ? 1 : 0)
            fired <= idx || continue
            local_fired = fired - base
            if local_fired == 1; n[k] += 1               # ∅→X
            elseif local_fired == 2; n[k] = max(0,n[k]-1) # X→∅
            elseif local_fired == 3; n[k] += 1            # 2X→3X
            elseif local_fired == 4; n[k] = max(0,n[k]-1) # 3X→2X
            elseif local_fired == 5
                nb = k > 1 ? k-1 : k+1
                n[k] = max(0,n[k]-1); n[nb] += 1          # diffusion left
            elseif local_fired == 6
                n[k] = max(0,n[k]-1); n[k+1] += 1         # diffusion right
            end
            break
        end
        while si <= length(save_times) && t >= save_times[si]
            push!(snaps, copy(n)); si += 1
        end
    end
    while length(snaps) <= length(save_times); push!(snaps, copy(n)); end
    snaps   # snaps[1]=t=0, snaps[2..end] = at save_times
end

using Random
const N_SSA    = e("N_SSA", 1000)
const ssa_save = collect(range(0.0, t_end; step=dt_fsp))

ssa_times = ssa_save[2:end]
println("Running $N_SSA SSA trajectories...")
n0_ssa = fill(n_uns, K)
@time ssa_trajs = [run_ssa(model, K, d, copy(n0_ssa), t_end, ssa_times; seed=i)
                   for i in 1:N_SSA]
println("SSA done.")

# Config probabilities from SSA: P(m voxels > n_uns) at each save time
function ssa_config_probs(trajs, ti, n_uns, K, N_traj)
    p_m = zeros(K+1)
    for traj in trajs
        snap = traj[ti]   # molecule counts at time index ti
        m = count(n -> n > n_uns, snap)
        p_m[m+1] += 1.0
    end
    p_m ./= N_traj
end

# SSA marginal histogram at voxel 1
function ssa_marginal(trajs, ti, n_max)
    p = zeros(n_max+1)
    for traj in trajs
        n = traj[ti][1]
        n <= n_max && (p[n+1] += 1.0)
    end
    p ./= length(trajs)
end

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(1350, 850))
Label(fig[0, 1:4],
    "Schlögl RDME — Symmetry Breaking from Unstable IC  |  K=$K voxels  D=$D  " *
    "n_uns=$n_uns  →  n_low=$n_low / n_high=$n_high  |  t_end=$t_end";
    fontsize=11, tellwidth=false)

vox_clrs = [:crimson, :darkorange, :steelblue4, :forestgreen][1:K]

# ── Row 1: Voxel means ────────────────────────────────────────────────────────
for (col, (title, t_hist, μ_hist)) in enumerate([
    ("1. Mean-field ODE\n(symmetric IC → stays at n_uns forever)", t_mf_hist, μ_mf_hist),
    ("2. Per-voxel FSP + mean-field coupling\n(bimodal P(n) but no spatial coordination)", t_pv_hist, μ_pv_hist),
    ("3. Multigrid FSP (V-cycle)\n(joint distribution: spatial patterns emerge)", t_mg_hist, [μ_mg_hist[i,:] for i in 1:length(t_mg_hist)]),
])
    ax = Axis(fig[1, col]; xlabel="time", ylabel="E[n]", title=title, titlesize=10)
    for k in 1:K
        μ_k = [μ_hist[i][k] for i in eachindex(t_hist)]
        lines!(ax, t_hist, μ_k; color=vox_clrs[k], linewidth=2,
               label="voxel $k")
    end
    hlines!(ax, [Float64(n_high)]; color=(:black,0.25), linestyle=:dash, linewidth=1)
    hlines!(ax, [Float64(n_uns)];  color=(:gray50,0.5), linestyle=:dot,  linewidth=1)
    hlines!(ax, [Float64(n_low)];  color=(:black,0.25), linestyle=:dash, linewidth=1)
    text!(ax, t_end*0.02, Float64(n_high)+3; text="n_high=$n_high", fontsize=8)
    text!(ax, t_end*0.02, Float64(n_uns)+3;  text="n_uns=$n_uns",   fontsize=8, color=:gray50)
    text!(ax, t_end*0.02, Float64(n_low)+3;  text="n_low=$n_low",   fontsize=8)
    col == 1 && axislegend(ax; position=:rc, labelsize=8, framevisible=false)
    ylims!(ax, -5, Float64(n_high)+20)
end

# ── Row 1 col 4: |S| ─────────────────────────────────────────────────────────
ax_ss = Axis(fig[1, 4]; xlabel="time", ylabel="|S|",
             title="Multigrid FSP state space\n(D=$D, prune=$prune_tol)", titlesize=10)
lines!(ax_ss, t_mg_hist, Float64.(sol_mg.state_space_sizes); color=:steelblue4, linewidth=2.5)

# ── Row 2: P(n) evolution at voxel 1 ─────────────────────────────────────────
n_show   = min(n_max, 210)
snap_clrs = Makie.resample_cmap(:plasma, length(snap_ts)+1)[1:end-1]
sorted_snap_steps = sort(collect(keys(snap_pvdists)))

# Per-voxel FSP P(n) at voxel 1 — skip t=0, auto-scale to see tails
ax_pv = Axis(fig[2, 1:2]; xlabel="n", ylabel="P(n)",
             title="Per-voxel FSP: P(n) at voxel 1  —  unimodal drift, no bistable splitting",
             titlesize=10)
non_zero_snaps = [ss for ss in sorted_snap_steps if ss > 0]
pv_clrs_nz = Makie.resample_cmap(:plasma, length(non_zero_snaps)+1)[1:end-1]
for (ci, ss) in enumerate(non_zero_snaps)
    t_s = t_pv_hist[ss+1]
    p   = snap_pvdists[ss][1]
    lines!(ax_pv, 0:n_show, p[1:n_show+1]; color=pv_clrs_nz[ci], linewidth=1.8,
           label="t=$(round(t_s,digits=1))")
end
vlines!(ax_pv, [Float64(n_low), Float64(n_high)]; color=(:black,0.2), linestyle=:dash, linewidth=1)
vlines!(ax_pv, [Float64(n_uns)];                  color=(:gray50,0.4), linestyle=:dot,  linewidth=1)
axislegend(ax_pv; position=:rt, labelsize=9, framevisible=false)
xlims!(ax_pv, 0, n_show)

# Multigrid FSP P(n) marginal at voxel 1
function mg_voxel_marginal(sol, ti, vox, n_max)
    states, probs = sol.snapshots[ti]
    p = zeros(n_max+1)
    for (s, pr) in zip(states, probs)
        n = Tuple(s)[vox]; n <= n_max && (p[n+1] += pr)
    end
    s = sum(p); s > 0 && (p ./= s); p
end

mg_non_zero_idxs = mg_snap_idxs[2:end]  # skip t=0
mg_clrs_nz = Makie.resample_cmap(:plasma, length(mg_non_zero_idxs)+1)[1:end-1]

ax_mg = Axis(fig[2, 3:4]; xlabel="n", ylabel="P(n)",
             title="Multigrid FSP: P(n) at voxel 1  —  bimodal splitting at n_low and n_high emerges",
             titlesize=10)
for (ci, (ti, ts)) in enumerate(zip(mg_non_zero_idxs, snap_ts[2:end]))
    p = mg_voxel_marginal(sol_mg, ti, 1, n_max)
    lines!(ax_mg, 0:n_show, p[1:n_show+1]; color=mg_clrs_nz[ci], linewidth=2,
           label="FSP t=$(round(t_mg_hist[ti],digits=0))")
    # SSA histogram at same snap time
    ssa_ti = argmin(abs.(ssa_times .- ts))
    p_ssa  = ssa_marginal(ssa_trajs, ssa_ti, n_show)
    stairs!(ax_mg, 0:n_show, p_ssa; color=(mg_clrs_nz[ci], 0.5), linewidth=1,
            linestyle=:dash)
end
text!(ax_mg, n_show*0.5, maximum(mg_voxel_marginal(sol_mg, mg_non_zero_idxs[end], 1, n_max))*0.8;
      text="── FSP   ╌╌ SSA", fontsize=9, color=:gray40)
vlines!(ax_mg, [Float64(n_low), Float64(n_high)]; color=(:black,0.2), linestyle=:dash, linewidth=1)
vlines!(ax_mg, [Float64(n_uns)];                  color=(:gray50,0.4), linestyle=:dot,  linewidth=1)
axislegend(ax_mg; position=:rt, labelsize=9, framevisible=false)
xlims!(ax_mg, 0, n_show)

# ── Row 3: Configuration probabilities over time (multigrid only) ─────────────
ax_cfg = Axis(fig[3, 1:4];
    xlabel="time", ylabel="probability",
    title="Multigrid FSP: P(exactly m voxels in HIGH state)  —  symmetry breaking visible as P(all-hi) and P(all-lo) grow equally",
    titlesize=10)

cfg_clrs   = [:crimson, :darkorange, :steelblue4, :forestgreen, :purple][1:K+1]
cfg_names  = ["all-lo", "1-hi", "2-hi", "3-hi", "all-hi"]

all_cfg = [config_probs(sol_mg, ti, n_uns, K) for ti in 1:length(t_mg_hist)]
for m in 0:K
    p_m   = [all_cfg[ti][m+1] for ti in eachindex(t_mg_hist)]
    lbl   = "m=$m ($(cfg_names[m+1]))"
    lines!(ax_cfg, t_mg_hist, p_m; color=cfg_clrs[m+1], linewidth=2.5, label=lbl)
end

# SSA comparison — same config probs as scatter/errorbars
ssa_times = ssa_save[2:end]
for m in 0:K
    p_ssa = [ssa_config_probs(ssa_trajs, ti, n_uns, K, N_SSA)[m+1]
             for ti in eachindex(ssa_times)]
    scatter!(ax_cfg, ssa_times[1:4:end], p_ssa[1:4:end];
             color=cfg_clrs[m+1], markersize=5, marker=:circle, strokewidth=0.5,
             strokecolor=:white)
end
text!(ax_cfg, t_end*0.6, 0.95; text="── multigrid FSP   ● SSA (N=$N_SSA)",
      fontsize=9, color=:gray30)
axislegend(ax_cfg; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_cfg, -0.02, 1.05)

rowgap!(fig.layout, 1, 10); rowgap!(fig.layout, 2, 10)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_schlogl_instability.pdf")
save(out, fig)
println("Saved → $out")
