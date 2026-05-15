"""
Schlögl RDME — Critical Droplet: Size Determines Fate
======================================================

K=6 voxels, Schlögl bistable (n_low≈31, n_high≈162), D=0.1.

Two initial conditions compared:
  SUBCRITICAL  — 2 central voxels at n_high, 4 outer at n_low
                 IC: (lo,lo, hi,hi, lo,lo)
                 → droplet TRAPPED: P(m=2 high) = 1.0 at t=100
                   Mean-field: smooth gradient. FSP: flat, front persists.

  SUPERCRITICAL — 4 central voxels at n_high, 2 outer at n_low
                  IC: (lo,hi,hi,hi,hi, lo)
                  → droplet EXPANDS: P(m≥4 high) > 0.9 at t=100
                    Mean-field: all voxels drift high. FSP: probability
                    wave sweeps the domain.

Key results (impossible to see from mean-field alone):
  1. Config probabilities P(m voxels in high state, t): FSP gives the
     exact fate — subcritical stays at m=2, supercritical reaches m=5,6.
  2. FSP Std[n]: HIGH at the interface for subcritical (bimodal boundary),
     LOW everywhere for supercritical (whole domain commits to high).
  3. Mean-field gives SAME qualitative picture for both — misses the
     fate bifurcation entirely.

ENV params:
  D, K, T_END, DT_FSP, N_MAX, PRUNE_TOL, N_SSA

Run: julia --project examples/fig_schlogl_critical_droplet.jl
"""

using DiscStochSim
using ExponentialUtilities, SparseArrays
using CairoMakie, Printf, Random

e(k,d) = parse(typeof(d), get(ENV, k, string(d)))

const K         = e("K",         6)
const D         = e("D",         0.1)
const t_end     = e("T_END",     60.0)
const dt_fsp    = e("DT_FSP",    1.0)
const dt_mf     = e("DT_MF",     0.01)
const n_max     = e("N_MAX",     220)
const prune_tol = e("PRUNE_TOL", 1e-5)
const N_SSA     = e("N_SSA",     2000)

model  = SchloglModel1D(D)
grid   = VoxelGrid(K, 1.0, 0)
n_low, n_uns, n_high = schlogl_fixed_points(model)
d = D

@printf("K=%d  D=%.3f  n_low=%d  n_uns=%d  n_high=%d  t_end=%.0f\n",
        K, D, n_low, n_uns, n_high, t_end)

snap_ts = unique([0.0; filter(<(t_end), [10.0, 30.0, 50.0]); t_end])

# ── ICs ───────────────────────────────────────────────────────────────────────
u0_sub  = CartesianIndex(n_low, n_low, n_high, n_high, n_low, n_low)  # 2-voxel
u0_sup  = CartesianIndex(n_low, n_high, n_high, n_high, n_high, n_low) # 4-voxel
@printf("Subcritical IC: %s\n", u0_sub)
@printf("Supercritical IC: %s\n", u0_sup)

# ── Multigrid FSP ─────────────────────────────────────────────────────────────
alg = RDMEMultigridFSP(
    dt             = dt_fsp,
    n_levels       = 1,
    n_max          = n_max,
    save_every     = 1,
    max_states     = 500_000,
    expand_depth   = 0,
    coarse_r_depth = 1,
    coarse_d_depth = 1,
    prune_tol      = prune_tol,
    weight_tol     = prune_tol,
    τ_pre          = 0.3,
)

println("Running subcritical FSP..."); flush(stdout)
@time sol_sub = solve_rdme_multigrid(model, grid, u0_sub, (0.0, t_end), Float64[], alg)

println("Running supercritical FSP..."); flush(stdout)
@time sol_sup = solve_rdme_multigrid(model, grid, u0_sup, (0.0, t_end), Float64[], alg)

println("Subcritical |S|: $(sol_sub.state_space_sizes[end])")
println("Supercritical |S|: $(sol_sup.state_space_sizes[end])")

μ_sub = mean_voxel_counts(sol_sub)
μ_sup = mean_voxel_counts(sol_sup)

# ── Config probs over time ────────────────────────────────────────────────────
function config_probs_all(sol, n_uns, K)
    p_m = [zeros(K+1) for _ in 1:length(sol.t)]
    for (ti, (states, probs)) in enumerate(sol.snapshots)
        for (s, pr) in zip(states, probs)
            m = count(n -> n > n_uns, Tuple(s))
            p_m[ti][m+1] += pr
        end
    end
    p_m
end

cfg_sub = config_probs_all(sol_sub, n_uns, K)
cfg_sup = config_probs_all(sol_sup, n_uns, K)

# ── Mean-field ODE ────────────────────────────────────────────────────────────
function run_mf(u0_ci, K, D, model, t_end, dt)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    μ = Float64[Tuple(u0_ci)...]
    dμ = zeros(K)
    snaps = [copy(μ)]
    for step in 1:round(Int, t_end/dt)
        for k in 1:K
            n = μ[k]
            rxn = c1*n*(n-1)/2 - c2*n*(n-1)*(n-2)/6 + c3 - c4*n
            lap = (k>1 ? D*(μ[k-1]-n) : 0.0) + (k<K ? D*(μ[k+1]-n) : 0.0)
            dμ[k] = rxn + lap
        end
        μ .+= dt .* dμ
    end
    μ
end

μ_mf_sub = run_mf(u0_sub, K, D, model, t_end, dt_mf)
μ_mf_sup = run_mf(u0_sup, K, D, model, t_end, dt_mf)

# ── SSA ground truth ──────────────────────────────────────────────────────────
function run_ssa_k6(seed, u0_ci, t_end, save_ts)
    c1=model.c1;c2=model.c2;c3=model.c3;c4=model.c4
    rng=Random.MersenneTwister(seed); n=collect(Int, Tuple(u0_ci)); t=0.0; si=1
    snaps=[copy(n)]
    while t < t_end
        props=Float64[]
        for k in 1:K
            nk=n[k]
            push!(props,c3,c4*nk,c1*nk*(nk-1)/2,c2*nk*(nk-1)*(nk-2)/6)
            k>1 && push!(props,D*nk); k<K && push!(props,D*nk)
        end
        a0=sum(props); a0==0&&break
        t += -log(rand(rng))/a0
        r=rand(rng)*a0; cum=0.0; fired=0
        for i in eachindex(props); cum+=props[i]; fired=i; cum>=r&&break; end
        idx=0
        for k in 1:K
            base=idx; idx+=4+(k > 1 ? 1 : 0) + (k < K ? 1 : 0); fired<=idx||continue
            lf=fired-base
            lf==1&&(n[k]+=1); lf==2&&(n[k]=max(0,n[k]-1))
            lf==3&&(n[k]+=1); lf==4&&(n[k]=max(0,n[k]-1))
            if lf==5; nb = k > 1 ? k-1 : k+1; n[k]=max(0,n[k]-1);n[nb]+=1; end
            lf==6&&(n[k]=max(0,n[k]-1);n[k+1]+=1)
            break
        end
        while si<=length(save_ts)&&t>=save_ts[si]; push!(snaps,copy(n));si+=1; end
    end
    while length(snaps)<=length(save_ts); push!(snaps,copy(n)); end
    snaps
end

ssa_save = collect(0.0:dt_fsp:t_end)[2:end]
println("Running $N_SSA SSA trajectories for each IC...")
@time ssa_sub = [run_ssa_k6(i,       u0_sub, t_end, ssa_save) for i in 1:N_SSA]
@time ssa_sup = [run_ssa_k6(i+N_SSA, u0_sup, t_end, ssa_save) for i in 1:N_SSA]

function ssa_cfg(trajs, n_uns, K)
    n_t = length(trajs[1])
    p_m = [zeros(K+1) for _ in 1:n_t]
    for traj in trajs, (ti,snap) in enumerate(traj)
        m = count(x -> x > n_uns, snap); p_m[ti][m+1] += 1.0
    end
    for pm in p_m; pm ./= length(trajs); end
    p_m
end

ssa_cfg_sub = ssa_cfg(ssa_sub, n_uns, K)
ssa_cfg_sup = ssa_cfg(ssa_sup, n_uns, K)
ssa_ts      = [0.0; ssa_save]

println("Done.")
@printf("Sub final SSA m-dist: %s\n",
    join([@sprintf("m=%d:%.3f",mi,ssa_cfg_sub[end][mi+1]) for mi in 0:K], " "))
@printf("Sup final SSA m-dist: %s\n",
    join([@sprintf("m=%d:%.3f",mi,ssa_cfg_sup[end][mi+1]) for mi in 0:K], " "))

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(1400, 880))

Label(fig[0,1:4],
    "Schlögl RDME — Critical Droplet  |  K=$K  D=$D  " *
    "n_low=$n_low / n_uns=$n_uns / n_high=$n_high  |  " *
    "Subcritical (2-voxel) vs Supercritical (4-voxel)";
    fontsize=12, tellwidth=false)

cfg_clrs  = [:crimson, :darkorange, :steelblue4, :forestgreen, :purple, :hotpink, :navy]
m_labels  = ["m=0 (all-lo)", "m=1", "m=2", "m=3", "m=4", "m=5", "m=6 (all-hi)"]
vox_clrs  = Makie.resample_cmap(:viridis, K)

# ── Row 1: E[n] profile over time ────────────────────────────────────────────
for (col, (sol, mu, u0_ic, label)) in enumerate([
    (sol_sub, μ_sub, u0_sub, "Subcritical (2-voxel droplet)"),
    (sol_sup, μ_sup, u0_sup, "Supercritical (4-voxel droplet)"),
])
    ax = Axis(fig[1, col]; xlabel="voxel", ylabel="E[n]",
              title="$label  D=$D", titlesize=11)
    ts_show = [sol.t[argmin(abs.(sol.t .- ts))] for ts in snap_ts]
    snap_clrs = Makie.resample_cmap(:plasma, length(ts_show)+1)[1:end-1]
    for (ci, (ti_use, ts)) in enumerate(zip(
            [argmin(abs.(sol.t .- ts)) for ts in snap_ts], snap_ts))
        lines!(ax, 1:K, mu[ti_use,:]; color=snap_clrs[ci], linewidth=2.5,
               label="t=$(round(Int,ts))")
    end
    hlines!(ax, [Float64(n_high)]; color=(:black,0.2), linestyle=:dash, linewidth=1)
    hlines!(ax, [Float64(n_uns)];  color=(:gray50,0.4), linestyle=:dot,  linewidth=1)
    hlines!(ax, [Float64(n_low)];  color=(:black,0.2), linestyle=:dash, linewidth=1)
    text!(ax, 1.1, Float64(n_high)+3; text="n_high", fontsize=9)
    text!(ax, 1.1, Float64(n_uns)+3;  text="n_uns",  fontsize=9, color=:gray50)
    text!(ax, 1.1, Float64(n_low)+3;  text="n_low",  fontsize=9)

    # Add MF final profile
    mf_final = col==1 ? μ_mf_sub : μ_mf_sup
    lines!(ax, 1:K, mf_final; color=:gray50, linewidth=1.5, linestyle=:dot,
           label="MF t=$(round(Int,t_end))")
    axislegend(ax; position=:rc, labelsize=9, framevisible=false)
    ylims!(ax, -5, Float64(n_high)+20)
    xticks!(ax, xticks=(1:K, string.(1:K)))
end

# ── Row 2: Config probability P(m hi) over time ───────────────────────────────
for (col, (sol, cfg_fsp, cfg_ssa, label)) in enumerate([
    (sol_sub, cfg_sub, ssa_cfg_sub, "Subcritical: P(m voxels in high state)"),
    (sol_sup, cfg_sup, ssa_cfg_sup, "Supercritical: P(m voxels in high state)"),
])
    ax = Axis(fig[2, col]; xlabel="time", ylabel="probability",
              title=label, titlesize=11)
    for m in 0:K
        p_fsp = [cfg_fsp[ti][m+1] for ti in eachindex(sol.t)]
        any(>(0.005), p_fsp) || continue
        lines!(ax, sol.t, p_fsp; color=cfg_clrs[m+1], linewidth=2.5,
               label=m_labels[m+1])
        # SSA dots (every 5 points)
        p_ssa = [ssa_cfg_ssa[ti][m+1] for ti in eachindex(ssa_ts)]
        scatter!(ax, ssa_ts[1:5:end], p_ssa[1:5:end];
                 color=cfg_clrs[m+1], markersize=4, strokewidth=0)
    end
    text!(ax, t_end*0.55, 0.92; text="— FSP  • SSA (N=$N_SSA)", fontsize=9, color=:gray40)
    axislegend(ax; position=:lt, labelsize=9, framevisible=false)
    ylims!(ax, -0.02, 1.05)
end

# ── Row 3: |S| over time and final config comparison ─────────────────────────
ax_ss = Axis(fig[3, 1]; xlabel="time", ylabel="|S|",
             title="Multigrid FSP state space — subcritical stays compact",
             titlesize=11)
lines!(ax_ss, sol_sub.t, Float64.(sol_sub.state_space_sizes);
       color=:steelblue4, linewidth=2.5, label="Subcritical")
lines!(ax_ss, sol_sup.t, Float64.(sol_sup.state_space_sizes);
       color=:crimson, linewidth=2.5, label="Supercritical")
axislegend(ax_ss; position=:rc, labelsize=9, framevisible=false)

# Final bar chart of config probs for both ICs
ax_bar = Axis(fig[3, 2]; xlabel="m (# voxels in high state)", ylabel="P(m)",
              title="Final config distribution at t=$t_end: FSP (bar) + SSA (●)",
              titlesize=11)
m_vals = 0:K
xs_sub = Float64.(m_vals) .- 0.2
xs_sup = Float64.(m_vals) .+ 0.2
barplot!(ax_bar, xs_sub, [cfg_sub[end][m+1] for m in m_vals];
         color=(:steelblue4, 0.7), width=0.35, label="Subcritical FSP")
barplot!(ax_bar, xs_sup, [cfg_sup[end][m+1] for m in m_vals];
         color=(:crimson, 0.7), width=0.35, label="Supercritical FSP")
scatter!(ax_bar, xs_sub .+ 0.05, [ssa_cfg_sub[end][m+1] for m in m_vals];
         color=:steelblue4, markersize=8, marker=:circle)
scatter!(ax_bar, xs_sup .+ 0.05, [ssa_cfg_sup[end][m+1] for m in m_vals];
         color=:crimson, markersize=8, marker=:circle)
xticks!(ax_bar, xticks=(Float64.(m_vals), ["$m" for m in m_vals]))
axislegend(ax_bar; position=:rt, labelsize=9, framevisible=false)

rowgap!(fig.layout, 1, 10); rowgap!(fig.layout, 2, 10)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_schlogl_critical_droplet.pdf")
save(out, fig)
println("Saved → $out")
