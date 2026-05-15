"""
Schlögl RDME — Bistable Front: Three-Method Comparison
=======================================================

K=4 voxels, 1D Schlögl (bistable: n_low≈31, n_high≈162).
IC: left 2 voxels in high state, right 2 in low state — a sharp front.

Demonstrates the failure modes of two older methods:

  1. Mean-field ODE
     - Solves dE[nₖ]/dt = f(E[nₖ]) + D·Lap E[nₖ]  (nonlinear PDE for means)
     - FAILURE: treats E[nₖ·(nₖ-1)] = E[nₖ]² — wrong for bistable distributions
     - Predicts front stays exactly fixed. No fluctuations, no interface motion.

  2. Per-voxel FSP with mean-field coupling
     - Each voxel solves full P(nₖ,t) under Schlögl kinetics
     - Neighboring coupling: inflow rate = D · E[nₖ_nb]  (mean-field)
     - FAILURE: front voxel's distribution is bimodal (P≈½ at 31, P≈½ at 162)
       so its mean ≈ 97 ≈ unstable fixed point. Mean-field coupling injects
       the UNSTABLE POINT into neighbors — completely wrong physics.

  3. Multigrid FSP (V-cycle, our method)
     - Pairs of adjacent voxels share a joint distribution P(n₁,n₂,t)
     - Captures which configuration each pair is in: (hi,hi),(hi,lo),(lo,hi),(lo,lo)
     - Front diffuses correctly: distribution of front position broadens over time.

Run: julia --project examples/fig_schlogl_front_comparison.jl
"""

using DiscStochSim
using ExponentialUtilities
using SparseArrays
using CairoMakie, Printf

# ── Parameters (all overridable via ENV) ──────────────────────────────────────
e(name, default) = parse(typeof(default), get(ENV, name, string(default)))

const K          = e("K",          4)
const D          = e("D",          0.05)   # slow diffusion → compact |S|
const h          = 1.0
const t_end      = e("T_END",      8.0)
const dt_mf      = e("DT_MF",      0.005)
const dt_fsp     = e("DT_FSP",     0.2)
const n_max      = e("N_MAX",      220)
const krylov_m   = e("KRYLOV_M",   30)
# Multigrid alg params — tuned for stable |S| ≈ 4000-5000
const n_levels      = e("N_LEVELS",   1)
const max_states    = e("MAX_STATES", 300_000)
const prune_tol     = e("PRUNE_TOL",  1e-4)   # aggressive: keeps 99.99% mass
const weight_tol    = e("WEIGHT_TOL", 1e-4)
const expand_depth  = e("EXP_DEPTH",  0)
const cr_depth      = e("CR_DEPTH",   1)
const cd_depth      = e("CD_DEPTH",   1)
const τ_pre         = e("TAU_PRE",    0.3)

model     = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)
n_low, n_uns, n_high = schlogl_fixed_points(model)
d = D / h^2

@printf("Schlögl: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
@printf("D=%.4f  d=%.4f  K=%d  t_end=%.1f  dt_fsp=%.2f\n", D, d, K, t_end, dt_fsp)
@printf("multigrid: n_levels=%d  prune_tol=%.0e  expand_depth=%d  cr_depth=%d  τ_pre=%.2f\n",
        n_levels, prune_tol, expand_depth, cr_depth, τ_pre)
const krylov_m = 30

# ── 1. Mean-field ODE ─────────────────────────────────────────────────────────
function schlogl_mf_rhs!(dμ, μ, model, K, d)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    for k in 1:K
        n = μ[k]
        rxn = c1*n*(n-1)/2 - c2*n*(n-1)*(n-2)/6 + c3 - c4*n
        lap = 0.0
        k > 1 && (lap += d*(μ[k-1] - n))
        k < K && (lap += d*(μ[k+1] - n))
        dμ[k] = rxn + lap
    end
end

μ_mf    = [Float64(n_high), Float64(n_high), Float64(n_low), Float64(n_low)]
dμ      = zeros(K)
t_mf_hist    = Float64[0.0]
μ_mf_hist    = [copy(μ_mf)]
n_steps_mf   = round(Int, t_end / dt_mf)
for _ in 1:n_steps_mf
    schlogl_mf_rhs!(dμ, μ_mf, model, K, d)
    μ_mf .+= dt_mf .* dμ
    push!(t_mf_hist, t_mf_hist[end] + dt_mf)
    push!(μ_mf_hist, copy(μ_mf))
end
println("Mean-field done.")

# ── 2. Per-voxel FSP with mean-field coupling ─────────────────────────────────
"""
Build 1D Schlögl CME generator for one voxel.
Inflow (from neighbours): Poisson source at rate μ_in = d·Σ E[n_nb]
Outflow: at rate d·n_nb·n (n = molecule count in this voxel)
"""
function build_schlogl_voxel_gen(model, d, n_nb, μ_in, n_max)
    c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4
    N = n_max + 1
    Is=Int[]; Js=Int[]; Vs=Float64[]
    addQ!(i,j,v) = (push!(Is,i); push!(Js,j); push!(Vs,v))

    for n in 0:n_max
        j = n + 1; diag = 0.0

        # Schlögl production: n → n+1
        r_up = c1*n*(n-1)/2 + c3
        if n < n_max && r_up > 0
            addQ!(j+1,j,r_up); diag -= r_up
        end
        # Schlögl degradation: n → n-1
        if n > 0
            r_dn = c2*n*(n-1)*(n-2)/6 + c4*n
            r_dn > 0 && (addQ!(j-1,j,r_dn); diag -= r_dn)
        end
        # Diffusion inflow (Poisson source)
        if n < n_max && μ_in > 0
            addQ!(j+1,j,μ_in); diag -= μ_in
        end
        # Diffusion outflow (per molecule)
        if n > 0
            r_out = d*n_nb*n
            addQ!(j-1,j,r_out); diag -= r_out
        end

        addQ!(j,j,diag)
    end
    sparse(Is,Js,Vs,N,N)
end

# Initialise: high state = Poisson(n_high) truncated; low = Poisson(n_low)
function poisson_vec(λ, n_max)
    p = zeros(n_max+1); p[1] = exp(-λ)
    for n in 1:n_max; p[n+1] = p[n]*λ/n; end
    s = sum(p); p ./= s; p
end

pvoxel_dists = [
    poisson_vec(Float64(n_high), n_max),
    poisson_vec(Float64(n_high), n_max),
    poisson_vec(Float64(n_low),  n_max),
    poisson_vec(Float64(n_low),  n_max),
]

pvoxel_means(p) = sum((i-1)*p[i] for i in eachindex(p))

t_pv_hist   = Float64[0.0]
μ_pv_hist   = [[pvoxel_means(p) for p in pvoxel_dists]]
n_steps_pv  = round(Int, t_end / dt_fsp)

snap_ts = [0.0, 2.0, 4.0, t_end]
snap_steps = Set(round(Int, ts/dt_fsp) for ts in snap_ts)
snap_pvdists = Dict{Int, Vector{Vector{Float64}}}()
snap_pvdists[0] = deepcopy(pvoxel_dists)

println("Running per-voxel FSP...")
for step in 1:n_steps_pv
    global pvoxel_dists
    μs = [pvoxel_means(p) for p in pvoxel_dists]
    new_dists = Vector{Float64}[]
    for k in 1:K
        n_nb = (k>1) + (k<K)
        μ_in = d * (k>1 ? μs[k-1] : 0.0) + d * (k<K ? μs[k+1] : 0.0)
        Q = build_schlogl_voxel_gen(model, d, n_nb, μ_in, n_max)
        p_new = expv(Float64(dt_fsp), Q, pvoxel_dists[k]; m=krylov_m)
        p_new .= max.(0.0, p_new); p_new ./= sum(p_new)
        push!(new_dists, p_new)
    end
    pvoxel_dists = new_dists
    push!(t_pv_hist, step*dt_fsp)
    push!(μ_pv_hist, [pvoxel_means(p) for p in pvoxel_dists])
    step in snap_steps && (snap_pvdists[step] = deepcopy(pvoxel_dists))
end
println("Per-voxel FSP done.")

# ── 3. Multigrid FSP (V-cycle) ────────────────────────────────────────────────
u0 = CartesianIndex(n_high, n_high, n_low, n_low)
alg = RDMEMultigridFSP(
    dt             = dt_fsp,
    n_levels       = n_levels,
    n_max          = n_max,
    save_every     = 1,
    max_states     = max_states,
    expand_depth   = expand_depth,
    coarse_r_depth = cr_depth,
    coarse_d_depth = cd_depth,
    prune_tol      = prune_tol,
    weight_tol     = weight_tol,
    τ_pre          = τ_pre,
)
println("Running multigrid FSP...")
@time sol_mg = solve_rdme_multigrid(model, fine_grid, u0, (0.0, t_end), Float64[], alg)
println("Multigrid done. $(length(sol_mg.t)) snapshots, |S|: $(sol_mg.state_space_sizes)")

μ_mg_hist  = mean_voxel_counts(sol_mg)      # n_t × K matrix
t_mg_hist  = sol_mg.t

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(1300, 820))

Label(fig[0, 1:4],
      "Schlögl RDME Bistable Front — $(K) voxels  D=$D  |  " *
      "n_low=$n_low  n_uns=$n_uns  n_high=$n_high  |  IC: (hi,hi,lo,lo)";
      fontsize=12, tellwidth=false)

colors = [:crimson, :darkorange, :steelblue4, :forestgreen]
vox_labels = ["voxel 1 (hi)", "voxel 2 (hi→?)", "voxel 3 (lo→?)", "voxel 4 (lo)"]

# ── Column 1: Mean-field voxel means ─────────────────────────────────────────
ax_mf = Axis(fig[1, 1]; xlabel="time", ylabel="E[n]",
             title="1. Mean-field ODE\n(E[n·(n-1)] ≈ E[n]²: wrong for bistable)",
             titlesize=10)
for k in 1:K
    μ_k = [μ_mf_hist[i][k] for i in eachindex(t_mf_hist)]
    lines!(ax_mf, t_mf_hist, μ_k; color=colors[k], linewidth=2, label=vox_labels[k])
end
hlines!(ax_mf, [Float64(n_high)]; color=(:black,0.3), linestyle=:dash, linewidth=1)
hlines!(ax_mf, [Float64(n_uns)];  color=(:gray50,0.5), linestyle=:dot,  linewidth=1)
hlines!(ax_mf, [Float64(n_low)];  color=(:black,0.3), linestyle=:dash, linewidth=1)
text!(ax_mf, t_end*0.02, Float64(n_high)+4; text="n_high=$n_high", fontsize=8)
text!(ax_mf, t_end*0.02, Float64(n_uns)+4;  text="n_uns=$n_uns",   fontsize=8, color=:gray50)
text!(ax_mf, t_end*0.02, Float64(n_low)+4;  text="n_low=$n_low",   fontsize=8)
axislegend(ax_mf; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_mf, -5, Float64(n_high)+20)

# ── Column 2: Per-voxel FSP means ─────────────────────────────────────────────
ax_pv = Axis(fig[1, 2]; xlabel="time", ylabel="E[n]",
             title="2. Per-voxel FSP + mean-field coupling\n(front mean ≈ n_uns → injects unstable point into neighbours)",
             titlesize=10)
for k in 1:K
    μ_k = [μ_pv_hist[i][k] for i in eachindex(t_pv_hist)]
    lines!(ax_pv, t_pv_hist, μ_k; color=colors[k], linewidth=2, label=vox_labels[k])
end
hlines!(ax_pv, [Float64(n_high)]; color=(:black,0.3), linestyle=:dash, linewidth=1)
hlines!(ax_pv, [Float64(n_uns)];  color=(:gray50,0.5), linestyle=:dot,  linewidth=1)
hlines!(ax_pv, [Float64(n_low)];  color=(:black,0.3), linestyle=:dash, linewidth=1)
axislegend(ax_pv; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_pv, -5, Float64(n_high)+20)

# ── Column 3: Multigrid FSP means ─────────────────────────────────────────────
ax_mg = Axis(fig[1, 3]; xlabel="time", ylabel="E[n]",
             title="3. Multigrid FSP (V-cycle)\n(joint P(n₁,n₂): correct front diffusion)",
             titlesize=10)
for k in 1:K
    μ_k = μ_mg_hist[:, k]
    lines!(ax_mg, t_mg_hist, μ_k; color=colors[k], linewidth=2, label=vox_labels[k])
end
hlines!(ax_mg, [Float64(n_high)]; color=(:black,0.3), linestyle=:dash, linewidth=1)
hlines!(ax_mg, [Float64(n_uns)];  color=(:gray50,0.5), linestyle=:dot,  linewidth=1)
hlines!(ax_mg, [Float64(n_low)];  color=(:black,0.3), linestyle=:dash, linewidth=1)
axislegend(ax_mg; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_mg, -5, Float64(n_high)+20)

# ── Column 4: State space size (multigrid) ─────────────────────────────────────
ax_ss = Axis(fig[1, 4]; xlabel="time", ylabel="|S|",
             title="Multigrid FSP state space\n(D=$(D), prune=$(prune_tol), expand=$(expand_depth), cr=$(cr_depth))",
             titlesize=10)
lines!(ax_ss, t_mg_hist, Float64.(sol_mg.state_space_sizes);
       color=:steelblue4, linewidth=2.5)

# ── Row 2: P(n) at front voxels (voxels 2 and 3) ─────────────────────────────
n_show = min(n_max, 220)

# Per-voxel FSP distributions at snap times
sorted_snap_steps = sort(collect(keys(snap_pvdists)))
snap_colors = Makie.resample_cmap(:plasma, length(sorted_snap_steps)+1)[1:end-1]

for (ki, vox) in enumerate([2, 3])
    ax = Axis(fig[2, ki]; xlabel="n", ylabel="P(n)",
              title="Per-voxel FSP: P(n) at voxel $vox ($(ki==1 ? "hi→?" : "lo→?"))",
              titlesize=10)
    for (ci, ss) in enumerate(sorted_snap_steps)
        t_s = t_pv_hist[ss + 1]
        p = snap_pvdists[ss][vox]
        lines!(ax, 0:n_show, p[1:n_show+1]; color=snap_colors[ci], linewidth=1.8,
               label="t=$(round(t_s,digits=1))")
    end
    vlines!(ax, [Float64(n_low), Float64(n_high)]; color=(:black,0.25), linestyle=:dash, linewidth=1)
    vlines!(ax, [Float64(n_uns)]; color=(:gray50,0.4), linestyle=:dot, linewidth=1)
    axislegend(ax; position=:rt, labelsize=8, framevisible=false)
    xlims!(ax, 0, n_show)
end

# Multigrid FSP marginals at snap times
function mg_voxel_marginal(sol, ti, vox, n_max)
    states, probs = sol.snapshots[ti]
    p = zeros(n_max+1)
    for (s, pr) in zip(states, probs)
        n = Tuple(s)[vox]; n <= n_max && (p[n+1] += pr)
    end
    s = sum(p); s > 0 && (p ./= s); p
end

mg_snap_ts = [0.0, 2.0, 4.0, t_end]
mg_snap_idxs = [argmin(abs.(t_mg_hist .- ts)) for ts in mg_snap_ts]

for (ki, vox) in enumerate([2, 3])
    ax = Axis(fig[2, ki+2]; xlabel="n", ylabel="P(n)",
              title="Multigrid FSP: P(n) at voxel $vox ($(ki==1 ? "hi→?" : "lo→?"))",
              titlesize=10)
    for (ci, ti) in enumerate(mg_snap_idxs)
        p = mg_voxel_marginal(sol_mg, ti, vox, n_max)
        lines!(ax, 0:n_show, p[1:n_show+1]; color=snap_colors[ci], linewidth=1.8,
               label="t=$(round(t_mg_hist[ti],digits=1))")
    end
    vlines!(ax, [Float64(n_low), Float64(n_high)]; color=(:black,0.25), linestyle=:dash, linewidth=1)
    vlines!(ax, [Float64(n_uns)]; color=(:gray50,0.4), linestyle=:dot, linewidth=1)
    axislegend(ax; position=:rt, labelsize=8, framevisible=false)
    xlims!(ax, 0, n_show)
end

rowgap!(fig.layout, 1, 12)
outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_schlogl_front_comparison.pdf")
save(out, fig)
println("Saved → $out")
