"""
Fine-grained FSP distribution in the active region — 2D BD wavefront.

What mean-field gives: E[n_k(t)]  (one number per voxel per time)
What FSP gives:        P(n_k=·, t) (full distribution at each active voxel)

Figure shows:
  Row 1 — E[n_k] heatmap at 5 times, with region borders (white=active, blue=equil)
  Row 2 — Std[n_k] heatmap at same times  (large in active front, shrinks at SS)
  Row 3 — Left:  E[n_k] ± Std along diagonal (the uncertainty envelope)
           Right: P(n) at the peak-front voxel, evolving from δ₀ → Poisson(μ_ss)

Run: julia --project examples/fig_bd2d_distribution.jl
"""

using DiscStochSim
using CairoMakie, Printf

# ─── Parameters ──────────────────────────────────────────────────────────────
const K     = parse(Int,     get(ENV, "SPATIAL_K",     "14"))
const k_b   = 1.0
const k_d   = 0.5         # μ_ss = 2.0
const D     = parse(Float64, get(ENV, "SPATIAL_D",     "1.0"))
const dx    = 1.0
const N0    = 0
const dt    = parse(Float64, get(ENV, "SPATIAL_DT",    "0.1"))
const t_end = parse(Float64, get(ENV, "SPATIAL_T_END", "30.0"))
const n_max = parse(Int,     get(ENV, "SPATIAL_N_MAX", "20"))

model = BirthDeathRDME(k_b, k_d, D, dx)
μ_ss  = ss_mean(model)    # 2.0

@printf("K=%d×%d  D=%.1f  μ_ss=%.1f  n_max=%d  dt=%.2f\n", K, K, D, μ_ss, n_max, dt)

# ─── Compute mean and std from a distribution vector ─────────────────────────
function dist_mean_std(p::Vector{Float64})
    μ  = sum((i-1)*p[i] for i in eachindex(p))
    σ2 = sum((i-1)^2 * p[i] for i in eachindex(p)) - μ^2
    μ, sqrt(max(0.0, σ2))
end

# ─── Build mean and std grids from SpatialFSP state ─────────────────────────
function mean_std_grids(s::SpatialFSP)
    mg = zeros(s.K_x, s.K_y)
    sg = zeros(s.K_x, s.K_y)
    for (idx, p) in s.dists
        μ, σ = dist_mean_std(p); mg[idx] = μ; sg[idx] = σ
    end
    for (idx, p) in s.equil_dists
        μ, σ = dist_mean_std(p); mg[idx] = μ; sg[idx] = σ
    end
    mg, sg
end

# ─── Voxel-border paths (batched) ────────────────────────────────────────────
function border_path(sg::Matrix{Int}, state::Int)
    xs = Float32[]; ys = Float32[]
    for idx in CartesianIndices(sg)
        sg[idx] == state || continue
        r, c = Tuple(idx)
        x0,x1 = c-0.5f0, c+0.5f0; y0,y1 = r-0.5f0, r+0.5f0
        append!(xs, (x0,x1,x1,x0,x0,NaN32)); append!(ys, (y0,y0,y1,y1,y0,NaN32))
    end
    xs, ys
end

# ─── Run SpatialFSP ──────────────────────────────────────────────────────────
const USE_FLUX_VCYCLE = get(ENV, "SPATIAL_SOLVER", "voxel") == "flux"

state  = SpatialFSP(model, K, K; n_max, ε_expand=0.3, ε_equil=0.08,
                    use_joint_active=false, use_block_vcycle=false)
center = CartesianIndex(K÷2+1, K÷2+1)
set_ic!(state, center, N0)

fg = USE_FLUX_VCYCLE ? build_flux_graph(CartesianIndex{2}[], K, K; α=0.7) : nothing
solver_label = USE_FLUX_VCYCLE ? "flux-vcycle" : "per-voxel"
println("Solver: $solver_label")

snap_ts    = unique([0.0; filter(<(t_end), [4.0, 9.0, 16.0]); t_end])
snap_steps = [max(0, round(Int, t/dt)) for t in snap_ts]

# Probe voxels along the diagonal — track their distributions over time
probe_offsets = [0, K÷6, K÷3]
probe_voxels  = [CartesianIndex(K÷2+1+o, K÷2+1+o) for o in probe_offsets]
probe_labels  = ["center (+$(probe_offsets[1]))", "+$(probe_offsets[2])", "+$(probe_offsets[3])"]

snap_mean    = Dict{Int, Matrix{Float64}}()
snap_std     = Dict{Int, Matrix{Float64}}()
snap_status  = Dict{Int, Matrix{Int}}()
probe_dists  = Dict{CartesianIndex{2}, Vector{Vector{Float64}}}(v => [] for v in probe_voxels)
probe_times  = Float64[]

t_hist       = Float64[]
n_act_hist   = Int[]
n_eq_hist    = Int[]
n_emp_hist   = Int[]

function record!(s, step)
    push!(t_hist, s.t); push!(n_act_hist, n_active(s))
    push!(n_eq_hist,  n_equilibrated(s)); push!(n_emp_hist, n_empty(s))
    if step ∈ snap_steps
        mg, sg = mean_std_grids(s)
        snap_mean[step] = mg; snap_std[step] = sg; snap_status[step] = status_grid(s)
    end
    # probe voxel distributions — save every 10 steps for smoothness
    if step % 10 == 0
        push!(probe_times, s.t)
        for v in probe_voxels
            p = haskey(s.dists, v)      ? copy(s.dists[v])       :
                haskey(s.equil_dists, v) ? copy(s.equil_dists[v]) :
                                           [1.0; zeros(n_max)]      # empty → δ_0
            # ensure length n_max+1
            if length(p) < n_max+1; append!(p, zeros(n_max+1-length(p))); end
            push!(probe_dists[v], p[1:n_max+1])
        end
    end
end

record!(state, 0)
n_steps = round(Int, t_end/dt)
println("Running $n_steps steps…")
for step in 1:n_steps
    USE_FLUX_VCYCLE ? step!(state, dt, fg) : step!(state, dt)
    record!(state, step)
    step % 50 == 0 && @printf("  t=%5.1f  active=%d  equil=%d\n",
                                state.t, n_active(state), n_equilibrated(state))
end
println("Done.")

# ─── Figure ──────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(1550, 1050))

μ_max_mean = μ_ss * 1.4
μ_max_std  = sqrt(μ_ss) * 1.5   # Poisson std at SS ≈ 1.41

# ── Rows 1–2: mean and std heatmaps ──────────────────────────────────────────
for (col, (step, t_s)) in enumerate(zip(snap_steps, snap_ts))
    mg = snap_mean[step]; sg = snap_std[step]; status = snap_status[step]

    # Row 1: mean
    ax_m = Axis(fig[1, col]; aspect=DataAspect(), title="t = $(round(t_s; digits=0))",
                titlesize=11, xticksvisible=false, yticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false)
    hm_m = heatmap!(ax_m, 1:K, 1:K, mg'; colormap=:inferno, colorrange=(0,μ_max_mean))
    xs_a, ys_a = border_path(status, 1)   # active
    xs_e, ys_e = border_path(status, 2)   # equil
    isempty(xs_a) || lines!(ax_m, xs_a, ys_a; color=:white,          linewidth=1.6)
    isempty(xs_e) || lines!(ax_m, xs_e, ys_e; color=:cornflowerblue, linewidth=0.8)
    col == length(snap_ts) && Colorbar(fig[1,col+1], hm_m; width=12, label="E[nₖ]",
        ticks=([0,μ_ss,μ_max_mean], ["0","μ_ss","$(round(μ_max_mean,digits=1))"]))

    # Row 2: std
    ax_s = Axis(fig[2, col]; aspect=DataAspect(),
                xticksvisible=false, yticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false)
    hm_s = heatmap!(ax_s, 1:K, 1:K, sg'; colormap=:magma, colorrange=(0,μ_max_std))
    isempty(xs_a) || lines!(ax_s, xs_a, ys_a; color=:white,          linewidth=1.6)
    isempty(xs_e) || lines!(ax_s, xs_e, ys_e; color=:cornflowerblue, linewidth=0.8)
    # Mark probe voxels
    for v in probe_voxels
        scatter!(ax_s, [Float32(Tuple(v)[2])], [Float32(Tuple(v)[1])];
                 color=:lime, marker=:cross, markersize=10, strokewidth=1)
    end
    col == length(snap_ts) && Colorbar(fig[2,col+1], hm_s; width=12, label="Std[nₖ]",
        ticks=([0,sqrt(μ_ss),μ_max_std], ["0","√μ_ss","$(round(μ_max_std,digits=1))"]))
end

# ── Row 3 left: E ± Std diagonal envelope over time ─────────────────────────
ax_diag = Axis(fig[3, 1:3]; xlabel="diagonal voxel k  (k,k)",
               ylabel="molecule count",
               title="E[nₖₖ] ± Std[nₖₖ]  —  FSP envelope along diagonal",
               titlesize=10)
clrs_d = [:steelblue4, :tomato, :darkorange, :purple, :forestgreen]
for (ci, (step, t_s)) in enumerate(zip(snap_steps, snap_ts))
    mg = snap_mean[step]; sg = snap_std[step]
    diag_k  = [k for k in 1:K]
    diag_mu = [mg[k,k] for k in 1:K]
    diag_sg = [sg[k,k] for k in 1:K]
    col = clrs_d[ci]
    band!(ax_diag, Float64.(diag_k), diag_mu .- diag_sg, diag_mu .+ diag_sg;
          color=(col, 0.18))
    lines!(ax_diag, diag_k, diag_mu; color=col, linewidth=2,
           label="t=$(round(Int,t_s))")
end
hlines!(ax_diag, [μ_ss]; color=(:gray50,0.5), linestyle=:dash, linewidth=1)
text!(ax_diag, K*0.65, μ_ss*1.05; text="μ_ss=$(μ_ss)", fontsize=9, color=:gray50)
axislegend(ax_diag; position=:rb, labelsize=9, framevisible=false)
ylims!(ax_diag, -0.15, μ_max_mean + 0.4)
# Add Std reference line at SS
hlines!(ax_diag, [μ_ss - sqrt(μ_ss), μ_ss + sqrt(μ_ss)];
        color=(:gray70,0.4), linestyle=:dot, linewidth=1)

# ── Row 3 right: P(n) at probe voxels over time ──────────────────────────────
probe_clrs = [:royalblue, :tomato, :forestgreen]
n_show = min(n_max, 12)

for (pi, v) in enumerate(probe_voxels)
    ax = Axis(fig[3, 3+pi]; xlabel="n",
              ylabel="P(nₖ = n)",
              title="P(n) at $(probe_labels[pi]) ($(Tuple(v)[1]),$(Tuple(v)[2]))",
              titlesize=10)

    # Select a few time indices to show
    n_probe = length(probe_times)
    show_idx = unique([1; [round(Int, f*n_probe) for f in [0.2, 0.45, 0.7, 1.0]]])
    show_idx = filter(i -> 1 ≤ i ≤ n_probe, show_idx)
    t_probe_clrs = Makie.resample_cmap(:plasma, length(show_idx)+1)[1:end-1]

    for (ci, idx) in enumerate(show_idx)
        p   = probe_dists[v][idx]
        t_l = probe_times[idx]
        lines!(ax, 0:n_show, p[1:n_show+1]; color=t_probe_clrs[ci], linewidth=1.6,
               label="t=$(round(t_l; digits=1))")
    end

    # Poisson(μ_ss) reference
    p_pois = [exp(-μ_ss) * μ_ss^n / factorial(n) for n in 0:n_show]
    lines!(ax, 0:n_show, p_pois; color=(:gray40,0.7), linestyle=:dash,
           linewidth=1.2, label="Poisson(μ_ss)")
    axislegend(ax; position=:rt, labelsize=8, framevisible=false)
    xlims!(ax, 0, n_show)
end

rowgap!(fig.layout, 1, 6); rowgap!(fig.layout, 2, 8)
colgap!(fig.layout, 5, 5)

Label(fig[0, 1:6],
      "2D BD-RDME wavefront  |  $(K)×$(K) grid  D=$D  μ_ss=$μ_ss  " *
      "|  FSP: full P(nₖ,t) per voxel  (not just mean-field E[nₖ])";
      fontsize=11, tellwidth=false)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_bd2d_distribution.pdf")
save(out, fig)
println("Saved → $out")
