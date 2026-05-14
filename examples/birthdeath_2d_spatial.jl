"""
2D Birth-Death Spatial FSP Demo — with voxel status visualization
==================================================================

Shows the full lifecycle of the spatially adaptive FSP:

  Phase 1 (SPREAD):    wavefront expands from center; active count 1 → K²
  Phase 2 (COLLAPSE):  equilibrated voxels merge back; active count K² → 0

Each voxel is in one of three states at any time:
  EMPTY       (dark)   — not yet reached by the wavefront
  ACTIVE      (orange) — currently tracked by the FSP
  EQUILIBRATED(blue)   — merged back; known to be Poisson(μ_ss)

Run: julia --project examples/birthdeath_2d_spatial.jl
"""

using DiscStochSim
using CairoMakie, Printf

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

const K     = 14
const k_b   = 1.0
const k_d   = 0.5     # μ_ss = k_b/k_d = 2
const D     = 1.0     # slow diffusion: wavefront spreads visibly
const dx    = 1.0
const N0    = 0       # start empty — birth drives the wavefront
const dt    = 0.1
const t_end = 30.0
const n_max = 30

model = BirthDeathRDME(k_b, k_d, D, dx)
μ_ss  = ss_mean(model)
@printf("SS mean μ_ss = %.1f   jump rate d = %.1f\n", μ_ss, jump_rate(model))
@printf("Grid %d×%d = %d voxels   IC: N0=%d at center\n", K, K, K^2, N0)
@printf("Spread timescale ~ %.1f  Equil timescale ~ %.1f\n",
        (K/2)^2/(2D), 1/k_d)

# ─────────────────────────────────────────────────────────────────────────────
# Initialize
# ─────────────────────────────────────────────────────────────────────────────

state  = SpatialFSP(model, K, K; n_max, ε_expand=0.3, ε_equil=0.08)
center = CartesianIndex(K÷2 + 1, K÷2 + 1)
set_ic!(state, center, N0)

# ─────────────────────────────────────────────────────────────────────────────
# Snap times — chosen to show all three phases
# ─────────────────────────────────────────────────────────────────────────────

snap_ts    = [0.0, 4.0, 9.0, 16.0, t_end]
snap_steps = [max(0, round(Int, t / dt)) for t in snap_ts]

snap_status = Dict{Int, Matrix{Int}}()      # 0=empty 1=active 2=equilibrated
snap_means  = Dict{Int, Matrix{Float64}}()  # E[n] per voxel
snap_center = Dict{Int, Vector{Float64}}()  # center voxel P(n)
center_mean_hist = Float64[]

t_hist         = Float64[]
n_act_hist     = Int[]
n_eq_hist      = Int[]
n_empty_hist   = Int[]
total_mean_hist = Float64[]

function record!(s, step)
    push!(t_hist,          s.t)
    push!(n_act_hist,      n_active(s))
    push!(n_eq_hist,       n_equilibrated(s))
    push!(n_empty_hist,    n_empty(s))
    push!(total_mean_hist,
          sum(voxel_mean(p) for (_, p) in s.dists; init=0.0) +
          sum(voxel_mean(p) for (_, p) in s.equil_dists; init=0.0))
    # Center mean: active → from dist, equilibrated → μ_ss, empty → 0
    push!(center_mean_hist,
          haskey(s.dists, center)      ? voxel_mean(s.dists[center]) :
          haskey(s.equil_dists, center) ? voxel_mean(s.equil_dists[center]) : 0.0)
    if step ∈ snap_steps
        snap_status[step] = status_grid(s)
        snap_means[step]  = mean_grid(s)   # includes equil means — used for heatmap
        p_c = get(s.dists, center, nothing)
        snap_center[step] = p_c !== nothing ? copy(p_c) : poisson_pmf_vec(μ_ss, n_max)
    end
end

record!(state, 0)

n_steps = round(Int, t_end / dt)
println("\nRunning $n_steps steps...")
for step in 1:n_steps
    step!(state, dt)
    record!(state, step)
    if step % 40 == 0
        @printf("  t=%5.1f  empty=%3d  active=%3d  equil=%3d  total_mean=%5.1f / %d\n",
                state.t, n_empty(state), n_active(state), n_equilibrated(state),
                total_mean_hist[end], K^2 * μ_ss)
    end
end
println("Done.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

# Color encoding:  empty = dark gray,  active = warm (by mean),  equilibrated = steel blue
# ─────────────────────────────────────────────────────────────────────────────
# Helpers for voxel-boundary visualization
# ─────────────────────────────────────────────────────────────────────────────

"""
Build a single batched path for all voxels of a given state.
Each voxel rectangle is appended as 5 points (closed loop) + NaN separator,
so one `lines!` call draws all borders for that state.
sg[row, col] = status; heatmap is plotted transposed so x=col, y=row.
"""
function voxel_border_path(sg::Matrix{Int}, state::Int)
    xs = Float32[];  ys = Float32[]
    for idx in CartesianIndices(sg)
        sg[idx] == state || continue
        row, col = Tuple(idx)
        x0, x1 = col - 0.5f0, col + 0.5f0
        y0, y1 = row - 0.5f0, row + 0.5f0
        append!(xs, [x0, x1, x1, x0, x0, NaN32])
        append!(ys, [y0, y0, y1, y1, y0, NaN32])
    end
    xs, ys
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

μ_max = 2.5μ_ss

fig = Figure(size = (1500, 1050))
Label(fig[0, 1:5],
      "2D Birth-Death Spatial FSP  |  $(K)×$(K) grid  D=$(D)  μ_ss=$(μ_ss)  " *
      "  borders: white=active, blue=equilibrated, none=empty";
      fontsize = 12)

# ── Row 1: Spatial snapshots with voxel borders ───────────────────────────────
for (col, (step, t)) in enumerate(zip(snap_steps, snap_ts))
    ax = Axis(fig[1, col];
              title     = @sprintf("t = %.0f", t),
              aspect    = DataAspect(),
              titlesize = 12,
              xticksvisible = false, yticksvisible = false,
              xticklabelsvisible = false, yticklabelsvisible = false)

    mg = snap_means[step]       # K×K matrix; mg[row,col] = E[n]
    sg = snap_status[step]

    # Heatmap: transposed so column index = x, row index = y
    hm = heatmap!(ax, 1:K, 1:K, mg';
                  colormap   = :inferno,
                  colorrange = (0.0, μ_max))

    # Active voxel borders — white, prominent
    xs_a, ys_a = voxel_border_path(sg, 1)
    isempty(xs_a) || lines!(ax, xs_a, ys_a; color = :white,          linewidth = 1.8)

    # Equilibrated voxel borders — cornflower blue, thin
    xs_e, ys_e = voxel_border_path(sg, 2)
    isempty(xs_e) || lines!(ax, xs_e, ys_e; color = :cornflowerblue, linewidth = 0.7)

    # Mark center voxel
    ci, cj = Tuple(center)
    scatter!(ax, [Float32(cj)], [Float32(ci)]; color = :red, markersize = 8, marker = :cross)

    na = count(==(1), sg);  ne = count(==(2), sg);  nu = count(==(0), sg)
    text!(ax, 0.5, 0.02;
          text = "∅=$nu  A=$na  ✓=$ne",
          space = :relative, align = (:center, :bottom),
          color = :white, fontsize = 9)

    col == length(snap_steps) &&
        Colorbar(fig[1, col+1], hm; label="E[n]", width=12)
end

# ── Row 2: Count arc + total mean ────────────────────────────────────────────
ax_counts = Axis(fig[2, 1:3];
                 xlabel = "time", ylabel = "voxel count",
                 title  = "Voxel states over time — active peaks at wavefront width, not K²",
                 titlesize = 11)
band!(ax_counts, t_hist, zeros(length(t_hist)), Float64.(n_act_hist); color=(:tomato,0.3))
band!(ax_counts, t_hist, Float64.(n_act_hist), Float64.(n_act_hist).+Float64.(n_eq_hist);
      color=(:cornflowerblue,0.3))
lines!(ax_counts, t_hist, Float64.(n_act_hist);  color=:tomato,          linewidth=2.5, label="Active (FSP)")
lines!(ax_counts, t_hist, Float64.(n_eq_hist);   color=:cornflowerblue,  linewidth=2.5, label="Equilibrated (merged)")
lines!(ax_counts, t_hist, Float64.(n_empty_hist); color=:gray,           linewidth=2,   linestyle=:dash, label="Empty")
hlines!(ax_counts, [Float64(K^2)]; color=:black, linewidth=1, linestyle=:dot)
text!(ax_counts, 1.0, K^2*0.97; text="K²=$(K^2)", fontsize=9, color=:black)
axislegend(ax_counts; position=:rc, labelsize=10)
ylims!(ax_counts, -5, K^2 + 20)

# Corrected ODE: birth only from activated voxels (not all K²)
# dE[N]/dt = (K² - n_empty(t)) · k_b  -  k_d · E[N]
# Numerically integrate using the recorded n_empty_hist.
N_corr = zeros(length(t_hist))
for i in 2:length(t_hist)
    dt_i = t_hist[i] - t_hist[i-1]
    n_producing = K^2 - n_empty_hist[i-1]   # active + equilibrated voxels
    dN = n_producing * k_b - k_d * N_corr[i-1]
    N_corr[i] = N_corr[i-1] + dN * dt_i
end

ax_mu = Axis(fig[2, 4:5];
             xlabel="time", ylabel="E[N_total]",
             title ="Why the naive ODE lags: empty voxels contribute zero birth rate",
             titlesize=11)
lines!(ax_mu, t_hist, total_mean_hist;
       color=:purple, linewidth=2.5, label="FSP (simulation)")
N_ode = K^2*μ_ss .+ (N0 - K^2*μ_ss) .* exp.(-k_d .* t_hist)
lines!(ax_mu, t_hist, N_ode;
       color=:black, linewidth=1.5, linestyle=:dash,
       label="Naive ODE (all K² voxels from t=0)")
lines!(ax_mu, t_hist, N_corr;
       color=:tomato, linewidth=2, linestyle=:dot,
       label="Corrected ODE: (K²-n_empty)·k_b - k_d·E[N]")
hlines!(ax_mu, [K^2*μ_ss]; color=:gray, linewidth=1, linestyle=:dot)
axislegend(ax_mu; position=:rb, labelsize=10)

# ── Row 3: Center voxel mean over time | Final total distribution ─────────────
ax_cm = Axis(fig[3, 1:2];
             xlabel="time", ylabel="E[n_center]",
             title ="Center voxel mean: fast rise → local SS → slow climb to global SS",
             titlesize=11)
lines!(ax_cm, t_hist, center_mean_hist; color=:darkorange, linewidth=2.5)
hlines!(ax_cm, [μ_ss]; color=:black, linewidth=1.5, linestyle=:dash)
# Mark snap times
for (step, t) in zip(snap_steps, snap_ts)
    vlines!(ax_cm, [t]; color=:gray, linewidth=1, linestyle=:dot)
end
text!(ax_cm, t_end*0.7, μ_ss*1.05; text="μ_ss=$(μ_ss)", fontsize=10, color=:black)

ax_td = Axis(fig[3, 3:5];
             xlabel="N_total", ylabel="P(N_total)",
             title ="Total count at t=$(t_end): collapsed to Poisson(K²·μ_ss=$(Int(K^2*μ_ss)))",
             titlesize=11)
p_final  = total_distribution(state)
μ_tot    = K^2 * μ_ss
N_show   = min(length(p_final)-1, round(Int, μ_tot + 5*sqrt(μ_tot)) + 10)
p_pois_T = poisson_pmf_vec(μ_tot, N_show)
lines!(ax_td, 0:N_show, p_final[1:N_show+1]; color=:steelblue, linewidth=2.5, label="FSP (product form)")
lines!(ax_td, 0:N_show, p_pois_T;            color=:crimson,   linewidth=2,   linestyle=:dash,
       label="Poisson($(Int(μ_tot)))")
axislegend(ax_td; position=:rt, labelsize=10)

outdir  = joinpath(@__DIR__, "output")
mkpath(outdir)
outpath = joinpath(outdir, "birthdeath_2d_spatial.png")
save(outpath, fig; px_per_unit = 2)
println("Saved: $outpath")
