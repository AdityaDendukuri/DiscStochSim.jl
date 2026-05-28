"""
fig_mg_validation.jl

Validate multigrid methods against per-voxel baseline and KPP wave speed.

For BranchingDeath/KPP waves: v_theory = 2√(D(k_b-k_d)) = 2√(0.4) ≈ 1.265 vox/s

Methods compared:
  - Per-voxel         exact mean-field RDME (baseline)
  - Block exact       exact Galerkin 2×2 joint CME — no total-count restriction
  - ADI pair vcycle   ADI pair joint CME with multi-level prolongation
  - MG L=2 (broken)  total-count coarsening — shows why restriction fails
"""

using DiscStochSim, CairoMakie, Printf, Statistics

const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4; const T_END=160.0; const N_FRM=100
const n_max=14
const v_theory = 2*sqrt(D*(k_b-k_d))
@printf("Theoretical KPP wave speed: %.4f vox/s\n\n", v_theory)

function run_sim(; label="", use_multigrid=false, mg_levels=0,
                   use_block_exact=false, use_block_vcycle=false,
                   use_block_exact_joint=false)
    @printf("Running %s ...\n", label); flush(stdout)
    model = BranchingDeathRDME(k_b, k_d, D, 1.0)
    state = SpatialFSP(model, K, K;
        n_max, ε_expand=0.02, ε_equil=0.15,
        use_multigrid, multigrid_levels=mg_levels,
        use_block_exact, use_block_vcycle, use_block_exact_joint)
    ci,cj = K÷2+1, K÷2+1
    for di in -4:4, dj in -4:4
        di^2+dj^2<=16 || continue
        1<=ci+di<=K && 1<=cj+dj<=K || continue
        set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
    end

    # Track outermost active voxel radius at each frame (true leading edge)
    MAX_R = min(ci-1, K-ci, cj-1, K-cj)
    t_frames  = collect(range(0.0, T_END; length=N_FRM))
    front_r   = zeros(N_FRM)   # max radius of any active voxel
    profiles  = Matrix{Float64}(undef, N_FRM, MAX_R)

    function radial_profile()
        prof = zeros(MAX_R); cnt = zeros(Int, MAX_R)
        mg = mean_grid(state)
        for i in 1:K, j in 1:K
            r = round(Int, sqrt((i-ci)^2+(j-cj)^2))
            1<=r<=MAX_R || continue
            prof[r] += mg[i,j]; cnt[r] += 1
        end
        [cnt[r]>0 ? prof[r]/cnt[r] : 0.0 for r in 1:MAX_R]
    end

    function max_active_radius()
        isempty(state.dists) && return 0.0
        maximum(sqrt((k[1]-ci)^2+(k[2]-cj)^2) for k in keys(state.dists))
    end

    profiles[1,:] = radial_profile()
    front_r[1]    = max_active_radius()

    t0 = time(); fi=2
    for step_i in 1:round(Int,T_END/dt)
        step!(state, dt)
        t = step_i*dt
        if fi<=N_FRM && t>=t_frames[fi]-1e-9
            profiles[fi,:] = radial_profile()
            front_r[fi]    = max_active_radius()
            fi+=1
        end
        fi>N_FRM && break
    end
    elapsed = time()-t0
    @printf("  done in %.1fs\n", elapsed); flush(stdout)
    profiles, t_frames, elapsed, front_r
end

# Each config: (label, kwargs...)
configs = [
    ("Per-voxel (baseline)",       (use_multigrid=false, mg_levels=0, use_block_exact=false, use_block_vcycle=false, use_block_exact_joint=false)),
    ("Pair exact (Galerkin)",      (use_multigrid=false, mg_levels=0, use_block_exact=true,  use_block_vcycle=false, use_block_exact_joint=false)),
    ("Block exact 2×2 (Galerkin)", (use_multigrid=false, mg_levels=0, use_block_exact=false, use_block_vcycle=false, use_block_exact_joint=true)),
    ("ADI pair vcycle (restrict)", (use_multigrid=false, mg_levels=0, use_block_exact=false, use_block_vcycle=true,  use_block_exact_joint=false)),
    ("MG L=2 (total-count)",       (use_multigrid=true,  mg_levels=2, use_block_exact=false, use_block_vcycle=false, use_block_exact_joint=false)),
]

results = [run_sim(; label=lbl, kw...) for (lbl,kw) in configs]

# Fit wave speed from outermost-active-radius vs time (leading edge)
function fit_speed(t_frames, front_r)
    # exclude early transient and boundary effects
    mask = (t_frames .> 15.0) .& (front_r .< 80.0) .& (front_r .> 0)
    ts = t_frames[mask]; rs = front_r[mask]
    length(ts) < 4 && return NaN
    t_mean = mean(ts); r_mean = mean(rs)
    num = sum((ts.-t_mean).*(rs.-r_mean))
    den = sum((ts.-t_mean).^2)
    den > 0 ? num/den : NaN
end

println("\nWave speed estimates (outermost active voxel):")
println("─"^60)
speeds = Float64[]
lbls   = [c[1] for c in configs]
for (r, (lbl,_)) in zip(results, configs)
    v = fit_speed(r[2], r[4])
    push!(speeds, v)
    err_pct = isnan(v) ? NaN : abs(v - v_theory)/v_theory*100
    @printf("  %-32s  v=%.3f  err=%.1f%%\n", lbl, v, err_pct)
end
@printf("  %-32s  v=%.3f  (reference)\n", "Theory (KPP)", v_theory)

outdir = joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)

# ── Figure 0: front radius vs time ───────────────────────────────────────────
cols_fig = [:gray50, :steelblue, :seagreen, :darkorange, :firebrick]
fig0 = Figure(size=(700,400))
ax0 = Axis(fig0[1,1];
    title="Leading-edge front radius vs time  (BranchingDeath RDME, K=$(K)×$(K))",
    xlabel="time", ylabel="max active radius (voxels)", titlesize=11)
for (r,(lbl,_),col) in zip(results,configs,cols_fig)
    lines!(ax0, r[2], r[4]; color=col, linewidth=2, label=lbl)
end
lines!(ax0, results[1][2], v_theory.*results[1][2]; color=:black,
       linewidth=1.5, linestyle=:dash, label="KPP theory (v=$(round(v_theory,digits=3)))")
axislegend(ax0; position=:lt, labelsize=9)
save(joinpath(outdir,"fig_mg_frontradius.png"), fig0; px_per_unit=3)
println("Saved → fig_mg_frontradius.png")

# ── Figure 1: wave speed bar chart ───────────────────────────────────────────
n_methods = length(configs)
fig1 = Figure(size=(800,420))
ax1 = Axis(fig1[1,1];
    title="Wave speed by method  (BranchingDeath RDME, K=$(K)×$(K))",
    ylabel="wave speed (vox/s)", titlesize=11,
    xticks=(1:n_methods, lbls), xticklabelsize=9)
barplot!(ax1, 1:n_methods, speeds; color=cols_fig[1:n_methods], strokewidth=0)
hlines!(ax1, [v_theory]; color=:black, linewidth=2, linestyle=:dash,
        label="KPP theory (v=$(round(v_theory,digits=3)))")
for (i,v) in enumerate(speeds)
    text!(ax1, i, v+0.01; text=@sprintf("%.3f", v),
          fontsize=9, align=(:center,:bottom))
end
axislegend(ax1; position=:rb, labelsize=9)
ylims!(ax1, 0, max(v_theory*1.3, maximum(filter(!isnan, speeds))*1.25))
save(joinpath(outdir,"fig_mg_speed.png"), fig1; px_per_unit=3)
println("\nSaved → fig_mg_speed.png")

# ── Figure 2: space-time heatmaps ────────────────────────────────────────────
fig2 = Figure(size=(n_methods*320+80, 500))
Label(fig2[0,1:n_methods],
    "Space-time diagram: wavefront propagation (BranchingDeath K=$(K)×$(K))",
    fontsize=11, font=:bold)

xt_max = maximum(maximum(r[1]) for r in results; init=1e-3)
t_ticks_vals = range(1,N_FRM,length=6)
t_ticks = (collect(t_ticks_vals),
           [@sprintf("%.0f", results[1][2][round(Int,f)]) for f in t_ticks_vals])

for (c,(r,(lbl,_))) in enumerate(zip(results,configs))
    ax = Axis(fig2[1,c]; title=lbl, titlesize=9,
              xlabel="radius", ylabel=c==1 ? "time" : "",
              xlabelsize=8, ylabelsize=8,
              xticklabelsize=7, yticklabelsize=7,
              yticks=t_ticks)
    heatmap!(ax, 1:size(r[1],2), 1:N_FRM, r[1]';
             colormap=:viridis, colorrange=(0,xt_max))
    r_theory = [v_theory * results[1][2][f] for f in 1:N_FRM]
    r_meas   = [isnan(speeds[c]) ? NaN : speeds[c]*results[1][2][f] for f in 1:N_FRM]
    lines!(ax, r_theory, 1:N_FRM; color=:white, linewidth=1.5, linestyle=:dash)
    isnan(speeds[c]) || lines!(ax, r_meas, 1:N_FRM; color=:red, linewidth=1.2,
           linestyle=:dot, label=@sprintf("v=%.3f",speeds[c]))
    c==1 && axislegend(ax; position=:lt, labelsize=7, framevisible=false)
end
Colorbar(fig2[1,n_methods+1]; colormap=:viridis, limits=(0,xt_max),
         label="⟨n⟩", labelsize=9, ticklabelsize=8, width=14)
colgap!(fig2.layout, 8); rowgap!(fig2.layout, 4)
save(joinpath(outdir,"fig_mg_validation.png"), fig2; px_per_unit=3)
println("Saved → fig_mg_validation.png")
