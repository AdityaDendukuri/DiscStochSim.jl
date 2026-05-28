"""
fig_ssa_comparison.jl

Validate FSP methods against truncated SSA ground truth for BranchingDeath RDME.

Truncated SSA: birth rate capped at n_max (same truncation as FSP), so the
comparison is apples-to-apples.  The KPP wave speed is the *deterministic PDE*
limit — not the right validation target.  SSA is.

Observable: wavefront radius r(t) = max radius of any voxel with n > 0.
            Wave speed: linear fit of r(t) in [t_fit_lo, t_fit_hi].

Parameters: K=80, T=60, n_max=14, N_TRAJ=20 trajectories.
  Wave radius reaches ~36 voxels at T=60 < K/2-2=38 (no boundary collision).
"""

using DiscStochSim, CairoMakie, Printf, Statistics, Random

const k_b    = 0.5;  const k_d = 0.1;  const D = 1.0
const K      = 80;   const N0  = 3;    const dt = 0.4
const T_END  = 60.0; const N_FRM = 80; const n_max = 14
const N_TRAJ = 20
const v_theory = 2 * sqrt(D * (k_b - k_d))

const ci = K÷2 + 1;  const cj = K÷2 + 1
const t_frames = collect(range(0.0, T_END; length=N_FRM))

# ── Truncated SSA ─────────────────────────────────────────────────────────────

function run_ssa(seed::Int)
    Random.seed!(seed)
    n = zeros(Int32, K, K)

    # Disk IC: same as FSP (radius 4, N0 particles)
    for di in -4:4, dj in -4:4
        di^2 + dj^2 <= 16 || continue
        1 <= ci+di <= K && 1 <= cj+dj <= K || continue
        n[ci+di, cj+dj] = N0
    end

    # Per-voxel propensity (birth capped at n_max)
    @inline function vrate(i, j)
        nv = n[i,j];  nv == 0 && return 0.0
        nb = Int(i>1) + Int(i<K) + Int(j>1) + Int(j<K)
        birth = nv < n_max ? k_b * nv : 0.0
        birth + (k_d + D * nb) * nv
    end

    rates  = [vrate(i,j) for i in 1:K, j in 1:K]
    R      = sum(rates)

    # Wavefront: max distance of any occupied voxel from centre
    function max_front()
        mx = 0.0
        @inbounds for i in 1:K, j in 1:K
            n[i,j] > 0 || continue
            mx = max(mx, sqrt(Float64(i-ci)^2 + Float64(j-cj)^2))
        end
        mx
    end

    front = zeros(N_FRM)
    front[1] = max_front()
    t = 0.0;  fi = 2

    while fi <= N_FRM && R > 1e-15
        # Next event time
        τ = -log(rand()) / R
        t += τ

        # Record all frames that fall before this event
        while fi <= N_FRM && t_frames[fi] <= t
            front[fi] = max_front();  fi += 1
        end
        fi > N_FRM && break

        # Select voxel proportional to rate (O(K²) linear scan)
        u   = rand() * R;  cum = 0.0
        ei  = ej = ci
        @inbounds for i in 1:K, j in 1:K
            cum += rates[i,j]
            if cum >= u;  ei = i;  ej = j;  break;  end
        end

        # Execute event at (ei, ej)
        nv = n[ei, ej];  nv == 0 && continue
        nb = Int(ei>1) + Int(ei<K) + Int(ej>1) + Int(ej<K)
        birth = nv < n_max ? k_b * nv : 0.0
        u2 = rand() * (birth + (k_d + D*nb) * nv)

        if u2 < birth
            n[ei, ej] += 1
        elseif u2 < birth + k_d * nv
            n[ei, ej] -= 1
        else
            # Diffusion: pick random neighbor
            dirs = Tuple{Int,Int}[]
            ei>1 && push!(dirs, (-1,0));  ei<K && push!(dirs, (1,0))
            ej>1 && push!(dirs, (0,-1));  ej<K && push!(dirs, (0,1))
            di, dj = rand(dirs)
            ni, nj = ei+di, ej+dj
            n[ei,ej] -= 1;  n[ni,nj] += 1
            R -= rates[ni,nj];  rates[ni,nj] = vrate(ni,nj);  R += rates[ni,nj]
        end

        # Update current voxel rate
        R -= rates[ei,ej];  rates[ei,ej] = vrate(ei,ej);  R += rates[ei,ej]
    end

    front
end

# ── Run SSA ───────────────────────────────────────────────────────────────────
@printf("Running %d SSA trajectories  (K=%d, T=%.0f, n_max=%d)...\n",
        N_TRAJ, K, T_END, n_max)
t0_ssa = time()
ssa_fronts = [run_ssa(s) for s in 1:N_TRAJ]
@printf("SSA done in %.1fs\n\n", time()-t0_ssa); flush(stdout)

ssa_mean = mean(ssa_fronts)
ssa_std  = std(ssa_fronts)

# ── FSP methods ───────────────────────────────────────────────────────────────
function run_fsp(; label="", kwargs...)
    @printf("FSP %-34s ...", label);  flush(stdout)
    model = BranchingDeathRDME(k_b, k_d, D, 1.0)
    state = SpatialFSP(model, K, K; n_max, ε_expand=0.02, ε_equil=0.15,
                       organic_activation=true, ε_organic=1e-6, kwargs...)

    for di in -4:4, dj in -4:4
        di^2+dj^2 <= 16 || continue
        1<=ci+di<=K && 1<=cj+dj<=K || continue
        set_ic!(state, CartesianIndex(ci+di, cj+dj), N0)
    end

    max_r(s) = isempty(s.dists) ? 0.0 :
               maximum(sqrt(Float64(k[1]-ci)^2 + Float64(k[2]-cj)^2)
                       for k in keys(s.dists))

    front = zeros(N_FRM);  front[1] = max_r(state)
    fi = 2;  t0 = time()
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t2 = step_i * dt
        if fi <= N_FRM && t2 >= t_frames[fi] - 1e-9
            front[fi] = max_r(state);  fi += 1
        end
        fi > N_FRM && break
    end
    @printf("  %.1fs\n", time()-t0);  flush(stdout)
    front
end

configs = [
    ("Per-voxel (baseline)",
     (use_block_vcycle=false, use_block_vcycle_4body=false, use_block_exact=false, use_multigrid=false)),
    ("Pair exact (Galerkin)",
     (use_block_vcycle=false, use_block_vcycle_4body=false, use_block_exact=true,  use_multigrid=false)),
    ("Multigrid L=1",
     (use_block_vcycle=false, use_block_vcycle_4body=false, use_block_exact=false, use_multigrid=true,
      multigrid_levels=1)),
    ("Block V-cycle 2×2",
     (use_block_vcycle=false, use_block_vcycle_4body=true,  use_block_exact=false, use_multigrid=false,
      block_vcycle_4body_max=3000, block_pre_frac=0.3)),
]
fsp_fronts = [run_fsp(; label=lbl, kw...) for (lbl, kw) in configs]

# ── Backbone: joint CME over all active voxels ────────────────────────────────
# IC: single center voxel (N0 particles) — IC-independent asymptotic wave speed
# Wavefront metric: E[max_r(t)] = Σ_s p[s]·max_{k: s[k]>0} r_k  (matches SSA mean max_r)
function run_backbone(; label="Backbone (joint CME)", ε_flux=1e-6, krylov_m=30)
    @printf("FSP %-34s ...", label);  flush(stdout)
    model = BranchingDeathRDME(k_b, k_d, D, 1.0)
    jf = JointRDMEFSP(model, K, K; n_max, ε_flux, krylov_m)

    set_ic!(jf, CartesianIndex(ci, cj), N0)   # single-voxel IC

    front = zeros(N_FRM);  front[1] = wavefront_radius_mean(jf, ci, cj)
    fi = 2;  t0 = time()
    n_steps = round(Int, T_END / dt)
    for step_i in 1:n_steps
        step!(jf, dt)
        t2 = step_i * dt
        while fi <= N_FRM && t_frames[fi] <= t2 + 1e-9
            front[fi] = wavefront_radius_mean(jf, ci, cj);  fi += 1
        end
        fi > N_FRM && break
        # progress line every ~5 time units
        step_i % round(Int, 5/dt) == 0 &&
            @printf("\r  [backbone t=%.1f  K_f=%d  |S|=%d]        ",
                    t2, n_fine(jf), n_joint_states(jf))
    end
    @printf("\r  %-34s  %.1fs  (K_f=%d  |S|=%d)\n",
            label, time()-t0, n_fine(jf), n_joint_states(jf))
    flush(stdout)
    front
end

backbone_front = run_backbone()

# ── Wave speed fitting ────────────────────────────────────────────────────────
const t_fit_lo = 10.0
const t_fit_hi = 50.0

function fit_speed(fr)
    mask = (t_frames .>= t_fit_lo) .& (fr .< K÷2 - 3) .& (fr .> 0) .& (t_frames .<= t_fit_hi)
    ts = t_frames[mask];  rs = fr[mask]
    length(ts) < 5 && return NaN
    tm = mean(ts);  rm = mean(rs)
    num = sum((ts.-tm).*(rs.-rm));  den = sum((ts.-tm).^2)
    den > 0 ? num/den : NaN
end

v_ssa_per_traj = [fit_speed(f) for f in ssa_fronts]
v_ssa_mean = mean(filter(!isnan, v_ssa_per_traj))
v_ssa_std  = std(filter(!isnan, v_ssa_per_traj))

println("\nWave speed (linear fit, r vs t):")
println("─"^62)
@printf("  %-34s  v=%.4f ± %.4f  (N=%d)\n",
        "SSA (truncated, ground truth)", v_ssa_mean, v_ssa_std, N_TRAJ)
for ((lbl,_), fr) in zip(configs, fsp_fronts)
    v = fit_speed(fr)
    δ = abs(v - v_ssa_mean) / v_ssa_mean * 100
    @printf("  %-34s  v=%.4f  err_SSA=%.1f%%\n", lbl, v, δ)
end
let v = fit_speed(backbone_front)
    δ = abs(v - v_ssa_mean) / v_ssa_mean * 100
    @printf("  %-34s  v=%.4f  err_SSA=%.1f%%\n", "Backbone (joint CME)", v, δ)
end
@printf("  %-34s  v=%.4f  (deterministic PDE limit)\n", "KPP theory", v_theory)

# ── Figure ────────────────────────────────────────────────────────────────────
outdir = joinpath(@__DIR__, "..", "paper", "figures");  mkpath(outdir)
cols   = [:steelblue, :darkorange, :seagreen, :firebrick, :mediumpurple]

fig = Figure(size=(760, 480))
ax  = Axis(fig[1,1];
    title = @sprintf("BranchingDeath RDME  K=%d×%d  n_max=%d  |  wavefront radius vs time",
                     K, K, n_max),
    xlabel = "time",  ylabel = "wavefront radius (voxels)",
    titlesize = 11, xlabelsize = 10, ylabelsize = 10)

# SSA ground truth band
band!(ax, t_frames, max.(0.0, ssa_mean .- ssa_std), ssa_mean .+ ssa_std;
      color = (:black, 0.13))
lines!(ax, t_frames, ssa_mean; color = :black, linewidth = 2.5,
       label = @sprintf("SSA ground truth (mean ± 1σ, N=%d)", N_TRAJ))

# FSP methods
for ((lbl,_), fr, col) in zip(configs, fsp_fronts, cols)
    v = fit_speed(fr)
    lines!(ax, t_frames, fr; color = col, linewidth = 2.0,
           label = @sprintf("%s  v=%.3f", lbl, v))
end

# Backbone
let v = fit_speed(backbone_front)
    lines!(ax, t_frames, backbone_front; color = :darkorchid, linewidth = 2.5,
           linestyle = :dash,
           label = @sprintf("Backbone (joint CME)  v=%.3f", v))
end

# KPP theory (dashed reference)
lines!(ax, t_frames, v_theory .* t_frames; color = :gray55, linewidth = 1.2,
       linestyle = :dot,
       label = @sprintf("KPP theory (v=%.3f, det. PDE limit)", v_theory))

xlims!(ax, 0, T_END);  ylims!(ax, 0, K÷2 + 2)
axislegend(ax; position = :lt, labelsize = 9, framevisible = true)

out = joinpath(outdir, "fig_ssa_comparison.png")
save(out, fig; px_per_unit = 3)
save(joinpath(outdir, "fig_ssa_comparison.pdf"), fig)
@printf("\nSaved → fig_ssa_comparison.{png,pdf}\n")
