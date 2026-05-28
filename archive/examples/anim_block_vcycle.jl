"""
anim_block_vcycle.jl

Animation of the 2×2 block V-cycle method for BranchingDeath RDME.

V-cycle structure per block:
  1. Pre-smooth:  per-voxel expv for τ = dt * 0.3
  2. Restrict:    top-3 states/voxel → product seed → rstep expand (depth=2)
  3. Coarse solve: full 4-body joint CME (all 4 diffusion directions)
  4. Prolongate:  multiplicative correction  p_new = p_pre * (marg_final/marg_init)

Visual style: plasma heatmap, orange active ring, blue equil interior,
quad-tree overlay (wavefront boundary + fine active + coarse equil grid).
"""

using DiscStochSim, CairoMakie, Printf, Statistics

const k_b = 0.5; const k_d = 0.1; const D = 1.0
const K = 200; const N0 = 3; const dt = 0.4
const T_END = 120.0; const N_FRM = 90
const n_max = 14
const v_theory = 2 * sqrt(D * (k_b - k_d))

# ── Colour helpers ────────────────────────────────────────────────────────────
const _PLASMA = CairoMakie.to_colormap(:plasma)
_plasma(t) = let c = _PLASMA[clamp(floor(Int, t * (length(_PLASMA)-1)) + 1,
                                    1, length(_PLASMA))];
             RGBf(c.r, c.g, c.b) end

const C_EMPTY = RGBf(0.04, 0.02, 0.06)
const C_ACT   = RGBf(0.88, 0.48, 0.22)
const C_EQ    = RGBf(0.27, 0.45, 0.72)

function make_heatmap(mg, sg, mu_max)
    img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        v  = clamp(sqrt(mg[i,j] / max(mu_max, 1e-3)), 0.0, 1.0)
        bc = _plasma(v)
        if sg[i,j] == 1
            α = 0.20          # lighter tint: let plasma gradient show through
            img[j,i] = RGBf((1-α)*bc.r + α*C_ACT.r,
                             (1-α)*bc.g + α*C_ACT.g,
                             (1-α)*bc.b + α*C_ACT.b)
        elseif sg[i,j] == 2
            α = 0.18
            img[j,i] = RGBf((1-α)*bc.r + α*C_EQ.r,
                             (1-α)*bc.g + α*C_EQ.g,
                             (1-α)*bc.b + α*C_EQ.b)
        else
            img[j,i] = bc
        end
    end
    img
end

function make_overlay(sg)
    bx = Float32[]; by = Float32[]
    fx = Float32[]; fy = Float32[]
    cx = Float32[]; cy = Float32[]
    hadd!(xs, ys, i, j) = (push!(xs, j-0.5f0, j+0.5f0, NaN32);
                             push!(ys, Float32(i)+0.5f0, Float32(i)+0.5f0, NaN32))
    vadd!(xs, ys, i, j) = (push!(xs, Float32(j)+0.5f0, Float32(j)+0.5f0, NaN32);
                             push!(ys, i-0.5f0, i+0.5f0, NaN32))
    for i in 1:K, j in 1:K
        s = sg[i,j]
        if i < K
            s2 = sg[i+1,j]
            s==1 && s2!=1 && hadd!(bx,by,i,j)
            s==1 && s2==1 && hadd!(fx,fy,i,j)
        end
        if j < K
            s2 = sg[i,j+1]
            s==1 && s2!=1 && vadd!(bx,by,i,j)
            s==1 && s2==1 && vadd!(fx,fy,i,j)
        end
        if s == 2
            i%8==0 && i<K && sg[i+1,j]==2 && hadd!(cx,cy,i,j)
            j%8==0 && j<K && sg[i,j+1]==2 && vadd!(cx,cy,i,j)
        end
    end
    (Point2f.(bx,by), Point2f.(fx,fy), Point2f.(cx,cy))
end

# ── Simulation ────────────────────────────────────────────────────────────────
ci, cj = K÷2+1, K÷2+1

model = BranchingDeathRDME(k_b, k_d, D, 1.0)
state = SpatialFSP(model, K, K;
    n_max, ε_expand=0.02, ε_equil=0.15,
    use_block_vcycle_4body = true,
    block_vcycle_4body_max = 3000,
    block_pre_frac         = 0.3,
    use_block_vcycle       = false)

for di in -4:4, dj in -4:4
    di^2 + dj^2 <= 16 || continue
    1 <= ci+di <= K && 1 <= cj+dj <= K || continue
    set_ic!(state, CartesianIndex(ci+di, cj+dj), N0)
end

t_frames = collect(range(0.0, T_END; length=N_FRM))
imgs     = Vector{Matrix{RGBf}}(undef, N_FRM)
ovls     = Vector{NTuple{3,Vector{Point2f}}}(undef, N_FRM)
front_r  = zeros(N_FRM)
mu_max   = Ref(1.0)

max_active_r(s) = isempty(s.dists) ? 0.0 :
    maximum(sqrt((k[1]-ci)^2 + (k[2]-cj)^2) for k in keys(s.dists))

mg0 = mean_grid(state); sg0 = status_grid(state)
mu_max[] = max(mu_max[], maximum(mg0; init=1e-3) * 1.02)
imgs[1]    = make_heatmap(mg0, sg0, mu_max[])
ovls[1]    = make_overlay(sg0)
front_r[1] = max_active_r(state)

@printf("Simulating block V-cycle  K=%d  T=%.0f  dt=%.2f  n_max=%d\n",
        K, T_END, dt, n_max); flush(stdout)

let fi = 2, t0 = time()
    for step_i in 1:round(Int, T_END / dt)
        step!(state, dt)
        t = step_i * dt
        if fi <= N_FRM && t >= t_frames[fi] - 1e-9
            mg = mean_grid(state); sg = status_grid(state)
            mu_max[] = max(mu_max[], maximum(mg; init=1e-3) * 1.02)
            imgs[fi]    = make_heatmap(mg, sg, mu_max[])
            ovls[fi]    = make_overlay(sg)
            front_r[fi] = max_active_r(state)
            if fi % 10 == 0
                @printf("  frame %3d/%d  t=%5.1f  elapsed=%.0fs\n",
                        fi, N_FRM, t, time()-t0); flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
    @printf("Done in %.1fs\n\n", time()-t0)
end

# Fit wave speed
function fit_speed(fr)
    mask = (t_frames .> 15.0) .& (fr .< 80.0) .& (fr .> 0)
    ts = t_frames[mask]; rs = fr[mask]
    length(ts) < 4 && return NaN
    tm = mean(ts); rm = mean(rs)
    num = sum((ts.-tm).*(rs.-rm)); den = sum((ts.-tm).^2)
    den > 0 ? num/den : NaN
end

v_meas = fit_speed(front_r)
@printf("Block V-cycle 2×2:  v=%.4f  err=%.1f%%\n", v_meas,
        abs(v_meas - v_theory)/v_theory*100)
@printf("KPP theory:         v=%.4f\n\n", v_theory)

# ── Build animation ───────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=10,
    Axis=(spinewidth=0.6, xgridvisible=false, ygridvisible=false)))

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)

ax_kw = (aspect=DataAspect(), yreversed=true,
          xaxisposition=:top, xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false)

title_str = @sprintf("Block V-cycle 2×2  v=%.3f  (KPP: %.3f, err %.0f%%)",
                     v_meas, v_theory, abs(v_meas-v_theory)/v_theory*100)

fig = Figure(size=(560, 440))
Label(fig[0, 1], title_str; fontsize=10, halign=:center)

ax      = Axis(fig[1,1]; ax_kw...)
img_obs = Observable(imgs[1])
bnd_obs = Observable(ovls[1][1])
fin_obs = Observable(ovls[1][2])
crs_obs = Observable(ovls[1][3])

image!(ax, 0.5..K+0.5, 0.5..K+0.5, img_obs)
linesegments!(ax, crs_obs; color=(:white,0.18), linewidth=0.4)
linesegments!(ax, fin_obs; color=(:white,0.12), linewidth=0.4)
linesegments!(ax, bnd_obs; color=(:white,0.90), linewidth=1.4)

t_lbl = Label(fig[2,1], "t = 0.0"; fontsize=10, color=:gray40)
Legend(fig[3,1],
    [LineElement(color=(:white,0.9), linewidth=2),
     LineElement(color=(:white,0.45), linewidth=1),
     PolyElement(color=RGBf(C_ACT.r, C_ACT.g, C_ACT.b)),
     PolyElement(color=RGBf(C_EQ.r,  C_EQ.g,  C_EQ.b))],
    ["Wave-front boundary", "Active fine grid", "Active ring", "Equil interior"],
    orientation=:horizontal, framevisible=false, labelsize=8, patchsize=(12,9))

rowgap!(fig.layout, 1, 4); rowgap!(fig.layout, 2, 6)

out = joinpath(outdir, "anim_block_vcycle.gif")
println("Rendering → $out"); flush(stdout)
@time record(fig, out, 1:N_FRM; framerate=12) do f
    img_obs[] = imgs[f]
    bnd_obs[] = ovls[f][1]
    fin_obs[] = ovls[f][2]
    crs_obs[] = ovls[f][3]
    t_lbl.text[] = @sprintf("t = %.1f", t_frames[f])
end
println("Saved → $out")

# ── Static snapshot at peak front ────────────────────────────────────────────
peak_f = argmax(front_r)
fig2   = Figure(size=(480, 460))
Label(fig2[0,1],
    @sprintf("Block V-cycle 2×2  K=%d  t=%.0f  |  v=%.3f  (KPP %.3f)",
             K, t_frames[peak_f], v_meas, v_theory);
    fontsize=10, font=:bold)
ax2 = Axis(fig2[1,1]; title="", ax_kw...)
image!(ax2, 0.5..K+0.5, 0.5..K+0.5, imgs[peak_f])
bnd, fin, crs = ovls[peak_f]
linesegments!(ax2, crs; color=(:white,0.18), linewidth=0.4)
linesegments!(ax2, fin; color=(:white,0.12), linewidth=0.4)
linesegments!(ax2, bnd; color=(:white,0.90), linewidth=1.4)
r_th = v_theory * t_frames[peak_f]
θ = range(0, 2π, 200)
lines!(ax2, cj .+ r_th.*cos.(θ), ci .+ r_th.*sin.(θ);
       color=(:yellow,0.7), linewidth=1.2, linestyle=:dash)
rowgap!(fig2.layout, 4)
out2 = joinpath(outdir, "fig_block_vcycle.png")
save(out2, fig2; px_per_unit=3)
println("Saved → $out2")
