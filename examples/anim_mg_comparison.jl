"""
anim_mg_comparison.jl

Side-by-side animation comparing multigrid methods for BranchingDeath RDME.
Four panels, same visual style as anim_branching_200x200.jl:
  plasma colormap, active ring orange tint, equil interior blue tint,
  quad-tree overlay (wave-front boundary + fine active grid + coarse equil grid).

Methods:
  1. Per-voxel (baseline)
  2. Pair exact — full 2-body joint CME, no restriction (exact Galerkin)
  3. ADI pair vcycle — intra-pair pre-smooth + total-count + multilevel prolong
  4. MG L=2 (total-count) — broken: restriction erases spatial gradient
"""

using DiscStochSim, CairoMakie, Printf, Statistics

const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4
const T_END=120.0; const N_FRM=90
const n_max=14
const v_theory = 2*sqrt(D*(k_b-k_d))

# ── Colourmap helpers (shared with anim_branching_200x200) ────────────────────
const _PLASMA = CairoMakie.to_colormap(:plasma)
_plasma(t) = let c=_PLASMA[clamp(floor(Int,t*(length(_PLASMA)-1))+1,1,length(_PLASMA))];
             RGBf(c.r,c.g,c.b) end

const C_EMPTY = RGBf(0.04, 0.02, 0.06)
const C_ACT   = RGBf(0.88, 0.48, 0.22)   # orange — active ring
const C_EQ    = RGBf(0.27, 0.45, 0.72)   # blue  — equil interior

function make_heatmap(mg, sg, mu_max)
    img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        v  = clamp(sqrt(mg[i,j] / max(mu_max, 1e-3)), 0.0, 1.0)
        bc = _plasma(v)
        if sg[i,j] == 1
            α = 0.35
            img[j,i] = RGBf((1-α)*bc.r+α*C_ACT.r, (1-α)*bc.g+α*C_ACT.g, (1-α)*bc.b+α*C_ACT.b)
        elseif sg[i,j] == 2
            α = 0.18
            img[j,i] = RGBf((1-α)*bc.r+α*C_EQ.r,  (1-α)*bc.g+α*C_EQ.g,  (1-α)*bc.b+α*C_EQ.b)
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
    function hadd!(xs,ys,i,j)
        push!(xs,j-0.5f0,j+0.5f0,NaN32); push!(ys,Float32(i)+0.5f0,Float32(i)+0.5f0,NaN32)
    end
    function vadd!(xs,ys,i,j)
        push!(xs,Float32(j)+0.5f0,Float32(j)+0.5f0,NaN32); push!(ys,i-0.5f0,i+0.5f0,NaN32)
    end
    for i in 1:K, j in 1:K
        s = sg[i,j]
        if i < K
            s2 = sg[i+1,j]
            if  s==1 && s2!=1;  hadd!(bx,by,i,j)
            elseif s==1 && s2==1; hadd!(fx,fy,i,j)
            end
        end
        if j < K
            s2 = sg[i,j+1]
            if  s==1 && s2!=1;  vadd!(bx,by,i,j)
            elseif s==1 && s2==1; vadd!(fx,fy,i,j)
            end
        end
        if s==2
            i%8==0 && i<K && sg[i+1,j]==2 && hadd!(cx,cy,i,j)
            j%8==0 && j<K && sg[i,j+1]==2 && vadd!(cx,cy,i,j)
        end
    end
    (Point2f.(bx,by), Point2f.(fx,fy), Point2f.(cx,cy))
end

# ── Simulation setup ─────────────────────────────────────────────────────────
ci, cj = K÷2+1, K÷2+1

methods = [
    ("Per-voxel (baseline)",     (use_multigrid=false, use_block_exact=false, use_block_vcycle=false, multigrid_levels=0)),
    ("Pair exact (Galerkin)",    (use_multigrid=false, use_block_exact=true,  use_block_vcycle=false, multigrid_levels=0)),
    ("ADI pair vcycle",          (use_multigrid=false, use_block_exact=false, use_block_vcycle=true,  multigrid_levels=0)),
    ("MG L=2 (total-count)",     (use_multigrid=true,  use_block_exact=false, use_block_vcycle=false, multigrid_levels=2)),
]
N_METHODS = length(methods)

function init_state(; kw...)
    model = BranchingDeathRDME(k_b, k_d, D, 1.0)
    s = SpatialFSP(model, K, K; n_max, ε_expand=0.02, ε_equil=0.15, kw...)
    for di in -4:4, dj in -4:4
        di^2+dj^2 <= 16 || continue
        1<=ci+di<=K && 1<=cj+dj<=K || continue
        set_ic!(s, CartesianIndex(ci+di,cj+dj), N0)
    end
    s
end

states = [init_state(; m[2]...) for m in methods]

t_frames = collect(range(0.0, T_END; length=N_FRM))
imgs   = [Vector{Matrix{RGBf}}(undef, N_FRM) for _ in 1:N_METHODS]
ovls   = [Vector{NTuple{3,Vector{Point2f}}}(undef, N_FRM) for _ in 1:N_METHODS]
front_r = [zeros(N_FRM) for _ in 1:N_METHODS]

mu_max = Ref(1.0)

function max_active_r(s)
    isempty(s.dists) && return 0.0
    maximum(sqrt((k[1]-ci)^2+(k[2]-cj)^2) for k in keys(s.dists))
end

# Initial frames
for m in 1:N_METHODS
    mg = mean_grid(states[m]); sg = status_grid(states[m])
    mu_max[] = max(mu_max[], maximum(mg; init=1e-3)*1.02)
    imgs[m][1]   = make_heatmap(mg, sg, mu_max[])
    ovls[m][1]   = make_overlay(sg)
    front_r[m][1] = max_active_r(states[m])
end

@printf("Simulating %d methods for T=%.0f (K=%d, dt=%.1f)...\n",
        N_METHODS, T_END, K, dt)

let fi = 2, t0 = time()
    for step_i in 1:round(Int, T_END/dt)
        for m in 1:N_METHODS
            step!(states[m], dt)
        end
        t = step_i * dt
        if fi <= N_FRM && t >= t_frames[fi]-1e-9
            for m in 1:N_METHODS
                mg = mean_grid(states[m]); sg = status_grid(states[m])
                mu_max[] = max(mu_max[], maximum(mg; init=1e-3)*1.02)
                imgs[m][fi]    = make_heatmap(mg, sg, mu_max[])
                ovls[m][fi]    = make_overlay(sg)
                front_r[m][fi] = max_active_r(states[m])
            end
            if fi % 10 == 0
                @printf("  frame %3d / %d  t=%5.1f  elapsed=%.0fs\n",
                        fi, N_FRM, t, time()-t0); flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
    @printf("Simulation done in %.1fs\n\n", time()-t0)
end

# Fit wave speeds
function fit_speed(fr)
    mask = (t_frames .> 15.0) .& (fr .< 80.0) .& (fr .> 0)
    ts = t_frames[mask]; rs = fr[mask]
    length(ts) < 4 && return NaN
    tm = mean(ts); rm = mean(rs)
    num = sum((ts.-tm).*(rs.-rm)); den = sum((ts.-tm).^2)
    den > 0 ? num/den : NaN
end

speeds = [fit_speed(front_r[m]) for m in 1:N_METHODS]
for (m, (lbl,_)) in enumerate(methods)
    @printf("  %-28s  v=%.3f  err=%.1f%%\n", lbl, speeds[m],
            abs(speeds[m]-v_theory)/v_theory*100)
end
@printf("  %-28s  v=%.3f  (KPP theory)\n\n", "", v_theory)

# ── Build animation ───────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=10,
    Axis=(spinewidth=0.6, xgridvisible=false, ygridvisible=false)))

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)

ax_kw = (aspect=DataAspect(), yreversed=true,
          xaxisposition=:top, xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false)

fig = Figure(size=(N_METHODS*280 + 60, 380))

speed_strs = [@sprintf("%s\nv=%.3f (err %.0f%%)",
                        methods[m][1], speeds[m],
                        abs(speeds[m]-v_theory)/v_theory*100)
              for m in 1:N_METHODS]
for m in 1:N_METHODS
    Label(fig[0,m], speed_strs[m]; fontsize=9, halign=:center)
end

axs     = [Axis(fig[1,m]; ax_kw...) for m in 1:N_METHODS]
img_obs = [Observable(imgs[m][1]) for m in 1:N_METHODS]
bnd_obs = [Observable(ovls[m][1][1]) for m in 1:N_METHODS]
fin_obs = [Observable(ovls[m][1][2]) for m in 1:N_METHODS]
crs_obs = [Observable(ovls[m][1][3]) for m in 1:N_METHODS]

for m in 1:N_METHODS
    image!(axs[m], 0.5..K+0.5, 0.5..K+0.5, img_obs[m])
    linesegments!(axs[m], crs_obs[m]; color=(:white,0.18), linewidth=0.4)
    linesegments!(axs[m], fin_obs[m]; color=(:white,0.40), linewidth=0.5)
    linesegments!(axs[m], bnd_obs[m]; color=(:white,0.88), linewidth=1.3)
    # Theory wavefront circle overlay
    θ = range(0, 2π, 200)
end

t_lbl = Label(fig[2,1:N_METHODS], "t = 0.0"; fontsize=10, color=:gray40)
Legend(fig[3,1:N_METHODS],
    [LineElement(color=(:white,0.9),linewidth=2),
     LineElement(color=(:white,0.45),linewidth=1),
     PolyElement(color=RGBf(C_ACT.r,C_ACT.g,C_ACT.b)),
     PolyElement(color=RGBf(C_EQ.r,C_EQ.g,C_EQ.b))],
    ["Wave-front boundary","Active fine grid","Active ring","Equil interior"],
    orientation=:horizontal, framevisible=false, labelsize=8, patchsize=(12,9))

colgap!(fig.layout, 6); rowgap!(fig.layout, 1, 4); rowgap!(fig.layout, 2, 6)

out = joinpath(outdir, "anim_mg_comparison.gif")
println("Rendering → $out"); flush(stdout)
@time record(fig, out, 1:N_FRM; framerate=12) do f
    for m in 1:N_METHODS
        img_obs[m][] = imgs[m][f]
        bnd_obs[m][] = ovls[m][f][1]
        fin_obs[m][] = ovls[m][f][2]
        crs_obs[m][] = ovls[m][f][3]
    end
    t_lbl.text[] = @sprintf("t = %.1f   (KPP theory: v=%.3f)", t_frames[f], v_theory)
end
println("Saved → $out")

# ── Static 4-panel at mid-simulation ──────────────────────────────────────────
mid_f = argmax(front_r[1])   # frame where per-voxel front is at max radius < 80
mid_f = clamp(mid_f, 1, N_FRM)
@printf("\nSaving static panel at frame %d (t=%.1f)\n", mid_f, t_frames[mid_f])

fig2 = Figure(size=(N_METHODS*280 + 60, 340))
Label(fig2[0,1:N_METHODS],
    @sprintf("BranchingDeath RDME  K=%d×%d  t=%.0f  |  KPP v=%.3f",
             K, K, t_frames[mid_f], v_theory);
    fontsize=11, font=:bold)

for m in 1:N_METHODS
    ax = Axis(fig2[1,m]; title=speed_strs[m], titlesize=9, ax_kw...)
    image!(ax, 0.5..K+0.5, 0.5..K+0.5, imgs[m][mid_f])
    bnd, fin, crs = ovls[m][mid_f]
    linesegments!(ax, crs; color=(:white,0.18), linewidth=0.4)
    linesegments!(ax, fin; color=(:white,0.40), linewidth=0.5)
    linesegments!(ax, bnd; color=(:white,0.88), linewidth=1.3)
    # KPP theory circle
    r_th = v_theory * t_frames[mid_f]
    θ = range(0, 2π, 200)
    lines!(ax, cj .+ r_th.*cos.(θ), ci .+ r_th.*sin.(θ);
           color=(:yellow,0.7), linewidth=1.2, linestyle=:dash)
end

colgap!(fig2.layout, 6); rowgap!(fig2.layout, 4)
out2 = joinpath(outdir, "fig_mg_comparison.png")
save(out2, fig2; px_per_unit=3)
println("Saved → $out2")
