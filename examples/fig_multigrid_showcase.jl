"""
fig_multigrid_showcase.jl

Generates four multigrid showcase figures from a single simulation run:

  GIF 1 — Correction field animation:
    Left:  plasma heatmap of ⟨n_k⟩ (wavefront)
    Right: correction magnitude ||p_new - p_old||₁ per voxel — shows WHERE
           the L2 quad level is doing work (peaks at wavefront ring, zero in interior)

  GIF 2 — Space-time xt animation (diagram builds up over time):
    Three panels: L0 (per-voxel), L1 (2×2 smooth), L2 (4×4 smooth) radial profiles
    x=radius, y=time, color=⟨n⟩ — each level sees the wavefront at a different scale

  PNG 1 — Correction field snapshot at peak-correction frame
  PNG 2 — Full space-time xt diagram (all frames visible)
"""

using DiscStochSim, CairoMakie, Printf, Statistics

# ── Simulation parameters ───────────────────────────────────────────────────
const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4; const T_END=200.0; const N_FRM=150
const n_max=14

model  = BranchingDeathRDME(k_b, k_d, D, 1.0)
v_wave = 2*sqrt(D*(k_b-k_d))
@printf("Multigrid showcase: v_wave≈%.2f  K=%d×%d  T=%.0f\n", v_wave, K, K, T_END)

state = SpatialFSP(model, K, K;
    n_max, ε_expand=0.02, ε_equil=0.15,
    use_multigrid=true, use_block_vcycle=false)

ci,cj = K÷2+1, K÷2+1
for di in -4:4, dj in -4:4
    di^2+dj^2<=16 || continue
    1<=ci+di<=K && 1<=cj+dj<=K || continue
    set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
end

# ── Colormaps ────────────────────────────────────────────────────────────────
const _PLASMA = CairoMakie.to_colormap(:plasma)
const _HOT    = CairoMakie.to_colormap(:hot)
const _VIRIDIS= CairoMakie.to_colormap(:viridis)

_cmap(t, cmap) = let c=cmap[clamp(floor(Int,t*(length(cmap)-1))+1,1,length(cmap))];
                 RGBf(c.r,c.g,c.b) end

# ── Helpers ──────────────────────────────────────────────────────────────────
const mu_max = Ref(1.0)
const corr_max = Ref(0.01)

function plasma_img(mg)
    cur = maximum(mg; init=1e-3)
    mu_max[] = max(mu_max[], cur*1.02)
    img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        img[j,i] = _cmap(clamp(sqrt(mg[i,j]/max(mu_max[],1e-3)),0,1), _PLASMA)
    end
    img
end

function corr_img(cg)
    cur = maximum(cg; init=1e-6)
    corr_max[] = max(corr_max[], cur)
    img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        t = clamp(cg[i,j]/max(corr_max[],1e-8), 0, 1)
        img[j,i] = _cmap(t, _HOT)
    end
    img
end

# Spatial box-smooth the mean grid (approximates what level-ℓ sees)
function box_smooth(mg, box)
    out = zeros(K, K)
    for i in 1:K, j in 1:K
        s=0.0; n=0
        for di in 0:box-1, dj in 0:box-1
            ii=i+di; jj=j+dj
            ii<=K && jj<=K && (s+=mg[ii,jj]; n+=1)
        end
        out[i,j] = n>0 ? s/n : 0.0
    end
    out
end

# Radial mean profile: average ⟨n⟩ in shells of radius r from centre
const MAX_R = min(ci-1, K-ci, cj-1, K-cj)
function radial_profile(mg)
    prof = zeros(MAX_R); cnt = zeros(Int, MAX_R)
    for i in 1:K, j in 1:K
        r = round(Int, sqrt((i-ci)^2+(j-cj)^2))
        1<=r<=MAX_R || continue
        prof[r] += mg[i,j]; cnt[r] += 1
    end
    [cnt[r]>0 ? prof[r]/cnt[r] : 0.0 for r in 1:MAX_R]
end

# ── Run simulation, collect frame data ───────────────────────────────────────
t_frames   = collect(range(0.0, T_END; length=N_FRM))
imgs_wave  = Vector{Matrix{RGBf}}(undef, N_FRM)    # plasma wavefront
imgs_corr  = Vector{Matrix{RGBf}}(undef, N_FRM)    # correction magnitude
profs_L0   = Matrix{Float64}(undef, N_FRM, MAX_R)  # radial profile L0
profs_L1   = Matrix{Float64}(undef, N_FRM, MAX_R)  # radial profile L1 (2×2 smooth)
profs_L2   = Matrix{Float64}(undef, N_FRM, MAX_R)  # radial profile L2 (4×4 smooth)
frame_times = zeros(N_FRM)

# IC frame
mg0 = mean_grid(state)
imgs_wave[1] = plasma_img(mg0)
imgs_corr[1] = zeros(RGBf, K, K)  # no correction at t=0
profs_L0[1,:] = radial_profile(mg0)
profs_L1[1,:] = radial_profile(box_smooth(mg0, 2))
profs_L2[1,:] = radial_profile(box_smooth(mg0, 4))
frame_times[1] = 0.0

println("Simulating $N_FRM frames…")
max_corr_per_frame = zeros(N_FRM)   # track peak correction value per frame
active_per_frame   = zeros(Int, N_FRM)  # track active voxel count per frame
let fi=2, corr_accum=zeros(K,K), n_accum=0
    for step_i in 1:round(Int, T_END/dt)
        p_before = Dict(k=>copy(v) for (k,v) in state.dists)
        step!(state, dt)
        for (k,p_new) in state.dists
            p_old = get(p_before, k, nothing)
            p_old === nothing && continue
            corr_accum[k[1],k[2]] += sum(abs, p_new .- p_old[1:length(p_new)]) / 2
        end
        n_accum += 1
        t = step_i*dt
        if fi<=N_FRM && t>=t_frames[fi]-1e-9
            mg = mean_grid(state)
            cg = corr_accum ./ max(n_accum,1)
            imgs_wave[fi] = plasma_img(mg)
            imgs_corr[fi] = corr_img(cg)
            max_corr_per_frame[fi] = maximum(cg)
            active_per_frame[fi]   = length(state.dists)
            profs_L0[fi,:] = radial_profile(mg)
            profs_L1[fi,:] = radial_profile(box_smooth(mg,2))
            profs_L2[fi,:] = radial_profile(box_smooth(mg,4))
            frame_times[fi] = t
            act=length(state.dists); eq=length(state.equil_dists)
            fi%20==0 && @printf("  frame %3d  t=%5.1f  act=%5d  eq=%5d\n",fi,t,act,eq)
            fi+=1
            fill!(corr_accum, 0.0); n_accum=0
        end
        fi>N_FRM && break
    end
end
println("Simulation done.")

# ── Global normalisation for xt diagrams ─────────────────────────────────────
xt_max = max(maximum(profs_L0), maximum(profs_L1), maximum(profs_L2), 1e-3)
r_axis = 1:MAX_R

# ── GIF 1: Correction field animation ────────────────────────────────────────
outdir = joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)

println("Rendering GIF 1 — correction field…")
fig1 = Figure(size=(900,460))
Label(fig1[0,1], "Wavefront  ⟨n_k⟩", fontsize=11, font=:bold)
Label(fig1[0,2], "Correction magnitude  ‖Δp_k‖₁  (L2 quad → L0 voxel)", fontsize=11, font=:bold)
ax_w = Axis(fig1[1,1]; aspect=DataAspect(), yreversed=true)
ax_c = Axis(fig1[1,2]; aspect=DataAspect(), yreversed=true)
hidedecorations!(ax_w); hidedecorations!(ax_c)
obs_w = Observable(imgs_wave[1]); obs_c = Observable(imgs_corr[1])
image!(ax_w, 0.5..K+0.5, 0.5..K+0.5, obs_w)
image!(ax_c, 0.5..K+0.5, 0.5..K+0.5, obs_c)
Colorbar(fig1[1,3]; colormap=:hot, limits=(0,1), label="‖Δp‖₁", labelsize=9,
         ticklabelsize=8, width=12, ticks=([0,1],["0","max"]))
t_lbl = Label(fig1[2,1:2], "t = 0.0", fontsize=10, color=:gray40)
colgap!(fig1.layout, 6); rowgap!(fig1.layout, 4)

@time record(fig1, joinpath(outdir,"anim_mg_correction.gif"), 1:N_FRM; framerate=12) do f
    obs_w[] = imgs_wave[f]; obs_c[] = imgs_corr[f]
    t_lbl.text[] = @sprintf("t = %.1f   active=%d", frame_times[f],
                             f<N_FRM ? length(state.dists) : 0)
end
println("Saved → anim_mg_correction.gif")

# ── GIF 2: Space-time xt animation (diagram builds up) ───────────────────────
println("Rendering GIF 2 — space-time xt…")
fig2 = Figure(size=(1200,500))
titles = ["L0  (per-voxel)", "L1  (2×2 smooth)", "L2  (4×4 smooth)"]
axs = [Axis(fig2[1,c]; title=titles[c], titlesize=11,
            xlabel="radius (voxels)", ylabel= c==1 ? "time" : "",
            xlabelsize=9, ylabelsize=9, xticklabelsize=7, yticklabelsize=7)
       for c in 1:3]

# Build full matrices upfront, mask non-yet-seen rows to black
profs_all = [profs_L0, profs_L1, profs_L2]
# Observable stores MAX_R × N_FRM (transposed) so heatmap(r_axis, frame_axis, data) works
hm_obs = [Observable(zeros(MAX_R, N_FRM)) for _ in 1:3]
for c in 1:3
    heatmap!(axs[c], r_axis, 1:N_FRM, hm_obs[c];
             colormap=:viridis, colorrange=(0, xt_max))
    axs[c].yticks = (range(1,N_FRM,length=6), [@sprintf("%.0f",frame_times[round(Int,f)])
                     for f in range(1,N_FRM,length=6)])
end
Colorbar(fig2[1,4]; colormap=:viridis, limits=(0,xt_max), label="⟨n⟩",
         labelsize=9, ticklabelsize=8, width=14)
Label(fig2[0,1:3], "Multigrid space-time: wavefront at each spatial scale  (radius vs time)",
      fontsize=12, font=:bold)
colgap!(fig2.layout, 8); rowgap!(fig2.layout, 4)

@time record(fig2, joinpath(outdir,"anim_mg_spacetime.gif"), 1:N_FRM; framerate=12) do f
    for c in 1:3
        M = zeros(MAX_R, N_FRM)
        M[:,1:f] = profs_all[c][1:f,:]'   # reveal columns up to current frame
        hm_obs[c][] = M
    end
end
println("Saved → anim_mg_spacetime.gif")

# ── PNG 1: Correction field snapshot at frame with peak correction ────────────
# Use frame with peak active voxels — that's when the wavefront ring is largest
# and the L2 correction is most physically interesting
peak_frame = argmax(active_per_frame)
@printf("Peak active frame: %d  t=%.1f  active=%d  max_corr=%.4f\n",
        peak_frame, frame_times[peak_frame],
        active_per_frame[peak_frame], max_corr_per_frame[peak_frame])

fig3 = Figure(size=(920,500))
Label(fig3[0,1:2],
    @sprintf("Multigrid correction field  ·  t = %.1f  ·  %d active voxels",
             frame_times[peak_frame], active_per_frame[peak_frame]),
    fontsize=11, font=:bold)
ax3w = Axis(fig3[1,1]; aspect=DataAspect(), yreversed=true,
            title="Wavefront  ⟨n_k⟩", titlesize=10)
ax3c = Axis(fig3[1,2]; aspect=DataAspect(), yreversed=true,
            title="Correction magnitude  ‖Δp_k‖₁  (L2→L0)", titlesize=10)
hidedecorations!(ax3w); hidedecorations!(ax3c)
image!(ax3w, 0.5..K+0.5, 0.5..K+0.5, imgs_wave[peak_frame])
image!(ax3c, 0.5..K+0.5, 0.5..K+0.5, imgs_corr[peak_frame])
Colorbar(fig3[1,3]; colormap=:hot, limits=(0,1), label="‖Δp‖₁",
         labelsize=9, ticklabelsize=8, width=14)
colsize!(fig3.layout, 1, Aspect(1, 1.0))
colsize!(fig3.layout, 2, Aspect(1, 1.0))
colgap!(fig3.layout, 8); rowgap!(fig3.layout, 4)
save(joinpath(outdir,"fig_mg_correction.png"), fig3; px_per_unit=3)
println("Saved → fig_mg_correction.png")

# ── PNG 2: Full space-time xt diagram ─────────────────────────────────────────
fig4 = Figure(size=(1300,480))
Label(fig4[0,1:3], "Multigrid space-time diagram: wavefront propagation at each spatial scale",
      fontsize=12, font=:bold)

t_ticks_vals = range(1,N_FRM,length=7)
t_ticks = (collect(t_ticks_vals),
           [@sprintf("%.0f", frame_times[round(Int,f)]) for f in t_ticks_vals])

for (c,(prof,ttl)) in enumerate(zip([profs_L0,profs_L1,profs_L2], titles))
    ax = Axis(fig4[1,c]; title=ttl, titlesize=11,
              xlabel="radius (voxels)", ylabel=c==1 ? "time" : "",
              xlabelsize=9, ylabelsize=9, xticklabelsize=7, yticklabelsize=7,
              yticks = t_ticks)
    heatmap!(ax, r_axis, 1:N_FRM, prof'; colormap=:viridis, colorrange=(0,xt_max))
    # Overlay wavefront velocity line
    v_r = v_wave   # voxels per time unit
    t_line = [frame_times[f] for f in 1:N_FRM]
    r_line = v_r .* t_line
    valid = r_line .≤ MAX_R
    lines!(ax, r_line[valid], (1:N_FRM)[valid]; color=:white, linewidth=1.5,
           linestyle=:dash, label="v_wave")
    c==1 && axislegend(ax; position=:lt, labelsize=8, framevisible=false)
end
Colorbar(fig4[1,4]; colormap=:viridis, limits=(0,xt_max), label="⟨n⟩",
         labelsize=9, ticklabelsize=8, width=14)
colgap!(fig4.layout, 8); rowgap!(fig4.layout, 4)
save(joinpath(outdir,"fig_mg_spacetime.png"), fig4; px_per_unit=3)
println("Saved → fig_mg_spacetime.png")
println("\nAll done.")
