using DiscStochSim, CairoMakie, Printf

# ── Model: Schlögl bistable RDME ───────────────────────────────────────────────
# Parameters: standard c3=21 (high state preferred, fast dynamics) with D=2.
# The P(hi) colormap makes the visualization physically honest:
#   - Equil-lo voxels → P(hi)=0 → SOLID RED   (no scatter)
#   - Equil-hi voxels → P(hi)=1 → SOLID BLUE  (no scatter)
#   - Active/transitioning voxels → P(hi)∈(0,1) → WHITE front
# Fast simulation (~40 min), clean figure due to P(hi) coloring.
const D   = 2.0   # diffusion
const K   = 60    # 60×60 grid

model = SchloglModel1D(D, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Schlögl fixed points: n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)
@printf("Grid: %d×%d  D=%.2f\n\n", K, K, D)

# ── Helper: initialise a circular droplet IC ───────────────────────────────────
function make_droplet_state(R_vox::Float64; n_max=220, ε_expand=0.15, ε_equil=0.12,
                             interface_width=3)
    state = SpatialFSP(model, K, K;
                       n_max, ε_expand, ε_equil)
    ci, cj = K÷2+1, K÷2+1
    n_max_eq = state.n_max_eq

    # Key optimization: only INTERFACE voxels are active (full CME).
    # Interior and exterior go directly to equil — huge speedup for 100×100 grid.
    # Interface = voxels within interface_width of the droplet boundary.
    for i in 1:K, j in 1:K
        r = sqrt((i-ci)^2 + (j-cj)^2)
        idx = CartesianIndex(i,j)
        if abs(r - R_vox) <= interface_width
            # Interface: put in active (full CME), delta IC at n_hi or n_lo
            n_ic = r <= R_vox ? n_hi : n_lo
            p = zeros(n_max + 1); p[n_ic + 1] = 1.0
            state.dists[idx] = p
        else
            # Interior/exterior: put directly in equil with delta IC
            n_ic = r <= R_vox ? n_hi : n_lo
            p = zeros(n_max_eq + 1); p[min(n_ic+1, n_max_eq+1)] = 1.0
            state.equil_dists[idx] = p
        end
    end
    state
end

# ── Colour helpers ─────────────────────────────────────────────────────────────
# Colour by P(hi) = probability mass in the high basin (n > n_un).
# This gives a physically clean picture:
#   P(hi)≈1 → deep blue   (equil-hi, firmly in high state)
#   P(hi)≈0 → deep red    (equil-lo, firmly in low state)
#   P(hi)∈(0,1) → white   (active transition zone / stochastic front)
# Much cleaner than colouring by mean, which makes transitioning voxels look
# "almost blue" and appear as scattered dots.
const _CMAP = CairoMakie.to_colormap(:RdBu)
_cmap(t) = let c=_CMAP[clamp(floor(Int,t*(length(_CMAP)-1))+1,1,length(_CMAP))];
           RGBf(c.r,c.g,c.b) end

function p_hi(p::Vector{Float64}, n_un_idx::Int)
    n_un_idx < length(p) ? sum(view(p, n_un_idx+1:length(p))) : 0.0
end

function make_heatmap(state)
    n_un_idx = n_un + 1   # 1-indexed for n = n_un
    img = zeros(RGBf, K, K)
    mg  = zeros(K, K)
    for (idx, p) in state.dists
        ph = p_hi(p, n_un_idx)
        img[idx[2], idx[1]] = _cmap(clamp(ph, 0.0, 1.0))
        mg[idx[1], idx[2]]  = voxel_mean(p)
    end
    for (idx, p) in state.equil_dists
        ph = p_hi(p, n_un_idx)
        img[idx[2], idx[1]] = _cmap(clamp(ph, 0.0, 1.0))
        mg[idx[1], idx[2]]  = voxel_mean(p)
    end
    img, mg
end

# ── Simulate one droplet and return history ─────────────────────────────────────
function run_droplet(R_vox; T_END=80.0, dt=0.5, N_FRM=40)
    state = make_droplet_state(R_vox)

    t_frames = collect(range(0.0, T_END; length=N_FRM))
    imgs     = Vector{Matrix{RGBf}}(undef, N_FRM)
    hi_frac  = zeros(N_FRM)   # fraction of voxels in high state

    # IC snapshot
    img0, mg0 = make_heatmap(state)
    imgs[1] = img0
    hi_frac[1] = count(img0[j,i].b > img0[j,i].r for i in 1:K, j in 1:K) / K^2

    @printf("  R=%.1f: frame 1 hi_frac=%.3f\n", R_vox, hi_frac[1])

    fi = 2
    for step in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step*dt
        if fi ≤ N_FRM && t ≥ t_frames[fi]-1e-9
            img, mg = make_heatmap(state)
            imgs[fi] = img
            hi_frac[fi] = count(img[j,i].b > img[j,i].r for i in 1:K, j in 1:K) / K^2
            if fi%5==0
                @printf("  R=%.1f  t=%.0f  hi_frac=%.3f  active=%d\n",
                        R_vox, t, hi_frac[fi], n_active(state))
                flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
    imgs, hi_frac, t_frames
end

# ── Run three radii: subcritical / near-critical / supercritical ───────────────
# D=2, c3=21: all three radii grow (high state preferred), but at different rates
radii = [2.0, 5.0, 12.0]
labels = ["R=2 (slow growth)", "R=5 (medium growth)", "R=12 (fast growth)"]

results = []
for R in radii
    @printf("Running R=%.1f...\n", R); flush(stdout)
    push!(results, run_droplet(R))
end

# ── Figure: 3×3 grid of snapshots (3 radii × early/mid/final) ─────────────────
set_theme!(Theme(fontsize=10,
    Axis=(spinewidth=0.6, xgridvisible=false, ygridvisible=false,
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false)))

fig = Figure(size=(1000, 850))

t_show = [1, 20, 40]   # frame indices: initial, mid, final
t_labs = ["t=0 (IC)", "t=40", "t=80 (final)"]

ax_kw = (aspect=DataAspect(), yreversed=true)

for (row, (R, label, (imgs, hi_frac, t_frames))) in enumerate(zip(radii, labels, results))
    for (col, (fi, tlab)) in enumerate(zip(t_show, t_labs))
        ax = Axis(fig[row, col]; title=col==2 ? label : tlab,
                  titlesize=9, ax_kw...)
        image!(ax, 0.5..K+0.5, 0.5..K+0.5, imgs[min(fi,length(imgs))])
        # Show droplet radius circle at t=0
        if fi == 1
            θ = LinRange(0, 2π, 200)
            lines!(ax, K÷2+1 .+ R.*cos.(θ), K÷2+1 .+ R.*sin.(θ);
                   color=:white, linewidth=1.5, linestyle=:dash)
        end
        # Hi-fraction annotation
        hf = hi_frac[min(fi,length(hi_frac))]
        text!(ax, 2, K-2; text="hi: $(round(100hf,digits=0))%",
              fontsize=7, color=:white, align=(:left,:bottom))
    end
end

# Right column: hi-fraction time series
ax_ts = Axis(fig[1:3, 4];
    xlabel="time", ylabel="fraction in n_hi state",
    title="Droplet fate: fraction of voxels at n_hi",
    titlesize=9, spinewidth=0.6)
cols_r = [:royalblue, :darkorange, :crimson]
for (R, label, (imgs, hi_frac, t_frames), col) in zip(radii, labels, results, cols_r)
    lines!(ax_ts, t_frames[1:length(hi_frac)], hi_frac;
           color=col, linewidth=2, label=label)
end
hlines!(ax_ts, [0.0, 1.0]; color=(:gray,0.4), linestyle=:dash, linewidth=1)
axislegend(ax_ts; position=:rc, labelsize=8, framevisible=false)

# Colorbar
Colorbar(fig[end+1, 1:3];
         colormap=:RdBu, limits=(0.0, 1.0),
         label="P(hi) = probability in high basin  (red=0, white=½, blue=1)",
         vertical=false, labelsize=9, ticklabelsize=8,
         ticks=([0.0, 0.5, 1.0], ["0 (lo)", "0.5 (front)", "1 (hi)"]))

rowgap!(fig.layout, 1, 6); rowgap!(fig.layout, 2, 6)
colgap!(fig.layout, 1, 6); colgap!(fig.layout, 2, 6); colgap!(fig.layout, 3, 12)

outdir = joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)
out = joinpath(outdir,"fig_schlogl_2d_droplet.png")
save(out, fig; px_per_unit=3)
println("\nSaved → $out")

# Also save animation for the near-critical case (R=6, index 2)
imgs7, hi_frac7, t_frames7 = results[2]
fig2 = Figure(size=(430,440))
ax2  = Axis(fig2[1,1]; title="Schlögl 2D droplet  R=5",
            ax_kw..., titlesize=10)
img_obs = Observable(imgs7[1])
image!(ax2, 0.5..K+0.5, 0.5..K+0.5, img_obs)
t_lab2 = Label(fig2[0,1],"t=0"; fontsize=11)
cb2 = Colorbar(fig2[1,2]; colormap=:RdBu, limits=(0.0,1.0),
               ticks=([0.0,0.5,1.0],["0","½","1"]),
               width=14, ticklabelsize=7, labelsize=7)
# Pin colorbar column to its natural width so the heatmap fills the rest
colsize!(fig2.layout, 1, Aspect(1, 1.0))  # heatmap column: square aspect
colsize!(fig2.layout, 2, Fixed(30))        # colorbar column: fixed 30px
out2 = joinpath(outdir,"anim_schlogl_2d_droplet.gif")
println("Rendering animation → $out2")
@time record(fig2, out2, 1:length(imgs7); framerate=10) do f
    img_obs[] = imgs7[f]
    t_lab2.text[] = "t = $(round(t_frames7[f],digits=1))"
end
println("Saved → $out2")
