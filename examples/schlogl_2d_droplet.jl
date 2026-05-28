"""
schlogl_2d_droplet.jl

Schlögl bistable RDME on 60×60 grid — multigrid FSP (L=3).
Shows critical-droplet dynamics: three initial radii, all growing (c3=21, high state preferred).
Output: paper/figures/fig_schlogl_2d_droplet.png
        paper/figures/anim_schlogl_2d_droplet.gif
"""

using DiscStochSim, CairoMakie, Printf, Statistics

# ── Model ──────────────────────────────────────────────────────────────────────
const D   = 2.0
const K   = 60

model = SchloglModel1D(D, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Schlögl fixed pts: n_lo=%d  n_un=%d  n_hi=%d  K=%d×%d  D=%.1f\n",
        n_lo, n_un, n_hi, K, K, D)

# ── State constructor with multigrid L=3 ──────────────────────────────────────
function make_droplet_state(R_vox::Float64; interface_width=3)
    state = SpatialFSP(model, K, K;
                       n_max       = 220,
                       n_max_eq    = 220,
                       ε_expand    = 0.15,           # flux-based activation for equil→active
                       ε_equil     = 0.12,
                       use_block_vcycle = false)
    ci, cj = K÷2+1, K÷2+1
    nm_eq  = state.n_max_eq
    for i in 1:K, j in 1:K
        r   = sqrt((i-ci)^2 + (j-cj)^2)
        idx = CartesianIndex(i,j)
        n_ic = r <= R_vox ? n_hi : n_lo
        if abs(r - R_vox) <= interface_width
            p = zeros(state.n_max + 1); p[n_ic+1] = 1.0
            state.dists[idx] = p
        else
            p = zeros(nm_eq + 1); p[min(n_ic+1, nm_eq+1)] = 1.0
            state.equil_dists[idx] = p
        end
    end
    state
end

# ── Colour: normalised mean (mean - n_lo)/(n_hi - n_lo) → RdBu ────────────────
# Using normalised mean rather than P(hi) because equil voxels are stored as
# delta functions (P(hi)=0 or 1 exactly), while active wavefront voxels sit at
# mean≈70 (≈25% of the way from n_lo to n_hi) — clearly visible as orange
# between the red exterior (equil-lo) and blue interior (equil-hi).
const _CMAP = CairoMakie.to_colormap(:RdBu)
_cmap(t) = let c = _CMAP[clamp(floor(Int, t*(length(_CMAP)-1))+1, 1, length(_CMAP))];
           RGBf(c.r, c.g, c.b) end

voxel_color_val(p) = clamp((sum((i-1)*p[i] for i in 1:length(p)) - n_lo) / (n_hi - n_lo),
                           0.0, 1.0)

function make_heatmap(state)
    img = zeros(RGBf, K, K)
    for (idx, p) in state.dists
        img[idx[2], idx[1]] = _cmap(voxel_color_val(p))
    end
    for (idx, p) in state.equil_dists
        img[idx[2], idx[1]] = _cmap(voxel_color_val(p))
    end
    img
end

# ── Run one droplet ─────────────────────────────────────────────────────────────
function run_droplet(R_vox; T_END=80.0, dt=0.5, N_FRM=40)
    state   = make_droplet_state(R_vox)
    t_frm   = collect(range(0.0, T_END; length=N_FRM))
    imgs    = Vector{Matrix{RGBf}}(undef, N_FRM)
    hi_frac = zeros(N_FRM)

    p0 = zeros(1); p0[1] = 1.0   # default for missing voxels (shouldn't occur for Schlogl)
    get_p(state, idx) = get(state.dists, idx, get(state.equil_dists, idx, p0))
    hi_frac_now(state) = mean(voxel_color_val(get_p(state, CartesianIndex(i,j)))
                              for i in 1:K, j in 1:K)

    imgs[1]    = make_heatmap(state)
    hi_frac[1] = hi_frac_now(state)

    fi = 2
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step_i * dt
        if fi ≤ N_FRM && t ≥ t_frm[fi] - 1e-9
            imgs[fi]    = make_heatmap(state)
            hi_frac[fi] = hi_frac_now(state)
            if fi % 10 == 0
                @printf("  R=%.0f  t=%.0f  hi=%.3f  active=%d  equil=%d\n",
                        R_vox, t, hi_frac[fi], n_active(state), n_equilibrated(state))
                flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
    imgs, hi_frac, t_frm
end

# ── Run all three radii ─────────────────────────────────────────────────────────
radii  = [2.0, 5.0, 12.0]
labels = ["R=2", "R=5", "R=12"]
@printf("\nRunning %d radii…\n\n", length(radii))
results = [begin @printf("R=%.0f…\n", R); flush(stdout); run_droplet(R) end for R in radii]

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false)))

fig      = Figure(size=(1060, 760))
t_show   = [1, 20, 40]
t_labs   = ["t = 0 (IC)", "t = 40", "t = 80"]
ax_kw    = (aspect=DataAspect(), yreversed=true)
row_lbls = ["R = 2", "R = 5", "R = 12"]

for (row, (R, lbl, (imgs, hi_frac, t_frm))) in enumerate(zip(radii, row_lbls, results))
    for (col, (fi, tl)) in enumerate(zip(t_show, t_labs))
        title_str = col == 2 ? lbl : (row == 1 ? tl : "")
        ax = Axis(fig[row, col]; title=title_str, titlesize=10, ax_kw...)
        image!(ax, 0.5..K+0.5, 0.5..K+0.5, imgs[min(fi, length(imgs))])
        if fi == 1
            θ = LinRange(0, 2π, 200)
            lines!(ax, K÷2+1 .+ R.*cos.(θ), K÷2+1 .+ R.*sin.(θ);
                   color=:white, linewidth=1.5, linestyle=:dash)
        end
        hf = hi_frac[min(fi, length(hi_frac))]
        text!(ax, 2, K-2; text="hi: $(round(Int, 100hf))%",
              fontsize=8, color=:white, align=(:left,:bottom))
    end
end

# Time-series panel
ax_ts = Axis(fig[1:3, 4];
    xlabel="time", ylabel="P(hi) — fraction of domain in high state",
    title="Droplet expansion (adaptive spatial FSP)",
    titlesize=10, spinewidth=0.7, xgridvisible=true, ygridvisible=true,
    xticklabelsvisible=true, yticklabelsvisible=true,
    xticksvisible=true, yticksvisible=true)
clrs = [:royalblue, :darkorange, :crimson]
for (R, lbl, (_, hi_frac, t_frm), c) in zip(radii, row_lbls, results, clrs)
    lines!(ax_ts, t_frm[1:length(hi_frac)], hi_frac; color=c, linewidth=2.2, label=lbl)
end
hlines!(ax_ts, [0.0, 1.0]; color=(:gray,0.4), linestyle=:dash, linewidth=1)
axislegend(ax_ts; position=:rc, labelsize=9, framevisible=false)

Colorbar(fig[end+1, 1:3];
         colormap=:RdBu, limits=(0.0, 1.0),
         label="(mean - n_lo) / (n_hi - n_lo)   (red = lo basin,  orange = wavefront,  blue = hi basin)",
         vertical=false, labelsize=9, ticklabelsize=8,
         ticks=([0.0, (n_un-n_lo)/(n_hi-n_lo), 1.0],
                ["n_lo=$n_lo", "n_un=$n_un", "n_hi=$n_hi"]))

rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 5)
colgap!(fig.layout, 1, 5); colgap!(fig.layout, 2, 5); colgap!(fig.layout, 3, 14)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out_png = joinpath(outdir, "fig_schlogl_2d_droplet.png")
save(out_png, fig; px_per_unit=3)
@printf("\nSaved → %s\n", out_png)

# ── Animation: R=5 droplet ──────────────────────────────────────────────────────
imgs5, _, t_frm5 = results[2]
fig2  = Figure(size=(440, 460))
ax2   = Axis(fig2[1,1]; title="Schlögl droplet  R=5  (adaptive FSP)", ax_kw..., titlesize=10)
img_o = Observable(imgs5[1])
image!(ax2, 0.5..K+0.5, 0.5..K+0.5, img_o)
t_lbl = Label(fig2[0,1], "t = 0"; fontsize=11)
Colorbar(fig2[1,2]; colormap=:RdBu, limits=(0.0,1.0),
         ticks=([0.0,0.5,1.0],["0","½","1"]),
         width=14, ticklabelsize=7, labelsize=7)
colsize!(fig2.layout, 1, Aspect(1, 1.0))
colsize!(fig2.layout, 2, Fixed(30))
out_gif = joinpath(outdir, "anim_schlogl_2d_droplet.gif")
@printf("Rendering animation → %s\n", out_gif)
@time record(fig2, out_gif, 1:length(imgs5); framerate=10) do f
    img_o[]     = imgs5[f]
    t_lbl.text[] = "t = $(round(t_frm5[f], digits=1))"
end
@printf("Saved → %s\n", out_gif)
