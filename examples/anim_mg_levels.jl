"""
anim_mg_levels.jl

Side-by-side comparison animation: multigrid with 2, 3, and 4 spatial levels.
Shows how additional coarsening levels affect the wavefront dynamics.

L=2: pairs only   (L0 voxels → L1 pairs)
L=3: quads        (L0 → L1 pairs → L2 2×2 quads)       [current default]
L=4: 8-blocks     (L0 → L1 → L2 → L3 2×4 blocks)
"""

using DiscStochSim, CairoMakie, Printf

const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4; const T_END=200.0; const N_FRM=120
const n_max=14

const _PLASMA = CairoMakie.to_colormap(:plasma)
_plasma(t) = let c=_PLASMA[clamp(floor(Int,t*(length(_PLASMA)-1))+1,1,length(_PLASMA))];
             RGBf(c.r,c.g,c.b) end

function run_level(n_lev)
    model  = BranchingDeathRDME(k_b, k_d, D, 1.0)
    state  = SpatialFSP(model, K, K;
        n_max, ε_expand=0.02, ε_equil=0.15,
        use_multigrid=true, multigrid_levels=n_lev, use_block_vcycle=false)
    ci,cj  = K÷2+1, K÷2+1
    for di in -4:4, dj in -4:4
        di^2+dj^2 <= 16 || continue
        1<=ci+di<=K && 1<=cj+dj<=K || continue
        set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
    end

    mu_max = Ref(1.0)
    function frame_img(s)
        mg = mean_grid(s)
        cur = maximum(mg; init=1e-3)
        mu_max[] = max(mu_max[], cur*1.02)
        img = zeros(RGBf, K, K)
        for j in 1:K, i in 1:K
            v = clamp(sqrt(mg[i,j]/max(mu_max[],1e-3)), 0, 1)
            img[j,i] = _plasma(v)
        end
        img
    end

    t_frames = collect(range(0.0, T_END; length=N_FRM))
    imgs = Vector{Matrix{RGBf}}(undef, N_FRM)
    imgs[1] = frame_img(state)
    acts = zeros(Int, N_FRM)

    fi = 2
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step_i * dt
        if fi <= N_FRM && t >= t_frames[fi]-1e-9
            imgs[fi] = frame_img(state)
            acts[fi] = length(state.dists)
            fi += 1
        end
        fi > N_FRM && break
    end
    imgs, acts, t_frames
end

println("Running L=2 …"); res2 = run_level(2)
println("Running L=3 …"); res3 = run_level(3)
println("Running L=4 …"); res4 = run_level(4)
println("Simulation done.")

# ── Animation ────────────────────────────────────────────────────────────────
outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)

fig = Figure(size=(1350, 500))
Label(fig[0,1:3],
    "Multigrid RDME: comparing spatial hierarchy depth  (200×200, A→2A)",
    fontsize=12, font=:bold)

lev_labels = ["L=2  (pairs only)", "L=3  (quads)", "L=4  (2×4 blocks)"]
results = [res2, res3, res4]
axs = [Axis(fig[1,c]; title=lev_labels[c], titlesize=11,
            aspect=DataAspect(), yreversed=true) for c in 1:3]
for ax in axs; hidedecorations!(ax); end

obs = [Observable(r[1][1]) for r in results]
for (c, (ax, ob)) in enumerate(zip(axs, obs))
    image!(ax, 0.5..K+0.5, 0.5..K+0.5, ob)
end

t_lbl = Label(fig[2,1:3], "t = 0.0", fontsize=10, color=:gray40)
colgap!(fig.layout, 8); rowgap!(fig.layout, 4)

out = joinpath(outdir, "anim_mg_levels.gif")
println("Rendering → $out")
@time record(fig, out, 1:N_FRM; framerate=12) do f
    for (c, r) in enumerate(results)
        obs[c][] = r[1][f]
    end
    acts3 = results[2][2][f]
    t_lbl.text[] = @sprintf("t = %.1f   |active| ≈ %d  (L=3)", results[1][3][f], acts3)
end
println("Saved → $out")

# ── Static comparison figure at peak wavefront ───────────────────────────────
peak_f = argmax(results[2][2])   # peak active count for L=3
@printf("Peak frame: %d  t=%.1f\n", peak_f, results[1][3][peak_f])

fig2 = Figure(size=(1350, 500))
Label(fig2[0,1:3],
    @sprintf("Multigrid wavefront at peak (t=%.1f): L=2 vs L=3 vs L=4",
             results[1][3][peak_f]),
    fontsize=12, font=:bold)

for (c, (r, lbl)) in enumerate(zip(results, lev_labels))
    ax = Axis(fig2[1,c]; title=lbl, titlesize=11,
              aspect=DataAspect(), yreversed=true)
    hidedecorations!(ax)
    image!(ax, 0.5..K+0.5, 0.5..K+0.5, r[1][peak_f])
    text!(ax, 5, K-5; text=@sprintf("active=%d", r[2][peak_f]),
          fontsize=9, color=:white, align=(:left,:bottom))
end
colgap!(fig2.layout, 8); rowgap!(fig2.layout, 4)
out2 = joinpath(outdir, "fig_mg_levels.png")
save(out2, fig2; px_per_unit=3)
println("Saved → $out2")
