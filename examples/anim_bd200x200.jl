"""
anim_bd200x200.jl — 200×200 BirthDeath RDME adaptive FSP animation.

Left panel : mean-field heatmap ⟨n_k⟩ (exact CME, V-cycle + Tridiagonal generator)
Right panel: adaptive state map  (black=empty, orange=active, blue=equil)

Run: julia --project examples/anim_bd200x200.jl
"""

using DiscStochSim, CairoMakie, Printf

# ── Parameters ────────────────────────────────────────────────────────────────
const K    = 200        # K×K grid
const k_b  = 1.0        # birth rate  →  μ_ss = 2 per voxel
const k_d  = 0.5
const D    = 1.0
const dx   = 1.0
const N0   = 1          # IC molecules at center
const dt   = 0.3
const T_END = 30.0
const N_FRM = 60        # 60 frames at 15fps = 4s animation

const n_max    = 16
const ε_expand = 0.15
const ε_equil  = 0.07

model  = BirthDeathRDME(k_b, k_d, D, dx)
μ_ss   = ss_mean(model)
center = CartesianIndex(K÷2+1, K÷2+1)

@printf("Grid %d×%d  μ_ss=%.1f  T=%.0f  N_frames=%d\n", K, K, μ_ss, T_END, N_FRM)

state = SpatialFSP(model, K, K;
                   n_max, ε_expand, ε_equil,
                   use_block_vcycle = true)
set_ic!(state, center, N0)

# ── Pre-compute frames ─────────────────────────────────────────────────────────
t_frames    = collect(range(0.0, T_END; length=N_FRM))
mean_frames = Vector{Matrix{Float64}}(undef, N_FRM)
stat_frames = Vector{Matrix{Int}}(undef, N_FRM)
kc_frames   = zeros(Int, N_FRM)

mean_frames[1] = mean_grid(state)
stat_frames[1] = status_grid(state)
kc_frames[1]   = n_active(state)

t0_sim = time()
println("Stepping FSP …")
let fi = 2
    for s in 1:round(Int, T_END/dt)
        step!(state, dt)
        t_now = s * dt
        if fi <= N_FRM && t_now >= t_frames[fi] - 1e-9
            mean_frames[fi] = mean_grid(state)
            stat_frames[fi] = status_grid(state)
            kc_frames[fi]   = n_active(state)
            if fi % 10 == 0
                @printf("  frame %3d/%d  t=%5.1f  active=%-5d  equil=%-5d  %.1fs\n",
                        fi, N_FRM, t_now, n_active(state),
                        n_equilibrated(state), time()-t0_sim)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
end
@printf("Simulation done in %.1f s.  Peak active: %d\n",
        time()-t0_sim, maximum(kc_frames))

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))

fig = Figure(size=(1000, 530))

MU_MAX = μ_ss * 3.0

mu_obs = Observable(mean_frames[1]')
sg_obs = Observable(Float32.(stat_frames[1]'))

ax_kw = (aspect=DataAspect(), yreversed=true, xaxisposition=:top, titlesize=11,
         xticksvisible=false, yticksvisible=false,
         xticklabelsvisible=false, yticklabelsvisible=false)

# Panel A: FSP mean heatmap
ax_mu = Axis(fig[1,1]; title="FSP  ⟨n_k⟩  (exact CME, V-cycle solver)", ax_kw...)
hm = heatmap!(ax_mu, 1:K, 1:K, mu_obs;
    colormap=:inferno, colorrange=(0, MU_MAX), highclip=:white, lowclip=:black)
scatter!(ax_mu, [center[2]], [center[1]]; marker=:cross, markersize=8,
    color=:cyan, strokewidth=0)
Colorbar(fig[1,1][1,2], hm; label="⟨n⟩", width=10, labelsize=8, ticklabelsize=7)

# Panel B: adaptive state map
PHASE_CMAP = cgrad([colorant"#1a1a2e", colorant"#e07b39", colorant"#4472c4"])
ax_ph = Axis(fig[1,2]; title="Adaptive state space  (orange=active, blue=equil)", ax_kw...)
heatmap!(ax_ph, 1:K, 1:K, sg_obs; colormap=PHASE_CMAP, colorrange=(0,2))
scatter!(ax_ph, [center[2]], [center[1]]; marker=:cross, markersize=8,
    color=:white, strokewidth=0)

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colgap!(fig.layout, 1, 10)

t_label  = Label(fig[0,1], "t = 0.0"; fontsize=13, halign=:left)
kc_label = Label(fig[0,2],
    "active: $(kc_frames[1])  equil: 0";
    fontsize=11, halign=:center)

Legend(fig[2,1:2],
    [PolyElement(color=colorant"#1a1a2e", strokecolor=:gray40, strokewidth=1),
     PolyElement(color=colorant"#e07b39"),
     PolyElement(color=colorant"#4472c4")],
    ["Empty (not yet reached by wave)", "Active (full CME, V-cycle)", "Equil (Poisson collapsed)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(16,11))

rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 4)

# ── Record ─────────────────────────────────────────────────────────────────────
outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out = joinpath(outdir, "anim_bd200x200.gif")
println("Rendering → $out")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]       = mean_frames[f]'
    sg_obs[]       = Float32.(stat_frames[f]')
    t_label.text[] = "t = $(round(t_frames[f], digits=1))"
    na = kc_frames[f]
    kc_label.text[] = "active: $na  equil: $(n_equilibrated(state))"
end
println("Saved → $out")
