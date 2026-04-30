# fig_talk_result.jl  --  RDME wavefront: probability propagating in space
#
# Birth-death RDME: k_b=5, k_d=0.5 (stationary = Poisson(10) per voxel)
# D=0.1, K=4 voxels
# IC: all molecules in v1 (n=80), rest empty
#
# Marginal heatmaps P(nk=n) at four times show the wavefront sweeping right:
#   - Probability mass in v1 decays (diffusion + death)
#   - v2, v3, v4 light up progressively as the wave arrives
#   - Each voxel's distribution spreads from a delta to Poisson(10)

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)
function save_fig(name, fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=3)
    println("  Saved: $name.{pdf,png}")
end

set_theme!(Theme(
    fontsize = 13,
    Axis = (spinewidth=0.8, xgridvisible=false, ygridvisible=false,
            xticklabelsize=11, yticklabelsize=11, xlabelsize=12, ylabelsize=12,
            titlesize=12),
))

# ── model ──────────────────────────────────────────────────────────────────────
D = 0.1; K = 4; h = 1.0/K; n_max = 100
k_b = 5.0; k_d = 0.5   # stationary mean = 10
model     = RDMEModel1D(D, k_b, k_d)
fine_grid = VoxelGrid(K, h, 0)
sys       = build_rdme_system(model, fine_grid)
bc        = s -> rdme_bc(s, n_max)
rates     = Float64[]

n_stat = round(Int, k_b/k_d)   # = 10
n_src  = 80

@printf("D=%.2f  d=%.1f  stationary mean=%d  K=%d\n", D, D/h^2, n_stat, K)

u0 = CartesianIndex(ntuple(k -> k==1 ? n_src : 0, K))
println("IC: v1=$(n_src), v2..v$(K)=0\n")

# ── simulation ─────────────────────────────────────────────────────────────────
T_SNAP = [0.0, 0.4, 1.2, 2.8]
dt = 0.05; T = maximum(T_SNAP)
n_steps = round(Int, T/dt)

sp = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp, u0, 1.0)

snap_states = Dict{Float64, Tuple{Vector,Vector}}()
snap_states[0.0] = ([u0], [1.0])

t_cur = 0.0
println("Running...")
for step in 1:n_steps
    global t_cur, sp
    dts = min(dt, T-t_cur)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t_cur)
    sp.probs .= max.(expv(dts, A, sp.probs; m=40), 0.0)
    t_cur += dts
    renormalize!(sp); prune_threshold!(sp, 1e-6)

    snap_t = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
    if abs(t_cur - snap_t) < dt/2 && !haskey(snap_states, snap_t)
        snap_states[snap_t] = (copy(sp.states), copy(sp.probs))
        mu = zeros(K)
        for (i,s) in enumerate(sp.states); tv=Tuple(s); p=sp.probs[i]
            for k in 1:K; mu[k]+=p*tv[k]; end; end
        @printf("  t=%.1f  |S|=%d  means: %s\n", t_cur, length(sp),
                join([@sprintf("%.1f",x) for x in mu], "  "))
    end
end
println("Done. Final |S|=$(length(sp))\n")

# ── marginal helper ─────────────────────────────────────────────────────────────
function voxel_marginals(states, probs, K, n_cap)
    M = zeros(K, n_cap+1)
    for (i,s) in enumerate(states)
        tv = Tuple(s); p = probs[i]
        for k in 1:K
            n = tv[k]; n <= n_cap && (M[k, n+1] += p)
        end
    end
    M  # K × (n_cap+1)
end

n_cap = n_src + 5   # only need to show up to ~85

global_max = maximum(
    maximum(voxel_marginals(snap_states[ts]..., K, n_cap))
    for ts in T_SNAP if haskey(snap_states, ts))

# ── figure ─────────────────────────────────────────────────────────────────────
fig = Figure(size=(980, 380))

panel_labels = ["(a)  t = 0", "(b)  t = 0.4", "(c)  t = 1.2", "(d)  t = 2.8"]
vox_labels   = ["v1","v2","v3","v4"]

for (pi, ts) in enumerate(T_SNAP)
    haskey(snap_states, ts) || continue
    sts, prs = snap_states[ts]
    M = voxel_marginals(sts, prs, K, n_cap)

    ax = Axis(fig[1, pi];
        title  = panel_labels[pi],
        xlabel = "voxel  k",
        ylabel = pi == 1 ? "molecules  n" : "",
        xticks = (1:K, vox_labels))

    heatmap!(ax, 0.5:(K+0.5), -0.5:(n_cap+0.5), M;
        colormap=:inferno,
        colorrange=(0, global_max*0.7))

    # Stationary mean reference line
    hlines!(ax, [Float64(n_stat)];
        color=(:white, 0.45), linewidth=1.0, linestyle=:dash)

    if pi == 1
        text!(ax, 4.3, Float64(n_stat);
              text="stationary\nmean = $n_stat",
              color=(:white,0.7), fontsize=8, align=(:right,:center))
        text!(ax, 1.0, n_src+1;
              text="initial\nsource", color=(:white,0.7),
              fontsize=8, align=(:center,:bottom))
    end

    ylims!(ax, 0, n_cap)
end

Colorbar(fig[1, K+1]; colormap=:inferno,
    colorrange=(0, global_max*0.7),
    label="P(nₖ = n)", width=14, labelsize=10, ticklabelsize=9)

Label(fig[2, 1:K];
    text="Probability wavefront sweeps left→right.  " *
         "v1 decays as molecules diffuse out and die (80→36).  " *
         "v2 rises, overshoots the stationary mean, then relaxes.  " *
         "v3, v4 activate progressively as the wave arrives.  " *
         "Diffusion and reactions both visible in the evolving distributions.",
    fontsize=9, color=(:black,0.6), tellwidth=false)

save_fig("fig_talk_result", fig)
println("Done.")
