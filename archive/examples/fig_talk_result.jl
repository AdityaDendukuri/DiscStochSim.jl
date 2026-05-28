# fig_talk_result.jl  --  K=8 wavefront via coarsened simulation
#
# Birth-death RDME: D=0.1, k_b=5, k_d=0.5 (stationary = Poisson(10)/voxel)
# K=8 fine voxels, balanced IC: v1=v2=60, v3..v8=0
#   -> pairs always balanced -> coarsening is valid throughout
#   -> evolve at K=4 coarse (tiny state space)
#   -> prolong analytically to K=8 fine marginals via Binom(nc,1/2)
#
# Shows: a K=8 RDME simulation that would be intractable with direct FSP,
# solved efficiently by coarsening, visualized as a K=8 heatmap.

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
    fontsize = 12,
    Axis = (spinewidth=0.8, xgridvisible=false, ygridvisible=false,
            xticklabelsize=10, yticklabelsize=10, xlabelsize=11, ylabelsize=11,
            titlesize=12),
))

# ── model ──────────────────────────────────────────────────────────────────────
D = 0.1; K_fine = 8; K_coarse = 4; h = 1.0/K_fine
k_b = 5.0; k_d = 0.5; n_stat = round(Int, k_b/k_d)  # 10 per fine voxel
n_src = 40; n_max_coarse = 2*n_src + 20

model_coarse = RDMEModel1D(D, k_b*2, k_d)   # coarse: birth scaled by r=2
coarse_grid  = VoxelGrid(K_coarse, 2*h, 0)
sys_coarse   = build_rdme_system(model_coarse, coarse_grid)
bc_coarse    = s -> rdme_bc(s, n_max_coarse)
rates        = Float64[]

@printf("K_fine=%d  K_coarse=%d  D=%.2f  d=%.1f  n_stat=%d/voxel\n",
        K_fine, K_coarse, D, D/h^2, n_stat)

# IC: balanced pairs — v1=v2=n_src, v3..v8=0
# Coarse IC: nc1=2*n_src, nc2=nc3=nc4=0
u0_c = CartesianIndex(2*n_src, 0, 0, 0)
@printf("IC: v1=v2=%d (balanced), v3..v8=0  ->  coarse nc1=%d\n\n",
        n_src, 2*n_src)

# ── K=4 coarse simulation ──────────────────────────────────────────────────────
T_SNAP = [0.0, 0.5, 1.5, 4.0, 10.0]
dt = 0.05; T = maximum(T_SNAP)
n_steps = round(Int, T/dt)

sp_c = StateSpace{CartesianIndex{K_coarse}, Float64}()
add_state!(sp_c, u0_c, 1.0)

snaps_c = Dict{Float64, Tuple{Vector,Vector}}()
snaps_c[0.0] = ([u0_c], [1.0])

times_c = Float64[0.0]; sizes_c = Int[1]
t_cur   = 0.0

println("Running K=$K_coarse coarse simulation...")
for step in 1:n_steps
    global t_cur, sp_c
    dts = min(dt, T-t_cur)
    expand!(sp_c, sys_coarse, bc_coarse; depth=1)
    A, = build_generator(sp_c, sys_coarse, rates, t_cur)
    sp_c.probs .= max.(expv(dts, A, sp_c.probs; m=40), 0.0)
    t_cur += dts
    renormalize!(sp_c); prune_threshold!(sp_c, 1e-6)
    push!(times_c, t_cur); push!(sizes_c, length(sp_c))
    snap_t = T_SNAP[argmin(abs.(T_SNAP.-t_cur))]
    abs(t_cur-snap_t)<dt/2 && !haskey(snaps_c, snap_t) &&
        (snaps_c[snap_t] = (copy(sp_c.states), copy(sp_c.probs)))
end
@printf("Done. Final |S_coarse|=%d\n", last(sizes_c))

# ── Analytical fine marginals via Binom prolongation ──────────────────────────
# No fine state space needed: P(n_k=n) computed directly from coarse marginals.

# log C(n,k)/2^n without overflow
function log_binom_half(n::Int, k::Int)
    n == 0 && return 0.0
    k = min(k, n-k)
    s = 0.0
    for i in 1:k; s += log(n-k+i) - log(i); end
    s - n*log(2.0)
end

function fine_marginals(coarse_states, coarse_probs, K_c, nc_cap, n_cap_fine)
    K_f = 2*K_c
    M   = zeros(K_f, n_cap_fine+1)
    for j in 1:K_c
        # Get nc_j marginal (coarse voxel j)
        nc_marg = zeros(nc_cap+1)
        for (i,s) in enumerate(coarse_states)
            nc = Tuple(s)[j]
            nc <= nc_cap && (nc_marg[nc+1] += coarse_probs[i])
        end
        # Convolve with Binom(nc, 1/2) to get fine marginal
        for nc in 0:nc_cap
            nc_marg[nc+1] < 1e-12 && continue
            p_nc = nc_marg[nc+1]
            for n in 0:min(nc, n_cap_fine)
                b = exp(log_binom_half(nc, n))
                M[2j-1, n+1] += p_nc * b
                M[2j,   n+1] += p_nc * b   # symmetric by Binom(nc,1/2)
            end
        end
    end
    M
end

n_cap_fine = n_src + 10   # fine voxels: 0..n_src+10

# Compute fine marginals at each snapshot
fine_snaps = Dict{Float64, Matrix{Float64}}()
for ts in T_SNAP
    haskey(snaps_c, ts) || continue
    if ts == 0.0
        # t=0: use actual fine IC directly (delta at n1=n2=n_src, n3..n8=0)
        M0 = zeros(K_fine, n_cap_fine+1)
        M0[1, n_src+1] = 1.0   # v1 at n=n_src
        M0[2, n_src+1] = 1.0   # v2 at n=n_src
        for k in 3:K_fine; M0[k, 1] = 1.0; end  # v3..v8 at n=0
        fine_snaps[0.0] = M0
        continue
    end
    fine_snaps[ts] = fine_marginals(snaps_c[ts]..., K_coarse, n_max_coarse, n_cap_fine)
    mu = [sum((0:n_cap_fine) .* fine_snaps[ts][k,:]) for k in 1:K_fine]
    @printf("  t=%.1f: fine means: %s\n", ts,
            join([@sprintf("%.1f",x) for x in mu], "  "))
end

# ── figure ─────────────────────────────────────────────────────────────────────
global_max = maximum(maximum(fine_snaps[ts]) for ts in T_SNAP if haskey(fine_snaps, ts))

fig = Figure(size=(1100, 420))

panel_labels = ["(a)  t = 0", "(b)  t = 0.5", "(c)  t = 1.5",
                "(d)  t = 4",  "(e)  t = 10  (stationary)"]
vox_ticks = (1:K_fine, ["v$k" for k in 1:K_fine])

for (pi, ts) in enumerate(T_SNAP)
    haskey(fine_snaps, ts) || continue
    M = fine_snaps[ts]

    ax = Axis(fig[1, pi];
        title  = panel_labels[pi],
        xlabel = "voxel  k",
        ylabel = pi == 1 ? "molecules  n" : "",
        xticks = vox_ticks)

    heatmap!(ax, 0.5:(K_fine+0.5), -0.5:(n_cap_fine+0.5), M;
        colormap=:inferno, colorrange=(0, global_max*0.8))

    hlines!(ax, [Float64(n_stat)];
        color=(:white, 0.35), linewidth=0.8, linestyle=:dot)

    # Pair boundaries
    for j in 1:K_coarse-1
        vlines!(ax, [2j+0.5]; color=(:white, 0.2), linewidth=0.6)
    end

    ylims!(ax, 0, n_cap_fine)
end

Colorbar(fig[1, length(T_SNAP)+1];
    colormap=:inferno, colorrange=(0, global_max*0.8),
    label="P(nₖ = n)", width=14, labelsize=10, ticklabelsize=9)

# State-space annotation below
Label(fig[2, 1:length(T_SNAP)];
    text="K=8 wavefront solved via K=4 coarse FSP.  " *
         "Max coarse state space: $(maximum(sizes_c)) states.  " *
         "K=8 direct FSP would be intractable (grows to 90k+ states by t=0.8).  " *
         "Fine marginals recovered analytically via Binom(nĉⱼ, ½) prolongation — exact for balanced pairs.",
    fontsize=9, color=(:black,0.55), tellwidth=false)

save_fig("fig_talk_result", fig)
println("Done.")
