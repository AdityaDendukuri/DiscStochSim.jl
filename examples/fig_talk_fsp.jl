# fig_talk_fsp.jl  --  Adaptive FSP visualization for the RDME
#
# Birth-death RDME, K=2 voxels, IC = (n=60, n=0)
# Shows the expand → solve → prune cycle at one time step,
# mirroring the reference figure in the paper.
#
# Top row: probability support (active states) at 4 time snapshots
# Bottom row: single time step zoom — current / after expand / after prune

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
    Axis = (spinewidth=0.9, xgridvisible=false, ygridvisible=false,
            xticklabelsize=10, yticklabelsize=10, xlabelsize=11, ylabelsize=11,
            titlesize=11),
))

# ── model ──────────────────────────────────────────────────────────────────────
D = 0.1; K = 2; h = 1.0/K; n_max = 80
k_b = 5.0; k_d = 0.5
model     = RDMEModel1D(D, k_b, k_d)
fine_grid = VoxelGrid(K, h, 0)
sys       = build_rdme_system(model, fine_grid)
bc        = s -> rdme_bc(s, n_max)
rates     = Float64[]

n_stat = round(Int, k_b/k_d)
n_src  = 60

@printf("K=%d  D=%.2f  d=%.1f  stationary mean=%d\n", K, D, D/h^2, n_stat)

u0 = CartesianIndex(n_src, 0)

# ── helper: 2D joint density matrix ───────────────────────────────────────────
function joint_matrix(states, probs, n_cap)
    M = zeros(n_cap+1, n_cap+1)
    for (i,s) in enumerate(states)
        tv = Tuple(s); p = probs[i]
        tv[1] <= n_cap && tv[2] <= n_cap && (M[tv[1]+1, tv[2]+1] += p)
    end
    M   # (n_cap+1) × (n_cap+1): rows=n1, cols=n2
end

function support_matrix(states, n_cap)
    # Binary matrix: 1 if state exists, 0 otherwise
    M = zeros(n_cap+1, n_cap+1)
    for s in states
        tv = Tuple(s)
        tv[1] <= n_cap && tv[2] <= n_cap && (M[tv[1]+1, tv[2]+1] = 1.0)
    end
    M
end

# ── run to t=0.8 saving snapshots ─────────────────────────────────────────────
T_SNAP = [0.0, 0.2, 0.5, 0.8]
dt = 0.05; T = 0.8; n_steps = round(Int, T/dt)

sp = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp, u0, 1.0)

snaps = Dict{Float64, Tuple{Vector,Vector}}()
snaps[0.0] = ([u0], [1.0])

t_cur = 0.0
for step in 1:n_steps
    global t_cur, sp
    dts   = min(dt, T-t_cur)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t_cur)
    sp.probs .= max.(expv(dts, A, sp.probs; m=40), 0.0)
    t_cur += dts
    renormalize!(sp); prune_threshold!(sp, 1e-6)
    snap_t = T_SNAP[argmin(abs.(T_SNAP .- t_cur))]
    abs(t_cur - snap_t) < dt/2 && !haskey(snaps, snap_t) &&
        (snaps[snap_t] = (copy(sp.states), copy(sp.probs)))
end
@printf("At t=0.8: |S|=%d\n", length(sp))

# ── capture expand → solve → prune at t=0.5 ───────────────────────────────────
sts_before, prs_before = snaps[0.5]

# Stage 1: before (current distribution)
sp1 = StateSpace{CartesianIndex{K}, Float64}()
for (i,s) in zip(eachindex(sts_before), sts_before)
    add_state!(sp1, s, prs_before[i])
end

# Stage 2: after expand (new states have p=0)
sp2 = deepcopy(sp1)
expand!(sp2, sys, bc; depth=1)

# Stage 3: after solve (p updated) then prune
sp3 = deepcopy(sp2)
A3, = build_generator(sp3, sys, rates, 0.5)
sp3.probs .= max.(expv(dt, A3, sp3.probs; m=40), 0.0)
renormalize!(sp3); prune_threshold!(sp3, 1e-6)

@printf("t=0.5: before=%d  after expand=%d  after prune=%d\n",
        length(sp1), length(sp2), length(sp3))

# ── figure ─────────────────────────────────────────────────────────────────────
n_cap = n_src + 5
fig   = Figure(size=(1000, 480))

# ── Top row: evolution snapshots (joint density P(n1,n2)) ─────────────────────
snap_titles = ["(a)  t = 0", "(b)  t = 0.2", "(c)  t = 0.5", "(d)  t = 0.8"]

for (pi, ts) in enumerate(T_SNAP)
    haskey(snaps, ts) || continue
    M = joint_matrix(snaps[ts]..., n_cap)
    panel_max = maximum(M)
    panel_max < 1e-10 && continue

    ax = Axis(fig[1, pi];
        title  = snap_titles[pi],
        xlabel = "n₁  (voxel 1)",
        ylabel = pi == 1 ? "n₂  (voxel 2)" : "")

    heatmap!(ax, 0:n_cap, 0:n_cap, M;
        colormap=:inferno, colorrange=(0, panel_max*0.9))

    # Stationary mean marker
    scatter!(ax, [n_stat], [n_stat]; color=:white, marker=:cross,
             markersize=10, strokewidth=2)
    if pi == 4
        text!(ax, n_stat+1, n_stat+2; text="stationary\nmean",
              color=:white, fontsize=8, align=(:left,:bottom))
    end
    xlims!(ax, 0, n_cap); ylims!(ax, 0, n_cap)
end

global_max = maximum(joint_matrix(snaps[0.5]..., n_cap))

# ── Bottom row: expand → solve → prune cycle ──────────────────────────────────
stage_titles = ["(e)  Current distribution",
                "(f)  After expand  (+boundary shell)",
                "(g)  After solve + prune"]

stages = [
    (sp1.states, sp1.probs, false),    # current
    (sp2.states, sp2.probs, true),     # after expand: show new states
    (sp3.states, sp3.probs, false),    # after prune
]

for (pi, (sts, prs, show_boundary)) in enumerate(stages)
    ax = Axis(fig[2, pi];
        title  = stage_titles[pi],
        xlabel = "n₁  (voxel 1)",
        ylabel = pi == 1 ? "n₂  (voxel 2)" : "")

    if show_boundary
        # Highlight new boundary states (p≈0) vs existing states
        S_before = Set(sts_before)
        new_sts = [s for s in sts if s ∉ S_before]
        old_sts_ids = [i for (i,s) in enumerate(sts) if s ∈ S_before]

        # Old states with probability
        M_old = joint_matrix(sts[old_sts_ids], prs[old_sts_ids], n_cap)
        heatmap!(ax, 0:n_cap, 0:n_cap, M_old;
            colormap=:inferno, colorrange=(0, global_max*0.8))

        # New boundary states (white dots)
        if !isempty(new_sts)
            n1s = [Tuple(s)[1] for s in new_sts]
            n2s = [Tuple(s)[2] for s in new_sts]
            scatter!(ax, Float64.(n1s), Float64.(n2s);
                color=(:white, 0.6), markersize=2, marker=:circle)
        end
        text!(ax, 0.05, 0.92; text="$(length(new_sts)) new\nboundary\nstates",
              color=:white, fontsize=8, space=:relative, align=(:left,:top))
    else
        M = joint_matrix(sts, prs, n_cap)
        heatmap!(ax, 0:n_cap, 0:n_cap, M;
            colormap=:inferno, colorrange=(0, global_max*0.8))
    end

    text!(ax, 0.98, 0.05; text="|S| = $(length(sts))",
          color=:white, fontsize=9, space=:relative, align=(:right,:bottom))
    xlims!(ax, 0, n_cap); ylims!(ax, 0, n_cap)
end

Colorbar(fig[2, 5]; colormap=:inferno,
    colorrange=(0, global_max*0.8),
    label="P(n₁, n₂)", width=14, labelsize=10, ticklabelsize=9)
Label(fig[1, 5]; text="per-panel\nscale", fontsize=8, color=:gray, tellheight=false)

save_fig("fig_talk_fsp", fig)
println("Done.")
