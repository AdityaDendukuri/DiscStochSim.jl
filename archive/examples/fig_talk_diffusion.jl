# fig_talk_diffusion.jl
#
# K=2 birth-death RDME, IC = (n_src, 0): all molecules in voxel 1
# Show joint distribution P(n1,n2) at 4 times as it spreads toward Binom equilibrium
#
# Key visual story:
#   - Diffusion moves probability ALONG constant-total diagonals (n1+n2=const)
#   - Reactions move probability ACROSS diagonals (change total)
#   - Equilibrium: Binom(nc, 1/2) distribution along each diagonal

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

# ── Row 1: birth-death (monostable) ───────────────────────────────────────────
D_bd = 0.5; K = 2; h = 1.0/K; n_max_bd = 80
k_b = 5.0; k_d = 0.5
n_stat = round(Int, k_b/k_d)   # = 10
n_src  = 60

model_bd  = RDMEModel1D(D_bd, k_b, k_d)
fine_grid = VoxelGrid(K, h, 0)
sys_bd    = build_rdme_system(model_bd, fine_grid)
bc_bd     = s -> rdme_bc(s, n_max_bd)
rates     = Float64[]

T_SNAP_BD = [0.0, 0.15, 0.4, 1.5]
dt_bd = 0.02; u0_bd = CartesianIndex(n_src, 0)
sp_bd = StateSpace{CartesianIndex{K}, Float64}(); add_state!(sp_bd, u0_bd, 1.0)
snaps_bd = Dict{Float64,Tuple{Vector,Vector}}(); snaps_bd[0.0] = ([u0_bd],[1.0])
t_cur = 0.0
for _ in 1:round(Int, maximum(T_SNAP_BD)/dt_bd)
    global t_cur, sp_bd; dts=min(dt_bd,maximum(T_SNAP_BD)-t_cur)
    expand!(sp_bd,sys_bd,bc_bd;depth=1); A,=build_generator(sp_bd,sys_bd,rates,t_cur)
    sp_bd.probs .= max.(expv(dts,A,sp_bd.probs;m=40),0.0); t_cur+=dts
    renormalize!(sp_bd); prune_threshold!(sp_bd,1e-8)
    snap_t=T_SNAP_BD[argmin(abs.(T_SNAP_BD.-t_cur))]
    abs(t_cur-snap_t)<dt_bd/2 && !haskey(snaps_bd,snap_t) &&
        (snaps_bd[snap_t]=(copy(sp_bd.states),copy(sp_bd.probs)))
end
println("Birth-death done. |S|=$(length(sp_bd))")

# ── Row 2: Schlögl (bistable front IC) ────────────────────────────────────────
D_sg = 0.005; n_max_sg = 200
model_sg  = SchloglModel1D(D_sg)
sys_sg    = build_schlogl_rdme_system(model_sg, fine_grid)
bc_sg     = s -> rdme_bc(s, n_max_sg)

fp = schlogl_fixed_points(model_sg; n_max=250)
n_lo, _, n_hi_sg = fp[1], fp[2], fp[3]
@printf("Schlogl: n_lo=%d  n_hi=%d\n", n_lo, n_hi_sg)

T_SNAP_SG = [0.0, 1.0, 5.0, 10.0]
dt_sg = 0.1; u0_sg = CartesianIndex(n_lo, n_hi_sg)
sp_sg = StateSpace{CartesianIndex{K}, Float64}(); add_state!(sp_sg, u0_sg, 1.0)
snaps_sg = Dict{Float64,Tuple{Vector,Vector}}(); snaps_sg[0.0] = ([u0_sg],[1.0])
t_cur = 0.0
for _ in 1:round(Int, maximum(T_SNAP_SG)/dt_sg)
    global t_cur, sp_sg; dts=min(dt_sg,maximum(T_SNAP_SG)-t_cur)
    expand!(sp_sg,sys_sg,bc_sg;depth=1); A_sg,=build_generator(sp_sg,sys_sg,rates,t_cur)
    sp_sg.probs .= max.(expv(dts,A_sg,sp_sg.probs;m=40),0.0); t_cur+=dts
    renormalize!(sp_sg); prune_threshold!(sp_sg,1e-10)
    snap_t=T_SNAP_SG[argmin(abs.(T_SNAP_SG.-t_cur))]
    abs(t_cur-snap_t)<dt_sg/2 && !haskey(snaps_sg,snap_t) &&
        (snaps_sg[snap_t]=(copy(sp_sg.states),copy(sp_sg.probs)))
end
println("Schlögl transient done. |S|=$(length(sp_sg))")

# Stationary distribution via: P(n1,n2) = P_coarse(nc) × π(n1|nc)
# Step 1: compute dynamic-π from K=2 stationary solve
println("Computing Schlögl stationary distribution...")
sp_stat = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:n_max_sg, n2 in 0:n_max_sg
    add_state!(sp_stat, CartesianIndex(n1,n2), 0.0)
end
sp_stat.probs[sp_stat.index[CartesianIndex(n_lo, n_hi_sg)]] = 1.0

# Stationary = symmetrize over the full state space.
# The model is symmetric n1<->n2, so π_stat(n1,n2) = π_stat(n2,n1).
# Build a lookup from the T=10 adaptive snapshot, then fill the FULL state space
# (sp_stat has all 201^2 states) with P(n1,n2) + P(n2,n1) from the snapshot.
sts_10, prs_10 = snaps_sg[10.0]
prob_lookup = Dict(sts_10[i] => prs_10[i] for i in eachindex(sts_10))
for (i,s) in enumerate(sp_stat.states)
    n1, n2 = Tuple(s)
    p_fwd  = get(prob_lookup, CartesianIndex(n1, n2), 0.0)
    p_mirr = get(prob_lookup, CartesianIndex(n2, n1), 0.0)
    sp_stat.probs[i] = p_fwd + p_mirr
end
sp_stat.probs ./= sum(sp_stat.probs)
snaps_sg[Inf] = (copy(sp_stat.states), copy(sp_stat.probs))
# Check bimodality
mid = (n_lo + n_hi_sg) ÷ 2
sts_s, prs_s = snaps_sg[Inf]
wlh = sum(prs_s[i] for (i,s) in enumerate(sts_s) if Tuple(s)[1]<mid && Tuple(s)[2]>mid; init=0.0)
whl = sum(prs_s[i] for (i,s) in enumerate(sts_s) if Tuple(s)[1]>mid && Tuple(s)[2]<mid; init=0.0)
@printf("  Stationary: P(lo,hi)=%.4f  P(hi,lo)=%.4f\n", wlh, whl)
println("Schlögl stationary done.")

# ── joint density matrix ───────────────────────────────────────────────────────
function joint_matrix(states, probs, n_cap)
    M = zeros(n_cap+1, n_cap+1)
    for (i,s) in enumerate(states)
        tv=Tuple(s); tv[1]<=n_cap && tv[2]<=n_cap && (M[tv[1]+1,tv[2]+1]+=probs[i])
    end
    M
end

function draw_row!(fig, row, snaps, T_SNAP, n_cap, row_label;
                   stat_pt=nothing, diag_ncs=Int[], annot_ic=nothing)
    labels_t = ["t = $(ts)" for ts in T_SNAP]
    for (pi, ts) in enumerate(T_SNAP)
        haskey(snaps, ts) || continue
        M = joint_matrix(snaps[ts]..., n_cap)
        panel_max = maximum(M); panel_max < 1e-10 && (panel_max = 1.0)

        ylabel_str = (pi == 1) ? "n₂  (voxel 2)\n$(row_label)" : ""
        ax = Axis(fig[row, pi];
            title  = pi==1 ? "t = $(ts)" : "t = $(ts)",
            xlabel = row==2 ? "n₁  (voxel 1)" : "",
            ylabel = ylabel_str)

        heatmap!(ax, 0:n_cap, 0:n_cap, M;
            colormap=:inferno, colorrange=(0, panel_max*0.9))

        for nc in diag_ncs
            nc>0 && nc<=2*n_cap && lines!(ax, [0,min(nc,n_cap)],[min(nc,n_cap),0];
                color=(:white,0.18), linewidth=0.8, linestyle=:dash)
        end

        if !isnothing(stat_pt)
            scatter!(ax, [Float64(stat_pt[1])],[Float64(stat_pt[2])];
                color=:white, marker=:cross, markersize=9, strokewidth=1.8)
        end

        if pi==1 && !isnothing(annot_ic)
            text!(ax, annot_ic[1], annot_ic[2];
                text=annot_ic[3], color=(:white,0.75), fontsize=8,
                align=(:left,:bottom))
        end
        xlims!(ax, 0, n_cap); ylims!(ax, 0, n_cap)
    end
end

# ── figure ─────────────────────────────────────────────────────────────────────
n_cap_bd = n_src + 5
n_cap_sg = n_hi_sg + 10

fig = Figure(size=(1100, 520))

# Row 1: birth-death (4 panels)
for (pi, ts) in enumerate(T_SNAP_BD)
    haskey(snaps_bd, ts) || continue
    M = joint_matrix(snaps_bd[ts]..., n_cap_bd)
    pm = maximum(M); pm < 1e-10 && (pm=1.0)
    ax = Axis(fig[1,pi]; title="t = $ts",
              xlabel=pi==4 ? "n₁" : "",
              ylabel=pi==1 ? "n₂\nbirth-death" : "")
    heatmap!(ax, 0:n_cap_bd, 0:n_cap_bd, M;
             colormap=:inferno, colorrange=(0,pm*0.9))
    for nc in [n_src, n_stat*2]
        nc>0 && lines!(ax,[0,min(nc,n_cap_bd)],[min(nc,n_cap_bd),0];
            color=(:white,0.18),linewidth=0.8,linestyle=:dash)
    end
    scatter!(ax,[Float64(n_stat)],[Float64(n_stat)];
        color=:white,marker=:cross,markersize=9,strokewidth=1.8)
    pi==1 && text!(ax, 2, n_cap_bd-6; text="IC: ($(n_src),0)",
                   color=(:white,0.8),fontsize=8)
    xlims!(ax,0,n_cap_bd); ylims!(ax,0,n_cap_bd)
end

# Row 2: Schlögl transient (4 panels) + stationary (1 panel)
sg_keys = [T_SNAP_SG..., Inf]
sg_titles = ["t = $(T_SNAP_SG[1])", "t = $(T_SNAP_SG[2])",
             "t = $(T_SNAP_SG[3])", "t = $(T_SNAP_SG[4])",
             "t → ∞\n(stationary)"]

for (pi, key) in enumerate(sg_keys)
    haskey(snaps_sg, key) || continue
    M = joint_matrix(snaps_sg[key]..., n_cap_sg)
    pm = maximum(M); pm < 1e-10 && (pm=1.0)
    ax = Axis(fig[2,pi]; title=sg_titles[pi],
              xlabel="n₁",
              ylabel=pi==1 ? "n₂\nSchlögl" : "")
    heatmap!(ax, 0:n_cap_sg, 0:n_cap_sg, M;
             colormap=:inferno, colorrange=(0,pm*0.9))
    for nc in [n_lo+n_hi_sg, 2*n_lo, 2*n_hi_sg]
        nc>0 && nc<=2*n_cap_sg && lines!(ax,[0,min(nc,n_cap_sg)],[min(nc,n_cap_sg),0];
            color=(:white,0.18),linewidth=0.8,linestyle=:dash)
    end
    if pi == 1
        text!(ax,2,n_cap_sg-15; text="IC:\n(n_lo,n_hi)",color=(:white,0.8),fontsize=8)
    end
    if pi == length(sg_keys)
        text!(ax, n_lo+3, n_hi_sg-8; text="(lo,hi)\nmode",
              color=(:white,0.8), fontsize=8, align=(:left,:top))
        text!(ax, n_hi_sg-3, n_lo+8; text="(hi,lo)\nmode",
              color=(:white,0.8), fontsize=8, align=(:right,:bottom))
        text!(ax, 0.05, 0.06;
              text="true stationary:\nboth orientations\nequally probable",
              color=(:white,0.55), fontsize=7, space=:relative, align=(:left,:bottom))
    end
    xlims!(ax,0,n_cap_sg); ylims!(ax,0,n_cap_sg)
end

save_fig("fig_talk_diffusion", fig)
println("Done.")
