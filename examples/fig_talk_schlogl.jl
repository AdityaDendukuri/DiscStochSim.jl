# fig_talk_schlogl.jl
#
# K=4 Schlögl: three-way comparison
#   1. Direct FSP (ground truth)
#   2. Coarse Binomial (K=2): fast but WRONG — front collapses
#   3. Adaptive V-cycle dynamic-π (K=4->K=2): fast AND correct
#
# Shows why Schlögl needs the adaptive approach, not just naive coarsening.

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
            titlesize=11),
))

const C_DIRECT = RGBf(0.208, 0.475, 0.694)
const C_COARSE = RGBf(0.835, 0.369, 0.000)
const C_VCYCLE = RGBf(0.000, 0.620, 0.451)
const C_GRAY   = RGBf(0.5, 0.5, 0.5)

D = 0.005; K = 4; h = 1.0/K; n_max = 180
model = SchloglModel1D(D); rates = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model; n_max=250)
n_lo, _, n_hi = fp[1], fp[2], fp[3]
@printf("n_lo=%d  n_hi=%d\n", n_lo, n_hi)

# IC: front — pair 1 (v1,v2) has front, pairs 2 is in hi state
u0 = CartesianIndex(n_lo, n_hi, n_hi, n_hi)
@printf("IC: (n_lo, n_hi, n_hi, n_hi) — front in pair 1\n\n")

dt = 0.2; T = 2.0; n_steps = round(Int, T/dt)
T_SNAP = 0.0:dt:T
T_HEAT = 1.0   # time for spatial heatmap snapshot

function vmeans(sp, K)
    mu = zeros(K)
    for (i,s) in enumerate(sp.states); tv=Tuple(s); p=sp.probs[i]
        for k in 1:K; mu[k]+=p*tv[k]; end; end
    mu
end

function fsp_step!(sp, sys, bc, t, dts; prune=1e-10)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t)
    sp.probs .= max.(expv(dts, A, sp.probs; m=40), 0.0)
    renormalize!(sp); prune_threshold!(sp, prune)
end

# ── 1. Direct K=4 FSP ────────────────────────────────────────────────────────
println("Running K=4 direct FSP...")
sys4 = build_schlogl_rdme_system(model, fine_grid)
bc4  = s -> rdme_bc(s, n_max)
sp4  = StateSpace{CartesianIndex{K}, Float64}(); add_state!(sp4, u0, 1.0)
times_d=[0.0]; means_d=[vmeans(sp4,K)]; t_cur=0.0
snap_d = nothing
for _ in 1:n_steps
    global t_cur, snap_d; dts=min(dt,T-t_cur); fsp_step!(sp4,sys4,bc4,t_cur,dts); t_cur+=dts
    push!(times_d,t_cur); push!(means_d,vmeans(sp4,K))
    abs(t_cur-T_HEAT)<dt/2 && isnothing(snap_d) && (snap_d=(copy(sp4.states),copy(sp4.probs)))
end
@printf("  Final |S|=%d\n", length(sp4))

# ── 2. Coarse Binomial (K=2) ──────────────────────────────────────────────────
println("Running K=2 coarse Binomial...")
sys_c = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
bc_c  = s -> rdme_bc(s, 2*n_max)
u0_c  = CartesianIndex(n_lo+n_hi, 2*n_hi)   # restrict IC to coarse
sp_c  = StateSpace{CartesianIndex{2}, Float64}(); add_state!(sp_c, u0_c, 1.0)
times_c=[0.0]
# Coarse means -> fine approx: nc/2 per voxel (Binomial split)
mu0_c = [Float64(n_lo+n_hi)/2, Float64(n_lo+n_hi)/2, Float64(n_hi), Float64(n_hi)]
means_c=[mu0_c]; t_cur=0.0; snap_c=nothing
for _ in 1:n_steps
    global t_cur, snap_c; dts=min(dt,T-t_cur); fsp_step!(sp_c,sys_c,bc_c,t_cur,dts); t_cur+=dts
    mu = vmeans(sp_c, 2)
    push!(times_c, t_cur)
    push!(means_c, [mu[1]/2, mu[1]/2, mu[2]/2, mu[2]/2])
    abs(t_cur-T_HEAT)<dt/2 && isnothing(snap_c) && (snap_c=(copy(sp_c.states),copy(sp_c.probs)))
end
@printf("  Final |S_coarse|=%d\n", length(sp_c))

# ── 3. V-cycle dynamic-π ──────────────────────────────────────────────────────
println("Computing dynamic-pi table...")
sp_pair = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max; add_state!(sp_pair,CartesianIndex(n1,n2),0.0); end
sp_pair.probs[sp_pair.index[CartesianIndex(n_lo,n_hi)]] = 1.0
A_p, = build_generator(sp_pair, build_schlogl_rdme_system(model,VoxelGrid(2,h,0)), rates, 0.0)
pi_table = compute_dynamic_pi(sp_pair, A_p; n_max=n_max)
println("  Done.")
println("Running K=4 V-cycle (dynamic-π)...")
sp_vc = StateSpace{CartesianIndex{K}, Float64}(); add_state!(sp_vc, u0, 1.0)
times_v=[0.0]; means_v=[vmeans(sp_vc,K)]; t_cur=0.0; snap_v=nothing
for _ in 1:n_steps
    global t_cur, sp_vc, snap_v; dts=min(dt,T-t_cur)
    sp_vc = two_level_vcycle_schlogl_injection(
        sp_vc, model, fine_grid, coarse_grid, pi_table, rates, t_cur, dts;
        coarse_n_max=2*n_max)
    t_cur+=dts; renormalize!(sp_vc); prune_threshold!(sp_vc,1e-10)
    push!(times_v,t_cur); push!(means_v,vmeans(sp_vc,K))
    abs(t_cur-T_HEAT)<dt/2 && isnothing(snap_v) && (snap_v=(copy(sp_vc.states),copy(sp_vc.probs)))
end
@printf("  Final |S|=%d\n", length(sp_vc))

# ── marginal helper ────────────────────────────────────────────────────────────
n_cap = n_hi + 20
function voxel_marginals(states, probs, K, n_cap)
    M = zeros(K, n_cap+1)
    for (i,s) in enumerate(states)
        tv=Tuple(s); p=probs[i]
        for k in 1:K; tv[k]<=n_cap && (M[k,tv[k]+1]+=p); end
    end
    M
end

# For coarse snap: prolong via Binom to get fine marginals
function coarse_marginals_binom(states_c, probs_c, K_c, n_cap)
    K_f = 2*K_c
    M   = zeros(K_f, n_cap+1)
    for (i,s) in enumerate(states_c)
        tv=Tuple(s); p=probs_c[i]
        for j in 1:K_c
            nc = tv[j]
            for n1 in 0:min(nc,n_cap)
                b = exp(sum(log(k) for k in (nc-n1+1):nc; init=0.0) -
                        sum(log(k) for k in 1:n1; init=0.0) - nc*log(2.0))
                M[2j-1,n1+1] += p*b
                M[2j,  n1+1] += p*b
            end
        end
    end
    M
end

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(860, 560))

vox_colors = [C_VCYCLE, C_VCYCLE, C_GRAY, C_GRAY]
vox_styles = [:solid, :dash, :dot, :dashdot]
vox_lw     = [2.5, 2.5, 1.5, 1.5]
vox_labels = ["v1 (front)","v2 (front)","v3","v4"]

function plot_means!(ax, times, means, title_str, show_legend)
    ax.title = title_str
    hlines!(ax, [Float64(n_lo), Float64(n_hi)]; color=(C_GRAY,0.3), linewidth=0.8, linestyle=:dot)
    text!(ax, 0.02, n_lo+4; text="n_lo", color=C_GRAY, fontsize=9)
    text!(ax, 0.02, n_hi+4; text="n_hi", color=C_GRAY, fontsize=9)
    for k in 1:K
        lines!(ax, times, [m[k] for m in means];
               color=vox_colors[k], linewidth=vox_lw[k],
               linestyle=vox_styles[k], label=vox_labels[k])
    end
    show_legend && axislegend(ax; position=:rc, framevisible=true, labelsize=9)
    ylims!(ax, 0, n_hi+25)
    ax.xlabel = "time  t"
    ax.ylabel = "E[nₖ]"
end

ax1 = Axis(fig[1,1]); plot_means!(ax1, times_d, means_d,
    "(a)  Direct K=4 FSP  (ground truth)", true)

ax2 = Axis(fig[1,2]); plot_means!(ax2, times_c, means_c,
    "(b)  Coarse Binomial K=2  [fast, but WRONG]", false)
# Highlight the wrong behavior: v1,v2 should stay at n_lo,n_hi but collapse to mean
text!(ax2, 0.35, 0.6; text="front collapses:\nv1=v2≈$(round(Int,(n_lo+n_hi)/2))\n(wrong!)",
      color=C_COARSE, fontsize=9, space=:relative, align=(:left,:center))

ax3 = Axis(fig[1,3]); plot_means!(ax3, times_v, means_v,
    "(c)  V-cycle dynamic-π  [fast, correct]", false)
text!(ax3, 0.35, 0.65; text="front preserved:\nv1≈n_lo, v2≈n_hi\n✓",
      color=C_VCYCLE, fontsize=9, space=:relative, align=(:left,:center))

# State space annotation
sz_d = length(sp4); sz_c = length(sp_c); sz_v = length(sp_vc)
text!(ax1, 0.98, 0.05; text="|S|=$(sz_d)", color=C_DIRECT, fontsize=9,
      space=:relative, align=(:right,:bottom))
text!(ax2, 0.98, 0.05; text="|S|=$(sz_c)", color=C_COARSE, fontsize=9,
      space=:relative, align=(:right,:bottom))
text!(ax3, 0.98, 0.05; text="|S|=$(sz_v)", color=C_VCYCLE, fontsize=9,
      space=:relative, align=(:right,:bottom))

# ── Row 2: spatial marginal heatmaps at t=T_HEAT ─────────────────────────────
M_d = !isnothing(snap_d) ? voxel_marginals(snap_d..., K, n_cap) : nothing
M_c = !isnothing(snap_c) ? coarse_marginals_binom(snap_c..., 2, n_cap) : nothing
M_v = !isnothing(snap_v) ? voxel_marginals(snap_v..., K, n_cap) : nothing

gmax = maximum(filter(!isnan, [
    isnothing(M_d) ? 0.0 : maximum(M_d),
    isnothing(M_v) ? 0.0 : maximum(M_v)]))

vox_ticks = (1:K, ["v$k" for k in 1:K])

for (col, M, ttl) in [
        (1, M_d, "(d)  Direct  [t=$(T_HEAT)]"),
        (2, M_c, "(e)  Coarse Binom  [t=$(T_HEAT)]"),
        (3, M_v, "(f)  V-cycle dyn-π  [t=$(T_HEAT)]")]
    isnothing(M) && continue
    ax = Axis(fig[2,col]; title=ttl,
              xlabel="voxel  k", ylabel=col==1 ? "molecules  n" : "",
              xticks=vox_ticks)
    cm = col==2 ? gmax : gmax   # same scale for all
    heatmap!(ax, 0.5:(K+0.5), -0.5:(n_cap+0.5), M;
             colormap=:inferno, colorrange=(0, max(gmax*0.8, 1e-10)))
    hlines!(ax, [Float64(n_lo), Float64(n_hi)];
            color=(:white,0.35), linewidth=0.8, linestyle=:dot)
    if col == 2
        text!(ax, 0.5, n_lo+4; text="← should be\nhere (n_lo)", color=C_COARSE,
              fontsize=8, align=(:center,:bottom))
        text!(ax, 1.5, n_hi-4; text="should be\nhere (n_hi) →", color=C_COARSE,
              fontsize=8, align=(:center,:top))
    end
    ylims!(ax, 0, n_cap)
end

save_fig("fig_talk_schlogl", fig)
println("Done.")
