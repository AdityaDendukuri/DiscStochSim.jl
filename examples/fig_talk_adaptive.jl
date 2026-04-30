"""
fig_talk_adaptive.jl  --  K=8 Schlogl adaptive mixed-resolution

IC = (n_lo x3, n_hi x5): front in pair 2 (voxels 3-4).
halo=0: only pair 2 stays fine; pairs 1,3,4 merged -> KM=5 vs K=8.

Three panels:
  (a) Voxel means over time
  (b) Pair-2 joint P(n3,n4) at t=0.5  -- bimodal, must keep fine
  (c) Pair-4 total P(nc4) at t=0.5    -- Binomial approx, safe to merge
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using Printf
using Distributions

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)
function save_fig(name, fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=3)
    println("  Saved: $name.{pdf,png}")
end

const C_FINE   = RGBf(0.835, 0.369, 0.000)
const C_COARSE = RGBf(0.208, 0.475, 0.694)
const C_GRAY   = RGBf(0.5, 0.5, 0.5)

set_theme!(Theme(
    fontsize = 12,
    Axis = (spinewidth=0.8, xgridvisible=false, ygridvisible=false,
            xticklabelsize=10, yticklabelsize=10, xlabelsize=11, ylabelsize=11,
            titlesize=11),
))

# model
D = 0.005; K = 8; h = 1.0/K; n_max = 180
model     = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)
rates     = Float64[]

fp = schlogl_fixed_points(model; n_max=250)
n_lo, n_un, n_hi = fp
@printf("Fixed points: n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

u0 = CartesianIndex(n_lo, n_lo, n_lo, n_hi, n_hi, n_hi, n_hi, n_hi)
println("IC: (n_lo x3, n_hi x5) -- front in pair 2 (voxels 3-4)")

rdme_sys      = build_schlogl_rdme_system(model, fine_grid)
bc_fine       = s -> rdme_bc(s, n_max)
rdme_diffonly = RDMEModel1D(model.D, 0.0, 0.0)

sp_fine = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp_fine, u0, 1.0)
expand!(sp_fine, rdme_sys, bc_fine; depth=1)

levels = select_coarsening_mask(sp_fine, rdme_diffonly, fine_grid; halo=0)
KM     = _mixed_dim(levels)
@printf("Initial levels: %s  (KM=%d vs K=%d)\n", string(levels), KM, K)

sp_mixed = partial_restrict(sp_fine, levels)
println("Initial |S_mixed|=$(length(sp_mixed.states))")

# time stepping in mixed space
dt = 0.2; T = 3.0; n_steps = round(Int, T/dt)
t_snaps = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

cache      = MixedSolverCache()
step_count = Ref(0)
wnd_center = Ref(-1.0)
mixed_snaps = Dict{Float64, Any}()

function compute_means_mixed(sp_m, lv, K)
    n_pairs = K ÷ 2
    mu  = zeros(K)
    pos = 1
    for j in 1:n_pairs
        for (i, s) in enumerate(sp_m.states)
            tv = Tuple(s); p = sp_m.probs[i]
            if lv[j] == 0
                mu[2j-1] += p * tv[pos]
                mu[2j]   += p * tv[pos+1]
            else
                mu[2j-1] += p * tv[pos] / 2
                mu[2j]   += p * tv[pos] / 2
            end
        end
        pos += lv[j] == 0 ? 2 : 1
    end
    mu
end

times     = Float64[0.0]
all_means = Vector{Float64}[compute_means_mixed(sp_mixed, levels, K)]

t_cur = 0.0
for step in 1:n_steps
    global sp_mixed, levels, t_cur
    dt_step = min(dt, T - t_cur)
    sp_mixed, levels = two_level_step_mixed(
        sp_mixed, levels, model, fine_grid, dt_step;
        cache=cache, mask_check_interval=1,
        step_count=step_count, halo=0,
        krylov_m=30, mixed_n_max=2*n_max,
        binom_tol=1e-6, prune_tol=1e-8, reexpand_depth=1,
        prev_window_center=wnd_center)
    t_cur += dt_step
    renormalize!(sp_mixed)

    push!(times, t_cur)
    push!(all_means, compute_means_mixed(sp_mixed, levels, K))

    snap_t = t_snaps[argmin(abs.(t_snaps .- t_cur))]
    if abs(t_cur - snap_t) < dt/2 && !haskey(mixed_snaps, snap_t)
        mixed_snaps[snap_t] = (copy(sp_mixed.states), copy(sp_mixed.probs), copy(levels))
        @printf("  t=%.1f  |S_mixed|=%d  levels=%s\n",
                t_cur, length(sp_mixed.states), string(levels))
    end
end
println("Done.")

# helpers
function mixed_offset(lv, pair_idx)
    pos = 1
    for j in 1:pair_idx-1; pos += lv[j] == 0 ? 2 : 1; end
    pos
end

function pair_joint_from_mixed(sts, prs, lv, j)
    @assert lv[j] == 0 "pair $j is not fine"
    op = mixed_offset(lv, j)
    d = Dict{Tuple{Int,Int},Float64}()
    for (i,s) in enumerate(sts)
        tv = Tuple(s)
        key = (tv[op], tv[op+1])
        d[key] = get(d, key, 0.0) + prs[i]
    end
    n1s = [k[1] for k in keys(d)]; n2s = [k[2] for k in keys(d)]
    lo1,hi1 = minimum(n1s),maximum(n1s); lo2,hi2 = minimum(n2s),maximum(n2s)
    M = zeros(hi1-lo1+1, hi2-lo2+1)
    for ((a,b),p) in d; M[a-lo1+1,b-lo2+1] += p; end
    M, lo1, lo2
end

function pair_total_from_mixed(sts, prs, lv, j)
    op = mixed_offset(lv, j)
    d = Dict{Int,Float64}()
    for (i,s) in enumerate(sts)
        tv = Tuple(s)
        nc = lv[j] == 0 ? tv[op]+tv[op+1] : tv[op]
        d[nc] = get(d, nc, 0.0) + prs[i]
    end
    d
end

# figure
t_ref = 1.0   # use t=1 for more spread
fig   = Figure(size=(950, 310))

ax1 = Axis(fig[1,1]; xlabel="time  t", ylabel="E[nk] (molecules)",
           title="(a)  Voxel means  [K=8 adaptive, halo=0]")
means_mat = hcat(all_means...)
n_pairs   = K ÷ 2
colors = [C_COARSE,C_COARSE, C_FINE,C_FINE, C_COARSE,C_COARSE, C_COARSE,C_COARSE]
styles = [:dash,:dash, :solid,:solid, :dot,:dot, :dashdot,:dashdot]
lws    = [1.5,1.5, 2.5,2.5, 1.5,1.5, 1.5,1.5]
lbls   = ["v1","v2","v3 (fine)","v4 (fine)","v5","v6","v7","v8"]
for k in 1:K
    lines!(ax1, times, means_mat[k,:]; color=colors[k],
           linewidth=lws[k], linestyle=styles[k], label=lbls[k])
end
hlines!(ax1, [Float64(n_lo), Float64(n_hi)]; color=(C_GRAY,0.4), linewidth=0.8, linestyle=:dot)
text!(ax1, 0.05, n_lo+4; text="n_lo", color=C_GRAY, fontsize=9)
text!(ax1, 0.05, n_hi+4; text="n_hi", color=C_GRAY, fontsize=9)
axislegend(ax1; position=:rc, framevisible=true, labelsize=8, nbanks=2)
ylims!(ax1, 0, n_hi+25)

# (b) pair-2 joint
if haskey(mixed_snaps, t_ref)
    sts_m, prs_m, lv_m = mixed_snaps[t_ref]
    if lv_m[2] == 0
        M, lo1, lo2 = pair_joint_from_mixed(sts_m, prs_m, lv_m, 2)
        ax2 = Axis(fig[1,2]; xlabel="n3 (voxel 3)", ylabel="n4 (voxel 4)",
                   title="(b)  Pair 2 joint at t=$(t_ref)  -- keep FINE")
        heatmap!(ax2, lo1:lo1+size(M,1)-1, lo2:lo2+size(M,2)-1, M;
                 colormap=:inferno, colorrange=(0, maximum(M)*0.95))
        scatter!(ax2, [n_lo, n_hi], [n_hi, n_lo]; color=:white,
                 marker=:circle, markersize=8, strokewidth=1.2, strokecolor=:white)
        text!(ax2, 0.05, 0.97; text="Bimodal -- cannot\nbe coarsened",
              color=:white, fontsize=9, space=:relative, align=(:left,:top))
    else
        Axis(fig[1,2]; title="(b) pair 2 is coarse")
    end
else
    Axis(fig[1,2]; title="(b) no snapshot at t=$(t_ref)")
end

# (c) pair-4 total
if haskey(mixed_snaps, t_ref)
    sts_m, prs_m, lv_m = mixed_snaps[t_ref]
    d4 = pair_total_from_mixed(sts_m, prs_m, lv_m, 4)
    nc_vals = sort(collect(keys(d4)))
    p_vals  = [d4[n] for n in nc_vals]
    ax3 = Axis(fig[1,3]; xlabel="nc4 = n7+n8", ylabel="P(nc4)",
               title="(c)  Pair 4 total at t=$(t_ref)  -- safe to MERGE")
    barplot!(ax3, nc_vals, p_vals; color=(C_COARSE,0.7), strokewidth=0, gap=0)
    nc_mode = nc_vals[argmax(p_vals)]
    if nc_mode > 0
        binom   = Binomial(nc_mode, 0.5)
        nr      = max(0,nc_mode-40):nc_mode+40
        b_probs = pdf.(binom, collect(nr))
        p_nc    = get(d4, nc_mode, 0.0)
        scale   = p_nc / max(maximum(b_probs), 1e-20)
        lines!(ax3, collect(nr), b_probs .* scale;
               color=(C_COARSE,0.5), linewidth=1.5, linestyle=:dash, label="Binom(nc,1/2)")
        axislegend(ax3; position=:rt, framevisible=false, labelsize=9)
    end
    text!(ax3, 0.05, 0.97; text="Approx Binomial --\nsafe to merge",
          color=C_COARSE, fontsize=9, space=:relative, align=(:left,:top))
else
    Axis(fig[1,3]; title="(c) no snapshot")
end

save_fig("fig_talk_adaptive", fig)
println("All done.")
