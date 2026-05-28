"""
K=8 Schlögl with adaptive coarsening K_eff=2.

Two super-voxels: left 4 voxels merged (hi state), right 4 merged (lo state).
The K=2 coarse FSP captures the bistable front dynamics with a tiny state space.
Per-voxel K=8 marginals are reconstructed analytically via Binomial splitting.

Key: corrected inter-patch diffusion uses D_c = M×D so that
     d_c = D_c/(M×h)² = M×D/(M²×h²) = D/(M×h²) = d/M  (boundary-molecule rate).

Run: julia --project examples/fig_k8_adaptive.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities

# ── model & coarsening ────────────────────────────────────────────────────────
K = 8;  D = 5.0;  h = 1.0
M = K ÷ 2        # patch size: each super-voxel contains M=4 fine voxels

model = SchloglModel1D(D)
n_low, n_uns, n_high = schlogl_fixed_points(model)
println("Fine: lo=$n_low  uns=$n_uns  hi=$n_high")

# Coarsened model for M=4 patch:
#   reactions: Binomial-corrected (coarsen_model by r=M)
#   diffusion: D_c = M×D  so that d_c = D_c/(Mh)² = D/(M h²) = d/M
model_rxn = coarsen_model(model, Float64(M))   # c1/M, c2/M², c3*M, c4, D (unchanged)
model_c   = SchloglModel1D(M * model_rxn.D,    # correct D_c = M×D
                            model_rxn.c1, model_rxn.c2,
                            model_rxn.c3, model_rxn.c4)
grid_c    = VoxelGrid(2, h * M, 0)             # K=2, dx=M×h

n_low_c, n_uns_c, n_high_c = schlogl_fixed_points(model_c; n_max = M * 400)
println("Coarse (M=$M): lo=$n_low_c  uns=$n_uns_c  hi=$n_high_c")
println("  ≈ M×fine:   lo=$(M*n_low)  uns=$(M*n_uns)  hi=$(M*n_high)")

# CoarseningMap: {v1..v4} → coarse 1,  {v5..v8} → coarse 2
cmap = CoarseningMapRegions([fill(1,M); fill(2,M)])

# ── K=2 coarse FSP ───────────────────────────────────────────────────────────
sys_c = build_schlogl_rdme_system(model_c, grid_c)
n_max_c = M * 300
bc_c    = s -> all(0 ≤ Tuple(s)[j] ≤ length(cmap.coarse_to_fine[j]) * n_max_fine
                   for j in 1:cmap.n_coarse)

# IC: left half all hi → n̄_L = M×n_hi,  right half all lo → n̄_R = M×n_lo
u0_c = CartesianIndex(M * n_high, M * n_low)

sp_c = StateSpace{CartesianIndex{2}, Float64}()
add_state!(sp_c, u0_c, 1.0)
expand!(sp_c, sys_c, bc_c; depth = 2)

dt        = 1.0
t_final   = 20.0
save_every = 4
prune_tol  = 1e-10

saved_t  = Float64[0.0]
saved_sp = [deepcopy(sp_c)]
ss_sizes = Int[length(sp_c)]

println("Running K=2 coarse FSP (IC: n̄_L=$(M*n_high), n̄_R=$(M*n_low)) ...")
A_c, = build_generator(sp_c, sys_c, Float64[], 0.0)

let t_cur = 0.0, step_cur = 0
    while t_cur < t_final
        dt_step = min(dt, t_final - t_cur)
        expand!(sp_c, sys_c, bc_c; depth = 1)
        A_c, = build_generator(sp_c, sys_c, Float64[], t_cur)
        sp_c.probs .= expv(Float64(dt_step), A_c, sp_c.probs; m = 30)
        t_cur += dt_step; step_cur += 1
        renormalize!(sp_c)
        prune_threshold!(sp_c, prune_tol)
        if step_cur % save_every == 0 || t_cur ≥ t_final
            push!(saved_t, t_cur)
            push!(saved_sp, deepcopy(sp_c))
            push!(ss_sizes, length(sp_c))
        end
    end
end

println("Done. $(length(saved_t)) snapshots  peak |S_coarse|=$(maximum(ss_sizes))")

# ── compute K=8 per-voxel marginals via Binomial reconstruction ───────────────
n_max_fine = 250
all_marg = [fine_marginals_from_coarse(sp, cmap, n_max_fine) for sp in saved_sp]

# Summary: basin masses for v1 (hi-side) and v5 (lo-side) at final time
function basin(marg, vox)
    P = marg[vox]
    sum(P[n_high-30:min(end,n_high+30)]), sum(P[n_low-20:min(end,n_low+20)])
end
hi1, lo1 = basin(all_marg[end], 1)
hi5, lo5 = basin(all_marg[end], 5)
println("v1 at t=$(round(Int,t_final)): hi-peak=$(round(hi1,digits=3))  lo-peak=$(round(lo1,digits=3))")
println("v5 at t=$(round(Int,t_final)): hi-peak=$(round(hi5,digits=3))  lo-peak=$(round(lo5,digits=3))")

# ── figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(960, 330))

# panel (a): coarse super-voxel marginals at final time
ax_a = Axis(fig[1,1];
    xlabel="total count  n̄",
    ylabel="P(n̄)",
    title="(a)  K=2 coarse marginals at t=$(round(Int,t_final))",
    titlesize=11)

sp_f = saved_sp[end]
P_L = zeros(n_max_c+1); P_R = zeros(n_max_c+1)
for (s,p) in zip(sp_f.states, sp_f.probs)
    nL, nR = Tuple(s)
    nL <= n_max_c && (P_L[nL+1] += p)
    nR <= n_max_c && (P_R[nR+1] += p)
end

# plot only the region with mass
rng = 50:min(n_max_c, M*220)
lines!(ax_a, collect(rng), P_L[rng.+1]; color=:firebrick4, linewidth=1.6, label="n̄_L (hi-side)")
lines!(ax_a, collect(rng), P_R[rng.+1]; color=:steelblue4, linewidth=1.6, label="n̄_R (lo-side)")
vlines!(ax_a, [n_high_c]; color=(:firebrick,0.5), linestyle=:dash, linewidth=1.0)
vlines!(ax_a, [n_low_c];  color=(:steelblue,0.5), linestyle=:dash, linewidth=1.0)
text!(ax_a, n_high_c+5, maximum(P_L[rng.+1])*0.95; text="n̄_hi≈$n_high_c", fontsize=8, color=:firebrick4)
text!(ax_a, n_low_c+5,  maximum(P_R[rng.+1])*0.95; text="n̄_lo≈$n_low_c",  fontsize=8, color=:steelblue4)
axislegend(ax_a; position=:rt, labelsize=9, framevisible=false)

# panel (b): all K=8 per-voxel marginals at final time
ax_b = Axis(fig[1,2];
    xlabel="molecules  n",
    ylabel="P(nₖ = n)",
    title="(b)  K=8 per-voxel marginals  (Binomial reconstruction)",
    titlesize=11)

colors_L = Makie.resample_cmap(Reverse(:Reds),   M+2)[2:end-1]
colors_R = Makie.resample_cmap(Reverse(:Blues),  M+2)[2:end-1]

for k in 1:M
    P = all_marg[end][k]
    lines!(ax_b, 0:n_max_fine, P; color=colors_L[k], linewidth=1.3,
           label= k==1 ? "v1–v4 (hi)" : "")
end
for k in (M+1):K
    P = all_marg[end][k]
    lines!(ax_b, 0:n_max_fine, P; color=colors_R[k-M], linewidth=1.3,
           label= k==M+1 ? "v5–v8 (lo)" : "")
end
vlines!(ax_b, [n_low];  color=(:steelblue,0.5), linestyle=:dash, linewidth=1.0)
vlines!(ax_b, [n_high]; color=(:firebrick,0.5), linestyle=:dash, linewidth=1.0)
text!(ax_b, n_low+2,  maximum(all_marg[end][1])*1.05; text="n_lo=$n_low",  fontsize=8, color=:steelblue4)
text!(ax_b, n_high+2, maximum(all_marg[end][1])*1.05; text="n_hi=$n_high", fontsize=8, color=:firebrick4)
xlims!(ax_b, 0, n_max_fine)
axislegend(ax_b; position=:rt, labelsize=9, framevisible=false)

# panel (c): state-space size over time
ax_c = Axis(fig[1,3];
    xlabel="time  t",
    ylabel="|S|  (coarse states)",
    title="(c)  K=2 coarse state-space  (K=$K voxels)",
    titlesize=11)

lines!(ax_c,  saved_t, Float64.(ss_sizes); color=:steelblue4, linewidth=1.6)
scatter!(ax_c, saved_t, Float64.(ss_sizes); color=:steelblue4, markersize=5)
text!(ax_c, saved_t[end]*0.05, maximum(ss_sizes)*1.1;
      text="peak: $(maximum(ss_sizes)) states\n(K=8 direct: ~250⁸≈4×10¹⁸)",
      fontsize=9, color=:steelblue4)

colgap!(fig.layout, 1, 20); colgap!(fig.layout, 2, 20)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_k8_adaptive.pdf")
save(out, fig)
println("\nSaved → $out")
fig
