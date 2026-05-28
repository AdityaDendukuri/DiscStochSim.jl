"""
K=32 Schlögl 1D persistent front with adaptive algebraic coarsening.

Adaptive CoarseningMap merges stable bulk regions into super-voxels,
keeping only the narrow front (2–4 voxels) at fine resolution.

The K_eff-voxel coarsened FSP uses build_schlogl_adaptive_system with
heterogeneous patch sizes and correct boundary-molecule diffusion rates.

Run: julia --project examples/fig_k32_adaptive.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities

# ── model & grid ──────────────────────────────────────────────────────────────
K   = 32
D   = 0.005
h   = 1.0 / 8          # same voxel size as the K=8 reference case
model     = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)

n_low, n_uns, n_high = schlogl_fixed_points(model)
println("K=$K  D=$D  h=$h")
println("Fixed points: lo=$n_low  uns=$n_uns  hi=$n_high")
println("Fine diffusion rate d = $(round(diffusion_rate(D, fine_grid), digits=3))")

# ── IC: left half hi, right half lo ──────────────────────────────────────────
half  = K ÷ 2
u0    = CartesianIndex(ntuple(k -> k ≤ half ? n_high : n_low, Val(K))...)
means_ic = Float64[k ≤ half ? n_high : n_low for k in 1:K]

# ── build adaptive CoarseningMap from IC means ────────────────────────────────
cmap = build_adaptive_cmap_1d(means_ic, n_low, n_high;
                               θ_lo = 25.0, θ_hi = 25.0, min_front = 1)

K_eff = cmap.n_coarse
println("\nAdaptive coarsening: K=$K → K_eff=$K_eff")
for j in 1:K_eff
    patch = cmap.coarse_to_fine[j]
    label = means_ic[patch[1]] > n_high - 30 ? "hi-bulk" :
            means_ic[patch[1]] < n_low  + 30 ? "lo-bulk" : "front"
    println("  region $j (M=$(length(patch))): voxels $(patch[1])–$(patch[end])  [$label]")
end

# ── build adaptive coarsened system ──────────────────────────────────────────
sys_c = build_schlogl_adaptive_system(model, fine_grid, cmap, Val(K_eff))

# BC: returns true for valid states (all components within range)
n_max_fine = 250
bc_c = let cmap=cmap, n_max_fine=n_max_fine
    s -> all(0 ≤ Tuple(s)[j] ≤ length(cmap.coarse_to_fine[j]) * n_max_fine
             for j in 1:cmap.n_coarse)
end

# ── restrict IC to coarse ─────────────────────────────────────────────────────
u0_tup = Tuple(u0)
u0_c_tup = ntuple(j -> sum(u0_tup[k] for k in cmap.coarse_to_fine[j]), Val(K_eff))
u0_c = CartesianIndex(u0_c_tup)

sp_c = StateSpace{CartesianIndex{K_eff}, Float64}()
add_state!(sp_c, u0_c, 1.0)
expand!(sp_c, sys_c, bc_c; depth = 2)

println("\nInitial coarse state: $u0_c_tup")
println("Initial |S_coarse| = $(length(sp_c))")

# ── run adaptive coarse FSP ───────────────────────────────────────────────────
dt        = 2.0
t_final   = 30.0
save_every = 3
prune_tol  = 1e-8

saved_t  = Float64[0.0]
saved_sp = [deepcopy(sp_c)]
ss_sizes = Int[length(sp_c)]

println("Running K_eff=$K_eff adaptive coarse FSP ...")
@time let t_cur = 0.0, step_cur = 0
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

# ── per-voxel marginals at final time ─────────────────────────────────────────
all_marg = [fine_marginals_from_coarse(sp, cmap, n_max_fine) for sp in saved_sp]

# basin masses for v1 (hi) and v$(half+1) (lo)
mf_v1  = all_marg[end][1]
mf_vlo = all_marg[end][half+1]
hi1  = sum(mf_v1[n_high-30:min(end,n_high+30)])
lo_m = sum(mf_vlo[n_low-15:min(end,n_low+15)])
println("v1  hi-peak = $(round(hi1,digits=3))")
println("v$(half+1) lo-peak = $(round(lo_m,digits=3))")

# ── figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(980, 320))

# panel (a): per-voxel mean profile at several times (spatial heatmap-like)
ax_a = Axis(fig[1,1];
    xlabel = "voxel  k",
    ylabel = "E[nₖ]",
    title  = "(a)  Per-voxel mean:  K=$K adaptive coarse",
    titlesize = 11,
    xticks = (1:4:K, string.(1:4:K)),
)

colors_t = Makie.resample_cmap(:Blues, length(saved_t)+1)[2:end]
for (ti, t_s) in enumerate(saved_t)
    means_t = [sum((0:n_max_fine) .* all_marg[ti][k]) for k in 1:K]
    lines!(ax_a, 1:K, means_t; color=colors_t[ti], linewidth=1.4,
           label="t=$(round(Int,t_s))")
end
hlines!(ax_a, [Float64(n_low)];  color=(:steelblue,0.5), linestyle=:dash, linewidth=1.0)
hlines!(ax_a, [Float64(n_high)]; color=(:firebrick,0.5), linestyle=:dash, linewidth=1.0)
vlines!(ax_a, [Float64(half)+0.5]; color=(:gray50,0.4), linestyle=:dot, linewidth=1.0)
text!(ax_a, 2.0,         n_high*0.95; text="n_hi=$n_high", fontsize=8, color=:firebrick4)
text!(ax_a, half+2.0,    n_low*1.15;  text="n_lo=$n_low",  fontsize=8, color=:steelblue4)
ylims!(ax_a, 0, n_high * 1.1)
axislegend(ax_a; position=:rc, labelsize=8, framevisible=false, nbanks=3)

# panel (b): state-space size over time
ax_b = Axis(fig[1,2];
    xlabel = "time  t",
    ylabel = "|S|  (coarse states)",
    title  = "(b)  State-space: K=$K, K_eff=$K_eff adaptive",
    titlesize = 11,
)

lines!(ax_b,  saved_t, Float64.(ss_sizes); color=:steelblue4, linewidth=1.6,
       label="K=$K (K_eff=$K_eff)")
scatter!(ax_b, saved_t, Float64.(ss_sizes); color=:steelblue4, markersize=5)

# panel (c): coarsening map visualisation
ax_c = Axis(fig[1,3];
    xlabel    = "voxel  k",
    ylabel    = "region  j",
    title     = "(c)  Adaptive CoarseningMap  (K_eff=$K_eff)",
    titlesize = 11,
    yticks    = (1:K_eff, ["r$j" for j in 1:K_eff]),
    xticks    = (1:4:K, string.(1:4:K)),
)

for j in 1:K_eff
    patch = cmap.coarse_to_fine[j]
    Mj    = length(patch)
    col   = Mj == 1 ? :orange : (means_ic[patch[1]] > n_high - 30 ? :firebrick : :steelblue)
    for k in patch
        poly!(ax_c, Rect(k-0.5, j-0.45, 1.0, 0.9); color=(col, 0.4), strokewidth=0.5)
    end
    text!(ax_c, Float64(patch[1] + (Mj-1)/2), Float64(j);
          text= Mj == 1 ? "v$( patch[1])" : "M=$Mj",
          fontsize=8, align=(:center,:center),
          color= Mj==1 ? :darkorange : :black)
end
xlims!(ax_c, 0.5, K+0.5)
ylims!(ax_c, 0.5, K_eff+0.5)

# ── K=64 run for comparison ────────────────────────────────────────────────────
println("\n── K=64 run ──")
K64 = 64; half64 = K64 ÷ 2
fine_grid64 = VoxelGrid(K64, h, 0)
u0_64       = CartesianIndex(ntuple(k -> k ≤ half64 ? n_high : n_low, Val(K64))...)
means_ic64  = Float64[k ≤ half64 ? n_high : n_low for k in 1:K64]

cmap64 = build_adaptive_cmap_1d(means_ic64, n_low, n_high;
                                 θ_lo = 25.0, θ_hi = 25.0, min_front = 1)
K_eff64 = cmap64.n_coarse
println("K=$K64 → K_eff=$K_eff64  patches: $(map(p->length(p), cmap64.coarse_to_fine))")

sys_c64 = build_schlogl_adaptive_system(model, fine_grid64, cmap64, Val(K_eff64))
bc_c64  = s -> all(0 ≤ Tuple(s)[j] ≤ length(cmap64.coarse_to_fine[j]) * n_max_fine
                   for j in 1:cmap64.n_coarse)

u0_c64_tup = ntuple(j -> sum(Tuple(u0_64)[k] for k in cmap64.coarse_to_fine[j]), Val(K_eff64))
u0_c64 = CartesianIndex(u0_c64_tup)
sp_c64 = StateSpace{CartesianIndex{K_eff64}, Float64}()
add_state!(sp_c64, u0_c64, 1.0)
expand!(sp_c64, sys_c64, bc_c64; depth = 2)

saved_t64  = Float64[0.0]
ss_sizes64 = Int[length(sp_c64)]

@time let t_cur = 0.0, step_cur = 0
    while t_cur < t_final
        dt_step = min(dt, t_final - t_cur)
        expand!(sp_c64, sys_c64, bc_c64; depth = 1)
        A_c64, = build_generator(sp_c64, sys_c64, Float64[], t_cur)
        sp_c64.probs .= expv(Float64(dt_step), A_c64, sp_c64.probs; m = 30)
        t_cur += dt_step; step_cur += 1
        renormalize!(sp_c64)
        prune_threshold!(sp_c64, prune_tol)
        if step_cur % save_every == 0 || t_cur ≥ t_final
            push!(saved_t64, t_cur)
            push!(ss_sizes64, length(sp_c64))
        end
    end
end
println("K=64 done. peak |S|=$(maximum(ss_sizes64))")

# ── update state-space panel with K=32 vs K=64 ────────────────────────────────
lines!(ax_b, saved_t64, Float64.(ss_sizes64);
       color=:firebrick4, linewidth=1.6, linestyle=:dash,
       label="K=$K64 (K_eff=$K_eff64)")
scatter!(ax_b, saved_t64, Float64.(ss_sizes64); color=:firebrick4, markersize=5)
axislegend(ax_b; position=:lt, labelsize=9, framevisible=false)
text!(ax_b, saved_t[end]*0.05, maximum(ss_sizes64)*1.08;
      text="K=$K64: $(maximum(ss_sizes64)) states", fontsize=9, color=:firebrick4)

# Re-draw K=32 label
text!(ax_b, saved_t[end]*0.55, maximum(ss_sizes)*1.12;
      text="K=$K: $(maximum(ss_sizes)) states", fontsize=9, color=:steelblue4)
axislegend(ax_b; position=:lt, labelsize=9, framevisible=false)

colgap!(fig.layout, 1, 20); colgap!(fig.layout, 2, 20)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_k32_adaptive.pdf")
save(out, fig)
println("\nSaved → $out")
fig
