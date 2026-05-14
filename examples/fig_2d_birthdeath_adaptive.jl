"""
2D birth-death RDME point-source wavefront via adaptive column-product coarsening.

IC: N molecules at corner voxel (1,1). Front spreads outward.
Column reduction: n̄_j = Σᵢ nᵢⱼ, K_y·λ production, same μ and D.
Adaptive coarsening: merge settled (n̄≈K_y·n_ss) and empty (n̄≈0) columns.
Key: K_eff grows as the wavefront spreads, then shrinks as columns settle.

θ_lo must be < N_ic to classify the IC column as active (not empty).

Run: julia --project examples/fig_2d_birthdeath_adaptive.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities

# ── parameters ────────────────────────────────────────────────────────────────
K_x = 20;  K_y = 20          # 20×20 square grid
D   = 1.0;  h = 1.0
λ   = 0.05; μ = 0.1        # n_ss = λ/μ = 0.5; n_ss_col = K_y*n_ss = 10 (small → tractable)
N_ic = 50                   # >> n_ss_col=10: col 1 starts settled; IC excess diffuses outward
t_final = 3 / μ             # 3 settling time constants = 30

n_ss     = λ / μ            # 0.5
n_ss_col = Float64(K_y) * n_ss   # 10

# Global birth: all voxels eventually settle to n_ss_col → K_eff→1 at SS.
# max_front caps active regions so K_eff ≤ max_front+2 during the transient.
θ_lo_abs  = 1.0             # empty: col sum < 1 molecule (absolute, independent of D/λ)
θ_hi_frac = 0.85            # settled: col sum > θ_hi_frac × n_ss_col
max_front = 2               # K_eff ≤ 4 throughout

println("K=$(K_x)×$(K_y), D=$D, n_ss=$n_ss, n_ss_col=$n_ss_col")
println("N_ic=$N_ic, t_final=$t_final, θ_lo_abs=$θ_lo_abs, θ_hi=$(round(θ_hi_frac*n_ss_col,digits=1))")

# ── 2D mean-field ODE — global birth everywhere ──────────────────────────────
function run_mf(K_x, K_y, λ, μ, D, h, N_ic, t_final, dt_mf, n_snap)
    μg = zeros(K_x, K_y);  μg[1,1] = Float64(N_ic)
    dμ = similar(μg);  d = D/h^2
    ts = range(0.0, t_final, length=n_snap)
    snaps = [copy(μg)];  si = 2;  t = 0.0
    while t < t_final
        fill!(dμ, 0.0)
        for i in 1:K_x, j in 1:K_y
            dμ[i,j] = λ - μ*μg[i,j]   # global birth everywhere
            i>1   && (dμ[i,j] += d*(μg[i-1,j]-μg[i,j]))
            i<K_x && (dμ[i,j] += d*(μg[i+1,j]-μg[i,j]))
            j>1   && (dμ[i,j] += d*(μg[i,j-1]-μg[i,j]))
            j<K_y && (dμ[i,j] += d*(μg[i,j+1]-μg[i,j]))
        end
        μg .+= dt_mf .* dμ;  t += dt_mf
        if si ≤ n_snap && t ≥ ts[si]-dt_mf/2
            push!(snaps, copy(μg)); si += 1
        end
    end
    collect(ts[1:length(snaps)]), snaps
end

println("2D mean-field ODE...")
@time t_mf, mu_mf = run_mf(K_x, K_y, λ, μ, D, h, N_ic, t_final, 0.01, 5)

# ── adaptive 2D FSP ───────────────────────────────────────────────────────────
# State: CartesianIndex{K_eff} where dim J = total molecules in 2D patch J.
# Coarsening: connected components of same-class 2D voxels (settled/active/empty).
fine_grid_2d = Grid2D(K_x, K_y, h, 0)
model_2d     = SchloglModel1D(D, 0.0, 0.0, λ, μ)   # c3=λ per fine voxel, c4=μ

# Mean-field state: K_x×K_y matrix (exact for linear reactions).
mf_state = zeros(Float64, K_x, K_y);  mf_state[1,1] = Float64(N_ic)
mf_d = D / h^2

# 2D coarsening: connected components of same-label voxels.
# θ_lo and θ_hi are per-fine-voxel thresholds (not column sums).
θ_lo_2d = θ_lo_abs / K_y        # per-voxel empty threshold
θ_hi_2d = θ_hi_frac * n_ss      # per-voxel settled threshold (n_ss = λ/μ)

cmap = build_adaptive_cmap_2d(mf_state, K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
K_eff = cmap.n_coarse
println("\nInitial CoarseningMap: $(K_x)×$(K_y) → K_eff=$K_eff")

sys_c = build_rdme_adaptive_2d(model_2d, fine_grid_2d, cmap, Val(K_eff))

# Per-region BC: 4× mean-field mean per patch (safety factor, no distributional assumption).
n_max_c_hard = round(Int, max(N_ic, n_ss_col) * 4)
function _region_bc(cm, mf, n_hard)
    means = [sum(mf[k] for k in cm.coarse_to_fine[r]) for r in 1:cm.n_coarse]
    caps  = [min(n_hard, max(5, round(Int, 4 * max(m, 1.0)))) for m in means]
    let c=cm, cp=caps; s -> all(0 ≤ Tuple(s)[r] ≤ cp[r] for r in 1:c.n_coarse); end
end

bc_c = _region_bc(cmap, mf_state, n_max_c_hard)

u0_tup = ntuple(j -> round(Int, sum(mf_state[k] for k in cmap.coarse_to_fine[j])), Val(K_eff))
u0     = CartesianIndex(u0_tup)
sp_c   = StateSpace{CartesianIndex{K_eff}, Float64}()
add_state!(sp_c, u0, 1.0)
expand!(sp_c, sys_c, bc_c; depth=4)
println("u0=$u0_tup  |S_init|=$(length(sp_c))")

dt_fsp = 0.05; prune_tol = 1e-9;  save_every = 40; rebuild_every = 20

function _mf_step!(μg, K_x, K_y, λ, μ, d, dt)
    dμ = similar(μg)
    fill!(dμ, 0.0)
    for i in 1:K_x, j in 1:K_y
        dμ[i,j] = λ - μ*μg[i,j]   # global birth
        i>1   && (dμ[i,j] += d*(μg[i-1,j]-μg[i,j]))
        i<K_x && (dμ[i,j] += d*(μg[i+1,j]-μg[i,j]))
        j>1   && (dμ[i,j] += d*(μg[i,j-1]-μg[i,j]))
        j<K_y && (dμ[i,j] += d*(μg[i,j+1]-μg[i,j]))
    end
    μg .+= dt .* dμ
end

saved_t      = Float64[0.0]
ss_sizes     = Int[length(sp_c)]
keff_hist    = Int[K_eff]
col_means_hist = [[sum(mf_state[i,j] for j in 1:K_y) for i in 1:K_x]]

# ── function barrier: forces Julia to specialise (and JIT once) per K_eff type ──
# depth=1 and max_states cap prevent |S| from growing unboundedly:
#   expand! is O(|S|×n_rxns) per step; at |S|≈5K that is ~10ms/step → 300 s total.
function _fsp_step!(sp::StateSpace{CartesianIndex{Ke}, T},
                     sys, bc, dt::Float64, prune_tol::Float64,
                     max_states::Int=3000) where {Ke, T}
    length(sp) < max_states && expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, Float64[], 0.0)
    sp.probs .= expv(dt, A, sp.probs; m=30)
    renormalize!(sp)
    prune_threshold!(sp, prune_tol)
    nothing
end

# ── pre-warm JIT for all likely K_eff types ──────────────────────────────────
print("Pre-compiling K_eff=2..$(min(K_x*K_y,8)): ")
for preKe in 2:min(K_x*K_y, 8)
    K_fine = K_x * K_y
    chunk = max(1, K_fine ÷ preKe)
    ctf_w = [collect((j-1)*chunk+1:min(j*chunk,K_fine)) for j in 1:preKe]
    ctf_w[end] = collect((preKe-1)*chunk+1:K_fine)
    ftc_w = zeros(Int, K_fine)
    for (j, ks) in enumerate(ctf_w), k in ks; ftc_w[k] = j; end
    cm_w  = CoarseningMap(K_fine, preKe, maximum(length, ctf_w), ftc_w, ctf_w)
    sys_w = build_rdme_adaptive_2d(model_2d, fine_grid_2d, cm_w, Val(preKe))
    bc_w  = s -> all(0 ≤ Tuple(s)[j] ≤ n_max_c_hard for j in 1:preKe)
    sp_w  = StateSpace{CartesianIndex{preKe}, Float64}()
    add_state!(sp_w, CartesianIndex(ntuple(_ -> round(Int, n_ss), Val(preKe))), 1.0)
    expand!(sp_w, sys_w, bc_w; depth=1)
    A_w, = build_generator(sp_w, sys_w, Float64[], 0.0)
    sp_w.probs .= expv(0.01, A_w, sp_w.probs; m=5)
    prune_threshold!(sp_w, 1e-20)
    print("$preKe ")
end
println("done.")

println("Running adaptive 2D FSP ...")
@time let t=0.0, s=0, sp=sp_c, cm=cmap, sys=sys_c, bc=bc_c, Ke=K_eff
    while t < t_final
        dt_s = min(dt_fsp, t_final-t)
        _fsp_step!(sp, sys, bc, Float64(dt_s), prune_tol)
        _mf_step!(mf_state, K_x, K_y, λ, μ, mf_d, dt_s)
        t += dt_s;  s += 1
        s % 40 == 0 && println("  step $s t=$(round(t,digits=1)) K_eff=$Ke |S|=$(length(sp))")

        if s % rebuild_every == 0
            new_cm = build_adaptive_cmap_2d(mf_state, K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
            new_Ke = new_cm.n_coarse
            new_bc = _region_bc(new_cm, mf_state, n_max_c_hard)
            if new_Ke != Ke
                new_sys   = build_rdme_adaptive_2d(model_2d, fine_grid_2d, new_cm, Val(new_Ke))
                new_means = [sum(mf_state[k] for k in new_cm.coarse_to_fine[j]) for j in 1:new_Ke]
                new_u0    = CartesianIndex(ntuple(j->round(Int,new_means[j]), Val(new_Ke)))
                new_sp    = StateSpace{CartesianIndex{new_Ke}, Float64}()
                add_state!(new_sp, new_u0, 1.0)
                expand!(new_sp, new_sys, new_bc; depth=4)
                println("  t=$(round(t,digits=1)) K_eff: $Ke→$new_Ke  |S|=$(length(new_sp))")
                sp=new_sp; cm=new_cm; sys=new_sys; bc=new_bc; Ke=new_Ke
            else
                bc=new_bc
            end
        end

        if s % save_every == 0 || t ≥ t_final
            push!(saved_t,  t)
            push!(ss_sizes, length(sp))
            push!(keff_hist, cm.n_coarse)
            push!(col_means_hist, [sum(mf_state[i,j] for j in 1:K_y) for i in 1:K_x])
            println("  t=$(round(t,digits=1)) K_eff=$Ke  |S|=$(length(sp))")
        end
    end
end

println("Done. peak |S|=$(maximum(ss_sizes))  K_eff: $(minimum(keff_hist))–$(maximum(keff_hist))")

# ── figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))
fig = Figure(size=(1100, 540))

# top row: 2D mean-field snapshots (circular wavefront from corner)
snap_idx = [1, 2, 3, 5]   # t = 0, 7.5, 15, 30 — straddles the K_eff collapse at t≈19
c_max = max(maximum(maximum.(mu_mf[2:end])), 1e-6)   # exclude t=0 IC spike from colorrange
c_max = min(c_max, n_ss * 3)                          # cap at 3× per-voxel SS for contrast
for (col, si) in enumerate(snap_idx)
    ax = Axis(fig[1,col]; title="t=$(round(t_mf[si],digits=1))", titlesize=10,
              aspect=DataAspect(),
              xlabel=col==2 ? "col j" : "", ylabel=col==1 ? "row i" : "")
    hm = heatmap!(ax, 1:K_x, 1:K_y, mu_mf[si];
                  colormap=:inferno, colorrange=(0, c_max), lowclip=:black)

    # All fine voxel grid lines (faint) — shows individual voxel structure
    for ii in 1:K_x-1
        hlines!(ax, [ii+0.5]; color=(:white, 0.15), linewidth=0.4)
    end
    for jj in 1:K_y-1
        vlines!(ax, [jj+0.5]; color=(:white, 0.15), linewidth=0.4)
    end

    # Region boundary edges (bright + thick) — shows which voxels are merged
    cmap_si = build_adaptive_cmap_2d(mu_mf[si], K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
    for i in 1:K_x, j in 1:K_y
        r = cmap_si.fine_to_coarse[(i-1)*K_y+j]
        if j < K_y && cmap_si.fine_to_coarse[(i-1)*K_y+j+1] != r
            linesegments!(ax, [Point2f(j+0.5,i-0.5), Point2f(j+0.5,i+0.5)];
                          color=(:white,0.95), linewidth=2.0)
        end
        if i < K_x && cmap_si.fine_to_coarse[i*K_y+j] != r
            linesegments!(ax, [Point2f(j-0.5,i+0.5), Point2f(j+0.5,i+0.5)];
                          color=(:white,0.95), linewidth=2.0)
        end
    end
    text!(ax, 0.5, K_y + 0.3; text="K_eff=$(cmap_si.n_coarse)",
          color=:white, fontsize=8, align=(:left, :bottom), space=:data)

    col==4 && Colorbar(fig[1,5], hm; label="E[nᵢⱼ]", width=12,
                        ticks=([0, c_max/2, c_max],
                               ["0", "$(round(c_max/2,digits=1))", "$(round(c_max,digits=1))"]),
                        labelsize=10, ticksize=3)
end
Label(fig[0,1:4],
      "2D birth-death $(K_x)×$(K_y): D=$D, n_ss=$n_ss  |  IC: N=$N_ic at corner (1,1)  " *
      "[mean-field = exact for linear reactions]",
      fontsize=10, tellwidth=false)

# bottom-left: column mean profile over time
ax_p = Axis(fig[2,1:2]; xlabel="column j", ylabel="E[n̄ⱼ]",
            title="Column means (adaptive FSP)", titlesize=10)
clrs = Makie.resample_cmap(:Blues, length(saved_t)+2)[2:end-1]
for (ti, t_s) in enumerate(saved_t)
    lines!(ax_p, 1:K_x, col_means_hist[ti]; color=clrs[ti], linewidth=1.3,
           label="t=$(round(Int,t_s))")
end
hlines!(ax_p,[n_ss_col]; color=(:gray50,0.35), linestyle=:dash, linewidth=0.9)
text!(ax_p, 1.2, n_ss_col*1.02; text="n_ss_col", fontsize=8, color=:gray50)
ylims!(ax_p, 0, n_ss_col*1.15)
axislegend(ax_p; position=:rc, labelsize=8, framevisible=false, nbanks=3)

# bottom-center: K_eff
ax_k = Axis(fig[2,3]; xlabel="time t", ylabel="K_eff",
            title="Active regions over time", titlesize=10)
lines!(ax_k,  saved_t, Float64.(keff_hist); color=:firebrick4, linewidth=1.6)
scatter!(ax_k, saved_t, Float64.(keff_hist); color=:firebrick4, markersize=6)
hlines!(ax_k,[Float64(K_x)]; color=(:gray50,0.35), linestyle=:dash, linewidth=0.9)
text!(ax_k, 0.3, K_x*0.96; text="K_x=$K_x", fontsize=8, color=:gray50)
ylims!(ax_k, 0, K_x+2)

# bottom-right: |S|
ax_s = Axis(fig[2,4]; xlabel="time t", ylabel="|S|",
            title="State-space size", titlesize=10)
lines!(ax_s,  saved_t, Float64.(ss_sizes); color=:steelblue4, linewidth=1.6)
scatter!(ax_s, saved_t, Float64.(ss_sizes); color=:steelblue4, markersize=6)
text!(ax_s, saved_t[end]*0.05, maximum(ss_sizes)*1.1;
      text="peak: $(maximum(ss_sizes))", fontsize=9, color=:steelblue4)

rowgap!(fig.layout, 1, 8); colgap!(fig.layout, 4, 6)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_2d_birthdeath_adaptive.pdf")
save(out, fig)
println("\nSaved → $out")
fig
