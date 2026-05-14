"""
2D birth-death RDME — large grid, fast diffusion, adaptive FSP with remap.

K=40×40 grid, D=10 (fast diffusion relative to kinetics).
IC: N=300 molecules at corner voxel (1,1). Wavefront spreads rapidly.

Key physics: d = D/h² = 10 >> λ=0.1, μ=0.5 → diffusion-dominated.
The wavefront reaches radius sqrt(2Dt)≈14 voxels at t=10,
filling the grid by t≈40. n_ss = λ/μ = 0.2 per voxel.

FSP: CoarseningMap keeps K_eff ≤ 4 (empty / 1-2 front / settled),
state space stays ≤ ~5k states despite 1600-voxel fine grid.
Topology changes use remap_sp_transition (no probability loss).

Run: julia --project examples/fig_2d_bd_large.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using SparseArrays
using Printf

# ── parameters ────────────────────────────────────────────────────────────────
# Tractability constraint: 2·d·N_ic·dt ≪ n_max_region
# so FSP can track one diffusion step.  Here: 2·2·100·0.1 = 40 ≪ 130.
const K_x = 30;  const K_y = 30      # 30×30 = 900 voxels
const D   = 2.0;  const h = 1.0      # 2× faster than D=1; d=2
const λ   = 0.1;  const μ = 0.2      # n_ss = 0.5 per voxel
const N_ic = 100                      # large IC: 200× above n_ss → visible pulse
const t_final = 40.0                  # ~8 settling time constants

n_ss     = λ / μ                   # 0.5
n_ss_col = Float64(K_y) * n_ss    # 15.0

println("Grid: $(K_x)×$(K_y)   D=$D  λ=$λ  μ=$μ")
println("n_ss=$n_ss  n_ss_col=$n_ss_col  d=$(D/h^2)")
println("IC: N=$N_ic at (1,1)   t_final=$t_final")
println("Spread timescale ~ $(round(K_x^2/(2D), digits=1))  equil ~ $(round(1/μ, digits=1))")
println("Max diffusion transfer/step ≈ $(round(2D*N_ic*0.1, digits=1)) molecules")
flush(stdout)

# ── mean-field ODE ────────────────────────────────────────────────────────────
function run_mf(K_x, K_y, λ, μ, D, h, N_ic, t_final, dt_mf, n_snap)
    μg = zeros(K_x, K_y);  μg[1,1] = Float64(N_ic)
    dμ = similar(μg);  d = D/h^2
    ts = range(0.0, t_final; length=n_snap)
    snaps = [copy(μg)];  si = 2;  t = 0.0
    while t < t_final
        fill!(dμ, 0.0)
        for i in 1:K_x, j in 1:K_y
            dμ[i,j] = λ - μ*μg[i,j]
            i>1   && (dμ[i,j] += d*(μg[i-1,j]-μg[i,j]))
            i<K_x && (dμ[i,j] += d*(μg[i+1,j]-μg[i,j]))
            j>1   && (dμ[i,j] += d*(μg[i,j-1]-μg[i,j]))
            j<K_y && (dμ[i,j] += d*(μg[i,j+1]-μg[i,j]))
        end
        μg .+= dt_mf .* dμ;  t += dt_mf
        if si ≤ n_snap && t ≥ ts[si] - dt_mf/2
            push!(snaps, copy(μg));  si += 1
        end
    end
    collect(ts[1:length(snaps)]), snaps
end

println("Mean-field ODE..."); flush(stdout)
@time t_mf, mu_mf = run_mf(K_x, K_y, λ, μ, D, h, N_ic, t_final, 0.005, 12)

# ── adaptive 2D FSP ───────────────────────────────────────────────────────────
fine_grid = Grid2D(K_x, K_y, h, 0)
# Schlögl1D used as 2D birth-death: c3=λ per voxel, c4=μ, D same
model_2d  = SchloglModel1D(D, 0.0, 0.0, λ, μ)

# Empty: below 70% of n_ss → background birth crosses this after ~6 time units,
# keeping the empty region visible during the spreading wavefront phase.
θ_lo_2d = 0.70 * n_ss             # per-voxel empty threshold (= 0.35)
θ_hi_2d = 0.90 * n_ss             # per-voxel settled threshold (= 0.45)

mf_state = zeros(Float64, K_x, K_y);  mf_state[1,1] = Float64(N_ic)
mf_d = D / h^2

cmap  = build_adaptive_cmap_2d(mf_state, K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
K_eff = cmap.n_coarse
println("\nInitial CoarseningMap: $(K_x)×$(K_y) → K_eff=$K_eff")

sys_c = build_rdme_adaptive_2d(model_2d, fine_grid, cmap, Val(K_eff))

# Per-region cap: generous enough to hold all N_ic molecules in one region.
# Constraint: 2·d·N_ic·dt = 2·2·100·0.1 = 40 ≪ n_max_region = 130.
n_max_region = round(Int, N_ic + n_ss_col + 15)   # = 130

function _region_bc(cm, mf, n_max_r)
    # Each region is capped at n_max_r regardless of its current mean.
    # This is conservative but safe: molecules can flow rapidly between regions.
    let c=cm, n=n_max_r; s -> all(0 ≤ Tuple(s)[r] ≤ n for r in 1:c.n_coarse); end
end

bc_c  = _region_bc(cmap, mf_state, n_max_region)
u0    = CartesianIndex(ntuple(j -> round(Int, sum(mf_state[k] for k in cmap.coarse_to_fine[j])), Val(K_eff)))
sp_c  = StateSpace{CartesianIndex{K_eff}, Float64}()
add_state!(sp_c, u0, 1.0)
expand!(sp_c, sys_c, bc_c; depth=3)
println("IC: $(Tuple(u0))  |S_init|=$(length(sp_c))")

# Pre-warm JIT for expected K_eff values (2..5)
print("Pre-compiling K_eff=2..5: "); flush(stdout)
for preKe in 2:5
    K_fine = K_x * K_y;  chunk = max(1, K_fine ÷ preKe)
    ctf_w = [collect((j-1)*chunk+1:min(j*chunk,K_fine)) for j in 1:preKe]
    ctf_w[end] = collect((preKe-1)*chunk+1:K_fine)
    ftc_w = zeros(Int, K_fine)
    for (j, ks) in enumerate(ctf_w), k in ks;  ftc_w[k] = j;  end
    cm_w  = CoarseningMap(K_fine, preKe, maximum(length, ctf_w), ftc_w, ctf_w)
    sys_w = build_rdme_adaptive_2d(model_2d, fine_grid, cm_w, Val(preKe))
    bc_w  = s -> all(0 ≤ Tuple(s)[j] ≤ n_max_region for j in 1:preKe)
    sp_w  = StateSpace{CartesianIndex{preKe}, Float64}()
    add_state!(sp_w, CartesianIndex(ntuple(_ -> round(Int, n_ss), Val(preKe))), 1.0)
    expand!(sp_w, sys_w, bc_w; depth=1)
    A_w, = build_generator(sp_w, sys_w, Float64[], 0.0)
    sp_w.probs .= expv(0.01, A_w, sp_w.probs; m=5)
    prune_threshold!(sp_w, 1e-20)
    print("$preKe ")
end
println("done.")

"""
Build the sparse diffusion+decay matrix A for mf_state ODE:
  dδ/dt = A·δ   where δ = μ - μ_ss, A = -μ·I + d·Lap₂D
"""
function build_mf_operator(K_x, K_y, d, μ_decay)
    K = K_x * K_y
    idx(i,j) = (i-1)*K_y + j
    I_list = Int[]; J_list = Int[]; V_list = Float64[]
    for i in 1:K_x, j in 1:K_y
        k  = idx(i,j)
        n_nb = (i>1)+(i<K_x)+(j>1)+(j<K_y)
        push!(I_list, k); push!(J_list, k); push!(V_list, -μ_decay - d*n_nb)
        i>1   && (push!(I_list,k); push!(J_list,idx(i-1,j)); push!(V_list,d))
        i<K_x && (push!(I_list,k); push!(J_list,idx(i+1,j)); push!(V_list,d))
        j>1   && (push!(I_list,k); push!(J_list,idx(i,j-1)); push!(V_list,d))
        j<K_y && (push!(I_list,k); push!(J_list,idx(i,j+1)); push!(V_list,d))
    end
    sparse(I_list, J_list, V_list, K, K)
end

function _fsp_step!(sp::StateSpace{CartesianIndex{Ke}, T},
                    sys, bc, dt::Float64, prune_tol::Float64,
                    max_states::Int=5000) where {Ke, T}
    length(sp) < max_states && expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, Float64[], 0.0)
    sp.probs .= expv(dt, A, sp.probs; m=30)
    renormalize!(sp)
    prune_threshold!(sp, prune_tol)
end

dt_fsp = 0.1;  prune_tol = 1e-9;  rebuild_every = 10;  save_every = 60

# Count fine voxels in each label at a given mf_state
function voxel_label_counts(mf, K_x, K_y, θ_lo, θ_hi)
    n_empty = n_front = n_settled = 0
    for i in 1:K_x, j in 1:K_y
        μ = mf[i,j]
        if μ ≤ θ_lo;   n_empty   += 1
        elseif μ ≥ θ_hi; n_settled += 1
        else;           n_front   += 1
        end
    end
    n_empty, n_front, n_settled
end

saved_t       = Float64[0.0]
ss_sizes      = Int[length(sp_c)]
keff_hist     = Int[K_eff]
col_mean_hist = [[sum(mf_state[i,j] for j in 1:K_y) for i in 1:K_x]]

# Fine voxel label counts over time (parallel to birthdeath_2d_spatial arcs)
let (ne, nf, ns) = voxel_label_counts(mf_state, K_x, K_y, θ_lo_2d, θ_hi_2d)
    global n_empty_hist   = Int[ne]
    global n_front_hist   = Int[nf]
    global n_settled_hist = Int[ns]
end

function run_adaptive_fsp!(sp_init, sys_init, bc_init, cm_init, Ke_init,
                            mf_state, A_mf, μ_ss,
                            saved_t, ss_sizes, keff_hist, col_mean_hist,
                            n_empty_hist, n_front_hist, n_settled_hist,
                            model_2d, fine_grid, dt_fsp, t_final, prune_tol,
                            rebuild_every, save_every, n_max_region, K_x, K_y,
                            θ_lo_2d, θ_hi_2d)
    sp = sp_init; sys = sys_init; bc = bc_init; cm = cm_init; Ke = Ke_init
    t = 0.0; s = 0

    while t < t_final
        dt_s = min(dt_fsp, t_final - t)
        _fsp_step!(sp, sys, bc, Float64(dt_s), prune_tol)

        # Exact mf update: advance deviation from SS via expv, no stability limit
        δ = vec(mf_state) .- μ_ss
        δ .= expv(Float64(dt_s), A_mf, δ; m = min(30, K_x*K_y))
        mf_state .= reshape(δ .+ μ_ss, K_x, K_y)
        t += dt_s;  s += 1

        if s % rebuild_every == 0
            new_cm = build_adaptive_cmap_2d(mf_state, K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
            new_Ke = new_cm.n_coarse
            new_bc = let n=n_max_region, c=new_cm
                s -> all(0 ≤ Tuple(s)[r] ≤ n for r in 1:c.n_coarse)
            end
            if new_Ke != Ke
                new_sys = build_rdme_adaptive_2d(model_2d, fine_grid, new_cm, Val(new_Ke))
                new_sp  = remap_sp_transition(sp, cm, new_cm, Val(new_Ke); weight_tol=prune_tol)
                if isempty(new_sp.states)
                    new_means = [sum(mf_state[k] for k in new_cm.coarse_to_fine[j]) for j in 1:new_Ke]
                    new_u0    = CartesianIndex(ntuple(j->round(Int,new_means[j]), Val(new_Ke)))
                    new_sp    = StateSpace{CartesianIndex{new_Ke}, Float64}()
                    add_state!(new_sp, new_u0, 1.0)
                end
                prune_threshold!(new_sp, prune_tol); renormalize!(new_sp)
                length(new_sp) < 4000 && expand!(new_sp, new_sys, new_bc; depth=1)
                @printf("  t=%5.1f  K_eff: %d→%d  |S|: %d→%d\n",
                        t, Ke, new_Ke, length(sp), length(new_sp))
                flush(stdout)
                sp=new_sp; cm=new_cm; sys=new_sys; bc=new_bc; Ke=new_Ke
            else
                sys = build_rdme_adaptive_2d(model_2d, fine_grid, new_cm, Val(Ke))
                bc  = new_bc; cm = new_cm
            end
        end

        if s % save_every == 0 || t ≥ t_final
            push!(saved_t, t); push!(ss_sizes, length(sp)); push!(keff_hist, Ke)
            push!(col_mean_hist, [sum(mf_state[i,j] for j in 1:K_y) for i in 1:K_x])
            ne, nf, ns = voxel_label_counts(mf_state, K_x, K_y, θ_lo_2d, θ_hi_2d)
            push!(n_empty_hist, ne); push!(n_front_hist, nf); push!(n_settled_hist, ns)
            @printf("  t=%5.1f  K_eff=%d  |S|=%d  empty=%d front=%d settled=%d\n",
                    t, Ke, length(sp), ne, nf, ns)
            flush(stdout)
        end
    end
end

# Build exact mf propagator (A_mf defined above)
A_mf    = build_mf_operator(K_x, K_y, mf_d, μ)
μ_ss_val = n_ss   # steady-state mean per voxel = λ/μ

println("\nRunning adaptive FSP on $(K_x)×$(K_y) grid..."); flush(stdout)
@time run_adaptive_fsp!(sp_c, sys_c, bc_c, cmap, K_eff,
                         mf_state, A_mf, μ_ss_val,
                         saved_t, ss_sizes, keff_hist, col_mean_hist,
                         n_empty_hist, n_front_hist, n_settled_hist,
                         model_2d, fine_grid, dt_fsp, t_final, prune_tol,
                         rebuild_every, save_every, n_max_region, K_x, K_y,
                         θ_lo_2d, θ_hi_2d)

println("Done.  peak |S|=$(maximum(ss_sizes))  K_eff range: $(minimum(keff_hist))–$(maximum(keff_hist))")

# ── figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

fig = Figure(size=(1400, 680))

# Helper: batch CoarseningMap region boundaries into one lines! call
function region_boundary_path(cm_snap, K_x, K_y)
    xs = Float32[]; ys = Float32[]
    for i in 1:K_x, j in 1:K_y
        r = cm_snap.fine_to_coarse[(i-1)*K_y + j]
        if j < K_y && cm_snap.fine_to_coarse[(i-1)*K_y + j+1] != r
            append!(xs, (j+0.5f0, j+0.5f0, NaN32)); append!(ys, (i-0.5f0, i+0.5f0, NaN32))
        end
        if i < K_x && cm_snap.fine_to_coarse[i*K_y + j] != r
            append!(xs, (j-0.5f0, j+0.5f0, NaN32)); append!(ys, (i+0.5f0, i+0.5f0, NaN32))
        end
    end
    xs, ys
end

# ── Row 1: spatial snapshots ──────────────────────────────────────────────────
snap_ts = [3.0, 5.0, 7.0, 15.0, t_final]
# Cap at 3×n_ss so settling gradient is clearly visible (IC pulse off-scale is OK)
c_max   = n_ss * 3.0

for (col, t_s) in enumerate(snap_ts)
    si  = argmin(abs.(t_mf .- t_s))
    ax  = Axis(fig[1, col]; aspect=DataAspect(),
               title="t = $(round(t_s; digits=0))", titlesize=11,
               xticksvisible=false, yticksvisible=false,
               xticklabelsvisible=false, yticklabelsvisible=false)

    hm = heatmap!(ax, 1:K_x, 1:K_y, mu_mf[si]';
                  colormap=:inferno, colorrange=(0.0, c_max), lowclip=:black)

    cm_snap = build_adaptive_cmap_2d(mu_mf[si], K_x, K_y; θ_lo=θ_lo_2d, θ_hi=θ_hi_2d)
    xs_b, ys_b = region_boundary_path(cm_snap, K_x, K_y)
    isempty(xs_b) || lines!(ax, xs_b, ys_b; color=(:white, 0.9), linewidth=2.0)

    text!(ax, 0.02, 0.96; text="K_eff=$(cm_snap.n_coarse)",
          color=:white, fontsize=9, align=(:left,:top), space=:relative,
          strokecolor=:black, strokewidth=2)

    col == length(snap_ts) &&
        Colorbar(fig[1, col+1], hm; label="E[nᵢⱼ]", width=12, labelsize=10,
                 ticks=([0, n_ss, c_max],
                        ["0", "n_ss=$(n_ss)", "$(round(c_max,digits=1))"]))
end

Label(fig[0, 1:5],
      "Adaptive FSP  |  $(K_x)×$(K_y) BD-RDME  D=$D  n_ss=$(n_ss)  IC: N=$(N_ic) at corner  " *
      "|  white lines = CoarseningMap super-region boundaries";
      fontsize=11, tellwidth=false)

# ── Row 2: diagnostics ────────────────────────────────────────────────────────
# Arc: fine voxel label counts (like birthdeath_2d_spatial)
ax_arc = Axis(fig[2, 1:2];
              xlabel="time", ylabel="# fine voxels  (out of $(K_x*K_y))",
              title="Fine-voxel states — $(K_x)×$(K_y) = $(K_x*K_y) voxels tracked via K_eff ≤ 3 super-regions",
              titlesize=10)
band!(ax_arc, saved_t, zeros(length(saved_t)), Float64.(n_settled_hist);
      color=(:royalblue, 0.25))
band!(ax_arc, saved_t, Float64.(n_settled_hist),
      Float64.(n_settled_hist) .+ Float64.(n_front_hist); color=(:tomato, 0.25))
# Compute label counts from the dense mf snapshots for smooth arc
mf_ne = Int[]; mf_nf = Int[]; mf_ns = Int[]
for snap in mu_mf
    ne, nf, ns = voxel_label_counts(snap, K_x, K_y, θ_lo_2d, θ_hi_2d)
    push!(mf_ne, ne); push!(mf_nf, nf); push!(mf_ns, ns)
end
band!(ax_arc, t_mf, zeros(length(t_mf)), Float64.(mf_ns); color=(:royalblue,0.25))
band!(ax_arc, t_mf, Float64.(mf_ns), Float64.(mf_ns).+Float64.(mf_nf); color=(:tomato,0.25))
lines!(ax_arc, t_mf, Float64.(mf_ns);
       color=:royalblue, linewidth=2.5, label="Settled (n̄ > θ_hi=$(round(θ_hi_2d,digits=2)))")
lines!(ax_arc, t_mf, Float64.(mf_nf);
       color=:tomato, linewidth=2.5, label="Front (θ_lo < n̄ < θ_hi)")
lines!(ax_arc, t_mf, Float64.(mf_ne);
       color=:gray50, linewidth=2, linestyle=:dash, label="Empty (n̄ < θ_lo=$(round(θ_lo_2d,digits=2)))")
hlines!(ax_arc, [Float64(K_x*K_y)]; color=(:black,0.3), linewidth=1, linestyle=:dot)
axislegend(ax_arc; position=:rc, labelsize=9, framevisible=false)
ylims!(ax_arc, -5, K_x*K_y + 20)

ax_s = Axis(fig[2, 3]; xlabel="time", ylabel="|S|",
            title="|S|: adaptive FSP vs full FSP scaling", titlesize=10)
lines!(ax_s, saved_t, Float64.(ss_sizes); color=:steelblue4, linewidth=2.5)
scatter!(ax_s, saved_t, Float64.(ss_sizes); color=:steelblue4, markersize=7)
text!(ax_s, saved_t[2], maximum(ss_sizes) * 1.08;
      text="peak |S| = $(maximum(ss_sizes))\n(full FSP: n_max^K — astronomical)",
      fontsize=9, color=:steelblue4)

# Total-count distribution at final time from FSP marginals
ax_p = Axis(fig[2, 4:5]; xlabel="diagonal voxel index",
            ylabel="E[nᵢᵢ]  (diagonal profile)",
            title="Mean-field diagonal profile at snapshots", titlesize=10)
diag_len = min(K_x, K_y)
clrs_p = Makie.resample_cmap(:plasma, length(snap_ts)+1)[1:end-1]
for (si_c, t_s) in enumerate(snap_ts)
    si = argmin(abs.(t_mf .- t_s))
    diag_vals = [mu_mf[si][k, k] for k in 1:diag_len]
    lines!(ax_p, 1:diag_len, diag_vals; color=clrs_p[si_c], linewidth=1.6,
           label="t=$(round(Int, t_s))")
end
hlines!(ax_p, [n_ss]; color=(:gray50, 0.5), linestyle=:dash, linewidth=1)
text!(ax_p, diag_len*0.6, n_ss*1.08; text="n_ss=$(n_ss)", fontsize=9, color=:gray50)
axislegend(ax_p; position=:rt, labelsize=9, framevisible=false)

rowgap!(fig.layout, 1, 8); colgap!(fig.layout, 5, 5)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_2d_bd_large.pdf")
save(out, fig)
println("Saved → $out")
