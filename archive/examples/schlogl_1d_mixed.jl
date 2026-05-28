"""
1D Schlögl Mixed-Resolution FSP
================================

Bistable Schlögl wavefront on a K-voxel 1D grid.
  n_hi-state occupies left half, n_lo-state occupies right half.
  n_hi is thermodynamically favored → front sweeps right.

Mixed-resolution strategy:
  - Active strip (K_act ≈ 2-4 voxels at the front): joint FSP, exact
  - Inactive voxels: mean-field ODE  dμ/dt = f_sch(μ) + D·Δμ

Key result: at the front, P(n) is BIMODAL — the voxel is quantum-superposed
between n_lo and n_hi states.  Mean-field sees only E[n], missing this entirely.
FSP, applied only to the 2-4 active voxels, captures it at negligible cost.

Run: julia --project examples/schlogl_1d_mixed.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using LinearAlgebra
using Printf
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Model parameters
# ─────────────────────────────────────────────────────────────────────────────

const C1 = 0.028
const C2 = 3.2e-4
const C3 = 19.5
const C4 = 1.0
const D  = 1.0      # D=1 → front width ≈ √(D/|f'_un|) ≈ 2 voxels: tight, clean

model = SchloglModel1D(D, C1, C2, C3, C4)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Fixed points:  n_lo=%d  n_unstable=%d  n_hi=%d\n", n_lo, n_un, n_hi)

const K        = 50      # total voxels (active front is independent of K)
const DX       = 1.0    # voxel size
const T_END    = 8.0
const DT       = 0.25
const N_MAX    = 200    # per-voxel FSP truncation (rstep bounds the actual support)
const MAX_KA   = 8      # hard cap on active voxels
const N_LEVELS = 1      # V-cycle levels on the active strip (0 = direct expv)

# ─────────────────────────────────────────────────────────────────────────────
# Mixed-resolution Schlögl runner
# ─────────────────────────────────────────────────────────────────────────────

# ── Type-stable inner evolve: dispatches on StateSpace{CartesianIndex{Ka}} ──────
# Function barrier: Julia compiles a specialization for each Ka value.
# V-cycle when Ka is even and ≥ 4 and n_levels ≥ 1; direct expv otherwise.
function _evolve_active_schlogl!(sp::StateSpace{CartesianIndex{Ka}, Float64},
                                  model, ga, lbc, rbc, μ, rates,
                                  t, dt, n_max, prune_tol, krylov_m,
                                  n_levels, τ_pre_frac) where Ka
    sys = build_active_schlogl_1d(model, ga, lbc, rbc, μ, Val(Ka))
    bc  = schlogl_mixed_bc(n_max)

    use_vcycle = n_levels >= 1 && iseven(Ka) && Ka >= 4

    if use_vcycle
        # Determine actual usable levels (each level halves Ka; stop if Ka/2 is 0)
        max_l = floor(Int, log2(Ka))
        nl    = min(n_levels, max_l)
        hier  = build_hierarchy(ga, nl)
        τ_pre = τ_pre_frac * dt

        # Expand before V-cycle so boundary states are in the space
        expand!(sp, sys, bc; depth = 1)

        sp = multi_level_vcycle_schlogl(sp, model, hier, lbc, rbc, μ, rates, t, dt;
                                         τ_pre, τ_post = τ_pre,
                                         krylov_m,
                                         weight_tol  = prune_tol,
                                         expand_coarse = true,
                                         coarse_r_depth = 2,
                                         coarse_d_depth = 1,
                                         coarse_n_max = 2 * n_max,
                                         max_states   = 2_000_000,
                                         prune_tol)
    else
        # Direct expv: small Ka or odd Ka
        expand!(sp, sys, bc; depth = 1)
        A, = build_generator(sp, sys, rates, t)
        sp.probs .= expv(dt, A, sp.probs; m = krylov_m)
    end

    renormalize!(sp)
    prune_threshold!(sp, prune_tol)
    sp
end

function run_schlogl_mixed(K, model, dx, n_max, T_END, DT, n_lo, n_hi;
                            max_k_act    = 4,
                            n_levels     = 1,        # V-cycle levels (0 = direct expv)
                            τ_pre_frac   = 0.1,      # pre/post smooth time as fraction of dt
                            promote_margin = 15.0,
                            demote_margin  = 8.0,
                            prune_tol    = 1e-12,
                            krylov_m     = 30)

    # IC: left half at n_hi, right half at n_lo
    μ = [i <= K÷2 ? Float64(n_hi) : Float64(n_lo) for i in 1:K]

    # Initial active set: the two voxels straddling the front
    k_lo0 = K÷2; k_hi0 = K÷2 + 1
    avs = ActiveVoxelSet1D(k_lo0, k_hi0, K)

    # Initial FSP: delta at (n_hi, n_lo)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    add_state!(sp, CartesianIndex(n_hi, n_lo), 1.0)

    # Seed expansion so boundary states exist in the state space
    grid_act = VoxelGrid(2, dx, 0)
    sys0 = build_active_schlogl_1d(model, grid_act, k_lo0-1, k_hi0+1, μ, Val(2))
    expand!(sp, sys0, schlogl_mixed_bc(n_max); depth = 4)

    n_steps  = round(Int, T_END / DT)
    rates    = Float64[]

    # Tracking
    times    = Float64[0.0]
    μ_hist   = [copy(μ)]
    sp_sizes = Int[length(sp)]
    ka_hist  = Int[n_active(avs)]
    alo_hist = Int[avs.k_lo]
    ahi_hist = Int[avs.k_hi]
    # Marginal P(n) snapshots at the middle of the active strip (for bimodality plot)
    marg_snaps = Tuple{Float64, Vector{Float64}, Vector{Float64}, Int, Int}[]

    snap_every = max(1, n_steps ÷ 8)

    for step in 1:n_steps
        t    = Float64((step - 1) * DT)
        Ka   = n_active(avs)
        lbc  = left_inactive(avs)
        rbc  = right_inactive(avs)
        ga   = VoxelGrid(Ka, dx, 0)

        # Evolve: V-cycle if Ka is even and ≥ 4, otherwise direct expv.
        # Function barrier (_evolve_active_schlogl!) makes inner loops type-stable
        # even though Ka (and thus the StateSpace type) changes between steps.
        sp = evolve_active_schlogl!(sp, model, ga, lbc, rbc, μ, rates, t, DT,
                                     n_max, prune_tol, krylov_m,
                                     n_levels, τ_pre_frac)

        # Update inactive means
        update_inactive_means_schlogl_1d!(μ, sp, avs, model, dx, DT)

        # Record marginal P(n) for both the leftmost and rightmost active voxels
        if step % snap_every == 0
            Ka = n_active(avs)
            P_left  = zeros(n_max + 1)   # local index 1 = k_lo (n_hi side)
            P_right = zeros(n_max + 1)   # local index Ka = k_hi (n_lo side, leading edge)
            for idx in eachindex(sp.states)
                p = sp.probs[idx]; p > 0.0 || continue
                tv = Tuple(sp.states[idx])
                P_left[tv[1]  + 1] += p
                P_right[tv[Ka] + 1] += p
            end
            push!(marg_snaps, (t + DT, P_left, P_right, avs.k_lo, avs.k_hi))
        end

        # Adapt active set
        sp, avs = adapt_schlogl_1d(sp, avs, μ, n_max, n_lo, n_hi;
                                    promote_margin, demote_margin,
                                    min_k_act = 2, max_k_act,
                                    weight_tol = prune_tol)

        push!(times, step * DT); push!(μ_hist, copy(μ))
        push!(sp_sizes, length(sp)); push!(ka_hist, n_active(avs))
        push!(alo_hist, avs.k_lo);  push!(ahi_hist, avs.k_hi)
    end

    times, μ_hist, sp_sizes, ka_hist, alo_hist, ahi_hist,
    marg_snaps  # Vector of (t, P_left, P_right, k_lo, k_hi)
end

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

println("Running 1D Schlögl mixed-resolution FSP  (K=$K, D=$D, n_levels=$N_LEVELS) ...")
@time times, μ_hist, sp_sizes, ka_hist, alo_hist, ahi_hist, marg_snaps =
    run_schlogl_mixed(K, model, DX, N_MAX, T_END, DT, n_lo, n_hi;
                      max_k_act = MAX_KA, n_levels = N_LEVELS,
                      τ_pre_frac = 0.1,
                      promote_margin = 15.0, demote_margin = 8.0)
println("Done. |S|_max=$(maximum(sp_sizes))  K_act_max=$(maximum(ka_hist))")

# ─────────────────────────────────────────────────────────────────────────────
# Mean-field reference (same ODE, no FSP)
# ─────────────────────────────────────────────────────────────────────────────

function run_schlogl_ode(K, D_val, T_END, DT)
    d = D_val   # dx=1
    f_sch(m) = C3 + C1*m*(m-1)/2 - C4*m - C2*m*(m-1)*(m-2)/6
    dt_stable = 1.0 / (2d + 10.0)
    n_sub  = max(1, ceil(Int, DT / dt_stable))
    dt_sub = DT / n_sub
    μ = [i <= K÷2 ? Float64(n_hi) : Float64(n_lo) for i in 1:K]
    n_steps = round(Int, T_END / DT)
    μ_hist_ode = [copy(μ)]
    Δ = zeros(K)
    for _ in 1:n_steps
        for _ in 1:n_sub
            fill!(Δ, 0.0)
            for k in 1:K
                Δ[k] = f_sch(μ[k])
                k > 1 && (Δ[k] += d*(μ[k-1]-μ[k]))
                k < K && (Δ[k] += d*(μ[k+1]-μ[k]))
            end
            μ .= max.(0.0, μ .+ Δ .* dt_sub)
        end
        push!(μ_hist_ode, copy(μ))
    end
    μ_hist_ode
end

println("Running mean-field ODE reference ...")
μ_hist_ode = run_schlogl_ode(K, D, T_END, DT)

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

# Choose 4 evenly-spaced snapshot indices
n_steps_total = length(times) - 1
snap_idx = [1, n_steps_total÷3+1, 2*n_steps_total÷3+1, n_steps_total+1]
snap_idx = clamp.(snap_idx, 1, length(times))
snap_colors = [:white, :cyan, :orange, :crimson]

# Choose 4 marginal snapshots
n_marg = length(marg_snaps)
marg_pick = n_marg >= 4 ? [1, n_marg÷3+1, 2*n_marg÷3+1, n_marg] : collect(1:n_marg)

fig = Figure(size = (1600, 900), backgroundcolor = :black)

# ─── Row 1: Mean-field E[n] profiles: FSP vs ODE ─────────────────────────────
ax_prof = Axis(fig[1, 1:3];
               title      = "Mean-field profile E[nₖ]  (FSP means vs ODE reference)",
               xlabel     = "voxel k",
               ylabel     = "E[nₖ]",
               titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
               xticklabelcolor = :white, yticklabelcolor = :white,
               backgroundcolor = :black,
               xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))

for (ci, si) in enumerate(snap_idx)
    t_s = times[si]
    lines!(ax_prof, 1:K, μ_hist[si];
           color = snap_colors[ci], linewidth = 3,
           label = @sprintf("FSP  t=%.1f", t_s))
    lines!(ax_prof, 1:K, μ_hist_ode[si];
           color = snap_colors[ci], linewidth = 1.5, linestyle = :dash)
end
hlines!(ax_prof, [Float64(n_lo)]; color = :royalblue, linewidth = 1,
        linestyle = :dash, label = "n_lo=$(n_lo)")
hlines!(ax_prof, [Float64(n_hi)]; color = :tomato,    linewidth = 1,
        linestyle = :dash, label = "n_hi=$(n_hi)")
axislegend(ax_prof; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 9, position = :lc)
ylims!(ax_prof, n_lo - 20, n_hi + 20)

# ─── Row 1 right: |S| and K_act over time ─────────────────────────────────────
ax_s = Axis(fig[1, 4];
            title      = "|S| active states",
            xlabel     = "time", ylabel = "state count",
            titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
            xticklabelcolor = :white, yticklabelcolor = :white,
            backgroundcolor = :black,
            xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))
lines!(ax_s, times, Float64.(sp_sizes); color = :white, linewidth = 2)

ax_k = Axis(fig[1, 5];
            title      = "K_act (active voxels)",
            xlabel     = "time", ylabel = "K_act",
            titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
            xticklabelcolor = :white, yticklabelcolor = :white,
            backgroundcolor = :black,
            xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))
lines!(ax_k, times, Float64.(ka_hist); color = :orange, linewidth = 2)
hlines!(ax_k, [Float64(MAX_KA)]; color = :gray, linewidth = 1, linestyle = :dash)
ylims!(ax_k, 0, MAX_KA + 2)

# ─── Row 2: Space-time heatmap ────────────────────────────────────────────────
ax_ht = Axis(fig[2, 1:3];
             title      = "E[nₖ] space-time  (white overlay = FSP active strip)",
             xlabel     = "voxel k", ylabel = "time",
             titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
             xticklabelcolor = :white, yticklabelcolor = :white,
             backgroundcolor = :black,
             xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1))

nT = length(times)
mat = zeros(nT, K)
for ti in 1:nT; mat[ti, :] .= μ_hist[ti]; end
heatmap!(ax_ht, 1:K, times, mat';
         colormap   = Reverse(:RdBu),
         colorrange = (Float32(n_lo - 15), Float32(n_hi + 15)))
# Active strip overlay
band!(ax_ht, times, Float64.(alo_hist), Float64.(ahi_hist);
      color = RGBAf(1,1,1,0.30))
lines!(ax_ht, times, Float64.(alo_hist); color = :white, linewidth = 1)
lines!(ax_ht, times, Float64.(ahi_hist); color = :white, linewidth = 1)

# ─── Row 2 right: Scaling comparison ─────────────────────────────────────────
# Full FSP cost ~ n_max^K; mixed-resolution cost ~ n_max^K_act (constant in K)
ax_scale = Axis(fig[2, 4];
                title      = "State-space cost vs grid size K",
                xlabel     = "K (total voxels)",
                ylabel     = "log₁₀ |S|",
                titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
                xticklabelcolor = :white, yticklabelcolor = :white,
                backgroundcolor = :black,
                xgridcolor = RGBAf(1,1,1,0.1), ygridcolor = RGBAf(1,1,1,0.1),
                titlesize  = 12)

Ks = [2, 4, 6, 8, 10, 12, 16]
n_eff = 180   # effective per-voxel state count (near n_lo or n_hi)
full_cost  = [Float64(K_i) * log10(n_eff) for K_i in Ks]   # n_eff^K
mixed_cost = fill(log10(Float64(maximum(sp_sizes))), length(Ks))  # actual peak |S|

lines!(ax_scale, Ks, full_cost;
       color = :tomato, linewidth = 3, label = "Full FSP: n_eff^K")
lines!(ax_scale, Ks, mixed_cost;
       color = :limegreen, linewidth = 2.5, linestyle = :dash,
       label = @sprintf("Mixed-res: |S|_peak=%d", maximum(sp_sizes)))
scatter!(ax_scale, [K], [log10(Float64(maximum(sp_sizes)))];
         color = :limegreen, markersize = 10)
axislegend(ax_scale; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 9, position = :lt)
ylims!(ax_scale, 0, Float64(maximum(Ks)) * log10(n_eff) * 1.05)

# Annotate intractable zone
hlines!(ax_scale, [6.0]; color = (:white, 0.3), linewidth = 1, linestyle = :dot)
text!(ax_scale, 2.5, 6.3; text = "10⁶ states", color = (:white, 0.5), fontsize = 9)

# ─── Row 2 right: Single-voxel bimodal P(n) ──────────────────────────────────
println("Computing single-voxel stationary distribution ...")
function schlogl_stationary_local(c1, c2, c3, c4; n_max=250, tol=1e-14)
    n_states = n_max + 1
    Q = zeros(n_states, n_states)
    for n in 0:n_max
        prod_rate = c3 + c1 * max(0, n*(n-1)) / 2
        deg_rate  = c4 * n + c2 * max(0, n*(n-1)*(n-2)) / 6
        if n < n_max
            Q[n+2, n+1] += prod_rate; Q[n+1, n+1] -= prod_rate
        end
        if n > 0
            Q[n,   n+1] += deg_rate;  Q[n+1, n+1] -= deg_rate
        end
    end
    vals, vecs = eigen(Q)
    idx = argmin(abs.(real.(vals)))
    π = real.(vecs[:, idx]); π .= abs.(π); π ./= sum(π)
    π[π .< tol] .= 0.0; π ./= sum(π)
    0:n_max, π
end
ns_stat, π_stat = schlogl_stationary_local(C1, C2, C3, C4; n_max=250)

ax_bi = Axis(fig[2, 5];
             title      = "Single-voxel P(n) at SS — FSP vs mean-field",
             xlabel     = "molecule count n",
             ylabel     = "probability",
             titlecolor = :white, xlabelcolor = :white, ylabelcolor = :white,
             xticklabelcolor = :white, yticklabelcolor = :white,
             backgroundcolor = :black,
             xgridcolor = RGBAf(1,1,1,0.08), ygridcolor = RGBAf(1,1,1,0.08),
             titlesize  = 12)

band!(ax_bi, Float64.(ns_stat), zeros(length(ns_stat)), Float64.(π_stat);
      color = (:white, 0.15))
lines!(ax_bi, Float64.(ns_stat), Float64.(π_stat);
       color = :white, linewidth = 2, label = "FSP: bimodal P(n)")

# Mean-field "prediction" for the same voxel at midpoint between fixed points
μ_mid = (n_lo + n_hi) / 2   # representative "front" mean
σ_mf  = sqrt(μ_mid)          # Poisson std (mean-field approximation)
n_mf  = 0:250
gauss = @. exp(-0.5 * ((n_mf - μ_mid) / σ_mf)^2) / (σ_mf * sqrt(2π))
gauss ./= sum(gauss)
lines!(ax_bi, Float64.(n_mf), gauss;
       color = :tomato, linewidth = 2, linestyle = :dash,
       label = @sprintf("ODE: single peak at E[n]=%.0f", μ_mid))

vlines!(ax_bi, [Float64(n_lo)]; color = :royalblue, linewidth = 1.5,
        linestyle = :dash, label = "n_lo=$(n_lo)")
vlines!(ax_bi, [Float64(n_hi)]; color = :crimson,   linewidth = 1.5,
        linestyle = :dash, label = "n_hi=$(n_hi)")
axislegend(ax_bi; labelcolor = :white, backgroundcolor = :transparent,
           framecolor = :transparent, labelsize = 8, position = :rt)
xlims!(ax_bi, 0, 220)

# ─── Title ────────────────────────────────────────────────────────────────────
Label(fig[0, 1:5],
      "1D Schlögl Mixed-Resolution FSP  |  K=$K voxels  D=$D  n_lo=$(n_lo)  n_hi=$(n_hi)  " *
      "White overlay = FSP active strip (K_act ≤ $(maximum(ka_hist)));  inactive voxels: mean-field ODE",
      color = :white, fontsize = 13)

Label(fig[3, 1:5],
      "FSP means (solid) match ODE reference (dashed) everywhere.  " *
      "Key: P(n) is BIMODAL near the front — ODE sees only E[n], missing the two-state structure.  " *
      "Mixed-resolution cost: |S| ≤ $(maximum(sp_sizes)) states (constant in K); full FSP scales as n_max^K.",
      color = :white, fontsize = 11)

png_path = joinpath(outdir, "schlogl_1d_mixed.png")
save(png_path, fig; px_per_unit = 3)
println("Saved: $png_path")
println("Done.")
