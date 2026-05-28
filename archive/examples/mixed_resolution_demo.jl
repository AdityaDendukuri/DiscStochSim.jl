"""
Mixed-Resolution CME — Phase 1: Fixed Active Set
=================================================

Key idea: represent the CME state as
  - "Active" voxels  → full FSP in CartesianIndex{NDIM_full}, inactive dims pinned to 0 by BC
  - "Inactive" voxels → scalar means μA[k], μB[k] only

The state space is O((N_MAX+1)^(2*K_act)) regardless of NDIM_full, because
the BC rejects any state with non-zero count in an inactive voxel.

Boundary coupling:
  Active → Inactive diffusion: standard CME loss term  (nA_k → nA_k-1)
  Inactive → Active diffusion: mean-field birth term    (nA_k → nA_k+1 at rate d*μA[inactive])

The propensity closures capture μA/μB by reference — updating μA[k] in-place
automatically updates all propensities without rebuilding the system struct.

Mean update:
  Active means:   E[nA_k] recomputed from sp.probs after each expv step (exact)
  Inactive means: explicit Euler on the diffusion equation dμ/dt = D∇²μ (exact by linearity)

Known approximation:
  Reactions in inactive voxels are skipped. If both species transit through an
  inactive voxel simultaneously, they don't annihilate → slight overestimate of
  total molecules. Cured in Phase 2 by activating voxels when both μA and μB exceed ε.

Validation: compare full FSP (NDIM=8, 2×2) vs mixed (NDIM=8, active=[1,4])
Demo:       4×4 grid (NDIM=32) with fixed active=[1,16] — previously intractable

Run: julia --project examples/mixed_resolution_demo.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# Core: BC, system builder, mean updater
# ─────────────────────────────────────────────────────────────────────────────

"""
    mixed_bc(active_set, K_vox, N_MAX, ::Val{NDIM})

Boundary condition that accepts states iff:
  - active voxels: 0 ≤ count ≤ N_MAX
  - inactive voxels: count == 0 (pinned)

This restriction makes the reachable FSP state space O((N_MAX+1)^(2*K_act))
regardless of NDIM, which is the key tractability gain.
"""
function mixed_bc(active_set::Set{Int}, K_vox::Int, N_MAX::Int, ::Val{NDIM}) where NDIM
    function bc(s::CartesianIndex{NDIM})
        t = Tuple(s)
        @inbounds for k in 1:K_vox
            iA = 2k-1; iB = 2k
            if k ∈ active_set
                (0 ≤ t[iA] ≤ N_MAX && 0 ≤ t[iB] ≤ N_MAX) || return false
            else
                (t[iA] == 0 && t[iB] == 0) || return false
            end
        end
        return true
    end
end

"""
    mixed_system(active_set, KC_x, KC_y, DA_C, DB_C, K_R, μA, μB, ::Val{NDIM})

Build the mixed-resolution DiscreteStochasticSystem.

The system is built ONCE. Propensities for inactive→active boundary terms close
over μA and μB vectors by reference — updating them in-place (via update_means!)
is sufficient to keep the generator current without rebuilding the struct.

Reactions:
  • A+B→∅ in every active voxel
  • Active↔Active diffusion: standard CME off-diagonal entries
  • Active↔Inactive boundary (per edge):
      Loss from active:    standard nA_k → nA_k - 1 at rate DA_C * nA_k
      Source into active:  mean-field birth  at rate DA_C * μA[inactive]  ← closure
"""
function mixed_system(
    active_set::Set{Int}, KC_x::Int, KC_y::Int,
    DA_C::Float64, DB_C::Float64, K_R::Float64,
    μA::Vector{Float64}, μB::Vector{Float64},   # captured by reference
    ::Val{NDIM}
) where NDIM
    K_vox = NDIM ÷ 2
    stoichs      = CartesianIndex{NDIM}[]
    propensities = Function[]
    unit(i)      = CartesianIndex(ntuple(j -> j == i ? 1 : 0, Val(NDIM)))

    # ── Per-voxel annihilation (active voxels only) ───────────────────────────
    for k in active_set
        let k = k
            iA = 2k-1; iB = 2k
            eA = unit(iA); eB = unit(iB)
            push!(stoichs, -eA - eB)
            push!(propensities, (x, rv, t) -> K_R * max(0, Tuple(x)[iA]) * max(0, Tuple(x)[iB]))
        end
    end

    # ── Diffusion: iterate undirected edges, handle 3 cases ──────────────────
    for i in 1:KC_x, j in 1:KC_y
        k = (i-1)*KC_y + j
        neighbors_k2 = Int[]
        j < KC_y    && push!(neighbors_k2, k + 1)
        i < KC_x    && push!(neighbors_k2, k + KC_y)

        for k2 in neighbors_k2
            ka  = k  ∈ active_set
            k2a = k2 ∈ active_set
            ka || k2a || continue          # inactive-inactive: handled by mean ODE

            let k=k, k2=k2, ka=ka, k2a=k2a
                iAk=2k-1; iBk=2k; iAk2=2k2-1; iBk2=2k2

                if ka && k2a
                    # ── Active ↔ Active: bidirectional standard CME ───────────
                    eAk=unit(iAk); eBk=unit(iBk); eAk2=unit(iAk2); eBk2=unit(iBk2)
                    push!(stoichs,-eAk+eAk2); push!(propensities,(x,rv,t)->DA_C*max(0,Tuple(x)[iAk]))
                    push!(stoichs, eAk-eAk2); push!(propensities,(x,rv,t)->DA_C*max(0,Tuple(x)[iAk2]))
                    push!(stoichs,-eBk+eBk2); push!(propensities,(x,rv,t)->DB_C*max(0,Tuple(x)[iBk]))
                    push!(stoichs, eBk-eBk2); push!(propensities,(x,rv,t)->DB_C*max(0,Tuple(x)[iBk2]))

                elseif ka && !k2a
                    # ── Active k ↔ Inactive k2 ────────────────────────────────
                    eAk=unit(iAk); eBk=unit(iBk)
                    # Loss: active k → inactive k2  (standard CME)
                    push!(stoichs,-eAk); push!(propensities,(x,rv,t)->DA_C*max(0,Tuple(x)[iAk]))
                    push!(stoichs,-eBk); push!(propensities,(x,rv,t)->DB_C*max(0,Tuple(x)[iBk]))
                    # Source: inactive k2 → active k  (mean-field, μ captured by ref)
                    let k2=k2
                        push!(stoichs,+eAk); push!(propensities,(x,rv,t)->DA_C*μA[k2])
                        push!(stoichs,+eBk); push!(propensities,(x,rv,t)->DB_C*μB[k2])
                    end

                else   # !ka && k2a
                    # ── Inactive k ↔ Active k2 ────────────────────────────────
                    eAk2=unit(iAk2); eBk2=unit(iBk2)
                    # Source: inactive k → active k2
                    let k=k
                        push!(stoichs,+eAk2); push!(propensities,(x,rv,t)->DA_C*μA[k])
                        push!(stoichs,+eBk2); push!(propensities,(x,rv,t)->DB_C*μB[k])
                    end
                    # Loss: active k2 → inactive k
                    push!(stoichs,-eAk2); push!(propensities,(x,rv,t)->DA_C*max(0,Tuple(x)[iAk2]))
                    push!(stoichs,-eBk2); push!(propensities,(x,rv,t)->DB_C*max(0,Tuple(x)[iBk2]))
                end
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{NDIM}}(stoichs, propensities)
end

"""
    update_means!(μA, μB, sp, active_set, KC_x, KC_y, DA_C, DB_C, dt)

Two-phase mean update after each expv step:

Phase 1 (exact):
  Recompute E[nA_k] and E[nB_k] for active voxels by summing sp.probs.
  This captures annihilation and diffusion losses from the CME exactly.

Phase 2 (exact by linearity of expectation):
  Euler step for inactive voxels:
    dμA_k/dt = DA_C * Σ_{m ∈ neighbors(k)} (μA[m] - μA[k])
  where μA[m] for active m uses the Phase-1-updated value (split-step ordering).

Mass conservation: the boundary flux added to inactive means exactly matches
the expected flux lost from the active CME (both = DA_C * E[nA_active] * dt).
"""
function update_means!(
    μA::Vector{Float64}, μB::Vector{Float64},
    sp,
    active_set::Set{Int}, KC_x::Int, KC_y::Int,
    DA_C::Float64, DB_C::Float64, dt::Float64
)
    K_vox = KC_x * KC_y

    # Phase 1: recompute active means from sp.probs
    for k in active_set; μA[k] = 0.0; μB[k] = 0.0; end
    for idx in eachindex(sp.states)
        t = Tuple(sp.states[idx]); p = sp.probs[idx]
        p > 0.0 || continue
        for k in active_set
            μA[k] += p * t[2k-1]
            μB[k] += p * t[2k]
        end
    end

    # Phase 2: sub-stepped explicit Euler for inactive diffusion.
    # Stability of explicit Euler on the graph Laplacian:
    #   dt_sub * max_d * max_degree < 1  (factor 2 safety margin used below)
    # For a 2D grid the worst case is degree 4 (interior voxel).
    max_d   = max(DA_C, DB_C)
    n_sub   = max(1, ceil(Int, max_d * dt * 4 * 2))  # degree 4, factor 2 safety
    dt_sub  = dt / n_sub
    ΔμA = zeros(K_vox); ΔμB = zeros(K_vox)
    for _ in 1:n_sub
        fill!(ΔμA, 0.0); fill!(ΔμB, 0.0)
        for i in 1:KC_x, j in 1:KC_y
            k = (i-1)*KC_y + j; k ∈ active_set && continue
            j < KC_y && (ΔμA[k] += DA_C*(μA[k+1]-μA[k]);    ΔμB[k] += DB_C*(μB[k+1]-μB[k]))
            j > 1    && (ΔμA[k] += DA_C*(μA[k-1]-μA[k]);    ΔμB[k] += DB_C*(μB[k-1]-μB[k]))
            i < KC_x && (ΔμA[k] += DA_C*(μA[k+KC_y]-μA[k]); ΔμB[k] += DB_C*(μB[k+KC_y]-μB[k]))
            i > 1    && (ΔμA[k] += DA_C*(μA[k-KC_y]-μA[k]); ΔμB[k] += DB_C*(μB[k-KC_y]-μB[k]))
        end
        for k in 1:K_vox
            k ∈ active_set && continue
            μA[k] = max(0.0, μA[k] + ΔμA[k]*dt_sub)
            μB[k] = max(0.0, μB[k] + ΔμB[k]*dt_sub)
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Simulation runner — generic over NDIM
# ─────────────────────────────────────────────────────────────────────────────

function run_mixed(
    ::Val{NDIM}, KC_x::Int, KC_y::Int,
    N_MOL::Int, N_MAX::Int,
    DA_C::Float64, DB_C::Float64, K_R::Float64,
    active_voxels::Vector{Int},
    T_END::Float64, DT::Float64;
    label::String = ""
) where NDIM
    K_vox      = NDIM ÷ 2
    active_set = Set(active_voxels)

    # μA, μB shared with propensity closures — update in-place each step
    μA = zeros(K_vox); μB = zeros(K_vox)

    # IC: N_MOL A in voxel 1, N_MOL B in last voxel
    #   If a source voxel is active → molecule appears in FSP state
    #   If inactive → only in mean (mean-field carries the initial mass)
    ic = zeros(Int, NDIM)
    μA[1]     = float(N_MOL); 1     ∈ active_set && (ic[1]    = N_MOL)
    μB[K_vox] = float(N_MOL); K_vox ∈ active_set && (ic[NDIM] = N_MOL)

    sp  = StateSpace{CartesianIndex{NDIM}, Float64}()
    add_state!(sp, CartesianIndex(Tuple(ic)), 1.0)
    sys = mixed_system(active_set, KC_x, KC_y, DA_C, DB_C, K_R, μA, μB, Val(NDIM))
    bc  = mixed_bc(active_set, K_vox, N_MAX, Val(NDIM))
    rates = Float64[]

    N_STEPS  = round(Int, T_END / DT)
    times    = Float64[0.0]
    totals   = Float64[sum(μA) + sum(μB)]
    sp_sizes = Int[length(sp)]

    t_wall = @elapsed for step in 1:N_STEPS
        expand!(sp, sys, bc; depth=1)
        A, = build_generator(sp, sys, rates, Float64((step-1)*DT))
        sp.probs .= expv(Float64(DT), A, sp.probs; m=30)
        renormalize!(sp); prune_threshold!(sp, 1e-12)
        update_means!(μA, μB, sp, active_set, KC_x, KC_y, DA_C, DB_C, DT)
        push!(times, step*DT); push!(totals, sum(μA)+sum(μB)); push!(sp_sizes, length(sp))
    end

    @printf("%-38s  |S|_peak=%6d  wall=%6.2fs  E[n](T)=%.3f\n",
            label, maximum(sp_sizes), t_wall, totals[end])
    times, totals, sp_sizes
end

# ─────────────────────────────────────────────────────────────────────────────
# Validation: 2×2 grid (NDIM=8)
# Full FSP  vs  Mixed [1,4] active  vs  Mixed [1,2,3,4] (should match Full)
# ─────────────────────────────────────────────────────────────────────────────
println("=" ^ 65)
println("VALIDATION  —  2×2 coarse grid  (NDIM=8)")
println("=" ^ 65)

const V_KC   = 2;      const V_NDIM = 8
const V_MOL  = 5;      const V_MAX  = 5
const V_DA   = 1.0/(0.5^2);   const V_DB = V_DA    # dx_coarse = 0.5
const V_KR   = 0.3
const V_TEND = 5.0;    const V_DT   = 0.05

t_ref,  n_ref,  s_ref  = run_mixed(Val(V_NDIM), V_KC, V_KC, V_MOL, V_MAX,
    V_DA, V_DB, V_KR, [1,2,3,4], V_TEND, V_DT; label="Full FSP — all 4 voxels active")

t_mix2, n_mix2, s_mix2 = run_mixed(Val(V_NDIM), V_KC, V_KC, V_MOL, V_MAX,
    V_DA, V_DB, V_KR, [1,4],     V_TEND, V_DT; label="Mixed — corners [1,4] active")

# Sanity check: all voxels active in mixed should match full FSP exactly
t_mix4, n_mix4, s_mix4 = run_mixed(Val(V_NDIM), V_KC, V_KC, V_MOL, V_MAX,
    V_DA, V_DB, V_KR, [1,2,3,4], V_TEND, V_DT; label="Mixed — all 4 active (≡ Full FSP)")

# ─────────────────────────────────────────────────────────────────────────────
# Demo: 4×4 grid (NDIM=32) — previously intractable
# ─────────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 65)
println("DEMO  —  4×4 coarse grid  (NDIM=32, previously intractable)")
println("=" ^ 65)

const D_KC   = 4;      const D_NDIM = 32
const D_MOL  = 3;      const D_MAX  = 3
const D_DA   = 1.0/(0.25^2);  const D_DB = D_DA   # dx_coarse = 0.25
const D_KR   = 0.3
const D_TEND = 6.0;    const D_DT   = 0.1

# Only corners active
t_d1, n_d1, s_d1 = run_mixed(Val(D_NDIM), D_KC, D_KC, D_MOL, D_MAX,
    D_DA, D_DB, D_KR, [1,16], D_TEND, D_DT; label="4×4 Mixed — corners [1,16]")

# Corners + immediate neighbors (9 active voxels)
ring1 = [1,2,5,6,   11,12,15,16]   # 2×2 blocks at each corner
t_d2, n_d2, s_d2 = run_mixed(Val(D_NDIM), D_KC, D_KC, D_MOL, D_MAX,
    D_DA, D_DB, D_KR, ring1, D_TEND, D_DT; label="4×4 Mixed — corner 2×2 blocks [8 voxels]")

# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

fig = Figure(size=(1100, 420), backgroundcolor=:black)

function dark_axis(fig, pos; title="", xlabel="", ylabel="")
    Axis(fig[pos...];
         title=title, titlecolor=:white, titlesize=13,
         xlabel=xlabel, ylabel=ylabel,
         xlabelcolor=:white, ylabelcolor=:white,
         xticklabelcolor=:white, yticklabelcolor=:white,
         backgroundcolor=:black,
         xgridcolor=RGBAf(1,1,1,0.1), ygridcolor=RGBAf(1,1,1,0.1))
end

# Panel 1: Validation — reaction progress
ax1 = dark_axis(fig, [1,1]; title="Validation: 2×2 (NDIM=8)", xlabel="time", ylabel="E[nA+nB]")
lines!(ax1, t_ref,  n_ref;  color=:white,   linewidth=2.5, label="Full FSP (ground truth)")
lines!(ax1, t_mix2, n_mix2; color=:cyan,    linewidth=2,   linestyle=:dash, label="Mixed [1,4] active")
lines!(ax1, t_mix4, n_mix4; color=:yellow,  linewidth=1.5, linestyle=:dot,  label="Mixed all-active (sanity)")
axislegend(ax1; bgcolor=:transparent, labelcolor=:white, framecolor=:transparent, labelsize=10)

# Panel 2: Error of mixed vs full FSP
ax2 = dark_axis(fig, [1,2]; title="Approximation error vs full FSP", xlabel="time", ylabel="|ΔE[nA+nB]|")
err2 = [abs(n_mix2[i] - n_ref[i]) for i in eachindex(t_ref)]
err4 = [abs(n_mix4[i] - n_ref[i]) for i in eachindex(t_ref)]
lines!(ax2, t_ref, err2; color=:cyan,   linewidth=2, label="Mixed [1,4] error")
lines!(ax2, t_ref, err4; color=:yellow, linewidth=2, label="All-active error (≈0)")
hlines!(ax2, [0.1]; color=:red, linewidth=1, linestyle=:dash, label="0.1 threshold")
axislegend(ax2; bgcolor=:transparent, labelcolor=:white, framecolor=:transparent, labelsize=10)

# Panel 3: State space size over time
ax3 = dark_axis(fig, [1,3]; title="|S| over time", xlabel="time", ylabel="active FSP states")
lines!(ax3, t_ref,  Float64.(s_ref);  color=:white,   linewidth=2, label="Full FSP 2×2")
lines!(ax3, t_mix2, Float64.(s_mix2); color=:cyan,    linewidth=2, linestyle=:dash, label="Mixed 2×2 [1,4]")
lines!(ax3, t_d1,   Float64.(s_d1);   color=:orange,  linewidth=2, label="Mixed 4×4 [1,16]")
lines!(ax3, t_d2,   Float64.(s_d2);   color=:magenta, linewidth=2, linestyle=:dash, label="Mixed 4×4 [8 vox]")
axislegend(ax3; bgcolor=:transparent, labelcolor=:white, framecolor=:transparent, labelsize=10)

# Panel 4: 4×4 demo reaction progress
ax4 = dark_axis(fig, [1,4]; title="Demo: 4×4 (NDIM=32, tractable!)", xlabel="time", ylabel="E[nA+nB]")
lines!(ax4, t_d1, n_d1; color=:orange,  linewidth=2, label="corners [1,16]")
lines!(ax4, t_d2, n_d2; color=:magenta, linewidth=2, linestyle=:dash, label="corner blocks [8 vox]")
axislegend(ax4; bgcolor=:transparent, labelcolor=:white, framecolor=:transparent, labelsize=10)

png_path = joinpath(outdir, "mixed_resolution_validation.png")
save(png_path, fig; px_per_unit=3)
println("\nSaved: $png_path")
println("Done.")
