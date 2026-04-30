"""
fig_criterion.jl

Two cases demonstrating the argmax-α merging criterion.

  (A) Birth-death K=8, IC=(5,0,0,0,0,0,0,0)  — 4 blocks
      Front in block 1; merged blocks 3-4 have 0 molecules → TV_adaptive≈0.

  (B) Schlögl K=4, IC=(n_lo,n_hi,n_hi,n_hi)  — 2 blocks
      With argmax+halo both blocks stay fine → TV_adaptive=0 exactly.

TV_naive  = TV(prolong(coarse_FSP), fine_FSP)  [from running coarse FSP]
TV_adaptive = sum of per-block TV for merged blocks only (fine-window blocks = 0)

Run:
  julia --project examples/fig_criterion.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

function save_fig(name, fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=2)
    println("  Saved: $name.{pdf,png}")
end

set_theme!(theme_light(); fontsize=11,
    Axis=(spinewidth=0.8, xtickwidth=0.7, ytickwidth=0.7,
          xgridcolor=(:black,0.06), ygridcolor=(:black,0.06)))

const BLUE   = Makie.wong_colors()[1]
const ORANGE = Makie.wong_colors()[2]
const GREEN  = Makie.wong_colors()[3]
const RED    = Makie.wong_colors()[6]
const GRAY   = RGBf(0.45,0.45,0.45)

const FLUX_TOL  = 1e-3
const PRUNE_TOL = 1e-6
const HALO      = 1

# ── utilities ─────────────────────────────────────────────────────────────────

function voxel_means(sp::StateSpace{CartesianIndex{K}}) where K
    μ = zeros(K)
    for (i,s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K; μ[k] += p*t[k]; end
    end
    μ
end

function pair_alpha(μ::Vector{Float64})
    n = length(μ)÷2
    [abs(μ[2j-1]-μ[2j]) / max(μ[2j-1]+μ[2j], 1e-12) for j in 1:n]
end

# TV between two state spaces on the same state type
function tv_between(sp_a::StateSpace{S}, sp_b::StateSpace{S}) where S
    total = 0.0
    for (i,s) in enumerate(sp_b.states)
        idx = get(sp_a.index, s, 0)
        pa  = idx > 0 ? sp_a.probs[idx] : 0.0
        total += abs(pa - sp_b.probs[i])
    end
    for (i,s) in enumerate(sp_a.states)
        get(sp_b.index, s, 0) == 0 && (total += abs(sp_a.probs[i]))
    end
    total / 2
end

# log C(n,k)/2^n without overflow
function log_binom_half(n::Int, k::Int)
    n == 0 && return 0.0
    k = min(k, n-k)
    s = 0.0
    for i in 1:k; s += log(n-k+i) - log(i); end
    s - n*log(2.0)
end

# 2D pair marginal from joint state
function pair_marginal(sp::StateSpace{CartesianIndex{K}}, j::Int) where K
    d = Dict{Tuple{Int,Int}, Float64}()
    for (i,s) in enumerate(sp.states)
        t = Tuple(s)
        key = (t[2j-1], t[2j])
        d[key] = get(d, key, 0.0) + sp.probs[i]
    end
    d
end

# TV(Bin(ñ,1/2), true 2D marginal) using 1 - overlap formula
function tv_binom_block(marg::Dict{Tuple{Int,Int},Float64})
    isempty(marg) && return 0.0
    ñ_prob = Dict{Int,Float64}()
    for ((n1,n2),p) in marg
        ñ = n1+n2; ñ_prob[ñ] = get(ñ_prob, ñ, 0.0)+p
    end
    overlap = 0.0
    for ((n1,n2),p) in marg
        ñ = n1+n2
        p_binom = ñ_prob[ñ] * exp(log_binom_half(ñ, n1))
        overlap += min(p, p_binom)
    end
    1.0 - overlap
end

# TV_adaptive: sum over merged blocks (fine-window blocks contribute 0)
function tv_adaptive_from_fine(sp::StateSpace{CartesianIndex{K}},
                                halo::Int=HALO) where K
    n_pairs = K÷2
    α       = pair_alpha(voxel_means(sp))
    j_star  = argmax(α)
    fl = max(1, j_star-halo); fh = min(n_pairs, j_star+halo)
    tv = 0.0
    for j in 1:n_pairs
        fl<=j<=fh && continue
        tv += tv_binom_block(pair_marginal(sp, j))
    end
    tv
end

# one FSP step: pre-expand so probability can flow, solve, prune
function fsp_step!(sp, sys, bc, rates, t, dt)
    expand!(sp, sys, bc; depth=1)
    A, = build_generator(sp, sys, rates, t)
    sp.probs .= max.(expv(dt, A, sp.probs; m=30), 0.0)
    renormalize!(sp); prune_threshold!(sp, PRUNE_TOL)
end

# ── main experiment runner ─────────────────────────────────────────────────────

struct RunResult
    times       :: Vector{Float64}
    alphas      :: Vector{Vector{Float64}}
    tv_naive    :: Vector{Float64}
    tv_adaptive :: Vector{Float64}
end

# TV_naive: one-shot restrict→prolong of sp_f itself.
# Measures "how much error does naive coarsening of the current distribution introduce?"
# No separate coarse FSP needed — avoids conflating merging error with FSP divergence.
function tv_naive_from_fine(sp::StateSpace{CartesianIndex{K}}) where K
    sp_c  = restrict(sp)
    sp_pr = prolong(sp_c, Val(K); binom_tol=1e-6); renormalize!(sp_pr)
    tv_between(sp_pr, sp)
end

function run_experiment(;
        fine_sys,
        sp_fine_ic :: StateSpace{CartesianIndex{K}},
        bc_fine, rates, dt, T, label) where K

    times    = Float64[0.0]
    α0       = pair_alpha(voxel_means(sp_fine_ic))
    tv_n0    = tv_naive_from_fine(sp_fine_ic)
    tv_a0    = tv_adaptive_from_fine(sp_fine_ic)
    alphas      = Vector{Float64}[α0]
    tv_naive    = Float64[tv_n0]
    tv_adaptive = Float64[tv_a0]

    @printf("  %s t=0  α=%s  TV_n=%.3f  TV_a=%.3f  |S|=%d\n",
            label, string(round.(α0;digits=2)), tv_n0, tv_a0,
            length(sp_fine_ic.states))

    sp_f = deepcopy(sp_fine_ic)
    t    = 0.0

    while t < T - 1e-10
        step = min(dt, T-t)
        fsp_step!(sp_f, fine_sys, bc_fine, rates, t, step)
        t += step

        α    = pair_alpha(voxel_means(sp_f))
        tv_n = tv_naive_from_fine(sp_f)
        tv_a = tv_adaptive_from_fine(sp_f)

        push!(times, t); push!(alphas, α)
        push!(tv_naive, tv_n); push!(tv_adaptive, tv_a)

        t % (5*dt) < dt+1e-9 &&
            @printf("  %s t=%.2f  α=%s  TV_n=%.3f  TV_a=%.3f  |Sf|=%d\n",
                    label, t, string(round.(α;digits=2)), tv_n, tv_a,
                    length(sp_f.states))
    end
    RunResult(times, alphas, tv_naive, tv_adaptive)
end

# ═══ Case A: Birth-death K=8 ══════════════════════════════════════════════════
println("\n━━━ Case A: Birth-Death K=8 ━━━")

BD_K=8; BD_D=1.0; BD_kb=2.0; BD_kd=1.0
BD_h=1/BD_K; BD_NMAX=18; BD_dt=0.05; BD_T=0.75

println("  building BD system..."); flush(stdout)
bd_model     = RDMEModel1D(BD_D, BD_kb, BD_kd)
bd_fine_grid = VoxelGrid(BD_K, BD_h, 0)
bd_fine_sys  = build_rdme_system(bd_model, bd_fine_grid)
println("  fine sys done"); flush(stdout)

sp_fine_bd = StateSpace{CartesianIndex{BD_K}, Float64}()
add_state!(sp_fine_bd, CartesianIndex(5,0,0,0,0,0,0,0), 1.0)

bd = run_experiment(;
    fine_sys   = bd_fine_sys,
    sp_fine_ic = sp_fine_bd,
    bc_fine    = s -> rdme_bc(s, BD_NMAX),
    rates=Float64[], dt=BD_dt, T=BD_T, label="BD")

# ═══ Case B: Schlögl K=4 ══════════════════════════════════════════════════════
println("\n━━━ Case B: Schlögl front K=4 ━━━")

SG_K=4; SG_D=0.005
SG_h=1/SG_K; SG_NMAX=200; SG_dt=0.2; SG_T=2.0

sg_model     = SchloglModel1D(SG_D)
sg_fine_grid = VoxelGrid(SG_K, SG_h, 0)
sg_fine_sys  = build_schlogl_rdme_system(sg_model, sg_fine_grid)

fp = schlogl_fixed_points(sg_model; n_max=SG_NMAX)
n_lo, _, n_hi = fp
@printf("  Fixed points: n_lo=%d  n_hi=%d\n", n_lo, n_hi)

sp_fine_sg = StateSpace{CartesianIndex{SG_K}, Float64}()
add_state!(sp_fine_sg, CartesianIndex(n_lo, n_hi, n_hi, n_hi), 1.0)

sg = run_experiment(;
    fine_sys   = sg_fine_sys,
    sp_fine_ic = sp_fine_sg,
    bc_fine    = s -> rdme_bc(s, SG_NMAX),
    rates=Float64[], dt=SG_dt, T=SG_T, label="SG")

# ═══ Figures ══════════════════════════════════════════════════════════════════
println("\n━━━ Plotting ━━━")

BLOCK_COLORS  = [BLUE, ORANGE, GREEN, RED]
BLOCK_STYLES  = [:solid, :dash, :dot, :dashdot]

function plot_alpha_panel!(ax, res; title="")
    ax.title = title
    ylims!(ax, -0.04, 1.08)
    ax.yticks = 0:0.25:1.0
    n = length(res.alphas[1])
    for j in 1:n
        αj = Float64[v[j] for v in res.alphas]
        lines!(ax, res.times, αj;
               color=BLOCK_COLORS[j], linewidth=2.0,
               linestyle=BLOCK_STYLES[j], label="block $j")
    end
    axislegend(ax; position=:rt, framevisible=false, labelsize=9, rowgap=2)
end

function plot_tv_panel!(ax, res; title="")
    ax.title = title
    ylims!(ax, -0.02, 1.08)
    lines!(ax, res.times, res.tv_naive;
           color=ORANGE, linewidth=2.0, label="naive  (merge all blocks)")
    lines!(ax, res.times, res.tv_adaptive;
           color=BLUE,   linewidth=2.0, label="adaptive  (argmax α_j ± 1)")
    axislegend(ax; position=:rt, framevisible=false, labelsize=9, rowgap=2)
end

fig1 = Figure(size=(860,340))
plot_alpha_panel!(Axis(fig1[1,1]; xlabel="time", ylabel="α_j  (imbalance)"), bd;
                  title="(a)  birth-death,  K=8,  IC=(5,0,…,0)")
plot_alpha_panel!(Axis(fig1[1,2]; xlabel="time", ylabel="α_j  (imbalance)"), sg;
                  title="(b)  Schlögl front,  K=4,  IC=(n_lo,n_hi,n_hi,n_hi)")
save_fig("fig_criterion", fig1)

fig2 = Figure(size=(860,340))
plot_tv_panel!(Axis(fig2[1,1]; xlabel="time", ylabel="TV vs full FSP"), bd;
               title="(a)  birth-death,  K=8,  IC=(5,0,…,0)")
plot_tv_panel!(Axis(fig2[1,2]; xlabel="time", ylabel="TV vs full FSP"), sg;
               title="(b)  Schlögl front,  K=4,  IC=(n_lo,n_hi,n_hi,n_hi)")
save_fig("fig_tv_criterion", fig2)

println("\nAll done.")
