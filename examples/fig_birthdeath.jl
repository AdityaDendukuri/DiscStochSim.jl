"""
Birth-death-diffusion figure for talk (slide: linear network benchmark).

Model: K=4 voxels, single species A
  ∅ →^{k_b} A  (per voxel)
  A →^{k_d} ∅  (per molecule)
  A_k ⇌ A_{k±1}  rate d = D/h²

Claim demonstrated:
  - Coarse generator (K=2) is EXACT for linear networks → accurate profiles
  - State space reduction: K=4 fine (28k) → K=2 coarse (625) = 46×

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_birthdeath.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0; const MED = 9f0; const LW = 1.4f0
const C_FINE   = RGBf(0.000, 0.447, 0.698)   # blue
const C_COARSE = RGBf(0.835, 0.369, 0.000)   # orange
const C_ANAL   = RGBf(0.4, 0.4, 0.4)         # gray

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

# ── Model ─────────────────────────────────────────────────────────────────────
K      = 4
D      = 1.0; k_b = 2.0; k_d = 1.0
h      = 1.0 / K
N_MAX  = 12           # fine per-voxel truncation: 13^4 = 28,561 states
T_END  = 2.0
dt     = 0.05

model       = RDMEModel1D(D, k_b, k_d)
rates       = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

d_fine = D / h^2
@printf("d_fine=%.1f  d_coarse=%.1f  k_b=%.1f  k_d=%.1f\n",
        d_fine, D/(2h)^2, k_b, k_d)

# ── Systems ───────────────────────────────────────────────────────────────────
fine_sys   = build_rdme_system(model, fine_grid)
coarse_sys = build_coarse_system(model, coarse_grid, fine_grid)
bc_fine    = s -> rdme_bc(s, N_MAX)
bc_coarse  = s -> rdme_bc(s, 2*N_MAX)

# ── Initial condition: point mass at (5,0,0,0) ────────────────────────────────
u0_fine   = CartesianIndex(5, 0, 0, 0)
u0_coarse = CartesianIndex(5, 0)

# ── Fine FSP (fixed full state space = ground truth) ──────────────────────────
sp_fine = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:N_MAX, n2 in 0:N_MAX, n3 in 0:N_MAX, n4 in 0:N_MAX
    add_state!(sp_fine, CartesianIndex(n1,n2,n3,n4), 0.0)
end
sp_fine.probs[sp_fine.index[u0_fine]] = 1.0

A_fine, = build_generator(sp_fine, fine_sys, rates, 0.0)
S_fine  = length(sp_fine)
@printf("Fine   |S| = %d\n", S_fine)

# ── Coarse FSP (fixed full state space) ───────────────────────────────────────
sp_coarse = StateSpace{CartesianIndex{2}, Float64}()
for nc1 in 0:2*N_MAX, nc2 in 0:2*N_MAX
    add_state!(sp_coarse, CartesianIndex(nc1, nc2), 0.0)
end
sp_coarse.probs[sp_coarse.index[u0_coarse]] = 1.0

A_coarse, = build_generator(sp_coarse, coarse_sys, rates, 0.0)
S_coarse  = length(sp_coarse)
@printf("Coarse |S| = %d  (%.0f× reduction)\n", S_coarse, S_fine / S_coarse)

# ── Helpers ───────────────────────────────────────────────────────────────────
function voxel_means(sp, K)
    mu = zeros(K)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

function tv_between(sp_a::StateSpace{CartesianIndex{Ka}},
                    sp_b::StateSpace{CartesianIndex{Kb}}) where {Ka, Kb}
    tv = 0.0
    for (i, s) in enumerate(sp_b.states)
        pa = get(sp_a.index, s, 0)
        pb = sp_b.probs[i]
        tv += abs((pa != 0 ? sp_a.probs[pa] : 0.0) - pb)
    end
    for (i, s) in enumerate(sp_a.states)
        get(sp_b.index, s, 0) == 0 && (tv += abs(sp_a.probs[i]))
    end
    tv / 2
end

function marginal_voxel1(sp::StateSpace{CartesianIndex{K}}, n_max) where K
    pm = zeros(n_max + 1)
    for (i, s) in enumerate(sp.states)
        n = Tuple(s)[1]
        n ≤ n_max && (pm[n+1] += sp.probs[i])
    end
    pm
end

# ── Evolve ────────────────────────────────────────────────────────────────────
n_steps = round(Int, T_END / dt)
t_vals  = Float64[]
mu_fine_t   = Vector{Vector{Float64}}()
mu_coarse_t = Vector{Vector{Float64}}()
tv_t        = Float64[]
marg_v1_fine   = Vector{Vector{Float64}}()
marg_v1_coarse = Vector{Vector{Float64}}()

p_fine   = copy(sp_fine.probs)
p_coarse = copy(sp_coarse.probs)

push!(t_vals, 0.0)
push!(mu_fine_t,   voxel_means(sp_fine,   K))
sp_prol0 = prolong(sp_coarse, Val(K); weight_tol=1e-14, binom_tol=1e-6)
push!(mu_coarse_t, voxel_means(sp_prol0, K))
push!(tv_t, tv_between(sp_prol0, sp_fine))
push!(marg_v1_fine,   marginal_voxel1(sp_fine, 10))
push!(marg_v1_coarse, marginal_voxel1(sp_prol0, 10))

println("Evolving..."); flush(stdout)
for step in 1:n_steps
    p_fine   .= max.(expv(dt, A_fine,   p_fine;   m=30), 0.0)
    p_coarse .= max.(expv(dt, A_coarse, p_coarse; m=30), 0.0)
    p_fine   ./= sum(p_fine)
    p_coarse ./= sum(p_coarse)
    sp_fine.probs   .= p_fine
    sp_coarse.probs .= p_coarse

    sp_prol = prolong(sp_coarse, Val(K); weight_tol=1e-14, binom_tol=1e-6)
    renormalize!(sp_prol)

    t = step * dt
    push!(t_vals, t)
    push!(mu_fine_t,   voxel_means(sp_fine, K))
    push!(mu_coarse_t, voxel_means(sp_prol, K))
    push!(tv_t,        tv_between(sp_prol, sp_fine))
    push!(marg_v1_fine,   marginal_voxel1(sp_fine, 10))
    push!(marg_v1_coarse, marginal_voxel1(sp_prol, 10))
end
println("Done."); flush(stdout)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(780, 280))

# Panel (a): voxel means over time
ax_a = Axis(fig[1,1];
    xlabel="time",
    ylabel="mean count ⟨nₖ⟩",
    title="(a)  Spatial profiles  [K=4 fine vs K=2 coarse+prolong]",
    titlesize=MED)

voxel_colors = [RGBf(0.0, 0.45, 0.70),
                RGBf(0.0, 0.62, 0.45),
                RGBf(0.80, 0.47, 0.65),
                RGBf(0.94, 0.89, 0.26)]

mu_star = 2.0    # k_b / k_d
hlines!(ax_a, [mu_star]; color=(:black, 0.25), linewidth=0.8, linestyle=:dot)

for k in 1:K
    mf = [mu[k] for mu in mu_fine_t]
    mc = [mu[k] for mu in mu_coarse_t]
    lines!(ax_a, t_vals, mf; color=voxel_colors[k], linewidth=LW,
           label="voxel $k")
    lines!(ax_a, t_vals, mc; color=voxel_colors[k], linewidth=LW,
           linestyle=:dash)
end

# Legend: solid=fine, dashed=coarse
elem_fine   = LineElement(color=:black, linewidth=LW, linestyle=:solid)
elem_coarse = LineElement(color=:black, linewidth=LW, linestyle=:dash)
elem_ss     = LineElement(color=(:black,0.35), linewidth=0.8, linestyle=:dot)
Legend(fig[1,1], [elem_fine, elem_coarse, elem_ss],
       ["fine FSP  (K=4, $(S_fine) states)",
        "coarse+prolong  (K=2, $(S_coarse) states)",
        "steady state  μ*=k_b/k_d"],
       framevisible=false, labelsize=SMALL-1,
       tellwidth=false, tellheight=false,
       halign=:right, valign=:top, margin=(4,4,4,4))

# Panel (b): TV distance vs time
ax_b = Axis(fig[1,2];
    xlabel="time",
    ylabel="TV(coarse+prolong, fine FSP)",
    title="(b)  Accuracy  [coarse generator exact for linear networks]",
    titlesize=MED,
    yscale=log10)

lines!(ax_b, t_vals[2:end], tv_t[2:end]; color=C_COARSE, linewidth=LW)
scatter!(ax_b, t_vals[2:end], tv_t[2:end]; color=C_COARSE, markersize=3)

text!(ax_b, t_vals[end]*0.05, maximum(tv_t[2:end])*1.5;
      text="$(S_fine)÷$(S_coarse) = $(round(Int, S_fine/S_coarse))× state-space reduction",
      fontsize=SMALL-1, color=:black)

colgap!(fig.layout, 12)
CairoMakie.save(joinpath(OUTDIR, "fig_birthdeath.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_birthdeath.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_birthdeath.{pdf,png}")
