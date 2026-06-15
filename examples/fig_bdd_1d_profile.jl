"""
fig_bdd_1d_profile.jl

The abstract's FIRST experiment: a 1-D birth–death–diffusion chain.

    voxel 1:  ∅ →(λ) X      (localized immigration source)
    every voxel:  X →(k_d n) ∅      (death)
    edges:  X_k ⇌(D) X_{k±1}        (diffusion)

The network is LINEAR, so spatial coarsening is EXACT (Kemeny–Snell-exact lumping,
no closure error) — this isolates the *reduction* from any nonlinear closure error.

Two panels:
  (a) spatial profile ⟨n_k⟩ vs k: the spatially-coarsened multigrid solution
      (merge adjacent voxel pairs → binomial split prolongation, A_c = R A P,
      coarse stationary, prolong back) overlaid on the direct fine FSP solution.
  (b) |S| vs K: the realised fine FSP support and the coarse DOF count, against
      the naive per-voxel product baseline (n_max+1)^K.

Output: paper/figures/fig_bdd_1d_profile.png
"""

using DiscStochSim, Printf, LinearAlgebra, SparseArrays, Plots
using DiscStochSim: StateSpace, build_generator

# ── linear birth–death–diffusion on a 1-D chain ─────────────────────────────────
const λ    = 3.0      # immigration rate at the source voxel
const k_d  = 1.0      # death rate coefficient
const D    = 1.5      # diffusion coefficient (dx = 1)
const NMAX = 16       # per-voxel count truncation

build_sys(K) = build_graph_rdme(
    K, Val(1), chain_edges(K),
    v -> v == 1 ?
        LocalReaction{1}[ LocalReaction{1}(( 1,), nv -> λ),
                          LocalReaction{1}((-1,), nv -> k_d * nv[1]) ] :
        LocalReaction{1}[ LocalReaction{1}((-1,), nv -> k_d * nv[1]) ],
    (D,))

chain_bc(K) = x -> all(0 <= n <= NMAX for n in Tuple(x))

# ── drive the genuine joint FSP to steady state (no mean-field) ──────────────────
# Backward Euler with a large dt is unconditionally stable, so the LINEAR chain
# relaxes in a handful of steps. We stop when the spatial profile ⟨n_k⟩ stops
# moving: the realised support keeps filling Poisson tails afterwards, but those
# states carry negligible mass and do not change the profile.
function solve_fine(K; dt = 2.0, maxsteps = 40, tol = 5e-4)
    f = UnifiedFSP(build_sys(K), chain_bc(K);
                   dt = dt, θ = 1.0, expand_depth = 1, prune_tol = 1e-8)
    μ_prev = fill(-1.0, K)
    for _ in 1:maxsteps
        step!(f)
        μ = mean_counts(f)
        maximum(abs.(μ .- μ_prev)) < tol && break
        μ_prev = μ
    end
    f
end

# ── spatial pair-merge transfer: aggregate = total count of an adjacent voxel pair
#    R : fine → coarse   sum fine mass over each aggregate (Kemeny–Snell lumping)
#    P : coarse → fine   redistribute over the within-pair stationary conditional
#                        (the binomial/multinomial split), taken from p_ref so the
#                        operator is the closed-form linear-exact prolongation.
#    R·P = I exactly;  A_c = R A P is the Galerkin coarse generator.
pair_image(s::CartesianIndex{K}) where {K} =
    CartesianIndex(ntuple(i -> s[2i-1] + s[2i], K ÷ 2))

function build_pair_transfer(sp::StateSpace{CartesianIndex{K}}, p_ref) where {K}
    n_f   = length(sp.states)
    cidx  = Dict{CartesianIndex{K ÷ 2}, Int}()
    cprob = Float64[]
    fc    = Vector{Int}(undef, n_f)
    for (i, s) in enumerate(sp.states)
        c = pair_image(s)
        j = get(cidx, c, 0)
        if j == 0
            push!(cprob, 0.0); j = length(cprob); cidx[c] = j
        end
        fc[i] = j
        cprob[j] += p_ref[i]
    end
    n_c = length(cprob)
    R   = sparse(fc, 1:n_f, ones(n_f), n_c, n_f)
    pv  = [cprob[fc[i]] > 0 ? p_ref[i] / cprob[fc[i]] : 0.0 for i in 1:n_f]
    P   = sparse(1:n_f, fc, pv, n_f, n_c)
    (R, P, n_c)
end

# stationary law of a (near-)generator: solve A π = 0 with a normalization row
function stationary(A)
    n = size(A, 1)
    M = Matrix(A); M[n, :] .= 1.0
    b = zeros(n);  b[n] = 1.0
    M \ b
end

profile(states, p, K) =
    [sum(Tuple(s)[k] * p[i] for (i, s) in enumerate(states)) for k in 1:K]

# ── (a) profile recovery on a storable chain ────────────────────────────────────
const K_PROF = 8
@printf("Fine FSP on K=%d chain (λ=%.0f, k_d=%.1f, D=%.1f)…\n", K_PROF, λ, k_d, D)
f = solve_fine(K_PROF)
π_f   = copy(f.sp.probs); π_f ./= sum(π_f)
μ_dir = profile(f.sp.states, π_f, K_PROF)
n_f   = length(f.sp.states)

A_f, = build_generator(f.sp, f.sys, Float64[], f.t)
R, P, n_c = build_pair_transfer(f.sp, π_f)
A_c   = R * A_f * P
π_c   = stationary(A_c)
π_mg  = P * π_c
μ_mg  = profile(f.sp.states, π_mg, K_PROF)

rp_err  = maximum(abs.(Matrix(R * P) - I))
prof_err = norm(μ_mg .- μ_dir) / norm(μ_dir)
@printf("  fine |S| = %d  →  coarse n_c = %d   (%.1f× fewer DOF)\n", n_f, n_c, n_f / n_c)
@printf("  max|R·P − I| = %.2e\n", rp_err)
@printf("  ‖profile_mg − profile_direct‖/‖direct‖ = %.2e\n", prof_err)

# ── (b) |S|-vs-K reduction sweep ────────────────────────────────────────────────
Ks = collect(4:2:12)
S_fine = Int[]; S_coarse = Int[]; S_naive = Float64[]
for K in Ks
    fK = solve_fine(K)
    πK = copy(fK.sp.probs); πK ./= sum(πK)
    _, _, ncK = build_pair_transfer(fK.sp, πK)
    push!(S_fine,   length(fK.sp.states))
    push!(S_coarse, ncK)
    push!(S_naive,  Float64(NMAX + 1)^K)
    @printf("  K=%2d   naive=(%d)^K=%.1e   fine|S|=%d   coarse n_c=%d\n",
            K, NMAX + 1, S_naive[end], S_fine[end], S_coarse[end])
end

# ── figure ──────────────────────────────────────────────────────────────────────
p1 = plot(1:K_PROF, μ_dir; lw = 2, color = :black, label = "direct FSP",
          xlabel = "voxel k", ylabel = "⟨n_k⟩",
          title = @sprintf("profile recovery  (rel err %.1e)", prof_err))
scatter!(p1, 1:K_PROF, μ_mg; marker = :star5, ms = 9, color = :crimson,
         label = @sprintf("multigrid (%d→%d DOF)", n_f, n_c))

p2 = plot(Ks, S_naive; yscale = :log10, lw = 2, ls = :dot, color = :gray,
          marker = :diamond, ms = 5, label = "per-voxel product (n_max+1)^K",
          xlabel = "chain length K", ylabel = "state count",
          title = "state-space reduction", legend = :topleft)
plot!(p2, Ks, S_fine;   lw = 2, color = :black,   marker = :circle, ms = 5,
      label = "fine adaptive FSP |S|")
plot!(p2, Ks, S_coarse; lw = 2, color = :crimson, marker = :star5,  ms = 7,
      label = "spatially-coarsened n_c")

plt = plot(p1, p2; layout = (1, 2), size = (1000, 400), left_margin = 5Plots.mm,
           bottom_margin = 5Plots.mm)
mkpath(joinpath(@__DIR__, "..", "paper", "figures"))
out = joinpath(@__DIR__, "..", "paper", "figures", "fig_bdd_1d_profile.png")
savefig(plt, out)
@printf("\nwrote %s\n", out)
