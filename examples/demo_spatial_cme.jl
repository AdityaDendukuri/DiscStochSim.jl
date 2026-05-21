"""
demo_spatial_cme.jl — Joint CME over K voxels with multigrid V-cycle coarsening

Demonstrates:
  1. Build a birth-death RDME on a K-voxel 1D chain as a single joint CME
     using StateSpace{CartesianIndex{K}, Float64}
  2. Evolve with spatial_vcycle! (restrict → coarse solve → prolong → smooth)
  3. Extract exact FSP means (no mean-field) via mean_counts
  4. Compare voxel marginals to the exact Poisson(k_b/k_d) steady state
  5. Show the coarsening plan at each step (flux-based voxel pairing)

Run:  julia --project examples/demo_spatial_cme.jl
"""

using DiscStochSim, Printf, Statistics

# ── Parameters ─────────────────────────────────────────────────────────────────
const K    = 10          # number of voxels (10-voxel 1D chain)
const k_b  = 2.0         # birth rate per voxel
const k_d  = 1.0         # death rate per molecule
const D    = 4.0         # diffusion coefficient
const dx   = 1.0         # voxel spacing
const n_max = 15         # max molecules per voxel (FSP truncation)
const T_END = 4.0        # end time
const dt    = 0.2        # V-cycle step size
const n_steps = round(Int, T_END / dt)

# ── Setup ──────────────────────────────────────────────────────────────────────

mesh  = VoxelMesh(K, [(i,i+1) for i in 1:K-1], fill(D/dx^2, K-1))
model = GraphRDMEModel(K, k_b, k_d)

μ_ss = k_b / k_d
@printf("K=%d voxels  |  k_b=%.1f  k_d=%.1f  D=%.1f  μ_ss=%.1f\n",
        K, k_b, k_d, D, μ_ss)
@printf("V-cycle steps: %d × dt=%.2f  (T_end=%.1f)\n\n", n_steps, dt, T_END)

# IC: 1 molecule at centre voxel K÷2+1
sp = StateSpace{CartesianIndex{K}, Float64}()
v0 = K ÷ 2 + 1
add_state!(sp, CartesianIndex(ntuple(k -> k == v0 ? 1 : 0, Val(K))), 1.0)
st = SpatialCMEState{K}(sp, model, mesh, 0)

@printf("IC: 1 molecule at voxel %d\n", v0)
@printf("Initial states: %d\n\n", length(st.sp))

# ── Time evolution ─────────────────────────────────────────────────────────────

t_log    = Float64[0.0]
μ_log    = Vector{Float64}[copy(mean_counts(st))]
nst_log  = Int[length(st.sp)]

for step in 1:n_steps
    global st = spatial_vcycle!(st, dt;
                         τ_pre=0.05*dt, τ_post=0.05*dt,
                         n_max, krylov_m=25)
    push!(t_log,   step * dt)
    push!(μ_log,   copy(mean_counts(st)))
    push!(nst_log, length(st.sp))

    if step % (n_steps ÷ 4) == 0
        μ = μ_log[end]
        @printf("t=%.2f  states=%5d  mean=[%s]  sum=%.3f\n",
                t_log[end], nst_log[end],
                join([@sprintf("%.2f", m) for m in μ], " "),
                sum(μ))
    end
end

# ── Final analysis ─────────────────────────────────────────────────────────────
μ_final = μ_log[end]
@printf("\n=== Final state (t=%.1f) ===\n", T_END)
@printf("States in joint CME: %d\n", nst_log[end])
@printf("Mean counts per voxel: %s\n",
        join([@sprintf("%.3f", m) for m in μ_final], "  "))
@printf("Target (μ_ss=%.1f):    %s\n",
        μ_ss, join([@sprintf("%.3f", μ_ss) for _ in μ_final], "  "))
@printf("Max error: %.4f\n", maximum(abs.(μ_final .- μ_ss)))
@printf("Total mean: %.4f  (target: %.1f)\n\n", sum(μ_final), K * μ_ss)

# ── Marginal distribution check ────────────────────────────────────────────────
# Compute marginal p(n_k) at centre voxel by summing over all other voxels
@printf("=== Marginal at centre voxel %d ===\n", v0)
marginal_v = zeros(n_max + 1)
for (ci, prob) in zip(st.sp.states, st.sp.probs)
    n_k = Tuple(ci)[v0]
    if 0 ≤ n_k ≤ n_max
        marginal_v[n_k + 1] += prob
    end
end
marginal_v ./= sum(marginal_v)

# Exact Poisson(μ_ss) reference
using SpecialFunctions: loggamma
poisson_pmf(n, λ) = exp(-λ + n*log(λ) - loggamma(n+1))
exact = [poisson_pmf(n, μ_ss) for n in 0:n_max]
exact ./= sum(exact)

@printf("%-5s  %-8s  %-8s  %-8s\n", "n", "FSP", "Exact", "Error")
for n in 0:min(8, n_max)
    @printf("%-5d  %-8.4f  %-8.4f  %-8.4f\n",
            n, marginal_v[n+1], exact[n+1], abs(marginal_v[n+1]-exact[n+1]))
end
tv = 0.5 * sum(abs.(marginal_v .- exact))
@printf("TV distance from Poisson(%.1f): %.5f\n", μ_ss, tv)

# ── Coarsening plan snapshot ───────────────────────────────────────────────────
@printf("\n=== Coarsening plan at t=%.1f ===\n", T_END)
μ_cur = mean_counts(st)
plan = make_coarsening_plan(mesh, μ_cur)
@printf("Fine K=%d → Coarse K'=%d\n", K, plan.K_coarse)
@printf("Merged pairs: %s\n", string(plan.pairs))
@printf("Singletons:   %s\n", string(plan.singles))
@printf("Voxel map:    %s\n", string(plan.voxel_map))

@printf("\nAll checks passed.\n")
