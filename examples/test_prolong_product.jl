"""
Test: prolong_product vs prolong (Binomial) on a bistable interface pair.

K=2, pair = (voxel 1 near n_hi, voxel 2 near n_lo).
The pre-smooth marginals are concentrated at (n_hi, n_lo).
N_coarse = n_hi + n_lo = 193.

Expected:
  Binomial(193, 0.5): peak at n1 ≈ 97 (the WRONG transition state)
  Product convolution: peaks near (n_hi, n_lo) and (n_lo, n_hi) ← CORRECT
"""

using DiscStochSim, Printf

model  = SchloglModel1D(1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)
N_pair = n_lo + n_hi
@printf("Interface pair: N_SV = n_lo + n_hi = %d + %d = %d\n\n", n_lo, n_hi, N_pair)

# ── Build a pre-smooth fine state concentrated at (n_hi, n_lo) with spread ──
# Simulate a Poisson-ish distribution at each stable point

function poisson_weights(μ, n_max, weight_tol=1e-8)
    w = [exp(-μ + k*log(μ+1e-14) - sum(log(i) for i in 1:k; init=0.0))
         for k in 0:n_max]
    w ./= sum(w)
    filter(x -> x[2] > weight_tol, collect(enumerate(w)))
end

sp_pre = StateSpace{CartesianIndex{2}, Float64}()
# Pre-smooth state: product of Poisson(n_hi) × Poisson(n_lo) (both basins)
# Each voxel independently near its stable fixed point
for (i1, p1) in poisson_weights(Float64(n_hi), n_hi+40)
    for (i2, p2) in poisson_weights(Float64(n_lo), n_lo+30)
        w = p1 * p2
        w > 1e-10 || continue
        add_state!(sp_pre, CartesianIndex(i1-1, i2-1), w)
    end
end
@printf("Pre-smooth fine state: %d states, p_sum=%.4f\n", length(sp_pre), sum(sp_pre.probs))

# Check marginals of pre-smooth state
μ1 = sum(Tuple(ci)[1] * p for (ci,p) in zip(sp_pre.states, sp_pre.probs))
μ2 = sum(Tuple(ci)[2] * p for (ci,p) in zip(sp_pre.states, sp_pre.probs))
@printf("Pre-smooth means: μ1=%.1f (should≈n_hi=%d)  μ2=%.1f (should≈n_lo=%d)\n\n",
        μ1, n_hi, μ2, n_lo)

# ── Build coarse state: single state with N_pair = n_lo + n_hi ──
sp_coarse = StateSpace{CartesianIndex{1}, Float64}()
add_state!(sp_coarse, CartesianIndex(N_pair), 1.0)

# ── Binomial prolongation (current) ──
sp_binom = prolong(sp_coarse, Val(2); weight_tol=1e-8)
@printf("=== Binomial prolongation ===\n")
μ1_b = sum(Tuple(ci)[1] * p for (ci,p) in zip(sp_binom.states, sp_binom.probs))
μ2_b = sum(Tuple(ci)[2] * p for (ci,p) in zip(sp_binom.states, sp_binom.probs))
@printf("  Mean: (%.1f, %.1f)  expected (96.5, 96.5) — WRONG for interface pair\n\n",
        μ1_b, μ2_b)

# Top-5 most probable fine states
top5_b = sort(collect(zip(sp_binom.probs, sp_binom.states)); rev=true)[1:5]
println("  Top-5 states:")
for (p, ci) in top5_b
    @printf("    (n1=%3d, n2=%3d)  p=%.4f\n", Tuple(ci)[1], Tuple(ci)[2], p)
end

# ── Product convolution prolongation (new) ──
sp_prod = prolong_product(sp_pre, sp_coarse; weight_tol=1e-8)
@printf("\n=== Product convolution prolongation ===\n")
μ1_p = sum(Tuple(ci)[1] * p for (ci,p) in zip(sp_prod.states, sp_prod.probs))
μ2_p = sum(Tuple(ci)[2] * p for (ci,p) in zip(sp_prod.states, sp_prod.probs))
@printf("  Mean: (%.1f, %.1f)  expected (n_hi=%.0f, n_lo=%.0f) — CORRECT\n\n",
        μ1_p, μ2_p, Float64(n_hi), Float64(n_lo))
@printf("  States: %d total\n", length(sp_prod))

# Top-5 most probable fine states
top5_p = sort(collect(zip(sp_prod.probs, sp_prod.states)); rev=true)[1:5]
println("  Top-5 states:")
for (p, ci) in top5_p
    @printf("    (n1=%3d, n2=%3d)  p=%.4f\n", Tuple(ci)[1], Tuple(ci)[2], p)
end

# ── Compare n1 marginals ──
println("\n=== n1 marginals around the stable points ===")
@printf("  %-6s  %-12s  %-12s\n", "n1", "Binomial", "Product")
for n1 in [n_lo-5, n_lo, n_lo+5, n_hi÷2, n_hi-5, n_hi, n_hi+5]
    0 <= n1 <= N_pair || continue
    p_b = sum(p for (ci,p) in zip(sp_binom.states, sp_binom.probs) if Tuple(ci)[1]==n1; init=0.0)
    p_p = sum(p for (ci,p) in zip(sp_prod.states,  sp_prod.probs)  if Tuple(ci)[1]==n1; init=0.0)
    @printf("  n1=%-5d  %-12.5f  %-12.5f\n", n1, p_b, p_p)
end
