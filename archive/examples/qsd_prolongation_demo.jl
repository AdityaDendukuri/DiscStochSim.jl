"""
QSD Prolongation Demo
=====================
Validates `build_qsd_table` / `prolong_qsd` against the Binomial baseline and
illustrates where the two diverge.

Three regimes are shown for a symmetric 2-voxel patch:

  1. Pure diffusion (no reactions)         → QSD ≡ Binomial(N, ½)
  2. Birth-death + fast diffusion          → QSD ≈ Binomial  (diffusion dominates)
  3. Schlögl (bistable) + slow diffusion   → QSD ≠ Binomial  (bimodal within pair)
"""

using DiscStochSim
using LinearAlgebra, Distributions

# ── helpers ───────────────────────────────────────────────────────────────────

binom_dist(N) = [pdf(Binomial(N, 0.5), n) for n in 0:N]

tv_distance(p, q) = 0.5 * sum(abs.(p .- q))

function print_comparison(label, N, qsd_col, binom_col)
    tv = tv_distance(qsd_col, binom_col)
    println("  $label  N=$N   TV(QSD, Binomial) = $(round(tv; digits=6))")
end

# ─────────────────────────────────────────────────────────────────────────────
# 1. Pure diffusion — must recover Binomial exactly
# ─────────────────────────────────────────────────────────────────────────────
println("=== 1. Pure diffusion (d=10.0) ===")
qt_pure = build_qsd_table(20, 10.0)

for N in [1, 5, 10, 20]
    qsd_col   = qt_pure.table[1:N+1, N+1]
    binom_col = binom_dist(N)
    print_comparison("pure diffusion ", N, qsd_col, binom_col)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Birth-death  ∅ →^λ A →^μ ∅  with fast diffusion
#    absorbing_pos(n) = λ,  absorbing_neg(n) = μ·n
#    Both channels are state-independent per voxel for birth; state-dependent for death.
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== 2. Birth-death (λ=2, μ=1) + fast diffusion (d=50) ===")
λ, μ = 2.0, 1.0
pp_bd = PatchPropensity(
    _ -> 0.0,       # no N-conserving forward
    _ -> 0.0,       # no N-conserving reverse
    _ -> λ,         # birth: abs rate per voxel (independent of n)
    n -> μ * n,     # death: abs rate per voxel (linear in n)
)
qt_bd = build_qsd_table(20, 50.0; propensity=pp_bd)

for N in [2, 5, 10]
    qsd_col   = qt_bd.table[1:N+1, N+1]
    binom_col = binom_dist(N)
    print_comparison("birth-death    ", N, qsd_col, binom_col)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. Schlögl model + slow diffusion (d=0.01)
#    Reactions per voxel:
#      ∅ →^c3  X                    (birth, absorbing_pos = c3)
#      X →^c4  ∅                    (death, absorbing_neg(n) = c4·n)
#      2X →^{c1/2} 3X               (autocatalysis, absorbing_neg(n) = c1·n(n-1)/2 ... net +1)
#      3X →^{c2/6} 2X               (suppression,   absorbing_pos(n) = c2·n(n-1)(n-2)/6 ... net -1)
#
# Note: autocatalysis increases N by 1 (absorbing_pos side),
#       suppression decreases N by 1 (absorbing_neg side).
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== 3. Schlögl bistable + slow diffusion (d=0.01, rescaled to n_max≈15) ===")

# Rescaled parameters so bistability occurs around n ≈ 5–10 (tractable for demo)
c1, c2, c3, c4 = 0.18, 2.5e-3, 2.2, 1.0
d_slow = 0.01

pp_schlogl = PatchPropensity(
    _ -> 0.0,
    _ -> 0.0,
    n -> c3 + c1 * n * max(n-1, 0) / 2,           # birth + autocatalysis (N increases)
    n -> c4 * n + c2 * n * max(n-1,0) * max(n-2,0) / 6,  # death + suppression (N decreases)
)
qt_schlogl = build_qsd_table(30, d_slow; propensity=pp_schlogl)

for N in [5, 10, 15, 20]
    qsd_col   = qt_schlogl.table[1:N+1, N+1]
    binom_col = binom_dist(N)
    tv = tv_distance(qsd_col, binom_col)
    println("  Schlögl slow diff  N=$N   TV(QSD, Binomial) = $(round(tv; digits=6))")
    # Show the QSD itself for N=10
    if N == 10
        println("    QSD  : $(round.(qsd_col; digits=3))")
        println("    Binom: $(round.(binom_col; digits=3))")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. Schlögl + fast diffusion — should recover Binomial
# ─────────────────────────────────────────────────────────────────────────────
println("\n=== 4. Schlögl bistable + FAST diffusion (d=100) — should → Binomial ===")
qt_schlogl_fast = build_qsd_table(30, 100.0; propensity=pp_schlogl)

for N in [5, 10, 15]
    qsd_col   = qt_schlogl_fast.table[1:N+1, N+1]
    binom_col = binom_dist(N)
    tv = tv_distance(qsd_col, binom_col)
    println("  Schlögl fast diff  N=$N   TV(QSD, Binomial) = $(round(tv; digits=6))")
end

println("\nDone.")
