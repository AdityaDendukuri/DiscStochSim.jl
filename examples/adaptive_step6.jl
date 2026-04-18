"""
Adaptive coarsening — Step 6: selector unit tests.

Two independent concerns are tested separately:

  Part 1 — Mask correctness (K=8, no V-cycle):
    Verify that the extended select_coarsening_mask correctly identifies
    fine pairs for each of the three front-placement scenarios.

    Benchmark A: front INSIDE pair 2
      IC = (n_low, n_low, n_low, n_high, ...)
      α_intra[2] ≈ 0.68  →  mask = [true, false, true, true]

    Benchmark B: front ON the pair 2|3 boundary
      IC = (n_low, n_low, n_low, n_low, n_high, ...)
      β_right[2] ≈ 0.68  →  mask = [true, false, false, true]

    Benchmark C: homogeneous (all-low)
      IC = (n_low, ×8)
      all α_eff ≈ 0  →  mask = [true, true, true, true]

  Part 2 — V-cycle probability conservation (K=4, capped at n_max=20):
    Run one adaptive V-cycle step with each mask scenario and verify Σp≈1.
    Small n_max keeps state spaces tractable regardless of model fixed points.
"""

using DiscStochSim
using Printf

D   = 0.005
model     = SchloglModel1D(D)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp

println("="^60)
@printf("n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)
println("="^60)

# ═══════════════════════════════════════════════════════════════
# Part 1: mask correctness  (K=8)
# ═══════════════════════════════════════════════════════════════

println("\n── Part 1: Mask correctness (K=8) ──")

K8  = 8
h8  = 1.0 / K8
g8  = VoxelGrid(K8, h8, 0)
n8  = n_high + 20
r8  = build_schlogl_rdme_system(model, g8)
rm8 = RDMEModel1D(model.D, 0.0, 0.0)

function make_sp8(u0; depth=1)
    sp = StateSpace{CartesianIndex{K8}, Float64}()
    add_state!(sp, u0, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n8, Tuple(s))
    expand!(sp, r8, bc; depth)
    sp
end

function check_mask(label, u0, expected)
    sp   = make_sp8(u0; depth=1)
    mask = select_coarsening_mask(sp, rm8, g8; α_lo=0.10, α_hi=0.25)
    KM   = K8 - sum(mask)

    α, φ = pair_admissibility(sp, rm8, g8)
    μ = zeros(K8)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K8; μ[k] += p * t[k]; end
    end
    β = [abs(μ[2j]-μ[2j+1])/max(μ[2j]+μ[2j+1],1e-14) for j in 1:K8÷2-1]

    @printf("  %-42s  mask=%s  KM=%d\n", label, string(mask), KM)
    @printf("    α_intra=[%s]\n",
            join([@sprintf("%.2f",x) for x in α], " "))
    @printf("    β_right=[%s]\n",
            join([@sprintf("%.2f",x) for x in β], " "))

    if mask == expected
        println("    PASS")
    else
        println("    FAIL — expected $expected")
        error("mask mismatch: $label")
    end
end

check_mask("A — front inside pair 2",
    CartesianIndex(n_low, n_low, n_low, n_high, n_high, n_high, n_high, n_high),
    BitVector([true, false, true, true]))

check_mask("B — front on pair 2|3 boundary",
    CartesianIndex(n_low, n_low, n_low, n_low, n_high, n_high, n_high, n_high),
    BitVector([true, false, false, true]))

check_mask("C — homogeneous (all-low)",
    CartesianIndex(n_low, n_low, n_low, n_low, n_low, n_low, n_low, n_low),
    BitVector([true, true, true, true]))

println("  All mask checks PASSED.")

# ═══════════════════════════════════════════════════════════════
# Part 2: V-cycle probability conservation  (K=4, n_max_vc=20)
# ═══════════════════════════════════════════════════════════════

println("\n── Part 2: V-cycle probability conservation (K=4, n_max_vc=20) ──")

K4        = 4
h4        = 1.0 / K4
g4        = VoxelGrid(K4, h4, 0)
n_max_vc  = 20          # small cap → tractable Binomial expansion
r4        = build_schlogl_rdme_system(model, g4)

# Synthetic low/high counts within [0, n_max_vc] for Part 2.
# These are not the model fixed points — they're chosen to create clear
# within-pair and inter-pair asymmetry while staying tractable.
lo4 = 3
hi4 = 15

function make_sp4(u0; depth=1)
    sp = StateSpace{CartesianIndex{K4}, Float64}()
    add_state!(sp, u0, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max_vc, Tuple(s))
    expand!(sp, r4, bc; depth)
    sp
end

function check_vcycle(label, u0, exp_mask)
    sp = make_sp4(u0; depth=1)
    sp_new, mask = two_level_vcycle_adaptive(sp, model, g4, Float64[], 0.0, 0.01;
                                              krylov_m=20,
                                              expand_mixed=true,
                                              mixed_expand_depth=1,
                                              mixed_n_max=2*n_max_vc,
                                              binom_tol=1e-4)
    KM    = K4 - sum(mask)
    total = sum(sp_new.probs)
    @printf("  %-32s  mask=%s  |S_out|=%d  Σp=%.7f\n",
            label, string(mask), length(sp_new.states), total)
    @assert abs(total - 1.0) < 5e-3  "probability not conserved: $total"
    println("    PASS")
end

# A: front inside pair 1 — α_intra[1] large
check_vcycle("A4 — front inside pair 1",
    CartesianIndex(lo4, hi4, hi4, hi4),
    BitVector([false, true]))

# B: front on pair 1|2 boundary — β_right[1] large → both pairs fine → direct solve
check_vcycle("B4 — front on pair 1|2 boundary",
    CartesianIndex(lo4, lo4, hi4, hi4),
    BitVector([false, false]))

# C: homogeneous — all coarsenable
check_vcycle("C4 — homogeneous (all-low)",
    CartesianIndex(lo4, lo4, lo4, lo4),
    BitVector([true, true]))

println("  All V-cycle probability checks PASSED.")

println("\n" * "="^60)
println("All adaptive_step6 checks PASSED.")
println("="^60)
