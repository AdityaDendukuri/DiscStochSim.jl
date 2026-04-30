"""
Adaptive coarsening — Step 2: partial_restrict.

Verifies that partial restriction produces the correct mixed state space
for three scenarios:

  all-coarse  mask=[true,true]    K=4 → K_mixed=2  (standard restrict)
  all-fine    mask=[false,false]  K=4 → K_mixed=4  (identity)
  mixed       mask=[false,true]   K=4 → K_mixed=3  (pair1 fine, pair2 coarse)

Checks:
  1. Output dimension K_mixed is correct
  2. Total probability mass is preserved
  3. Mixed case: fine pair {1,2} dimensions survive intact;
     coarse pair {3,4} is summed to nc=n3+n4
"""

using DiscStochSim
using Printf

# Build a small K=4 state space with known states
function make_sp4()
    sp = StateSpace{CartesianIndex{4}, Float64}()
    # States: (1,2,3,4), (2,1,3,4), (1,2,4,3), (2,1,4,3) — equal weight
    for s in [(1,2,3,4),(2,1,3,4),(1,2,4,3),(2,1,4,3)]
        add_state!(sp, CartesianIndex(s), 0.25)
    end
    sp
end

sp4 = make_sp4()
println("Fine K=4 states:")
for i in eachindex(sp4.states)
    @printf("  %s  p=%.4f\n", Tuple(sp4.states[i]), sp4.probs[i])
end
println()

# ── test 1: all-coarse ─────────────────────────────────────────────────────
mask_all = BitVector([true, true])
sp_c = partial_restrict(sp4, mask_all)
println("All-coarse mask=[1,1]  →  K_mixed=$(length(Tuple(sp_c.states[1])))")
for i in eachindex(sp_c.states)
    @printf("  %s  p=%.4f\n", Tuple(sp_c.states[i]), sp_c.probs[i])
end
@assert length(Tuple(sp_c.states[1])) == 2  "expected K_mixed=2"
@assert abs(sum(sp_c.probs) - 1.0) < 1e-12  "mass not conserved"
# All states map to nc1=3, nc2=7
@assert length(sp_c.states) == 1  "all four states collapse to one coarse state"
@assert Tuple(sp_c.states[1]) == (3, 7)  "expected coarse state (3,7)"
println("  PASS\n")

# ── test 2: all-fine ──────────────────────────────────────────────────────
mask_none = BitVector([false, false])
sp_f = partial_restrict(sp4, mask_none)
println("All-fine  mask=[0,0]   →  K_mixed=$(length(Tuple(sp_f.states[1])))")
for i in eachindex(sp_f.states)
    @printf("  %s  p=%.4f\n", Tuple(sp_f.states[i]), sp_f.probs[i])
end
@assert length(Tuple(sp_f.states[1])) == 4  "expected K_mixed=4"
@assert abs(sum(sp_f.probs) - 1.0) < 1e-12  "mass not conserved"
@assert length(sp_f.states) == 4  "all four states should be distinct"
println("  PASS\n")

# ── test 3: mixed — pair1 fine, pair2 coarse ─────────────────────────────
mask_mix = BitVector([false, true])
sp_m = partial_restrict(sp4, mask_mix)
println("Mixed     mask=[0,1]   →  K_mixed=$(length(Tuple(sp_m.states[1])))")
for i in eachindex(sp_m.states)
    @printf("  %s  p=%.4f\n", Tuple(sp_m.states[i]), sp_m.probs[i])
end
@assert length(Tuple(sp_m.states[1])) == 3  "expected K_mixed=3"
@assert abs(sum(sp_m.probs) - 1.0) < 1e-12  "mass not conserved"
# States (1,2,3,4) and (1,2,4,3) both → (1,2,7);  (2,1,3,4) and (2,1,4,3) → (2,1,7)
@assert length(sp_m.states) == 2  "expected 2 mixed states"
mixed_tuples = sort([Tuple(s) for s in sp_m.states])
@assert mixed_tuples == [(1,2,7),(2,1,7)]  "expected mixed states (1,2,7) and (2,1,7)"
println("  PASS\n")

println("All partial_restrict tests passed.")
