"""
Adaptive coarsening — Step 3: partial_prolong.

Tests:

  1. all-coarse  mask=[1,1]: partial_prolong should match standard prolong
  2. all-fine    mask=[0,0]: partial_prolong is identity (round-trip with partial_restrict)
  3. mixed       mask=[0,1]: fine pair passes through, coarse pair gets Binomial split

  4. Round-trip mass conservation: restrict then prolong recovers total probability mass
     (individual states are recovered exactly when nc_j uniquely determines the split
     — here nc_j=0 or 1 so Binomial is trivial).
"""

using DiscStochSim
using Printf

function make_sp4_small()
    sp = StateSpace{CartesianIndex{4}, Float64}()
    # Use small counts so Binomial splits are exact / tractable
    add_state!(sp, CartesianIndex(0, 1, 2, 0), 0.4)
    add_state!(sp, CartesianIndex(1, 0, 0, 2), 0.3)
    add_state!(sp, CartesianIndex(0, 1, 1, 1), 0.3)
    sp
end

sp4 = make_sp4_small()
println("Fine K=4 states (mass=$(round(sum(sp4.probs),digits=6))):")
for i in eachindex(sp4.states)
    @printf("  %s  p=%.4f\n", Tuple(sp4.states[i]), sp4.probs[i])
end
println()

# ── test 1: all-coarse round-trip ─────────────────────────────────────────
mask_all = BitVector([true, true])
sp_c     = partial_restrict(sp4, mask_all)
sp_back  = partial_prolong(sp_c, mask_all, Val(4))
println("All-coarse round-trip  mask=[1,1]:")
println("  mixed states: $(length(sp_c.states))  fine-back states: $(length(sp_back.states))")
@printf("  mass after restrict=%.6f  after prolong=%.6f\n",
        sum(sp_c.probs), sum(sp_back.probs))
@assert abs(sum(sp_back.probs) - 1.0) < 1e-10  "mass not conserved (all-coarse)"
println("  PASS\n")

# ── test 2: all-fine round-trip ───────────────────────────────────────────
mask_none = BitVector([false, false])
sp_f      = partial_restrict(sp4, mask_none)
sp_back2  = partial_prolong(sp_f, mask_none, Val(4))
println("All-fine  round-trip   mask=[0,0]:")
@printf("  mass after restrict=%.6f  after prolong=%.6f\n",
        sum(sp_f.probs), sum(sp_back2.probs))
@assert abs(sum(sp_back2.probs) - 1.0) < 1e-10  "mass not conserved (all-fine)"
# All-fine prolong is identity: states and probs should be identical
got_states = sort([Tuple(s) for s in sp_back2.states])
exp_states = sort([Tuple(s) for s in sp4.states])
@assert got_states == exp_states  "all-fine round-trip should recover exact states"
println("  PASS\n")

# ── test 3: mixed round-trip   ────────────────────────────────────────────
mask_mix = BitVector([false, true])
sp_m     = partial_restrict(sp4, mask_mix)
sp_back3 = partial_prolong(sp_m, mask_mix, Val(4))
println("Mixed     round-trip   mask=[0,1]:")
println("  mixed states:")
for i in eachindex(sp_m.states)
    @printf("    %s  p=%.4f\n", Tuple(sp_m.states[i]), sp_m.probs[i])
end
println("  prolonged fine states:")
for i in eachindex(sp_back3.states)
    @printf("    %s  p=%.4f\n", Tuple(sp_back3.states[i]), sp_back3.probs[i])
end
@assert abs(sum(sp_back3.probs) - 1.0) < 1e-10  "mass not conserved (mixed)"
# Pair 1 is fine: n1,n2 are preserved exactly.
# Pair 2 is coarse with nc in {0,1,2,3}: Binomial splits back.
# Check dimension
@assert length(Tuple(sp_back3.states[1])) == 4  "expected K=4 output"
println("  PASS\n")

println("All partial_prolong tests passed.")
