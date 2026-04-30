"""
Adaptive coarsening — Step 5: two_level_vcycle_adaptive.

Tests on K=4 (2 pairs) with tight truncation to keep state spaces small.
Verifies:
  1. Returns (sp_new, mask) with probability conserved
  2. Output state space has CartesianIndex{K} type
  3. Hysteresis: mask is stable across multiple steps
"""

using DiscStochSim
using Printf

D   = 0.005
K   = 4
h   = 1.0 / K
model = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)

# Initial distribution near the low fixed point
n0 = 5
sp = StateSpace{CartesianIndex{K}, Float64}()
add_state!(sp, CartesianIndex(ntuple(_ -> n0, K)), 1.0)

rates_fine = build_schlogl_rdme_system(model, fine_grid)
bc = s -> all(c -> 0 ≤ c ≤ 30, Tuple(s))
expand!(sp, rates_fine, bc; depth = 2)
println("Initial fine state space size: $(length(sp.states))")

dt = 0.01

# ── Single step ───────────────────────────────────────────────────────────────
sp_new, mask = two_level_vcycle_adaptive(sp, model, fine_grid, Float64[], 0.0, dt;
                                          krylov_m = 30,
                                          expand_mixed = true,
                                          mixed_expand_depth = 1,
                                          mixed_n_max = 60,
                                          binom_tol = 1e-4)

KM = K - sum(mask)
total_prob = sum(sp_new.probs)
@printf("Mask: %s  →  KM=%d\n", string(mask), KM)
@printf("Output fine state space size: %d\n", length(sp_new.states))
@printf("Total probability after step: %.8f\n", total_prob)

@assert abs(total_prob - 1.0) < 1e-5   "probability not conserved: $total_prob"
@assert length(sp_new.states) > 0       "empty output state space"
@assert eltype(sp_new.states) == CartesianIndex{K}   "wrong output dimension"
println("PASS: single adaptive V-cycle step")

# ── Hysteresis loop ───────────────────────────────────────────────────────────
println("\nRunning 5 V-cycle steps with mask hysteresis...")
prev_mask = mask
sp_cur = sp_new
for step in 1:5
    global sp_cur, prev_mask
    sp_cur, prev_mask = two_level_vcycle_adaptive(sp_cur, model, fine_grid, Float64[], step*dt, dt;
                                                    prev_mask,
                                                    expand_mixed = true,
                                                    mixed_expand_depth = 1,
                                                    mixed_n_max = 60,
                                                    binom_tol = 1e-4)
    @printf("  step %d  |S|=%4d  mask=%s  Σp=%.6f\n",
            step, length(sp_cur.states), string(prev_mask), sum(sp_cur.probs))
    @assert abs(sum(sp_cur.probs) - 1.0) < 1e-4  "step $step: probability not conserved"
end
println("\nAll adaptive V-cycle steps passed.")
