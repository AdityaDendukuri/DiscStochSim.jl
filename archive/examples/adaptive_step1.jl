"""
Adaptive coarsening — Step 1: pair admissibility and hysteretic mask.

Tests on three Schlögl ICs whose expected behaviour is known:

  IC-A  K=2  (31, 162)              → pair 1 asymmetric → keep fine
  IC-B  K=4  (31,162,162,162)       → pair 1 front, pair 2 homogeneous
  IC-C  K=4  (162,162,162,162)      → both pairs homogeneous → both coarsenable

Also tests hysteresis: a pair in the band [α_lo, α_hi] should keep
its previous mask entry rather than flip.
"""

using DiscStochSim
using Printf

D   = 0.005
K4  = 4
h4  = 1.0 / K4
K2  = 2
h2  = 1.0 / K2

model4 = SchloglModel1D(D)
grid4  = VoxelGrid(K4, h4, 0)
model2 = SchloglModel1D(D)
grid2  = VoxelGrid(K2, h2, 0)

fp = schlogl_fixed_points(model4)
n_low, _, n_high = fp

function point_sp(u0)
    K = length(Tuple(u0))
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    sp
end

function show_diagnostics(label, sp, model, grid; prev_mask=nothing)
    K      = length(Tuple(sp.states[1]))
    α, φn  = pair_admissibility(sp, model, grid)
    mask   = select_coarsening_mask(sp, model, grid, prev_mask)
    println(label)
    for j in 1:K÷2
        k1, k2 = 2j-1, 2j
        status = mask[j] ? "COARSEN" : "keep FINE"
        prev   = isnothing(prev_mask) ? "—" : (prev_mask[j] ? "coarsen" : "fine")
        @printf("  pair {%d,%d}  α=%.3f  φ_norm=%.2f  prev=%-7s → %s\n",
                k1, k2, α[j], φn[j], prev, status)
    end
    println()
    mask
end

println("=" ^ 60)
println("Pair admissibility diagnostics  (δ defaults: α_lo=0.10, α_hi=0.25)")
println("=" ^ 60)
println()

show_diagnostics(
    "IC-A  K=2  (n_low=$(n_low), n_high=$(n_high)) — expect: keep fine",
    point_sp(CartesianIndex(n_low, n_high)), model2, grid2)

show_diagnostics(
    "IC-B  K=4  (n_low, n_high, n_high, n_high) — expect: pair 1 fine, pair 2 coarsen",
    point_sp(CartesianIndex(n_low, n_high, n_high, n_high)), model4, grid4)

show_diagnostics(
    "IC-C  K=4  (n_high, n_high, n_high, n_high) — expect: both coarsen",
    point_sp(CartesianIndex(n_high, n_high, n_high, n_high)), model4, grid4)

# Hysteresis test: pair with α in [α_lo, α_hi] should stick to prev_mask
println("─" ^ 60)
println("Hysteresis test: α in the band keeps previous decision")
println("─" ^ 60)
# Construct a pair with α ≈ 0.15 (between 0.10 and 0.25)
# E[n1]≈60, E[n2]≈100 → total=160, α=40/160=0.25 (edge of band)
# Use a point IC: (60, 100) for K=2
sp_mid = point_sp(CartesianIndex(60, 100))
println("  α ≈ $(round(abs(60-100)/(60+100), digits=3)) for (60,100)")
mask_prev_fine   = BitVector([false])
mask_prev_coarse = BitVector([true])
show_diagnostics("  prev=fine   → should stay fine",   sp_mid, model2, grid2; prev_mask=mask_prev_fine)
show_diagnostics("  prev=coarse → should stay coarse", sp_mid, model2, grid2; prev_mask=mask_prev_coarse)
