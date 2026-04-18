"""
Adaptive coarsening — Step 4: mixed-grid generator.

Checks build_schlogl_mixed_system on K=4 with three masks:
  all-fine    mask=[0,0]  → K_mixed=4, should match build_schlogl_rdme_system
  all-coarse  mask=[1,1]  → K_mixed=2, should match build_schlogl_coarse_system
  mixed       mask=[0,1]  → K_mixed=3

For each case verifies:
  1. Generator builds without error
  2. Column sums are zero (probability conservation)
  3. Off-diagonal entries are non-negative
"""

using DiscStochSim
using SparseArrays
using LinearAlgebra
using Printf

D   = 0.005
K   = 4
h   = 1.0 / K
model = SchloglModel1D(D)
fine_grid  = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

# Build a small bounded state space for testing
n_max = 20
bc(s) = all(c -> 0 ≤ c ≤ n_max, Tuple(s))

function make_sp(KM)
    sp = StateSpace{CartesianIndex{KM}, Float64}()
    # seed with a few states near the low fixed point
    if KM == 4
        for a in 0:5, b in 0:5
            add_state!(sp, CartesianIndex(a, b, a, b), 1.0/(6^2))
        end
    elseif KM == 2
        for a in 0:10
            add_state!(sp, CartesianIndex(a, a), 1.0/11)
        end
    else  # KM == 3
        for a in 0:5, b in 0:10
            add_state!(sp, CartesianIndex(a, b, a+b), 1.0/(6*11))
        end
    end
    sp
end

function check_generator(label, sys, sp)
    A, = build_generator(sp, sys, Float64[], 0.0)
    @assert issparse(A)
    n = size(A, 1)
    col_sums = vec(sum(A, dims=1))
    max_col_err = maximum(abs.(col_sums))
    # Check off-diagonal entries are non-negative
    rows, cols, vals = findnz(A)
    min_offdiag = minimum(vals[rows .!= cols]; init=0.0)
    @printf("  %s  |S|=%d  max|col_sum|=%.2e  min_offdiag=%.2e\n",
            label, length(sp.states), max_col_err, min_offdiag)
    @assert max_col_err < 1e-10  "column sums not zero for $label"
    @assert min_offdiag ≥ -1e-14  "negative off-diagonal for $label"
    println("  PASS")
end

println("Mixed-grid generator checks  (K=$K, n_max=$n_max)")
println("─" ^ 55)

# all-fine
rates_fine  = build_schlogl_rdme_system(model, fine_grid)
rates_mixed0 = build_schlogl_mixed_system(model, fine_grid, BitVector([false,false]))
sp4 = make_sp(4)
expand!(sp4, rates_fine, bc; depth=1)
check_generator("Fine   (ref)", rates_fine,   sp4)
check_generator("Mixed [0,0]  ", rates_mixed0, sp4)

# all-coarse
rates_coarse  = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
rates_mixed1  = build_schlogl_mixed_system(model, fine_grid, BitVector([true,true]))
sp2 = make_sp(2)
expand!(sp2, rates_coarse, s->all(c->0≤c≤2n_max, Tuple(s)); depth=1)
check_generator("Coarse (ref)", rates_coarse, sp2)
check_generator("Mixed [1,1]  ", rates_mixed1, sp2)

# mixed
rates_mix = build_schlogl_mixed_system(model, fine_grid, BitVector([false,true]))
sp3 = make_sp(3)
expand!(sp3, rates_mix, s->all(c->0≤c≤2n_max, Tuple(s)); depth=1)
check_generator("Mixed [0,1]  ", rates_mix, sp3)

println("\nAll mixed-grid generator checks passed.")
