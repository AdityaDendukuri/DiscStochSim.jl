# K=2 Schlögl full FSP: run to long time to see bimodal development
# Full state space (all n1,n2 from 0 to n_max) — no adaptive truncation

using DiscStochSim
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Printf

D = 0.005; K = 2; h = 1.0/K; n_max = 200
model = SchloglModel1D(D)
fine_grid = VoxelGrid(K, h, 0)
sys   = build_schlogl_rdme_system(model, fine_grid)
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=300)
n_lo, n_un, n_hi = fp[1], fp[2], fp[3]
@printf("n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

# Full state space
println("Building full state space $(n_max+1)^2 = $((n_max+1)^2) states...")
sp = StateSpace{CartesianIndex{K}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp, CartesianIndex(n1,n2), 0.0)
end

# IC: front state (n_lo, n_hi)
sp.probs[sp.index[CartesianIndex(n_lo, n_hi)]] = 1.0

t_gen = @elapsed begin; result = build_generator(sp, sys, rates, 0.0); end
A = result[1]
@printf("Generator built in %.1fs  (nnz=%d)\n", t_gen, nnz(A))

# Evolve in steps, check bimodality at each
mid = (n_lo + n_hi) ÷ 2
function bimodal_weights(states, probs)
    w_lohi = sum(probs[i] for (i,s) in enumerate(states)
                 if Tuple(s)[1] < mid && Tuple(s)[2] > mid; init=0.0)
    w_hilo = sum(probs[i] for (i,s) in enumerate(states)
                 if Tuple(s)[1] > mid && Tuple(s)[2] < mid; init=0.0)
    w_lohi, w_hilo
end

p = copy(sp.probs)
t_prev = 0.0
println("\n  t        P(lo,hi)  P(hi,lo)  time/step")
for t in [1.0, 5.0, 10.0, 50.0, 200.0, 500.0, 2000.0, 10000.0]
    global p, t_prev
    dt = t - t_prev
    elapsed = @elapsed p = expv(dt, A, p; m=80)
    p ./= sum(p)
    t_prev = t
    sp.probs .= p
    wlh, whl = bimodal_weights(sp.states, sp.probs)
    @printf("  %7.1f   %.4f    %.4f    %.2fs\n", t, wlh, whl, elapsed)
end
println("Done.")
