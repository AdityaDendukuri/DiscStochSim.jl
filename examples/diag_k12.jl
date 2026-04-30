using DiscStochSim
using Printf

D = 0.005; K = 12; h = 1.0/K
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_buf = 30; n̄_low = 2*n_low; n̄_high = 2*n_high

levels = [1, 1, 0, 1, 1, 1]
off = DiscStochSim._mixed_offsets(levels)
op3 = off[3]
KM  = _mixed_dim(levels)
@printf("KM=%d  op3=%d\n", KM, op3); flush(stdout)

u0a = CartesianIndex(n̄_low, n̄_low, n_low,  n_high, n̄_high, n̄_high, n̄_high)
u0b = CartesianIndex(n̄_low, n̄_low, n_high, n_low,  n̄_high, n̄_high, n̄_high)
sp = StateSpace{CartesianIndex{KM}, Float64}()
add_state!(sp, u0a, 0.5); add_state!(sp, u0b, 0.5)

@printf("building mixed system...\n"); flush(stdout)
@time mxsys = build_schlogl_mixed_system(model, g, levels)
@printf("done. expanding...\n"); flush(stdout)

n_max_coarse = 2*(n_high + n_buf); n_max_fine = n_high + n_buf
function bc_mix(s)
    t = Tuple(s)
    t[1] ≤ n_max_coarse && t[2] ≤ n_max_coarse &&
    t[op3]   ≤ n_max_fine && t[op3+1] ≤ n_max_fine &&
    t[5] ≤ n_max_coarse && t[6] ≤ n_max_coarse && t[7] ≤ n_max_coarse
end
@time expand!(sp, mxsys, bc_mix; depth=1)
@printf("|S|=%d\n", length(sp.states)); flush(stdout)

@printf("building generator...\n"); flush(stdout)
@time A, = build_generator(sp, mxsys, Float64[], 0.0)
@printf("generator: %dx%d\n", size(A)...); flush(stdout)

using SparseArrays, ExponentialUtilities
@printf("running expv...\n"); flush(stdout)
@time v = expv(0.1, A, sp.probs; m=40)
@printf("expv done. sum(v)=%.6f\n", sum(v)); flush(stdout)

# Now simulate what two_level_step_mixed does: permissive BC
sp2 = StateSpace{CartesianIndex{KM}, Float64}()
add_state!(sp2, u0a, 0.5); add_state!(sp2, u0b, 0.5)
expand!(sp2, mxsys, bc_mix; depth=1)   # initial tight expand → 54 states
@printf("\nPermissive expand (all dims ≤ 384)...\n"); flush(stdout)
bc_perm = s -> all(c -> 0 ≤ c ≤ 384, Tuple(s))
@time expand!(sp2, mxsys, bc_perm; depth=1)
@printf("|S| after permissive expand=%d\n", length(sp2.states)); flush(stdout)

@printf("building generator for 618-state sp...\n"); flush(stdout)
@time A2, = build_generator(sp2, mxsys, Float64[], 0.0)
@printf("generator: %dx%d\n", size(A2)...); flush(stdout)
@printf("running expv on 618 states...\n"); flush(stdout)
p2 = zeros(length(sp2.states)); p2[1] = 0.5; p2[2] = 0.5
@time v2 = expv(0.1, A2, p2; m=40)
@printf("expv done. sum(v2)=%.6f\n", sum(v2)); flush(stdout)
