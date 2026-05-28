using DiscStochSim, Printf

# Scan c3 values to find the Maxwell point (v_flat = 0) for the Schlogl model.
# We use a 1D strip (1×K grid) with a flat front IC: left half = equil-hi, right half = equil-lo.
# Measure how far the front has moved after T time units; positive = hi expands, negative = lo expands.

K     = 80
T     = 10.0
dt    = 0.5
D     = 2.0

# fixed Schlogl params, scan c3
k_b   = 0.028
k_d   = 3.2e-4
c4    = 1.0

function flat_front_displacement(c3)
    model = SchloglModel1D(D, k_b, k_d, c3, c4)
    n_lo, n_un, n_hi = schlogl_fixed_points(model)
    n_lo > 0 && n_un > n_lo && n_hi > n_un || return NaN, n_lo, n_un, n_hi

    s = SpatialFSP(model, K, 1;
                   n_max=n_hi+30, n_max_eq=n_hi+30,
                   ε_expand=0.15, ε_equil=0.12, use_block_vcycle=false)

    # flat IC: left half = equil-hi, right half = equil-lo
    for i in 1:K
        p = zeros(n_hi+31)
        p[(i <= K÷2 ? n_hi : n_lo) + 1] = 1.0
        s.equil_dists[CartesianIndex(i,1)] = p
    end

    voxel_mean_f(p) = sum((n-1)*p[n] for n in 1:length(p))
    front_pos(s) = begin
        means = [voxel_mean_f(get(s.dists, CartesianIndex(i,1),
                  get(s.equil_dists, CartesianIndex(i,1), [1.0]))) for i in 1:K]
        # front = first voxel from right where mean > n_un
        idxs = findall(m -> m > n_un, means)
        isempty(idxs) ? K÷2 : last(idxs)
    end

    pos0 = K÷2
    for _ in 1:round(Int, T/dt); step!(s, dt); end
    pos1 = front_pos(s)
    displacement = pos1 - pos0
    return Float64(displacement), n_lo, n_un, n_hi
end

@printf("Scanning c3 for Maxwell point  (K=%d, T=%.0f, D=%.1f)\n\n", K, T, D)
@printf("  %6s  %5s  %5s  %5s  %8s\n", "c3", "n_lo", "n_un", "n_hi", "Δfront")

for c3 in [14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0]
    disp, n_lo, n_un, n_hi = flat_front_displacement(c3)
    direction = disp > 1 ? "hi expands" : disp < -1 ? "lo expands" : "near zero ✓"
    @printf("  %6.1f  %5d  %5d  %5d  %8.1f   %s\n", c3, n_lo, n_un, n_hi, disp, direction)
    flush(stdout)
end
