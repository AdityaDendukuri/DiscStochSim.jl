"""
Test: multi_level_vcycle_schlogl on wide-front Schlogl (large D)

With D=5: front width ≈ √(D/|f'(n_un)|) ≈ √(5/0.195) ≈ 5 voxels → K_act=4-6
With D=1: front width ≈ 2-3 voxels → K_act=2-3 (V-cycle not triggered)

Compare: direct expv (n_levels=0) vs V-cycle (n_levels=1)
"""

using DiscStochSim, Printf

model = SchloglModel1D(1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Fixed points: n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

# Use K=20, IC with 4 voxels in transition zone to force K_act=4 (V-cycle triggers)
K = 20; DX = 1.0; N_MAX = 200; DT = 0.5; T_END = 4.0; MAX_KA = 6

# Transition IC: middle 4 voxels at n_un (unstable fixed point) → bistable zone
ic_fn(i) = (i == 1 || i == K)   ? Float64(n_lo) :
            (i >= K÷2-1 && i <= K÷2+2) ? Float64(n_un) :
            i < K÷2 ? Float64(n_hi) : Float64(n_lo)

function run_test(n_levels)
    μ = [i <= K÷2 ? Float64(n_hi) : Float64(n_lo) for i in 1:K]
    # 4-voxel active strip straddling the interface: (hi,hi,lo,lo)
    # V-cycle triggers because Ka=4 is even and >= 4
    k_lo0 = K÷2 - 1; k_hi0 = K÷2 + 2
    avs = ActiveVoxelSet1D(k_lo0, k_hi0, K)
    sp  = StateSpace{CartesianIndex{4}, Float64}()
    add_state!(sp, CartesianIndex(n_hi, n_hi, n_lo, n_lo), 1.0)
    grid0 = VoxelGrid(4, DX, 0)
    sys0  = build_active_schlogl_1d(model, grid0, k_lo0-1, k_hi0+1, μ, Val(4))
    expand!(sp, sys0, schlogl_mixed_bc(N_MAX); depth=1)

    rates = Float64[]
    n_steps = round(Int, T_END / DT)
    max_states = Int[]; ka_vals = Int[]

    t0 = time()
    for step in 1:n_steps
        t   = Float64((step-1) * DT)
        Ka  = n_active(avs)
        lbc = left_inactive(avs); rbc = right_inactive(avs)
        ga  = VoxelGrid(Ka, DX, 0)

        sp = evolve_active_schlogl!(sp, model, ga, lbc, rbc, μ, rates, t, DT,
                                     N_MAX, 1e-3, 30, n_levels, 0.1)

        update_inactive_means_schlogl_1d!(μ, sp, avs, model, DX, DT)
        sp, avs = adapt_schlogl_1d(sp, avs, μ, N_MAX, n_lo, n_hi;
                                    promote_margin=15.0, demote_margin=8.0,
                                    min_k_act=2, max_k_act=MAX_KA)
        push!(max_states, length(sp)); push!(ka_vals, n_active(avs))
    end
    elapsed = time() - t0

    @printf("  n_levels=%d  time=%.2fs  |S|_max=%d  K_act_max=%d  front_pos=[%d,%d]\n",
            n_levels, elapsed, maximum(max_states), maximum(ka_vals),
            avs.k_lo, avs.k_hi)
    μ, avs.k_lo, avs.k_hi
end

println("\nDirect expv (n_levels=0):")
μ0, lo0, hi0 = run_test(0)

println("V-cycle     (n_levels=1):")
μ1, lo1, hi1 = run_test(1)

println("\nFront position: direct=$(lo0)-$(hi0)  vcycle=$(lo1)-$(hi1)")
println("Mean profile match: max|Δμ| = $(maximum(abs.(μ0 .- μ1)))")
