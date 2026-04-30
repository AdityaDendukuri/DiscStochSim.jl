"""
Mixed-FSP diagnostic: localise the 10% mean mismatch between FSP and SSA.

Tests (in order):
  1. Rate-table comparison at u0: print all 16 generator propensities vs SSA rates.
  2. Analytic drift d/dt E[nc₁,n₃,n₄] at u0 from both sides.
  3. One-step comparison at T=0.1 — vary expand_depth (1..16) and Krylov m.
  4. Short-time SSA reference (N=100k) for comparison.

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/diag_mixed_fsp.jl
"""

using DiscStochSim
using ExponentialUtilities
using Printf, Random, Statistics, LinearAlgebra

# ── Model ─────────────────────────────────────────────────────────────────────
D = 0.005; K = 4; h = 1.0/8
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)
rates_vec = Float64[]
d     = D / h^2

c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
@printf("Model: c1=%.4f  c2=%.2e  c3=%.2f  c4=%.2f  d=D/h²=%.4f\n",
        c1, c2, c3, c4, d)

mask = BitVector([true, false])   # pair 1 coarse, pair 2 fine
mxsys = build_schlogl_mixed_system(model, g, mask)
KM    = K - sum(mask)             # = 3

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
u0  = CartesianIndex(2*n_low, n_low, n_high)  # (62, 31, 162)
nc1, n3, n4 = Tuple(u0)
@printf("u0 = (nc₁=%d, n₃=%d, n₄=%d)\n\n", nc1, n3, n4)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Rate-table comparison
# ═══════════════════════════════════════════════════════════════════════════════
println("="^70)
println("1. RATE TABLE AT u0")
println("="^70)

# Compute all propensities from the generator system at u0
println("\n  Generator reactions:")
println("  idx  stoich      propensity(u0)  description")
println("  " * "-"^60)
for (i, (stoich, prop)) in enumerate(zip(mxsys.stoichvecs, mxsys.propensities))
    r = prop(u0, rates_vec, 0.0)
    s = Tuple(stoich)
    @printf("  %2d   (%+d,%+d,%+d)  %12.6f\n", i, s[1], s[2], s[3], r)
end

# Analytic SSA rates at u0 = (nc1=62, n3=31, n4=162)
ssa_rates = [
    2*c3,                                            # r1:  ∅→A in pair1 (coarse)
    c4*nc1,                                          # r2:  A→∅ in pair1
    c1*max(0,nc1*(nc1-1))/4,                        # r3:  2A→3A pair1 (Bin moment)
    c2*max(0,nc1*(nc1-1)*(nc1-2))/24,               # r4:  3A→2A pair1 (Bin moment)
    c3,                                              # r5:  ∅→A in voxel3 (fine)
    c4*n3,                                           # r6:  A→∅ in voxel3
    c1*max(0,n3*(n3-1))/2,                          # r7:  2A→3A voxel3
    c2*max(0,n3*(n3-1)*(n3-2))/6,                   # r8:  3A→2A voxel3
    c3,                                              # r9:  ∅→A in voxel4 (fine)
    c4*n4,                                           # r10: A→∅ in voxel4
    c1*max(0,n4*(n4-1))/2,                          # r11: 2A→3A voxel4
    c2*max(0,n4*(n4-1)*(n4-2))/6,                   # r12: 3A→2A voxel4
    d*n3,                                            # r13: voxel3→voxel4 (intra-pair2)
    d*n4,                                            # r14: voxel4→voxel3 (intra-pair2)
    d*nc1/2,                                         # r15: pair1→voxel3 (cross-pair, coarse side)
    d*n3,                                            # r16: voxel3→pair1 (cross-pair, fine side)
]

stoich_nc1 = [+1,-1,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,+1]
stoich_n3  = [ 0, 0, 0, 0,+1,-1,+1,-1, 0, 0, 0, 0,-1,+1,+1,-1]
stoich_n4  = [ 0, 0, 0, 0, 0, 0, 0, 0,+1,-1,+1,-1,+1,-1, 0, 0]

println("\n  SSA rates vs generator rates:")
println("  idx   SSA rate     gen rate     match?")
println("  " * "-"^50)

n_gen = length(mxsys.stoichvecs)
let all_match = true
    for i in 1:16
        ssa_r = ssa_rates[i]
        gen_r = i ≤ n_gen ? mxsys.propensities[i](u0, rates_vec, 0.0) : NaN
        ok = abs(ssa_r - gen_r) < 1e-10 * max(ssa_r, 1e-12)
        all_match &= ok
        @printf("  %2d    %10.5f   %10.5f   %s\n", i, ssa_r, gen_r, ok ? "✓" : "✗ ← MISMATCH")
    end
    println("\n  All rates match: $(all_match ? "YES" : "NO")")
end

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Analytic drift at u0
# ═══════════════════════════════════════════════════════════════════════════════
println()
println("="^70)
println("2. ANALYTIC DRIFT d/dt E[nc₁,n₃,n₄] AT u0")
println("="^70)

drift_nc1 = sum(stoich_nc1[i] * ssa_rates[i] for i in 1:16)
drift_n3  = sum(stoich_n3[i]  * ssa_rates[i] for i in 1:16)
drift_n4  = sum(stoich_n4[i]  * ssa_rates[i] for i in 1:16)

@printf("\n  d/dt E[nc₁] = %+.4f  (reaction: %+.4f  diffusion: %+.4f)\n",
    drift_nc1,
    sum(stoich_nc1[i]*ssa_rates[i] for i in 1:12),
    sum(stoich_nc1[i]*ssa_rates[i] for i in 13:16))
@printf("  d/dt E[n₃]  = %+.4f  (reaction: %+.4f  diffusion: %+.4f)\n",
    drift_n3,
    sum(stoich_n3[i]*ssa_rates[i] for i in 1:12),
    sum(stoich_n3[i]*ssa_rates[i] for i in 13:16))
@printf("  d/dt E[n₄]  = %+.4f  (reaction: %+.4f  diffusion: %+.4f)\n",
    drift_n4,
    sum(stoich_n4[i]*ssa_rates[i] for i in 1:12),
    sum(stoich_n4[i]*ssa_rates[i] for i in 13:16))

@printf("\n  Expected change per step (dt=0.1): Δnc₁≈%.2f  Δn₃≈%.2f  Δn₄≈%.2f\n",
    drift_nc1*0.1, drift_n3*0.1, drift_n4*0.1)

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Short-time SSA reference (T=0.1, N=100k)
# ═══════════════════════════════════════════════════════════════════════════════
println()
println("="^70)
println("3. SSA REFERENCE AT T=0.1 (N=100k independent trajectories)")
println("="^70)

function mixed_ssa_step!(s, t_final, rng)
    r = Vector{Float64}(undef, 16)
    function fill!()
        nc1, n3, n4 = s[1], s[2], s[3]
        r[1]  = 2*c3;               r[2]  = c4*nc1
        r[3]  = c1*max(0,nc1*(nc1-1))/4
        r[4]  = c2*max(0,nc1*(nc1-1)*(nc1-2))/24
        r[5]  = c3;                 r[6]  = c4*n3
        r[7]  = c1*max(0,n3*(n3-1))/2
        r[8]  = c2*max(0,n3*(n3-1)*(n3-2))/6
        r[9]  = c3;                 r[10] = c4*n4
        r[11] = c1*max(0,n4*(n4-1))/2
        r[12] = c2*max(0,n4*(n4-1)*(n4-2))/6
        r[13] = d*n3;               r[14] = d*n4
        r[15] = d*nc1/2;            r[16] = d*n3
    end
    function apply!(idx)
        if     idx==1||idx==3;  s[1]+=1
        elseif idx==2||idx==4;  s[1]-=1
        elseif idx==5||idx==7;  s[2]+=1
        elseif idx==6||idx==8;  s[2]-=1
        elseif idx==9||idx==11; s[3]+=1
        elseif idx==10||idx==12;s[3]-=1
        elseif idx==13; s[2]-=1; s[3]+=1
        elseif idx==14; s[3]-=1; s[2]+=1
        elseif idx==15; s[1]-=1; s[2]+=1
        elseif idx==16; s[2]-=1; s[1]+=1
        end
    end
    t = 0.0; fill!()
    while t < t_final
        a0 = sum(r); a0 ≤ 0.0 && break
        t += -log(rand(rng))/a0; t > t_final && break
        rv = rand(rng)*a0; cs = 0.0; idx = 16
        for i in eachindex(r); cs += r[i]; if cs ≥ rv; idx=i; break; end; end
        apply!(idx); fill!()
    end
end

u0_ssa = [nc1, n3, n4]
N_ssa  = 100_000
rng    = MersenneTwister(42)
samps  = [let s=copy(u0_ssa); mixed_ssa_step!(s,0.1,rng); copy(s) end for _ in 1:N_ssa]

ssa_mean_nc1 = mean(s[1] for s in samps)
ssa_mean_n3  = mean(s[2] for s in samps)
ssa_mean_n4  = mean(s[3] for s in samps)
@printf("\n  SSA (N=%dk):  nc₁=%.3f  n₃=%.3f  n₄=%.3f\n",
        N_ssa÷1000, ssa_mean_nc1, ssa_mean_n3, ssa_mean_n4)
@printf("  SSA deltas:   Δnc₁=%+.3f  Δn₃=%+.3f  Δn₄=%+.3f\n",
        ssa_mean_nc1-nc1, ssa_mean_n3-n3, ssa_mean_n4-n4)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. FSP at T=0.1 — sweep expand_depth and Krylov m
# ═══════════════════════════════════════════════════════════════════════════════
println()
println("="^70)
println("4. FSP AT T=0.1 — SWEEP expand_depth AND KRYLOV m")
println("="^70)

n_max_c = 2*(n_high + 60)
bc      = s -> all(c -> 0 ≤ c ≤ n_max_c, Tuple(s))
dt_test = 0.1

function run_fsp_one_step(depth, krylov_m)
    sp = StateSpace{CartesianIndex{KM}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, mxsys, bc; depth=depth)
    A, = build_generator(sp, mxsys, rates_vec, 0.0)
    p  = copy(sp.probs)
    p  = max.(expv(dt_test, A, p; m=krylov_m), 0.0)
    s  = sum(p); p ./= s   # renormalize (s ≈ 1 if no truncation)
    mass_before_clamp = sum(max.(expv(dt_test, A, sp.probs; m=krylov_m), 0.0))

    # compute means
    mnc1 = sum(Tuple(sp.states[i])[1] * p[i] for i in eachindex(sp.states))
    mn3  = sum(Tuple(sp.states[i])[2] * p[i] for i in eachindex(sp.states))
    mn4  = sum(Tuple(sp.states[i])[3] * p[i] for i in eachindex(sp.states))

    mnc1, mn3, mn4, length(sp), mass_before_clamp
end

println()
@printf("  %-8s  %-6s  %-8s  %-8s  %-8s  %-6s  %-10s\n",
        "depth", "m", "E[nc₁]", "E[n₃]", "E[n₄]", "|S|", "Σp (pre-clamp)")
println("  " * "-"^70)

for depth in [1, 2, 4, 8, 16, 32]
    for m in [30, 60, 120]
        try
            mnc1, mn3, mn4, sz, mass = run_fsp_one_step(depth, m)
            @printf("  %-8d  %-6d  %-8.3f  %-8.3f  %-8.3f  %-6d  %.6f\n",
                    depth, m, mnc1, mn3, mn4, sz, mass)
        catch e
            @printf("  %-8d  %-6d  ERROR: %s\n", depth, m, string(e)[1:40])
        end
        flush(stdout)
    end
end

println()
@printf("  SSA reference:  nc₁=%.3f  n₃=%.3f  n₄=%.3f\n",
        ssa_mean_nc1, ssa_mean_n3, ssa_mean_n4)
@printf("  Analytic (1st-order): nc₁≈%.3f  n₃≈%.3f  n₄≈%.3f\n",
        nc1+drift_nc1*dt_test, n3+drift_n3*dt_test, n4+drift_n4*dt_test)

println()
println("="^70)
println("INTERPRETATION GUIDE")
println("="^70)
println("""
  If depth-1 FSP already matches SSA: expansion is not the bottleneck.
  If means increase monotonically with depth: truncation is the bug.
  If means match at large depth but m matters: Krylov convergence is the bug.
  If rates don't match (section 1): generator vs SSA mismatch.
""")
