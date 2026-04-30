"""
Level-2 error isolation: K=4, 2 pairs, tractable fine reference.

Three solvers at the same physical model:
  fine    : P(n₁,n₂,n₃,n₄)           — ground truth
  level-1 : P(n̄₁, n̄₂)  levels=[1,1]  — pair-coarse
  level-2 : P(n̄₁₂)      levels=[2,2]  — quad-coarse

IC: both pairs in the high state (bulk, far from interface) — this is
    exactly when level-2 is triggered in the K=12 adaptive run.

We also compare with a "front" IC: pair 1 low, pair 2 high — to show
that level-2 is only safe in the bulk.

Metrics (all at t=T_END):
  TV(P_fine→(n̄₁,n̄₂),  P_level1)   : level-1 coarsening error
  TV(P_fine→n̄₁₂,        P_level2)   : level-2 coarsening error
  TV(P_level1→n̄₁₂,      P_level2)   : level-1 vs level-2 directly

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/isolate_level2_error.jl
"""

using DiscStochSim
using ExponentialUtilities
using Printf

D = 0.005
K = 4
h = 1.0 / 12   # match K=12 voxel spacing so d = D/h² matches K=12 scenario
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)
rates = Float64[]

d_fine = diffusion_rate(model.D, g)
fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
@printf("n_low=%d  n_uns=%d  n_high=%d  d_fine=%.4f\n", n_low, n_uns, n_high, d_fine)

n_buf  = 40
n_max_fine   = n_high + n_buf   # per fine voxel
n_max_coarse = 2*(n_high + n_buf)  # per pair voxel (level-1)
n_max_quad   = 4*(n_high + n_buf)  # quad voxel (level-2)

T_END  = 0.5
dt     = 0.1
n_steps = round(Int, T_END / dt)

# ─────────────────────────────────────────────────────────────────────────────
# Build systems
# ─────────────────────────────────────────────────────────────────────────────
fine_sys = build_schlogl_rdme_system(model, g)             # P(n₁,n₂,n₃,n₄)
sys_l1   = build_schlogl_mixed_system(model, g, [1,1])     # P(n̄₁, n̄₂)
sys_l2   = build_schlogl_mixed_system(model, g, [2,2])     # P(n̄₁₂)

bc_fine   = s -> all(c -> 0 <= c <= n_max_fine,   Tuple(s))
bc_l1     = s -> all(c -> 0 <= c <= n_max_coarse, Tuple(s))
bc_l2     = s -> all(c -> 0 <= c <= n_max_quad,   Tuple(s))

# ─────────────────────────────────────────────────────────────────────────────
# Generic adaptive runner
# ─────────────────────────────────────────────────────────────────────────────
function run_direct(sys, u0, bc, n_steps, dt; label="", krylov_m=40, prune_tol=1e-12)
    sp = StateSpace{typeof(u0), Float64}()
    add_state!(sp, u0, 1.0)
    for step in 1:n_steps
        expand!(sp, sys, bc; depth=2)
        A, = build_generator(sp, sys, rates, 0.0)
        sp.probs .= max.(expv(dt, A, sp.probs; m=krylov_m), 0.0)
        sp.probs ./= sum(sp.probs)
        prune_threshold!(sp, prune_tol)
    end
    isempty(label) || @printf("  %-12s  |S|=%d\n", label, length(sp))
    sp
end

# ─────────────────────────────────────────────────────────────────────────────
# Marginal projections
# ─────────────────────────────────────────────────────────────────────────────

# Fine → (n̄₁, n̄₂) : sum n₁+n₂ and n₃+n₄
function project_fine_to_l1(sp::StateSpace)
    mg = Dict{Tuple{Int,Int}, Float64}()
    for (i,s) in enumerate(sp.states)
        n1,n2,n3,n4 = Tuple(s)
        k = (n1+n2, n3+n4)
        mg[k] = get(mg, k, 0.0) + sp.probs[i]
    end
    mg
end

# Fine → n̄₁₂ : sum all four
function project_fine_to_l2(sp::StateSpace)
    mg = Dict{Int, Float64}()
    for (i,s) in enumerate(sp.states)
        k = sum(Tuple(s))
        mg[k] = get(mg, k, 0.0) + sp.probs[i]
    end
    mg
end

# Level-1 → n̄₁₂ : sum both pair coords
function project_l1_to_l2(sp::StateSpace)
    mg = Dict{Int, Float64}()
    for (i,s) in enumerate(sp.states)
        nc1, nc2 = Tuple(s)
        mg[nc1+nc2] = get(mg, nc1+nc2, 0.0) + sp.probs[i]
    end
    mg
end

# Level-2 → Dict{Int,Float64}
function l2_dist(sp::StateSpace)
    mg = Dict{Int, Float64}()
    for (i,s) in enumerate(sp.states)
        k = Tuple(s)[1]
        mg[k] = get(mg, k, 0.0) + sp.probs[i]
    end
    mg
end

# Level-1 → Dict{Tuple{Int,Int},Float64}
function l1_dist(sp::StateSpace)
    mg = Dict{Tuple{Int,Int}, Float64}()
    for (i,s) in enumerate(sp.states)
        mg[Tuple(s)] = get(mg, Tuple(s), 0.0) + sp.probs[i]
    end
    mg
end

tv(a::Dict{K,V}, b::Dict{K,V}) where {K,V} =
    0.5 * sum(abs(get(a,k,0.0) - get(b,k,0.0)) for k in union(keys(a),keys(b)))

# ─────────────────────────────────────────────────────────────────────────────
# Scenario A: both pairs in the HIGH state (bulk scenario)
# ─────────────────────────────────────────────────────────────────────────────
println("\n══ Scenario A: both pairs HIGH (bulk; level-2 should be safe) ══")
u0_fine_A = CartesianIndex(n_high, n_high, n_high, n_high)
u0_l1_A   = CartesianIndex(2*n_high, 2*n_high)
u0_l2_A   = CartesianIndex(4*n_high)

@printf("Running fine (K=4)...\n"); flush(stdout)
@time sp_fine_A = run_direct(fine_sys, u0_fine_A, bc_fine, n_steps, dt; label="fine")

@printf("Running level-1...\n"); flush(stdout)
@time sp_l1_A = run_direct(sys_l1, u0_l1_A, bc_l1, n_steps, dt; label="level-1")

@printf("Running level-2...\n"); flush(stdout)
@time sp_l2_A = run_direct(sys_l2, u0_l2_A, bc_l2, n_steps, dt; label="level-2")

fine_l1_A = project_fine_to_l1(sp_fine_A)
fine_l2_A = project_fine_to_l2(sp_fine_A)
l1_l2_A   = project_l1_to_l2(sp_l1_A)
l2d_A     = l2_dist(sp_l2_A)
l1d_A     = l1_dist(sp_l1_A)

@printf("\n[Scenario A — both HIGH]\n")
@printf("  TV(fine→(n̄₁,n̄₂),  level-1)   = %.6f   (level-1 error vs fine)\n",
        tv(fine_l1_A, l1d_A))
@printf("  TV(fine→n̄₁₂,        level-2)   = %.6f   (level-2 error vs fine)\n",
        tv(fine_l2_A, l2d_A))
@printf("  TV(level-1→n̄₁₂,      level-2)   = %.6f   (level-1 vs level-2)\n",
        tv(l1_l2_A, l2d_A))

# ─────────────────────────────────────────────────────────────────────────────
# Scenario B: both pairs in the LOW state (bulk scenario)
# ─────────────────────────────────────────────────────────────────────────────
println("\n══ Scenario B: both pairs LOW (bulk) ══")
u0_fine_B = CartesianIndex(n_low, n_low, n_low, n_low)
u0_l1_B   = CartesianIndex(2*n_low, 2*n_low)
u0_l2_B   = CartesianIndex(4*n_low)

@printf("Running fine...\n"); flush(stdout)
@time sp_fine_B = run_direct(fine_sys, u0_fine_B, bc_fine, n_steps, dt; label="fine")

@printf("Running level-1...\n"); flush(stdout)
@time sp_l1_B = run_direct(sys_l1, u0_l1_B, bc_l1, n_steps, dt; label="level-1")

@printf("Running level-2...\n"); flush(stdout)
@time sp_l2_B = run_direct(sys_l2, u0_l2_B, bc_l2, n_steps, dt; label="level-2")

fine_l1_B = project_fine_to_l1(sp_fine_B)
fine_l2_B = project_fine_to_l2(sp_fine_B)
l1_l2_B   = project_l1_to_l2(sp_l1_B)
l2d_B     = l2_dist(sp_l2_B)
l1d_B     = l1_dist(sp_l1_B)

@printf("\n[Scenario B — both LOW]\n")
@printf("  TV(fine→(n̄₁,n̄₂),  level-1)   = %.6f   (level-1 error vs fine)\n",
        tv(fine_l1_B, l1d_B))
@printf("  TV(fine→n̄₁₂,        level-2)   = %.6f   (level-2 error vs fine)\n",
        tv(fine_l2_B, l2d_B))
@printf("  TV(level-1→n̄₁₂,      level-2)   = %.6f   (level-1 vs level-2)\n",
        tv(l1_l2_B, l2d_B))

# ─────────────────────────────────────────────────────────────────────────────
# Scenario C: pair 1 low, pair 2 high (front scenario — level-2 should be bad)
# ─────────────────────────────────────────────────────────────────────────────
println("\n══ Scenario C: pair 1 LOW, pair 2 HIGH (front; level-2 should be risky) ══")
u0_fine_C = CartesianIndex(n_low, n_low, n_high, n_high)
u0_l1_C   = CartesianIndex(2*n_low, 2*n_high)
u0_l2_C   = CartesianIndex(2*n_low + 2*n_high)

@printf("Running fine...\n"); flush(stdout)
@time sp_fine_C = run_direct(fine_sys, u0_fine_C, bc_fine, n_steps, dt; label="fine")

@printf("Running level-1...\n"); flush(stdout)
@time sp_l1_C = run_direct(sys_l1, u0_l1_C, bc_l1, n_steps, dt; label="level-1")

@printf("Running level-2...\n"); flush(stdout)
@time sp_l2_C = run_direct(sys_l2, u0_l2_C, bc_l2, n_steps, dt; label="level-2")

fine_l1_C = project_fine_to_l1(sp_fine_C)
fine_l2_C = project_fine_to_l2(sp_fine_C)
l1_l2_C   = project_l1_to_l2(sp_l1_C)
l2d_C     = l2_dist(sp_l2_C)
l1d_C     = l1_dist(sp_l1_C)

@printf("\n[Scenario C — pair 1 LOW, pair 2 HIGH]\n")
@printf("  TV(fine→(n̄₁,n̄₂),  level-1)   = %.6f   (level-1 error vs fine)\n",
        tv(fine_l1_C, l1d_C))
@printf("  TV(fine→n̄₁₂,        level-2)   = %.6f   (level-2 error vs fine)\n",
        tv(fine_l2_C, l2d_C))
@printf("  TV(level-1→n̄₁₂,      level-2)   = %.6f   (level-1 vs level-2)\n",
        tv(l1_l2_C, l2d_C))

println("\nDone.")
