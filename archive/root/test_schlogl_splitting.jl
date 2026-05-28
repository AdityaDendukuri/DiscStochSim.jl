using DiscStochSim
using LinearAlgebra

# ── Schlögl unstable-state splitting demo ────────────────────────────────────
# D=5: intra-pair equilibration τ_eq ≈ 1/(2D) = 0.1  <<  reaction τ_rx ≈ 5
# → V-cycle approximation (fast intra-pair mixing) is well-justified here.
# IC: all voxels at the unstable Schlögl fixed point.
# Expected: probability mass splits toward low (~31) and high (~162) stable states.

K     = 4
D     = 5.0
model = SchloglModel1D(D)
fine_grid = VoxelGrid(K, 1.0, 0)
rates = Float64[]

n_low, n_unstable, n_high = schlogl_fixed_points(model)
println("Schlögl fixed points: low=$n_low  unstable=$n_unstable  high=$n_high")
println("D=$D  →  intra-pair eq time τ_eq≈$(round(1/(2D),digits=2))  reaction timescale τ_rx≈5")
println()

u0    = CartesianIndex(ntuple(_ -> n_unstable, Val(K)))
rates = Float64[]
tspan = (0.0, 12.0)

alg = RDMEMultigridFSP(
    dt                = 0.1,
    n_levels          = 1,
    n_max             = 300,
    save_every        = 10,       # snapshot every 1 time unit
    max_states        = 500_000,
    expand_depth      = 0,
    prune_tol         = 1e-10,
    coarse_expand_depth = 2,
    τ_pre             = 0.5,      # ≫ τ_eq: fully equilibrates within pairs
)

println("IC: all voxels at n=$n_unstable (unstable fixed point)")
println("Running t=0 to $(tspan[2])...")
sol = solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)
println("Done.  $(length(sol.t)) snapshots.\n")

# ── basin boundaries ──────────────────────────────────────────────────────────
b_lo = (n_low + n_unstable) ÷ 2      # midpoint low–unstable
b_hi = (n_unstable + n_high) ÷ 2     # midpoint unstable–high

function regime_prob(sol, ti, k)
    states, probs = sol.snapshots[ti]
    pl = pm = ph = 0.0
    for (s, p) in zip(states, probs)
        n = Tuple(s)[k]
        if n < b_lo;      pl += p
        elseif n > b_hi;  ph += p
        else;             pm += p
        end
    end
    pl, pm, ph
end

# ── summary table ─────────────────────────────────────────────────────────────
mu = mean_voxel_counts(sol)

println("Basin boundaries: low < $(b_lo)  |  unstable $(b_lo)-$(b_hi)  |  high > $(b_hi)")
println()
println("  t       |S|     μ[v1..v4]                  → P(low) / P(unstable) / P(high) per voxel")
println("-"^90)

for (ti, t) in enumerate(sol.t)
    ns = sol.state_space_sizes[ti]
    μ  = round.(mu[ti,:], digits=1)
    rps = [regime_prob(sol, ti, k) for k in 1:K]
    regime_str = join(["$(round(100r[1],digits=0))%/$(round(100r[3],digits=0))%" for r in rps], "  ")
    println("  t=$(lpad(round(t,digits=1),4))  |S|=$(lpad(ns,7))  μ=$μ    $regime_str")
end

println()
println("(columns show L%/H% for each voxel — watch high% grow as splitting happens)")
println()

# ── marginal at final time ────────────────────────────────────────────────────
println("Marginal distributions at t=$(round(sol.t[end],digits=1)):")
states_f, probs_f = sol.snapshots[end]
for k in 1:K
    dist = Dict{Int,Float64}()
    for (s, p) in zip(states_f, probs_f)
        n = Tuple(s)[k]; dist[n] = get(dist, n, 0.0) + p
    end
    ns_sorted = sort!(collect(keys(dist)))
    p_lo = sum((get(dist,n,0.0) for n in ns_sorted if n < b_lo);  init=0.0)
    p_mi = sum((get(dist,n,0.0) for n in ns_sorted if b_lo<=n<=b_hi); init=0.0)
    p_hi = sum((get(dist,n,0.0) for n in ns_sorted if n > b_hi);  init=0.0)
    println("  voxel $k: P(low<$b_lo)=$(round(p_lo,digits=3))  P(unstable $b_lo-$b_hi)=$(round(p_mi,digits=3))  P(high>$b_hi)=$(round(p_hi,digits=3))")
    lo_ns = [n for n in ns_sorted if n < b_lo]
    hi_ns = [n for n in ns_sorted if n > b_hi]
    if !isempty(lo_ns)
        println("    low-basin peak:  n=$(lo_ns[argmax([dist[n] for n in lo_ns])])  (stable low: $n_low)")
    end
    if !isempty(hi_ns)
        println("    high-basin peak: n=$(hi_ns[argmax([dist[n] for n in hi_ns])])  (stable high: $n_high)")
    end
end
