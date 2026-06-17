# Stochastic Toggle Switch — FP-FSP simulation (3 configurations)
# Saves 2D probability grids at t≈2, t≈15, t≈30 plus size trajectories.
#
# Config 1: prob_quantile=0.3, flux_tolerance=1.0  (quantile only)
# Config 2: prob_quantile=0.3, flux_tolerance=1e-7 (quantile + flux filter)
# Config 3: prob_quantile=0.001, flux_tolerance=1.0 (tight quantile, no flux)
#
# Run: julia --project examples/toggle_switch.jl  (~10 min total)
using DiscStochSim, Catalyst, Serialization

rn = @reaction_network begin
    @species U(t) V(t)
    @parameters η α1 β1 K1 d1 dU α2 β2 K2 d2
    η*(α1 + β1*K1^3/(K1^3 + V^3)),   0 --> U
    d1 + dU,                           U --> 0
    η*(α2 + β2*K2^3/(K2^3 + U^3)),   0 --> V
    d2,                                V --> 0
end

# Parameters from Tian & Burrage (2006) / Dinh (2020), scale=100
rates = [
    1.0,          # η
    20.0,         # α1
    400.0,        # β1
    100.0,        # K1
    1.0,          # d1
    0.1/1.1,      # dU = s*γ/(1+s), s=0.1, γ=1
    20.0,         # α2
    400.0,        # β2
    100.0,        # K2
    1.0,          # d2
]

u0   = CartesianIndex(85, 5)
tf   = 30.0
NGRID = 300  # grid size for 2D probability arrays

prob = FSPProblem(rn, u0, (0.0, tf), rates; bounds=(0, NGRID - 1))

# Target snapshot times for the figure
TARGET_TIMES = [2.0, 15.0, 30.0]

function snap_at_times(sol, targets)
    # For each target time, find the nearest snapshot index
    results = []
    for t_target in targets
        idx = argmin(abs.(sol.t .- t_target))
        states, probs = sol.snapshots[idx]
        # Build 2D grid: P_binary (1 if active), P_prob (probability value)
        P_binary = zeros(Float32, NGRID, NGRID)
        P_prob   = zeros(Float32, NGRID, NGRID)
        for (s, p) in zip(states, probs)
            u, v = s[1]+1, s[2]+1   # 1-indexed
            1 <= u <= NGRID && 1 <= v <= NGRID || continue
            P_binary[u, v] = 1.0f0
            P_prob[u, v]   = Float32(p)
        end
        push!(results, (t=sol.t[idx], P_binary=P_binary, P_prob=P_prob))
    end
    results
end

configs = [
    (label="quantile_only",  prob_quantile=0.3,   flux_tolerance=1.0),
    (label="flux_filter",    prob_quantile=0.3,   flux_tolerance=1e-7),
    (label="tight_quantile", prob_quantile=0.001, flux_tolerance=1.0),
]

results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)

for (i, cfg) in enumerate(configs)
    println("\n=== Config $i: $(cfg.label) ===")
    cb = terminal_progress(["U", "V"]; tf=tf, title="Toggle $(cfg.label)")
    t_start = time()
    sol, diag = solve(prob, AdaptiveFSP(
        ε_dt           = 1.0,
        prob_quantile  = cfg.prob_quantile,
        flux_tolerance = cfg.flux_tolerance,
        save_interval  = 150,
        expand_method  = :stoich,
        progress_callback = cb,
    ))
    t_wall = time() - t_start

    traj = mean_trajectory(sol)
    println("Done: $(diag.total_iters) iters, $(round(t_wall; digits=1))s, " *
            "$(length(sol)) snapshots, final|S|=$(sol.state_space_sizes[end])")

    grids = snap_at_times(sol, TARGET_TIMES)

    out = Dict(
        "snap_t"       => sol.t,
        "snap_sizes"   => sol.state_space_sizes,
        "snap_traj"    => traj,
        "t_log"        => diag.t_log,
        "size_log"     => diag.size_log,
        "wall_time"    => t_wall,
        "total_iters"  => diag.total_iters,
        "label"        => cfg.label,
        "prob_quantile" => cfg.prob_quantile,
        "flux_tolerance"=> cfg.flux_tolerance,
        "target_times" => TARGET_TIMES,
        "grids"        => grids,   # Vector of (t, P_binary, P_prob) NamedTuples
        "NGRID"        => NGRID,
    )
    serialize(joinpath(results_dir, "toggle_switch_$(i).jls"), out)
    println("Saved → examples/results/toggle_switch_$(i).jls")
end
println("\nAll 3 toggle switch configs complete.")
