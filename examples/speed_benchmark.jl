# Per-step wall-clock comparison: forward enumeration vs standard construction (speed.pdf)
# Runs both solver modes on Oregonator, Robertson, and Toggle Switch.
# Saves timing data to results/speed_*.jls, then generates paper/speed.pdf
#
# Run: julia --project examples/speed_benchmark.jl  (~hours total)
using DiscStochSim, Catalyst, Serialization

# ── Oregonator ───────────────────────────────────────────────────────────────
rn_orag = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end
let y1s=500.0, y2s=1000.0, y3s=2000.0, mu1s=2000.0, mu2s=50000.0
    global rates_orag = [mu1s/y2s, mu2s/(y1s*y2s), (mu1s+mu2s)/y1s,
                         2*mu1s/y1s^2, (mu1s+mu2s)/y3s]
end
prob_orag = FSPProblem(rn_orag, CartesianIndex(500, 1000, 2000), (0.0, 0.11), rates_orag;
                       bounds=(0, 50_000))

# ── Robertson ────────────────────────────────────────────────────────────────
rn_rober = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
rates_rober = [0.04, 3e6, 1.0]
prob_rober  = FSPProblem(rn_rober, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates_rober;
                         bounds=(0, 100_001))

# ── Toggle switch ────────────────────────────────────────────────────────────
rn_tog = @reaction_network begin
    @species U(t) V(t)
    @parameters η α1 β1 K1 d1 dU α2 β2 K2 d2
    η*(α1 + β1*K1^3/(K1^3 + V^3)),   0 --> U
    d1 + dU,                           U --> 0
    η*(α2 + β2*K2^3/(K2^3 + U^3)),   0 --> V
    d2,                                V --> 0
end
rates_tog = [1.0, 20.0, 400.0, 100.0, 1.0, 0.1/1.1, 20.0, 400.0, 100.0, 1.0]
prob_tog  = FSPProblem(rn_tog, CartesianIndex(85, 5), (0.0, 30.0), rates_tog;
                       bounds=(0, 299))

results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)

systems = [
    (name="oregonator", prob=prob_orag,
     alg_kw=(ε_dt=1.0, prob_quantile=0.1, flux_tolerance=1.0,
             expansion_depth=1, save_interval=50_000, expand_method=:ssa)),
    (name="robertson",  prob=prob_rober,
     alg_kw=(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
             expansion_depth=1, save_interval=500_000, expand_method=:stoich)),
    (name="toggle",     prob=prob_tog,
     alg_kw=(ε_dt=1.0, prob_quantile=0.3, flux_tolerance=1.0,
             expansion_depth=1, save_interval=200_000, expand_method=:stoich)),
]

for sys in systems
    for (mode, label) in [(true, "forward"), (false, "standard")]
        fname = "speed_$(sys.name)_$(label).jls"
        fpath = joinpath(results_dir, fname)
        if isfile(fpath)
            println("Skipping $(fname) (already exists)")
            continue
        end

        println("\n=== $(sys.name) / $(label) ===")
        alg = AdaptiveFSP(; sys.alg_kw..., use_incremental=mode, record_step_time=true)
        t_wall_start = time()
        sol, diag = solve(sys.prob, alg)
        t_wall = time() - t_wall_start
        println("Done: $(diag.total_iters) iters, $(round(t_wall; digits=1))s")

        # Subsample to ~2000 points for plotting
        n_save = 2000
        t_full = diag.t_log
        t_grid = range(t_full[1], t_full[end]; length=n_save)
        raw_idx = searchsortedfirst.(Ref(t_full), t_grid)
        idx = unique(clamp.(raw_idx, 1, length(t_full)))

        out = Dict(
            "t_log"          => diag.t_log[idx],
            "step_time_log"  => diag.step_time_log[idx],
            "size_log"       => diag.size_log[idx],
            "wall_time"      => t_wall,
            "total_iters"    => diag.total_iters,
            "label"          => label,
            "system"         => sys.name,
        )
        serialize(fpath, out)
        println("Saved → examples/results/$(fname)")
    end
end

println("\nAll speed benchmark runs complete. Run speed_plot.jl to generate paper/speed.pdf")
