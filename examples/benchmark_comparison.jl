# Generator construction benchmark (benchmark_comparison.pdf)
# Compares build_generator (standard) vs reconstruct_generator (forward enumeration)
# on the Robertson system at a representative time step (|S|≈1100→1430).
#
# Run: julia --project examples/benchmark_comparison.jl  (~5 min for BenchmarkTools)
using DiscStochSim, Catalyst, Expokit,
      BenchmarkTools, Plots, Serialization, LinearAlgebra, SparseArrays

# ── Robertson reaction network ──────────────────────────────────────────────
rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
rates = [0.04, 3e6, 1.0]
prob  = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates;
                   bounds=(0, 100_001))

# ── Run solver until |S| first exceeds 1100 (≈ t=170), capture state ────────
println("Running Robertson to capture representative state space (|S|≈1100)...")

model = prob.model
bc    = prob.bc
t0, tf = prob.tspan

sp = StateSpace{CartesianIndex{3}, Float64}()
add_state!(sp, prob.u0, 1.0)

A_old    = spzeros(Float64, 0, 0)
gids_old = Int[]
t_cur    = Float64(t0)
iter     = 0

captured_sp    = nothing
captured_A_old = nothing
captured_gids_old = nothing
captured_t     = 0.0

alg = AdaptiveFSP(ε_dt=0.01, prob_quantile=0.4, flux_tolerance=1e-9,
                  expand_method=:stoich)

while t_cur < tf && iter < 200_000
    iter += 1

    expand!(sp, model, bc; depth=alg.expansion_depth)

    # Capture just before we cross 1100 states
    if length(sp) >= 1100 && captured_sp === nothing
        captured_sp       = deepcopy(sp)
        captured_A_old    = copy(A_old)
        captured_gids_old = copy(gids_old)
        captured_t        = t_cur
        println("Captured at iter=$iter, t=$(round(t_cur; sigdigits=4)), |S|=$(length(sp))")
        break
    end

    A, _, _, _ = reconstruct_generator(sp, model, rates, t_cur, A_old, gids_old)
    A_old    = A
    gids_old = get_global_ids(sp)

    Φ = dot(-LinearAlgebra.diag(A), sp.probs)
    dt = Φ > 0 ? alg.ε_dt / Φ : (tf - t_cur)
    dt = min(dt, tf - t_cur)

    sp.probs = expv(dt, A, sp.probs)
    compress!(sp, model, rates, t_cur, alg.prob_quantile;
              flux_tolerance=alg.flux_tolerance)
    renormalize!(sp)
    t_cur += dt
end

if captured_sp === nothing
    error("Did not reach |S|=1100 within 200,000 iterations")
end

sp_bench    = captured_sp
A_old_bench = captured_A_old
gids_bench  = captured_gids_old
t_bench     = captured_t

println("Benchmarking at t≈$(round(t_bench; sigdigits=4)), |S|=$(length(sp_bench))")

# ── Micro-benchmarks ────────────────────────────────────────────────────────
println("Running BenchmarkTools (10^4 evals each) — may take a few minutes...")

b_standard = @benchmark build_generator(
    $sp_bench, $model, $rates, $t_bench) samples=10000 evals=1

b_forward  = @benchmark reconstruct_generator(
    $sp_bench, $model, $rates, $t_bench, $A_old_bench, $gids_bench) samples=10000 evals=1

t_std  = b_standard.times .* 1e-6   # nanoseconds → milliseconds
t_fwd  = b_forward.times  .* 1e-6

med_std = median(t_std)
med_fwd = median(t_fwd)
speedup = med_std / med_fwd

println("Standard:   median=$(round(med_std; digits=2)) ms, " *
        "mean=$(round(mean(t_std); digits=2)) ms, " *
        "std=$(round(std(t_std); digits=2)) ms")
println("Forward:    median=$(round(med_fwd; digits=2)) ms, " *
        "mean=$(round(mean(t_fwd); digits=2)) ms, " *
        "std=$(round(std(t_fwd); digits=2)) ms")
println("Speedup: $(round(speedup; digits=1))×")

# ── PDF and CDF figure ───────────────────────────────────────────────────────
gr()
default(
    fontfamily    = "Computer Modern",
    titlefontsize  = 11,
    guidefontsize  = 10,
    tickfontsize   = 9,
    legendfontsize = 9,
    linewidth = 1.5,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

# Log-spaced histogram bins spanning both distributions
all_t = vcat(t_std, t_fwd)
tmin, tmax = minimum(all_t), maximum(all_t)
edges = 10 .^ range(log10(max(tmin, 1e-3)), log10(tmax * 1.1); length=60)

function log_hist_density(vals, edges)
    counts = zeros(length(edges) - 1)
    for v in vals
        idx = searchsortedfirst(edges, v) - 1
        1 <= idx <= length(counts) && (counts[idx] += 1)
    end
    widths  = diff(edges)
    density = counts ./ (sum(counts) .* widths)
    mids    = sqrt.(edges[1:end-1] .* edges[2:end])
    mids, density
end

mids_std, dens_std = log_hist_density(t_std, edges)
mids_fwd, dens_fwd = log_hist_density(t_fwd, edges)

p1 = plot(xlabel="Time (ms)", ylabel="Density", title="(a) PDF",
          xscale=:log10, yscale=:log10)
plot!(p1, mids_std, dens_std; label="Standard", color=:black, linestyle=:solid)
plot!(p1, mids_fwd, dens_fwd; label="Forward enum.", color=:black, linestyle=:dash)
vline!(p1, [med_std]; label="median std $(round(med_std; digits=2)) ms",
       color=:gray, linestyle=:dot)
vline!(p1, [med_fwd]; label="median fwd $(round(med_fwd; digits=2)) ms",
       color=:darkgray, linestyle=:dot)

# CDF
sort_std = sort(t_std)
sort_fwd = sort(t_fwd)
cdf_std  = (1:length(sort_std)) ./ length(sort_std)
cdf_fwd  = (1:length(sort_fwd)) ./ length(sort_fwd)

p2 = plot(xlabel="Time (ms)", ylabel="Cumulative fraction", title="(b) CDF",
          xscale=:log10)
plot!(p2, sort_std, cdf_std; label="Standard", color=:black, linestyle=:solid)
plot!(p2, sort_fwd, cdf_fwd; label="Forward enum.", color=:black, linestyle=:dash)

fig = plot(p1, p2; layout=(1, 2), size=(900, 380), margin=5Plots.mm)

paper_dir = joinpath(@__DIR__, "..", "paper")
savefig(fig, joinpath(paper_dir, "benchmark_comparison.pdf"))
savefig(fig, joinpath(paper_dir, "benchmark_comparison.png"))
println("Saved benchmark_comparison.{pdf,png} → paper/")

# Save raw timing data for reference
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
serialize(joinpath(results_dir, "benchmark_comparison.jls"), Dict(
    "t_standard_ms" => t_std,
    "t_forward_ms"  => t_fwd,
    "median_std"    => med_std,
    "median_fwd"    => med_fwd,
    "speedup"       => speedup,
    "state_size"    => length(sp_bench),
    "t_captured"    => t_bench,
))
println("Saved raw timing → examples/results/benchmark_comparison.jls")
