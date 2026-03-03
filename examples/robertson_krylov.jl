# Robertson Autocatalytic System — KrylovFSP baseline
# Runs KrylovFSP with the correct absorbing-form generator (Algorithm 2,
# Vo & Sidje 2016) on the stiff Robertson system over [0, 1e4].
# Saves diagnostics to examples/results/robertson_krylov.jls for plotting.
using DiscStochSim, Catalyst, Serialization

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end

rates = [0.04, 3e6, 1.0]
prob = FSPProblem(rn, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates;
                  bounds=(0, 100_001))

println("Running KrylovFSP (absorbing generator, max_iter=500)...")
t_wall_start = time()
sol, diag = solve(prob, KrylovFSP(
    ε=1e-2, τ_init=0.1, r=1, save_interval=10, max_iter=500,
))
t_wall = time() - t_wall_start

t_reached = sol.t[end]
println("Done: $(diag.total_iters) iters, $(diag.reject_count) rejections, " *
        "$(round(t_wall; digits=2))s")
println("t_reached = $(round(t_reached; sigdigits=4)) / 1e4 " *
        "($(round(100*t_reached/1e4; digits=3))%)")
if !isempty(diag.τ_log)
    valid = filter(x -> x > 0, diag.τ_log)
    println("τ range: $(minimum(valid)) → $(maximum(valid))")
end

# Save diagnostics for plotting
results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
out = Dict(
    "τ_log"       => diag.τ_log,
    "time_log"    => diag.time_log,
    "total_iters" => diag.total_iters,
    "reject_count"=> diag.reject_count,
    "wall_time"   => t_wall,
    "t_reached"   => t_reached,
    "S_final"     => isempty(diag.size_log) ? 0 : diag.size_log[end],
)
serialize(joinpath(results_dir, "robertson_krylov.jls"), out)
println("Saved diagnostics → examples/results/robertson_krylov.jls")
