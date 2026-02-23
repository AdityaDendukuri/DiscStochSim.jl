# Quick test of gradual bottleneck figure
# δt=1000, depth=100, flux pruning, tf=3e6
using DiscStochSim, Catalyst, Random
using Expokit: expmv
using Plots; gr()
default(fontfamily="Computer Modern", linewidth=2, framestyle=:box, grid=false, dpi=300)

const OUTDIR = joinpath(@__DIR__, "output", "comparison")
mkpath(OUTDIR)

rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

k1, k2 = 1e-6, 1e-1
rates = [k1, k2]
u0    = CartesianIndex(1, 0, 0)
tf    = 3e6
prob  = FSPProblem(rn, u0, (0.0, tf), rates; bounds=(0, Int(4e8)))

E_C_exact = t -> (k2/k1) * (k1*t - 1 + exp(-k1*t))
println("Analytical E[C](tf=$tf) = $(round(E_C_exact(tf); sigdigits=6))")

# ── FP-FSP: δt=1000, depth=100, flux pruning ─────────────────────────
function fpfsp_fixed(prob, δt, depth, qile, ftol, save_interval)
    model = prob.model; rates_ = prob.rates; bc = prob.bc
    t0, tf_ = prob.tspan; n_species = length(Tuple(prob.u0))
    sp = StateSpace{typeof(prob.u0), Float64}()
    add_state!(sp, prob.u0, 1.0)
    sol_t = Float64[t0]
    sol_snaps = [(copy(sp.states), copy(sp.probs))]
    sol_sizes = [length(sp)]
    t = Float64(t0); iter = 0
    while t < tf_
        iter += 1
        expand!(sp, model, bc; depth=depth)
        A, _, _ = build_generator(sp, model, rates_, t)
        adt = min(δt, tf_ - t)
        sp.probs = expmv(adt, A, sp.probs)
        compress!(sp, model, rates_, t, qile; flux_tolerance=ftol)
        renormalize!(sp)
        t += adt
        if iter % save_interval == 0 || t >= tf_
            push!(sol_t, t); push!(sol_snaps, (copy(sp.states), copy(sp.probs))); push!(sol_sizes, length(sp))
        end
    end
    FSPSolution{typeof(prob.u0)}(sol_t, sol_snaps, sol_sizes, n_species)
end

println("\nRunning FP-FSP (δt=1000, depth=100, tf=$tf)...")
t0 = time()
sol_afsp = fpfsp_fixed(prob, 1000.0, 100, 0.9, 1e-6, 10)
t_afsp = time() - t0
traj_afsp = mean_trajectory(sol_afsp)
println("  Done: $(round(t_afsp; digits=1))s, max|S|=$(maximum(sol_afsp.state_space_sizes))")
println("  Final ⟨C⟩=$(round(traj_afsp[end,3]; sigdigits=5)), analytical=$(round(E_C_exact(tf); sigdigits=5))")

# ── KrylovFSP ────────────────────────────────────────────────────────
println("\nRunning KrylovFSP (ε=0.01, τ_init=1e4, max_iter=500)...")
Random.seed!(1)
t0 = time()
sol_kfsp, diag_kfsp = solve(prob, KrylovFSP(ε=0.01, τ_init=1e4, ℓ=5, r=1, save_interval=10, max_iter=500))
t_kfsp = time() - t0
traj_kfsp = mean_trajectory(sol_kfsp)
println("  Done: $(round(t_kfsp; digits=1))s, $(diag_kfsp.total_iters) iters, t_reached=$(round(sol_kfsp.t[end]; sigdigits=3))")
if size(traj_kfsp, 1) > 0
    println("  Final ⟨C⟩=$(round(traj_kfsp[end,3]; sigdigits=5))")
end

# ── Plot ─────────────────────────────────────────────────────────────
t_ref = range(0.0, tf; length=300)
p1 = plot(title="(a) Mean ⟨C⟩", xlabel="Time", ylabel="⟨C⟩")
plot!(p1, t_ref, E_C_exact.(t_ref), label="Analytical", color=:black, linestyle=:dot)
plot!(p1, sol_afsp.t, traj_afsp[:, 3], label="FP-FSP", color=:blue, linestyle=:solid)
length(sol_kfsp.t) > 1 && plot!(p1, sol_kfsp.t, traj_kfsp[:, 3], label="Krylov-FSP-SSA", color=:red, linestyle=:dash)

p2 = plot(title="(b) Mean ⟨B⟩", xlabel="Time", ylabel="⟨B⟩")
plot!(p2, sol_afsp.t, traj_afsp[:, 2], label="FP-FSP", color=:blue, linestyle=:solid)
length(sol_kfsp.t) > 1 && plot!(p2, sol_kfsp.t, traj_kfsp[:, 2], label="Krylov-FSP-SSA", color=:red, linestyle=:dash)

fig = plot(p1, p2; layout=(1,2), size=(900,350), margin=5Plots.mm)
savefig(fig, joinpath(OUTDIR, "bottleneck_comparison.png"))
savefig(fig, joinpath(OUTDIR, "bottleneck_comparison.pdf"))
println("\nSaved bottleneck_comparison.{png,pdf}")
