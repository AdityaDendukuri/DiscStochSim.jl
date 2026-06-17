# Bottleneck Paper Figures: flux_vs_naive_comparison.png + error_accumulation_analysis.png
#
# System: A →(k1) B, B →(k2) B+C
# k1 = 1e-6, k2 = 0.1, initial state (1,0,0)
#
# Run: julia --project examples/bottleneck_paper_plots.jl

using DiscStochSim, Catalyst, Plots, Expokit, SparseArrays, LinearAlgebra

gr()
default(
    fontfamily    = "Computer Modern",
    titlefontsize  = 10,
    guidefontsize  = 9,
    tickfontsize   = 8,
    legendfontsize = 8,
    linewidth = 1.5,
    framestyle = :box,
    grid = false,
    dpi = 300,
)

const N_SPECIES = 3

rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

const MODEL = DiscreteStochasticSystem(rn)
const RATES = [1e-6, 1e-1]
const BC    = x -> RectLatticeBoundaryCondition(x, (0, Int(1e11)))

function compute_mean(sp::StateSpace)
    μ = zeros(N_SPECIES)
    for (i, s) in enumerate(sp.states)
        for d in 1:N_SPECIES
            μ[d] += s[d] * sp.probs[i]
        end
    end
    μ
end

# ── Simulation runners ───────────────────────────────────────────────────────
function run_fsp(; δt, tf, prob_quantile, flux_tolerance, save_every, log_every)
    sp = StateSpace{CartesianIndex{3}, Float64}()
    add_state!(sp, CartesianIndex(1, 0, 0), 1.0)

    n_steps = round(Int, tf / δt)
    times = Float64[0.0]
    means = zeros(1, N_SPECIES)
    means[1, :] = compute_mean(sp)
    sizes = Int[length(sp)]

    gids_old = Int[]
    A_old = spzeros(Float64, 0, 0)

    for iter in 1:n_steps
        t = iter * δt

        expand!(sp, MODEL, BC; depth=1)

        A, _, _, _ = isempty(gids_old) ?
            build_generator(sp, MODEL, RATES, t) :
            reconstruct_generator(sp, MODEL, RATES, t, A_old, gids_old)

        sp.probs = expmv(δt, A, sp.probs)

        compress!(sp, MODEL, RATES, t, prob_quantile; flux_tolerance=flux_tolerance)
        renormalize!(sp)

        gids_old = get_global_ids(sp)
        A_old = A

        if iter % save_every == 0 || iter == n_steps
            push!(times, t)
            means = vcat(means, compute_mean(sp)')
            push!(sizes, length(sp))
        end

        if iter % log_every == 0
            μ = compute_mean(sp)
            println("  t=$t, |S|=$(length(sp)), <A>=$(round(μ[1];sigdigits=3)), <C>=$(round(μ[3];sigdigits=4))")
        end
    end

    (times=times, means=means, sizes=sizes)
end

function run_exact_cme()
    sp = StateSpace{CartesianIndex{3}, Float64}()
    for i in 0:2, j in 0:2, k in 0:10000
        p0 = (i == 1 && j == 0 && k == 0) ? 1.0 : 0.0
        add_state!(sp, CartesianIndex(i, j, k), p0)
    end
    println("  State space: $(length(sp)) states")

    A, _, _, _ = build_generator(sp, MODEL, RATES, 0.0)
    println("  Generator built, nnz=$(nnz(A))")

    exact_times = collect(0.0:1e4:1e5)
    exact_means = zeros(length(exact_times), N_SPECIES)
    p_init = copy(sp.probs)

    for (idx, t) in enumerate(exact_times)
        p_t = t == 0.0 ? p_init : expmv(t, A, p_init)
        for (i, s) in enumerate(sp.states)
            for d in 1:N_SPECIES
                exact_means[idx, d] += s[d] * p_t[i]
            end
        end
        println("  t=$t: <A>=$(round(exact_means[idx,1];sigdigits=4)), <C>=$(round(exact_means[idx,3];sigdigits=4))")
    end

    (times=exact_times, means=exact_means)
end

# ── Run simulations ──────────────────────────────────────────────────────────
println("=" ^ 60)
println("Simulation 1: Flux-aware pruning (δt=1000, α=0.9, ε_flux=1e-12)")
println("=" ^ 60)
res_flux = run_fsp(δt=1000.0, tf=3e6, prob_quantile=0.9, flux_tolerance=1e-12,
                   save_every=100, log_every=500)

println("\n" * "=" ^ 60)
println("Simulation 2: No-flux pruning (δt=1000, α=0.9, ε_flux=0)")
println("=" ^ 60)
res_noflux = run_fsp(δt=1000.0, tf=3e6, prob_quantile=0.9, flux_tolerance=0.0,
                     save_every=100, log_every=500)

println("\n" * "=" ^ 60)
println("Simulation 3: Exact CME (tf=1e5)")
println("=" ^ 60)
res_exact = run_exact_cme()

println("\n" * "=" ^ 60)
println("Simulation 4: Fine adaptive (δt=10, ε_flux=1e-6, tf=1e5)")
println("=" ^ 60)
res_fine = run_fsp(δt=10.0, tf=1e5, prob_quantile=0.9, flux_tolerance=1e-6,
                   save_every=100, log_every=2000)

paper_dir = joinpath(@__DIR__, "..", "paper")

# ── Figure 1: flux_vs_naive_comparison.png ───────────────────────────────────
println("\nGenerating flux_vs_naive_comparison.png...")

mf  = res_flux.times   .> 0
mnf = res_noflux.times .> 0

C_flux = res_flux.means[mf, 3]
C_flux_norm = C_flux ./ max(maximum(C_flux), 1.0)
p1a = plot(res_flux.times[mf], res_flux.means[mf, 1];
           label="⟨A⟩", color=:black, linestyle=:solid,
           xlabel="Time (s)", ylabel="Mean population / max",
           title="(a) Flux-based pruning (α=0.9, ε_flux=10⁻¹²)")
plot!(p1a, res_flux.times[mf], res_flux.means[mf, 2]; label="⟨B⟩", color=:black, linestyle=:dash)
plot!(p1a, res_flux.times[mf], C_flux_norm; label="⟨C⟩/max(⟨C⟩)", color=:black, linestyle=:dot)

p1b = plot(res_flux.times[mf], res_flux.sizes[mf];
           label="", color=:black,
           xlabel="Time (s)", ylabel="|S|",
           title="(b) State space size (flux-aware)")

C_noflux = res_noflux.means[mnf, 3]
C_noflux_norm = C_noflux ./ max(maximum(C_noflux), 1.0)
p1c = plot(res_noflux.times[mnf], res_noflux.means[mnf, 1];
           label="⟨A⟩", color=:black, linestyle=:solid,
           xlabel="Time (s)", ylabel="Mean population / max",
           title="(c) Probability-only pruning (α=0.9, no flux)")
plot!(p1c, res_noflux.times[mnf], res_noflux.means[mnf, 2]; label="⟨B⟩", color=:black, linestyle=:dash)
plot!(p1c, res_noflux.times[mnf], C_noflux_norm; label="⟨C⟩/max(⟨C⟩)", color=:black, linestyle=:dot)

p1d = plot(res_noflux.times[mnf], res_noflux.sizes[mnf];
           label="", color=:black,
           xlabel="Time (s)", ylabel="|S|",
           title="(d) State space size (probability-only)")

fig1 = plot(p1a, p1b, p1c, p1d; layout=(2,2), size=(1000, 700), margin=5Plots.mm)
savefig(fig1, joinpath(paper_dir, "flux_vs_naive_comparison.png"))
println("  Saved flux_vs_naive_comparison.png")

# ── Figure 2: error_accumulation_analysis.png ────────────────────────────────
println("Generating error_accumulation_analysis.png...")

adaptive_at_exact = zeros(length(res_exact.times), N_SPECIES)
for (i, t_ex) in enumerate(res_exact.times)
    _, idx = findmin(abs.(res_fine.times .- t_ex))
    adaptive_at_exact[i, :] = res_fine.means[idx, :]
end

abs_error_C = abs.(adaptive_at_exact[:, 3] .- res_exact.means[:, 3])
rel_error_C = abs_error_C ./ max.(res_exact.means[:, 3], 1e-15)
error_rate  = diff(abs_error_C) ./ diff(res_exact.times)
t_mid       = (res_exact.times[1:end-1] .+ res_exact.times[2:end]) ./ 2

valid = res_exact.means[:, 3] .> 0.1

p2a = plot(res_exact.times[valid], rel_error_C[valid];
           marker=:circle, markersize=5, color=:black, label="",
           xlabel="Time (s)", ylabel="Relative error",
           title="(a) Relative error in ⟨C⟩")

# Drop last point: finite difference at simulation endpoint is unreliable
# (fine run's final snapshot may not align exactly with t=tf)
er_interior = error_rate[1:end-1]
tm_interior = t_mid[1:end-1]
vr = abs.(er_interior) .> 0
p2b = plot(tm_interior[vr], abs.(er_interior[vr]);
           marker=:circle, markersize=5, color=:black, label="",
           xlabel="Time (s)", ylabel="Error rate (s⁻¹)",
           title="(b) Instantaneous error rate")

p2c = plot(res_exact.times, res_exact.means[:, 3];
           color=:black, linestyle=:solid, linewidth=2, label="⟨C⟩ exact",
           xlabel="Time (s)", ylabel="Mean population",
           title="(c) Mean trajectories")
scatter!(p2c, res_exact.times, adaptive_at_exact[:, 3];
         marker=:circle, markersize=6, color=:black, label="⟨C⟩ adaptive")
plot!(p2c, res_exact.times, res_exact.means[:, 1];
      color=:black, linestyle=:dash, label="⟨A⟩ exact")

vd = res_exact.means[:, 3] .> 0.1
p2d = scatter(res_exact.means[vd, 3], abs_error_C[vd];
              color=:black, markersize=6, label="data",
              xlabel="⟨C⟩ (exact)", ylabel="Absolute error in ⟨C⟩",
              title="(d) Absolute error scaling")
if sum(vd) >= 2
    x_fit = res_exact.means[vd, 3]
    slope = sum(x_fit .* abs_error_C[vd]) / sum(x_fit .^ 2)
    x_line = range(0, maximum(x_fit); length=100)
    plot!(p2d, collect(x_line), slope .* collect(x_line);
          color=:black, linestyle=:dash,
          label="slope=$(round(slope; sigdigits=3))")
end

fig2 = plot(p2a, p2b, p2c, p2d; layout=(2,2), size=(1000, 700), margin=5Plots.mm)
savefig(fig2, joinpath(paper_dir, "error_accumulation_analysis.png"))
println("  Saved error_accumulation_analysis.png")

println("\nFlux-aware: final <A>=$(round(res_flux.means[end,1];sigdigits=4)), <C>=$(round(res_flux.means[end,3];sigdigits=4))")
println("No-flux:    final <A>=$(round(res_noflux.means[end,1];sigdigits=4)), <C>=$(round(res_noflux.means[end,3];sigdigits=4))")
println("Exact(1e5): <A>=$(round(res_exact.means[end,1];sigdigits=4)), <C>=$(round(res_exact.means[end,3];sigdigits=4))")
println("Adaptive:   <C>=$(round(adaptive_at_exact[end,3];sigdigits=4))")
println("Done!")
