# Bottleneck Paper Figures: flux_vs_naive_comparison.png + error_accumulation_analysis.png
#
# System: A →(k1) B, B →(k2) B+C
# k1 = 1e-6, k2 = 0.1, initial state (1,0,0)
#
# Uses explicit FSP loops (not solve()) because the bottleneck system needs
# small fixed time steps to track the wave front with incremental expansion.

using DiscStochSim, Catalyst, CairoMakie, Expokit, SparseArrays, LinearAlgebra

const N_SPECIES = 3

# ============================================================================
# Model definition
# ============================================================================
rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

const MODEL = DiscreteStochasticSystem(rn)
const RATES = [1e-6, 1e-1]
const BC = x -> RectLatticeBoundaryCondition(x, (0, Int(1e11)))

function compute_mean(sp::StateSpace)
    μ = zeros(N_SPECIES)
    for (i, s) in enumerate(sp.states)
        for d in 1:N_SPECIES
            μ[d] += s[d] * sp.probs[i]
        end
    end
    μ
end

# ============================================================================
# Simulation runners
# ============================================================================
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

        if isempty(gids_old)
            A, _, _, _ = build_generator(sp, MODEL, RATES, t)
        else
            A, _, _, _ = reconstruct_generator(sp, MODEL, RATES, t, A_old, gids_old)
        end

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

# ============================================================================
# Run simulations
# ============================================================================
println("=" ^ 60)
println("Simulation 1: Flux-aware pruning (δt=10, α=0.9, ε_flux=1e-12)")
println("=" ^ 60)
res_flux = run_fsp(δt=10.0, tf=3e6, prob_quantile=0.9, flux_tolerance=1e-12,
                   save_every=1000, log_every=50000)

println("\n" * "=" ^ 60)
println("Simulation 2: No-flux pruning (δt=100, α=0.9, ε_flux=0)")
println("=" ^ 60)
res_noflux = run_fsp(δt=100.0, tf=3e6, prob_quantile=0.9, flux_tolerance=0.0,
                     save_every=1000, log_every=5000)

println("\n" * "=" ^ 60)
println("Simulation 3: Exact CME")
println("=" ^ 60)
res_exact = run_exact_cme()

println("\n" * "=" ^ 60)
println("Simulation 4: Fine adaptive (δt=10, ε_flux=1e-6, for error analysis)")
println("=" ^ 60)
res_fine = run_fsp(δt=10.0, tf=1e5, prob_quantile=0.9, flux_tolerance=1e-6,
                   save_every=100, log_every=2000)

# ============================================================================
# Figure 1: flux_vs_naive_comparison.png
# ============================================================================
println("\nGenerating flux_vs_naive_comparison.png...")

fig1 = Figure(size=(1200, 900), fontsize=14)

# Panel (a): Flux-based means
ax1a = Axis(fig1[1,1], xlabel="Time (s)", ylabel="Mean population",
            title="(a) Flux-based pruning (α=0.9, ε_flux=10⁻¹²)")
mf = res_flux.times .> 0
lines!(ax1a, res_flux.times[mf], res_flux.means[mf, 1], label="⟨A⟩", linewidth=2, color=:blue)
lines!(ax1a, res_flux.times[mf], res_flux.means[mf, 2], label="⟨B⟩", linewidth=2, color=:orange, linestyle=:dash)
lines!(ax1a, res_flux.times[mf], res_flux.means[mf, 3], label="⟨C⟩", linewidth=2, color=:green, linestyle=:dot)
axislegend(ax1a, position=:lt)

# Panel (b): State space size
ax1b = Axis(fig1[1,2], xlabel="Time (s)", ylabel="State space |S|",
            title="(b) State space growth (flux-aware)")
lines!(ax1b, res_flux.times[mf], res_flux.sizes[mf], color=:blue, linewidth=2)

# Panel (c): No-flux means
ax1c = Axis(fig1[2,1], xlabel="Time (s)", ylabel="Mean population",
            title="(c) Probability-only pruning (α=0.9, no flux)")
mnf = res_noflux.times .> 0
lines!(ax1c, res_noflux.times[mnf], res_noflux.means[mnf, 1], label="⟨A⟩", linewidth=2, color=:blue)
lines!(ax1c, res_noflux.times[mnf], res_noflux.means[mnf, 2], label="⟨B⟩", linewidth=2, color=:orange, linestyle=:dash)
lines!(ax1c, res_noflux.times[mnf], res_noflux.means[mnf, 3], label="⟨C⟩", linewidth=2, color=:green, linestyle=:dot)
axislegend(ax1c, position=:rt)

# Panel (d): No-flux state space size
ax1d = Axis(fig1[2,2], xlabel="Time (s)", ylabel="State space |S|",
            title="(d) State space (probability-only)")
lines!(ax1d, res_noflux.times[mnf], res_noflux.sizes[mnf], color=:blue, linewidth=2)

save("/home/adityad/DiscStochSim.jl/paper/flux_vs_naive_comparison.png", fig1, px_per_unit=3)
println("  Saved flux_vs_naive_comparison.png")

# ============================================================================
# Figure 2: error_accumulation_analysis.png
# ============================================================================
println("Generating error_accumulation_analysis.png...")

# Match fine adaptive to exact timepoints
adaptive_at_exact = zeros(length(res_exact.times), N_SPECIES)
for (i, t_ex) in enumerate(res_exact.times)
    _, idx = findmin(abs.(res_fine.times .- t_ex))
    adaptive_at_exact[i, :] = res_fine.means[idx, :]
end

abs_error_C = abs.(adaptive_at_exact[:, 3] .- res_exact.means[:, 3])
rel_error_C = abs_error_C ./ max.(res_exact.means[:, 3], 1e-15)
error_rate = diff(abs_error_C) ./ diff(res_exact.times)

fig2 = Figure(size=(1200, 900), fontsize=14)

# Panel (a): Relative error
ax2a = Axis(fig2[1,1], xlabel="Time (s)", ylabel="Relative error",
            title="(a) Relative error in ⟨C⟩")
valid = res_exact.means[:, 3] .> 0.1
if any(valid)
    scatter!(ax2a, res_exact.times[valid], rel_error_C[valid], markersize=10, color=:blue)
    lines!(ax2a, res_exact.times[valid], rel_error_C[valid], color=:blue, linewidth=2)
end

# Panel (b): Error rate
ax2b = Axis(fig2[1,2], xlabel="Time (s)", ylabel="Error rate (s⁻¹)",
            title="(b) Instantaneous error rate")
t_mid = (res_exact.times[1:end-1] .+ res_exact.times[2:end]) ./ 2
vr = abs.(error_rate) .> 0
if any(vr)
    scatter!(ax2b, t_mid[vr], abs.(error_rate[vr]), markersize=10, color=:red)
    lines!(ax2b, t_mid[vr], abs.(error_rate[vr]), color=:red, linewidth=2)
end

# Panel (c): Mean trajectories
ax2c = Axis(fig2[2,1], xlabel="Time (s)", ylabel="Mean population",
            title="(c) Mean trajectories: exact vs adaptive")
lines!(ax2c, res_exact.times, res_exact.means[:, 3], label="⟨C⟩ exact", linewidth=3, color=:black)
scatter!(ax2c, res_exact.times, adaptive_at_exact[:, 3], label="⟨C⟩ adaptive",
         markersize=12, color=:red)
lines!(ax2c, res_exact.times, res_exact.means[:, 1], label="⟨A⟩ exact", linewidth=2,
       color=:gray, linestyle=:dash)
axislegend(ax2c, position=:rt)

# Panel (d): Error scaling
ax2d = Axis(fig2[2,2], xlabel="⟨C⟩ (exact)", ylabel="Absolute error in ⟨C⟩",
            title="(d) Absolute error scaling")
vd = res_exact.means[:, 3] .> 0.1
if any(vd)
    scatter!(ax2d, res_exact.means[vd, 3], abs_error_C[vd], markersize=10, color=:blue)
    x_fit = res_exact.means[vd, 3]
    y_fit = abs_error_C[vd]
    if length(x_fit) >= 2
        slope = sum(x_fit .* y_fit) / sum(x_fit .^ 2)
        x_line = range(0, maximum(x_fit), length=100)
        lines!(ax2d, collect(x_line), slope .* collect(x_line), color=:red, linewidth=2,
               linestyle=:dash, label="slope=$(round(slope; sigdigits=3))")
        axislegend(ax2d, position=:lt)
    end
end

save("/home/adityad/DiscStochSim.jl/paper/error_accumulation_analysis.png", fig2, px_per_unit=3)
println("  Saved error_accumulation_analysis.png")

# ============================================================================
# Summary
# ============================================================================
println("\n" * "=" ^ 60)
println("Summary")
println("=" ^ 60)
println("Flux-aware: final <A>=$(round(res_flux.means[end,1];sigdigits=4)), <C>=$(round(res_flux.means[end,3];sigdigits=4))")
println("No-flux:    final <A>=$(round(res_noflux.means[end,1];sigdigits=4)), <C>=$(round(res_noflux.means[end,3];sigdigits=4))")
println("Exact(1e5): final <A>=$(round(res_exact.means[end,1];sigdigits=4)), <C>=$(round(res_exact.means[end,3];sigdigits=4))")

println("\nError at t=1e5:")
println("  Adaptive <C> = $(round(adaptive_at_exact[end,3]; sigdigits=4))")
println("  Exact    <C> = $(round(res_exact.means[end,3]; sigdigits=4))")
println("  Abs error    = $(round(abs_error_C[end]; sigdigits=3))")
println("  Rel error    = $(round(rel_error_C[end]*100; sigdigits=3))%")
println("\nDone!")
