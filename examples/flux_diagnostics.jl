# Flux diagnostics: log Φ_total, Φ_max, |J|, dt at every step
# for Robertson (tf=50) and Oregonator (tf=0.02, ~half-oscillation)
# Goal: understand the relationship between these quantities to design
# a principled unified time-stepping rule.

using DiscStochSim, Catalyst, SparseArrays, ExponentialUtilities, Statistics

# ─── Generic diagnostic solver ───────────────────────────────────────────────
function run_flux_diagnostics(prob; ε_dt, prob_quantile, flux_tolerance,
                               expand_method, max_iter, tf_override=nothing)
    model = prob.model; rates = prob.rates; bc = prob.bc
    t0, tf = prob.tspan
    tf = something(tf_override, tf)

    sp = DiscStochSim.StateSpace{typeof(prob.u0), Float64}()
    DiscStochSim.add_state!(sp, prob.u0, 1.0)

    A_old = spzeros(Float64, 0, 0)
    gids_old = Int[]
    t = Float64(t0); iter = 0

    # Log: t, |J|, Φ_total, Φ_max, dt
    log_t      = Float64[]
    log_J      = Int[]
    log_Φtot   = Float64[]
    log_Φmax   = Float64[]
    log_dt     = Float64[]

    while t < tf && iter < max_iter
        iter += 1

        if expand_method === :ssa
            DiscStochSim.expand_ssa!(sp, model, rates, t, bc; depth=1)
        else
            DiscStochSim.expand!(sp, model, bc; depth=1)
        end

        A, in_flow, out_flow = DiscStochSim.reconstruct_generator(
            sp, model, rates, t, A_old, gids_old)
        A_old = A; gids_old = DiscStochSim.get_global_ids(sp)

        out_flux  = Vector(-out_flow * sp.probs)
        Φ_total   = sum(out_flux)
        Φ_max     = maximum(out_flux)
        J_size    = length(sp)

        # Record BEFORE taking step (so we know what governed the dt choice)
        push!(log_t,    t)
        push!(log_J,    J_size)
        push!(log_Φtot, Φ_total)
        push!(log_Φmax, Φ_max)

        # Use :total stepping for Oregonator; record both quantities for both
        dt_total = Φ_total > 0 ? ε_dt / Φ_total : (tf - t)
        push!(log_dt, min(dt_total, tf - t))

        dt = min(dt_total, tf - t)

        sp.probs = expv(dt, A, sp.probs)
        DiscStochSim.compress!(sp, model, rates, t, prob_quantile;
                               flux_tolerance=flux_tolerance)
        DiscStochSim.renormalize!(sp)
        t += dt
    end

    (log_t=log_t, log_J=log_J, log_Φtot=log_Φtot, log_Φmax=log_Φmax,
     log_dt=log_dt, t_final=t, total_iters=iter)
end

# ─── Robertson (tf=50, :stoich) ───────────────────────────────────────────────
println("=== Robertson (tf=50, ε_dt=0.01) ===")
rn_rober = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
rates_rober = [0.04, 3e6, 1.0]
prob_rober = FSPProblem(rn_rober, CartesianIndex(10_000, 0, 0), (0.0, 50.0), rates_rober;
                        bounds=(0, 100_001))
t0 = time()
dr = run_flux_diagnostics(prob_rober; ε_dt=0.01, prob_quantile=0.4,
                           flux_tolerance=1e-9, expand_method=:stoich,
                           max_iter=500_000)
println("Done: $(dr.total_iters) iters, $(round(time()-t0;digits=1))s, t_final=$(dr.t_final)")

# Key statistics
ratio_r = dr.log_Φtot ./ dr.log_Φmax   # = Φ_total / Φ_max ≈ |J| if uniform?
println("\nRobertson step-level statistics:")
println("  |J|:           min=$(minimum(dr.log_J)), median=$(round(Int,median(dr.log_J))), max=$(maximum(dr.log_J))")
println("  Φ_total:       min=$(round(minimum(dr.log_Φtot);sigdigits=3)), max=$(round(maximum(dr.log_Φtot);sigdigits=3))")
println("  Φ_max:         min=$(round(minimum(dr.log_Φmax);sigdigits=3)), max=$(round(maximum(dr.log_Φmax);sigdigits=3))")
println("  Φ_tot/Φ_max:   min=$(round(minimum(ratio_r);sigdigits=3)), median=$(round(median(ratio_r);sigdigits=3)), max=$(round(maximum(ratio_r);sigdigits=3))")
println("  Φ_tot/(|J|·Φ_max): $(round(median(ratio_r ./ dr.log_J);sigdigits=3)) (median; 1.0 = perfectly uniform)")

# What dt(max) gives vs dt(total):
dt_from_max   = 0.01 ./ dr.log_Φmax
dt_from_total = dr.log_dt
speedup = dt_from_max ./ dt_from_total
println("  dt(max)/dt(total): min=$(round(minimum(speedup);sigdigits=3)), median=$(round(median(speedup);sigdigits=3)), max=$(round(maximum(speedup);sigdigits=3))")

# ─── Oregonator (tf=0.02, :ssa) ───────────────────────────────────────────────
println("\n=== Oregonator (tf=0.02, ε_dt=1.0) ===")
rn_oreg = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end
y1s, y2s, y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s = 2000.0, 50000.0
rates_oreg = [mu1s/y2s, mu2s/(y1s*y2s), (mu1s+mu2s)/y1s,
              2*mu1s/y1s^2, (mu1s+mu2s)/y3s]
prob_oreg = FSPProblem(rn_oreg, CartesianIndex(500, 1000, 2000), (0.0, 0.02), rates_oreg;
                       bounds=(0, 50_000))
t0 = time()
do_ = run_flux_diagnostics(prob_oreg; ε_dt=1.0, prob_quantile=0.1,
                            flux_tolerance=1.0, expand_method=:ssa,
                            max_iter=50_000)
println("Done: $(do_.total_iters) iters, $(round(time()-t0;digits=1))s, t_final=$(do_.t_final)")

ratio_o = do_.log_Φtot ./ do_.log_Φmax
println("\nOregonator step-level statistics:")
println("  |J|:           min=$(minimum(do_.log_J)), median=$(round(Int,median(do_.log_J))), max=$(maximum(do_.log_J))")
println("  Φ_total:       min=$(round(minimum(do_.log_Φtot);sigdigits=3)), max=$(round(maximum(do_.log_Φtot);sigdigits=3))")
println("  Φ_max:         min=$(round(minimum(do_.log_Φmax);sigdigits=3)), max=$(round(maximum(do_.log_Φmax);sigdigits=3))")
println("  Φ_tot/Φ_max:   min=$(round(minimum(ratio_o);sigdigits=3)), median=$(round(median(ratio_o);sigdigits=3)), max=$(round(maximum(ratio_o);sigdigits=3))")
println("  Φ_tot/(|J|·Φ_max): $(round(median(ratio_o ./ do_.log_J);sigdigits=3)) (median; 1.0 = perfectly uniform)")

dt_from_max_o   = 1.0 ./ do_.log_Φmax
dt_from_total_o = do_.log_dt
speedup_o = dt_from_max_o ./ dt_from_total_o
println("  dt(max)/dt(total): min=$(round(minimum(speedup_o);sigdigits=3)), median=$(round(median(speedup_o);sigdigits=3)), max=$(round(maximum(speedup_o);sigdigits=3))")

# ─── Summary table ────────────────────────────────────────────────────────────
println("\n=== Summary: what governs the difference? ===")
println("System       | median Φ_tot/Φ_max | median |J| | Φ_tot/(|J|·Φ_max)")
println("Robertson    | $(round(median(ratio_r);sigdigits=3))               | $(round(Int,median(dr.log_J)))          | $(round(median(ratio_r./dr.log_J);sigdigits=3))")
println("Oregonator   | $(round(median(ratio_o);sigdigits=3))               | $(round(Int,median(do_.log_J)))          | $(round(median(ratio_o./do_.log_J);sigdigits=3))")
println("\nIf Φ_tot/(|J|·Φ_max) ≈ 1: distribution is uniform across J (Φ_max ≈ Φ_tot/|J|)")
println("If Φ_tot/(|J|·Φ_max) >> 1: one high-flux state dominates (Φ_max >> Φ_tot/|J|)")
