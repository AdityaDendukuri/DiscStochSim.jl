# Verify FP-FSP E[C] at tf=1.1e6 for the bottleneck system
# A → B (k1=1e-6),  B → B+C (k2=0.1)
# Analytical: E[C(tf)] = (k2/k1)*(k1*tf - 1 + exp(-k1*tf))
#
# Strategy: fixed δt=1000, depth=100 per step (no pruning).
# After 1100 steps the state space has ~110,001 states (chain A + C-chain),
# but the generator is banded so expmv stays fast.

using DiscStochSim, Catalyst, Expokit

rn = @reaction_network begin
    k1, A --> B
    k2, B --> B + C
end

k1, k2 = 1e-6, 0.1
tf     = 1.1e6
rates  = [k1, k2]

E_C_exact = k2/k1 * (k1*tf - 1 + exp(-k1*tf))
println("Analytical E[C](tf=$(tf)) = $(round(E_C_exact; sigdigits=6))")

prob = FSPProblem(rn, CartesianIndex(1, 0, 0), (0.0, tf), rates; bounds=(0, Int(2e6)))

function fpfsp_nopruning(prob::FSPProblem{E,T};
                          δt::Float64 = 1000.0,
                          expand_depth::Int = 100,
                          report_every::Int = 100) where {E,T}
    model     = prob.model
    rates_    = prob.rates
    bc        = prob.bc
    t0, tf_   = prob.tspan
    n_species = length(Tuple(prob.u0))

    sp = StateSpace{E, Float64}()
    add_state!(sp, prob.u0, 1.0)

    t    = Float64(t0)
    iter = 0

    while t < tf_
        iter += 1
        expand!(sp, model, bc; depth=expand_depth)
        A_mat, _, _ = build_generator(sp, model, rates_, t)
        actual_dt   = min(δt, tf_ - t)
        sp.probs    = expmv(actual_dt, A_mat, sp.probs)
        @inbounds for i in eachindex(sp.probs)
            if sp.probs[i] < 0; sp.probs[i] = 0.0; end
        end
        renormalize!(sp)
        t += actual_dt

        if iter % report_every == 0 || t >= tf_
            # quick E[C] from current distribution
            E_C = sum(sp.probs[sp.index[s]] * s[3] for s in sp.states)
            println("  iter=$(lpad(iter,5)), t=$(rpad(round(t; sigdigits=4),10)), " *
                    "|S|=$(lpad(length(sp),7)),  ⟨C⟩=$(round(E_C; sigdigits=5))")
        end
    end

    # Final E[C]
    E_C_final = sum(sp.probs[sp.index[s]] * s[3] for s in sp.states)
    return E_C_final, length(sp)
end

println("\nRunning FP-FSP (δt=1000, depth=3, no pruning)...")
t_start = time()
E_C_fsp, final_size = fpfsp_nopruning(prob; δt=1000.0, expand_depth=3, report_every=110)
elapsed = time() - t_start

println("\n── Results ──────────────────────────────────")
println("  FSP E[C]        = $(round(E_C_fsp;   sigdigits=6))")
println("  Analytical E[C] = $(round(E_C_exact; sigdigits=6))")
rel_err = 100*abs(E_C_fsp - E_C_exact)/E_C_exact
println("  Relative error  = $(round(rel_err; sigdigits=3))%")
println("  Wall time       = $(round(elapsed; digits=1))s")
println("  Final |S|       = $final_size")
