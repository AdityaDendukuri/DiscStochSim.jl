# Oregonator — KrylovFSP baseline (non-absorbing generator, original algorithm)
# Replicates the comparison from comparison.jl (commit cf90ad7):
#   - Non-absorbing (conservative) generator
#   - expmv from Expokit (τ-lock via NaN/overflow, not FSP mass-loss)
#   - No τ growth: τ stays at τ_init=1e-6 throughout
#   - Lookback window (ℓ=5) for probability-max pruning
#   - 3000 iterations → t≈0.003 (2.7% of tf=0.11), |S| grows to ~13,000
#
# Saves diagnostics to examples/results/oregonator_krylov.jls for plotting.
using DiscStochSim, Catalyst, Serialization
using Expokit: expmv

rn = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end

y1s, y2s, y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s    = 2000.0, 50000.0
rates = [
    mu1s / y2s,
    mu2s / (y1s * y2s),
    (mu1s + mu2s) / y1s,
    2 * mu1s / y1s^2,
    (mu1s + mu2s) / y3s,
]

tf   = 0.11
prob = FSPProblem(rn, CartesianIndex(500, 1000, 2000), (0.0, tf), rates;
                  bounds=(0, 50_000))

# ── Inline old KrylovFSP algorithm (wrapped in function for Julia scope) ───
function run_krylov_oregonator(prob, tf;
    ε=1e-2, τ_init=1e-6, ℓ=5, r=1, max_iter=3000, save_interval=1000)

    E = CartesianIndex{3}
    model  = prob.model
    rates_ = prob.rates
    bc     = prob.bc

    sp = StateSpace{E, Float64}()
    add_state!(sp, prob.u0, 1.0)

    τ_log        = Float64[]
    mass_log     = Float64[]
    time_log     = Float64[]
    walltime_log = Float64[]
    size_log     = Int[]
    reject_count = 0

    prob_history = Vector{Dict{E,Float64}}()

    t    = 0.0
    τ    = τ_init
    iter = 0
    wall_start = time()

    while t < tf && iter < max_iter
        iter += 1

        # Threshold pruning with lookback protection
        threshold = ε / max(length(sp), 1)
        current_probs = Dict{E,Float64}(
            s => sp.probs[sp.index[s]] for s in sp.states)
        push!(prob_history, current_probs)
        if length(prob_history) > ℓ
            popfirst!(prob_history)
        end
        max_probs = Dict{E,Float64}()
        for hist in prob_history
            for (s, p) in hist
                max_probs[s] = max(get(max_probs, s, 0.0), p)
            end
        end
        prune_threshold!(sp, threshold; max_probs=max_probs)

        # SSA-driven expansion + r-step reachability
        expand_ssa!(sp, model, rates_, t, bc; depth=1)
        if r > 0
            expand!(sp, model, bc; depth=r)
        end

        # Non-absorbing generator (rows sum to 0)
        A, _, _ = build_generator(sp, model, rates_, t)

        τ_step = min(τ, tf - t)

        # Attempt matrix exponential; τ-lock if NaN/Inf (numerical overflow)
        local p_new
        expmv_ok = true
        try
            p_new = expmv(τ_step, A, sp.probs)
            if any(isnan, p_new) || any(isinf, p_new)
                expmv_ok = false
            end
        catch
            expmv_ok = false
        end

        if !expmv_ok
            τ /= 2
            reject_count += 1
            push!(τ_log, τ_step)
            push!(mass_log, NaN)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))
            continue
        end

        # Clamp small negatives
        @inbounds for i in eachindex(p_new)
            if p_new[i] < 0; p_new[i] = 0.0; end
        end

        mass = sum(p_new)

        if mass >= 1.0 - ε
            # Accept (τ does NOT grow — this is the τ-lock mechanism)
            sp.probs = p_new
            renormalize!(sp)
            t += τ_step

            push!(τ_log, τ_step)
            push!(mass_log, mass)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))

            if iter % save_interval == 0 || t >= tf
                @info "KrylovFSP: t=$(round(t;sigdigits=4)), |S|=$(length(sp)), " *
                      "τ=$(round(τ_step;sigdigits=3)), mass=$(round(mass;sigdigits=6))"
            end
        else
            # Reject: halve τ
            τ /= 2
            reject_count += 1
            push!(τ_log, τ_step)
            push!(mass_log, mass)
            push!(time_log, t)
            push!(walltime_log, time() - wall_start)
            push!(size_log, length(sp))
        end
    end

    t_wall = time() - wall_start
    return (τ_log=τ_log, mass_log=mass_log, time_log=time_log,
            walltime_log=walltime_log, size_log=size_log,
            reject_count=reject_count, total_iters=iter,
            t_reached=t, t_wall=t_wall)
end

println("Running KrylovFSP Oregonator (non-absorbing generator, max_iter=3000)...")
res = run_krylov_oregonator(prob, tf)

println("Done: $(res.total_iters) iters, $(res.reject_count) rejections, " *
        "$(round(res.t_wall; digits=1))s")
println("t_reached = $(round(res.t_reached; sigdigits=4)) / $tf " *
        "($(round(100*res.t_reached/tf; digits=2))%)")
if !isempty(res.size_log)
    println("max|S|=$(maximum(res.size_log)), final|S|=$(res.size_log[end])")
end
valid_τ = filter(x -> !isnan(x) && x > 0, res.τ_log)
if !isempty(valid_τ)
    println("τ range: $(minimum(valid_τ)) → $(maximum(valid_τ))")
end

results_dir = joinpath(@__DIR__, "results")
mkpath(results_dir)
out = Dict(
    "τ_log"       => res.τ_log,
    "time_log"    => res.time_log,
    "size_log"    => res.size_log,
    "total_iters" => res.total_iters,
    "reject_count"=> res.reject_count,
    "wall_time"   => res.t_wall,
    "t_reached"   => res.t_reached,
    "S_max"       => isempty(res.size_log) ? 0 : maximum(res.size_log),
    "S_final"     => isempty(res.size_log) ? 0 : res.size_log[end],
)
serialize(joinpath(results_dir, "oregonator_krylov.jls"), out)
println("Saved diagnostics → examples/results/oregonator_krylov.jls")
