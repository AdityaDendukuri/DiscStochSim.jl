#!/usr/bin/env julia
# Compare old Pluto notebook logic vs current refactored code
# Run from project root: julia --project test_old_vs_new.jl

using Pkg; Pkg.activate(".")
using Catalyst, SparseArrays, LinearAlgebra, ExponentialUtilities, Expokit

# ─── Load the refactored module ───
using DiscStochSim

rn = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
model = DiscreteStochasticSystem(rn)
rates = [0.04, 3e6, 1.0]
U₀ = CartesianIndex(10_000, 0, 0)
bounds = (0, 100_001)
bc_fn(x) = DiscStochSim.RectLatticeBoundaryCondition(x, bounds)

# ═══════════════════════════════════════════════════════
# OLD PLUTO STYLE (Set-based, expand!, maximum, expv)
# ═══════════════════════════════════════════════════════

# We need the old expand_forward and expand1! from the notebook
function expand_forward_old(x, model, bc_fn)
    states = map(model.stoichvecs) do S
        x_new = x + S
        bc_fn(x_new) ? x_new : nothing
    end
    return Set(filter(!isnothing, states))
end

function expand1_old!(X::Set, model, bc_fn)
    for x in copy(X)
        union!(X, expand_forward_old(x, model, bc_fn))
    end
    return X
end

function expand_old!(X::Set, pₜ::Vector, model, bc_fn, N::Int)
    X_prev = collect(X)
    for _ in 1:N
        expand1_old!(X, model, bc_fn)
    end
    X_vec = collect(X)
    idxs = [findfirst(==(x), X_vec) for x in X_prev]
    qₜ = zeros(length(X_vec))
    qₜ[idxs] = pₜ
    return X, qₜ
end

# Old MasterEquation (from commit 9d63463)
function MasterEquation_old(S::Set{E}, model, rates, bc_fn, t) where E
    S_vec = collect(S)
    n = length(S_vec)
    id = Dict{E,Int}(s => i for (i,s) in enumerate(S_vec))

    I_idx = Int[]; J_idx = Int[]; V_val = Float64[]
    for (j, x) in enumerate(S_vec)
        col_sum = 0.0
        for (k, ν) in enumerate(model.stoichvecs)
            y = x + ν
            i = get(id, y, 0)
            if i != 0
                α = model.propensities[k](x, rates, t)
                if α > 0
                    push!(I_idx, i); push!(J_idx, j); push!(V_val, α)
                    col_sum += α
                end
            end
        end
        if col_sum > 0
            push!(I_idx, j); push!(J_idx, j); push!(V_val, -col_sum)
        end
    end
    A = sparse(I_idx, J_idx, V_val, n, n)
    out_flow = spdiagm(0 => diag(A))
    in_flow = A - out_flow
    return A, in_flow, out_flow
end

# Old robust_purge! (from notebook)
function robust_purge_old!(X::Set, p::Vector, model, rates, t, prob_quantile; flux_tolerance=1e-9)
    X_vec = collect(X)
    n = length(p)

    idx = sortperm(p; rev=false)
    k = round(Int, prob_quantile * length(idx))
    candidate_idxs = k > 0 ? idx[1:k] : Int[]

    flux_vector = zeros(Float64, n)
    for i in eachindex(X_vec)
        weight = sum(prop(X_vec[i], rates, t) for prop in model.propensities)
        flux_vector[i] = p[i] * weight
    end
    total_flux = sum(flux_vector)
    flux_threshold = total_flux * flux_tolerance

    final_idxs = Set{Int}()
    for idx in candidate_idxs
        if flux_vector[idx] < flux_threshold
            push!(final_idxs, idx)
        end
    end

    new_p = [p[i] for i in eachindex(p) if i ∉ final_idxs]
    states_to_remove = Set(X_vec[collect(final_idxs)])
    new_X = setdiff(X, states_to_remove)
    return new_X, new_p, total_flux
end

# ═══════════════════════════════════════════════════════
# NEW REFACTORED STYLE (StateSpace, expand_ssa!, sum, expmv)
# ═══════════════════════════════════════════════════════
# (uses the module's solve, but we'll step manually for diagnostics)

println("=" ^ 80)
println("OLD PLUTO STYLE: expand!(all neighbors) + maximum(out_flux) + expv")
println("=" ^ 80)

let
    𝒮ₜ = Set([U₀])
    pₜ = [1.0]
    t = 0.0
    tf = 1e4
    ϵ_dt = 0.01

    for iter in 1:500
        # 1) Expand
        𝒮ₜ, pₜ = expand_old!(𝒮ₜ, pₜ, model, bc_fn, 1)

        # 2) Build generator
        A, in_flow, out_flow = MasterEquation_old(𝒮ₜ, model, rates, bc_fn, t)

        # 3) Time step
        out_flux = Vector(-out_flow * pₜ)
        max_flux = maximum(out_flux)
        δt = max_flux > 0 ? ϵ_dt / max_flux : (tf - t)
        δt = min(δt, tf - t)

        # 4) Propagate
        pₜ = expv(δt, A, pₜ)

        # 5) Prune + renormalize
        𝒮ₜ, pₜ, total_flux = robust_purge_old!(𝒮ₜ, pₜ, model, rates, t, 0.4; flux_tolerance=1e-9)
        s = sum(pₜ)
        if s > 0; pₜ ./= s; end

        t += δt

        if iter <= 10 || iter % 100 == 0
            println("iter=$iter: t=$(round(t;sigdigits=4)), |S|=$(length(𝒮ₜ)), dt=$(round(δt;sigdigits=3)), max_flux=$(round(max_flux;sigdigits=4)), sum_p=$(round(sum(pₜ);sigdigits=6))")
        end
    end
end

println()
println("=" ^ 80)
println("NEW REFACTORED: expand_ssa! + maximum(out_flux) + expmv")
println("=" ^ 80)

let
    sp = StateSpace{CartesianIndex, Float64}()
    add_state!(sp, U₀, 1.0)
    t = 0.0
    tf = 1e4
    ϵ_dt = 0.01

    for iter in 1:500
        # 1) Expand (SSA style)
        expand_ssa!(sp, model, rates, t, bc_fn; depth=1)

        # 2) Build generator
        A, in_flow, out_flow = build_generator(sp, model, rates, t)

        # 3) Time step
        out_flux = Vector(-out_flow * sp.probs)
        max_flux = maximum(out_flux)
        δt = max_flux > 0 ? ϵ_dt / max_flux : (tf - t)
        δt = min(δt, tf - t)

        # 4) Propagate
        sp.probs = expmv(δt, A, sp.probs)

        # 5) Prune + renormalize
        compress!(sp, model, rates, t, 0.4; flux_tolerance=1e-9)
        renormalize!(sp)

        t += δt

        if iter <= 10 || iter % 100 == 0
            println("iter=$iter: t=$(round(t;sigdigits=4)), |S|=$(length(sp)), dt=$(round(δt;sigdigits=3)), max_flux=$(round(max_flux;sigdigits=4)), sum_p=$(round(sum(sp.probs);sigdigits=6))")
        end
    end
end

println()
println("=" ^ 80)
println("NEW REFACTORED: expand_ssa! + sum(out_flux) + expmv (working commit)")
println("=" ^ 80)

let
    sp = StateSpace{CartesianIndex, Float64}()
    add_state!(sp, U₀, 1.0)
    t = 0.0
    tf = 1e4
    ϵ_dt = 0.01

    for iter in 1:500
        expand_ssa!(sp, model, rates, t, bc_fn; depth=1)
        A, in_flow, out_flow = build_generator(sp, model, rates, t)

        out_flux = Vector(-out_flow * sp.probs)
        total_flux = sum(out_flux)
        δt = total_flux > 0 ? ϵ_dt / total_flux : (tf - t)
        δt = min(δt, tf - t)

        sp.probs = expmv(δt, A, sp.probs)
        compress!(sp, model, rates, t, 0.4; flux_tolerance=1e-9)
        renormalize!(sp)
        t += δt

        if iter <= 10 || iter % 100 == 0
            println("iter=$iter: t=$(round(t;sigdigits=4)), |S|=$(length(sp)), dt=$(round(δt;sigdigits=3)), sum_flux=$(round(total_flux;sigdigits=4)), sum_p=$(round(sum(sp.probs);sigdigits=6))")
        end
    end
end

println()
println("=" ^ 80)
println("OLD PLUTO STYLE with expand! + sum(out_flux)")
println("=" ^ 80)

let
    𝒮ₜ = Set([U₀])
    pₜ = [1.0]
    t = 0.0
    tf = 1e4
    ϵ_dt = 0.01

    for iter in 1:500
        𝒮ₜ, pₜ = expand_old!(𝒮ₜ, pₜ, model, bc_fn, 1)
        A, in_flow, out_flow = MasterEquation_old(𝒮ₜ, model, rates, bc_fn, t)

        out_flux = Vector(-out_flow * pₜ)
        total_flux = sum(out_flux)
        δt = total_flux > 0 ? ϵ_dt / total_flux : (tf - t)
        δt = min(δt, tf - t)

        pₜ = expv(δt, A, pₜ)
        𝒮ₜ, pₜ, _ = robust_purge_old!(𝒮ₜ, pₜ, model, rates, t, 0.4; flux_tolerance=1e-9)
        s = sum(pₜ); if s > 0; pₜ ./= s; end
        t += δt

        if iter <= 10 || iter % 100 == 0
            println("iter=$iter: t=$(round(t;sigdigits=4)), |S|=$(length(𝒮ₜ)), dt=$(round(δt;sigdigits=3)), sum_flux=$(round(total_flux;sigdigits=4)), sum_p=$(round(sum(pₜ);sigdigits=6))")
        end
    end
end
