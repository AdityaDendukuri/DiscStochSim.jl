"""
    expand1!(X, model, boundary_condition)

In-place expansion of the set `X` by applying `ExpandForward` to each element.
Each element `x` is replaced by the union `x ∪ ExpandForward(x, model, boundary_condition)`.
"""
function expand1!(X::Set{Element}, model::Model, boundary_condition::Function) where {Element,Model}
    # Iterate over a copy so we can safely modify X during iteration.
    for x in copy(X)
        union!(X, expand_forward(x, model, boundary_condition))
    end
    return X
end

"""
    expand!(X, model, boundary_condition, N)

Repeatedly applies `expand1!` to `X` for `N` iterations.
"""
function expand!(X::Set{Element}, model::Model, boundary_condition::Function, N::Int) where {Element,Model}
    for _ in 1:N
        expand1!(X, model, boundary_condition)
    end
    return X
end

"""
    expand!(X, pₜ, model, boundary_condition, N)

Expands `X` (using `expand1!`) for `N` iterations and updates the probability vector.
The vector `pₜ` corresponds to the probabilities of the elements of `X` *before* expansion.
After expansion, the probabilities are stored in `qₜ` at the positions corresponding
to the original elements.
Returns the expanded set and the new probability vector.
"""
function expand!(X::Set{Element}, pₜ::Vector, model::Model, boundary_condition::Function, N::Int) where {Element,Model}
    X_prev = collect(X)
    for _ in 1:N
        expand1!(X, model, boundary_condition)
    end
    X_vec = collect(X)
    idxs = [findfirst(==(x), X_vec) for x in X_prev]
    qₜ = zeros(length(X_vec))
    qₜ[idxs] = pₜ
    return X, qₜ
end

"""
    purge!(X, p, percentage)

Removes from `X` the elements with the lowest probability values.
The indices to remove are determined by `findLowestValuesPercent_naive(p, percentage)`.
Returns the purged set and the updated probability vector.
"""
function purge!(X::Set{Element}, p::Vector{T}, percentage::Number) where {Element,T}
    X_vec = collect(X)
    idxs = findLowestValuesPercent_naive(p, percentage)
    new_p = [p[i] for i in eachindex(p) if i ∉ idxs]
    new_X = setdiff(X, Set(X_vec[idxs]))
    return new_X, new_p
end

"""
    purge!(X, p, percentage)

Removes from `X` the elements with the lowest probability values.
The indices to remove are determined by `findLowestValuesPercent_naive(p, percentage)`.
Returns the purged set and the updated probability vector.
"""
function purge!(X::Set{Element}, p::Vector{T}, percentage::Number) where {Element,T}
    X_vec = collect(X)
    idxs = findLowestValuesPercent_naive(p, percentage)
    new_p = [p[i] for i in eachindex(p) if i ∉ idxs]
    new_X = setdiff(X, Set(X_vec[idxs]))
    return new_X, new_p
end

function purge!(
    X::Set{Element},
    p::Vector{T},
    flux_vector::AbstractVector,   # a0(x) per state (sum of propensities)
    prob_quantile::Real;           # in [0,1], fraction of states (by count) to drop by φ
    renormalize::Bool = false      # if true, renormalize remaining mass
) where {Element,T<:Real}

    @assert length(p) == length(flux_vector) "p and flux_vector size mismatch"
    @assert 0 ≤ prob_quantile ≤ 1 "prob_quantile must be in [0,1]"

    X_vec = collect(X)

    # φ_i = p_i * a0(x_i)  (outgoing probability flux from state i)
    ϕ = p .* flux_vector

    # Rank states by ϕ ascending; choose bottom fraction to remove
    idx = sortperm(ϕ; rev=false)
    k = round(Int, prob_quantile * length(idx))            # number of states to remove
    drop_idxs = k > 0 ? idx[1:k] : Int[]

    # Build new containers
    keep_mask = trues(length(p))
    keep_mask[drop_idxs] .= false
    new_p = p[keep_mask]
    new_X = Set(X_vec[findall(keep_mask)])

    # Optional renormalization (typical FSP uses a sink instead; renorm biases)
    if renormalize
        s = sum(new_p)
        if s > 0
            new_p ./= s
        end
    end

    # Diagnostics (optional/handy)
    total_out_flux = sum(ϕ)
    removed_out_flux = sum(ϕ[drop_idxs])
    kept_out_flux = total_out_flux - removed_out_flux

    return new_X, new_p, removed_out_flux, kept_out_flux
end


function purge2!(
    X::Set{Element}, 
    p::Vector{T}, 
    flux_vector::Vector,
    model::Model,
    rates,
    t,
    prob_quantile::Number;
    flux_tolerance::Number = 1e-9 
) where {Element, T, Model}
    
    X_prev = collect(X)
    
    candidate_idxs = findLowestValuesPercent_naive(p, prob_quantile)
    
    total_flux = sum(flux_vector)
    flux_threshold = total_flux * flux_tolerance
    
    final_idxs = Int[]  
    for idx in candidate_idxs  
        if flux_vector[idx] < flux_threshold
            push!(final_idxs, idx)  
        end
    end
    
    
    new_p = [p[i] for i in eachindex(p) if i ∉ final_idxs]
    
    states_to_remove = Set(X_prev[i] for i in final_idxs)
    new_X = setdiff(X, states_to_remove)
    
    return new_X, new_p, flux_threshold, total_flux
end

function purge1!(
    X::Set{Element}, 
    p::Vector{T}, 
    model::Model,
    rates,
    t,
    prob_quantile::Number;
    flux_tolerance::Number = 1e-9 
) where {Element, T, Model}
    
    X_prev = collect(X)

    #  Candidate Selection: Find states with low probability mass
    candidate_idxs = findLowestValuesPercent_naive(p, prob_quantile)

    # compute total flux  
    flux_vector = zeros(T, length(p))
    for i in eachindex(X_prev)
        current_state = X_prev[i]
        weight = sum(prop(current_state, rates, t) for prop in model.propensities)
        flux_vector[i] = p[i] * weight
    end
    total_flux = sum(flux_vector)
    flux_threshold = total_flux * flux_tolerance

    low_flux_states = findall(x->x<flux_threshold, flux_vector)
    if length(low_flux_states) != 0
        final_idxs = candidate_idxs[low_flux_states[1]]
    else
        final_idxs = candidate_idxs
    end
    new_p = [p[i] for i in eachindex(p) if i ∉ final_idxs]
    new_X = setdiff!(X, Set(X_prev[final_idxs]))
    
    return new_X, new_p, flux_threshold, total_flux
end




# -------------------------------------------------------------------
# Overloads that use additional parameters (rates and time `t`)
# -------------------------------------------------------------------


"""
    expand1!(X, model, rates, t, boundary_condition)

Performs an in-place expansion of `X` using `SSA_STEP` with the provided `rates` and time `t`.
For each element `x` in `X`, computes `new_x = SSA_STEP(x, model, rates, t)`.
If `boundary_condition(new_x)` is `true`, `new_x` is added to `X`; otherwise `x` is retained.
"""
function expand1!(X::Set{Element}, model::Model, rates::AbstractArray, t::Number, boundary_condition::Function) where {Element,Model}
    for x in copy(X)
        new_x = ssa_step(x, model, rates, t)
        # Add new_x if it meets the boundary condition; otherwise, re-add x (which is redundant if x ∈ X)
        union!(X, boundary_condition(new_x) ? Set([new_x]) : Set([x]))
    end
    return X
end

"""
    expand!(X, model, rates, t, boundary_condition, N)

Repeatedly applies the above `expand1!` (with `rates` and `t`) for `N` iterations.
    println(X)
"""
function expand!(X::Set{Element}, model::Model, rates::AbstractArray, t::Number, boundary_condition::Function, N::Int) where {Element,Model}
    for _ in 1:N
        expand1!(X, model, rates, t, boundary_condition)
    end
    return X
end

"""
    expand!(X, pₜ, model, rates, t, boundary_condition, N)

Expands `X` using `SSA_STEP` (with `rates` and `t`) for `N` iterations and updates the probability vector.
Returns the expanded set and the new probability vector.
"""
function expand!(X::Set{Element}, pₜ::Vector, model::Model, rates::AbstractArray, t::Number, boundary_condition::Function, N::Int) where {Element,Model}
    X_prev = collect(X)
    for _ in 1:N
        expand1!(X, model, rates, t, boundary_condition)
    end
    X_vec = collect(X)
    idxs = [findfirst(==(x), X_vec) for x in X_prev]
    qₜ = zeros(length(X_vec))
    qₜ[idxs] = pₜ
    return X, qₜ
end

# Export the public functions.
export expand1!, expand!, purge!

