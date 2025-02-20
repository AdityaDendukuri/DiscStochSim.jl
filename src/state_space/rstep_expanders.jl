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

