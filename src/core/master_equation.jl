using SparseArrays

"""
    struct MasterEquationData{T, ElementType, ModelType}

Container for data needed to assemble the master equation operator.

# Fields
- `active_states_vec::Vector{ElementType}`: Currently active states.
- `state_id_map::Dict{ElementType, Int}`: Mapping from each state to its integer index.
- `bound_cond::Function`: Boundary condition function.
- `model::ModelType`: The underlying model of the system.
- `rates::Vector{T}`: Reaction channel rates.
- `t::T`: The current time.
"""
struct MasterEquationData{T, ElementType, ModelType}
    active_states_vec::Vector{ElementType}
    state_id_map::Dict{ElementType, Int}
    bound_cond::Function
    model::ModelType
    rates::Vector{T}
    t::T
end

"""
    MasterEquationNonDiagonal(data::MasterEquationData{T, ElementType, ModelType}) where {T, ElementType, ModelType}

Constructs the sparse non-diagonal part of the master equation operator.

# Returns
A sparse matrix representing the non-diagonal contributions.
"""
function MasterEquationNonDiagonal(data::MasterEquationData{T, ElementType, ModelType}) where {T, ElementType, ModelType}
    # Unpack values
    active_states = data.active_states_vec
    state_id_map = data.state_id_map
    boundary_condition = data.bound_cond
    model = data.model
    rates = data.rates
    t = data.t
    n_active = length(active_states)

    # Get source state for each active state.
    source_states = map(state -> begin
        sources = expand_backward(state, model, boundary_condition)
        filter(x -> x ∈ active_states, sources)
    end, active_states)

    # Construct non-diagonal entries.
    non_diag = [begin
                   j = state_id_map[xj]  # Get connected state index.
                   xᵢ = active_states[i]
                   S = xᵢ - xj         # Reaction vector.
                   k = FindElement(S, model.stoichvecs)
                   α = model.propensities[k](xᵢ, rates, t)
                   (i, j, α)
               end for (i, sources) in enumerate(source_states) for xj in sources]

    # Extract row indices, column indices, and values.
    I = @inbounds [entry[1] for entry in non_diag]
    J = @inbounds [entry[2] for entry in non_diag]
    K = @inbounds [entry[3] for entry in non_diag]

    # Construct and return the sparse matrix.
    sparse(I, J, K, n_active, n_active)
end

"""
    MasterEquationDiagonal(data::MasterEquationData{T, ElementType, ModelType}) where {T, ElementType, ModelType}

Constructs the sparse diagonal part of the master equation operator.

# Returns
A sparse matrix representing the diagonal contributions.
"""
function MasterEquationDiagonal(data::MasterEquationData{T, ElementType, ModelType}) where {T, ElementType, ModelType}
    # Unpack values.
    active_states = data.active_states_vec
    model = data.model
    rates = data.rates
    t = data.t
    n_active = length(active_states)

    # Compute diagonal entries.
    diagonal_entries = [begin
                            αs = [ ( (xᵢ + S) in active_states ) ? model.propensities[k](xᵢ + S, rates, t) : 0.0
                                   for (k, S) in enumerate(model.stoichvecs) ]
                            (i, -sum(αs))
                        end for (i, xᵢ) in enumerate(active_states)]

    I = @inbounds [entry[1] for entry in diagonal_entries]
    K = @inbounds [entry[2] for entry in diagonal_entries]

    # Construct and return the sparse diagonal matrix.
    sparse(I, I, K, n_active, n_active)
end

"""
    MasterEquation(active_states::Set{ElementType}, model::ModelType, rates::Vector{T},
                   boundary_condition::Function, t::T) where {T, ElementType, ModelType}

Assembles the complete master equation operator by summing its non-diagonal and diagonal parts.

# Returns
A sparse matrix representing the full master equation operator.
"""
function MasterEquation(active_states::Set{ElementType}, model::ModelType, rates::Vector{T},
                        boundary_condition::Function, t::T) where {T, ElementType, ModelType}
    # Prepare the active states and mapping.
    active_states_vec = collect(active_states)
    state_id_map = Dict(s => i for (i, s) in enumerate(active_states_vec))

    # Create the data container.
    data = MasterEquationData(active_states_vec, state_id_map, boundary_condition, model, rates, t)

    # Assemble the non-diagonal and diagonal parts.
    S₁ = MasterEquationNonDiagonal(data)
    S₂ = MasterEquationDiagonal(data)

    # Return the complete operator.
    S₁ + S₂
end

# Export the MasterEquation function for external use.
export MasterEquation

