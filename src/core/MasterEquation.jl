using SparseArrays


"""
    MasterEquationBuilder{T}(space, model, rates)

Container for building (and reconstructing) sparse CME generators
restricted to the current active `space`.
- `model.stoichvecs :: Vector{E}`
- `model.propensities :: Vector{(state, rates, t) -> T}`
"""
struct MasterEquationBuilder{T}
    space::StateSpace
    model::Any
    rates::Vector{T}
end

"""
    build_sparse_matrix(builder, t) -> SparseMatrixCSC{T,Int}

Build the restricted generator from scratch at time `t`.
Columns correspond to sources; off-diagonals are inflow; diagonal is
minus total outflow that **stays** in the active set.
"""
function build_sparse_matrix(builder::MasterEquationBuilder{T}, t::Real) where {T}
    space = builder.space
    model = builder.model
    rates = builder.rates

    n = length(space.state_vec)
    # Map current active states to compact indices
    state_idx = Dict(s => i for (i, s) in enumerate(space.state_vec))

    I = Int[]; J = Int[]; V = T[]

    @inbounds for (j, x) in enumerate(space.state_vec)
        colsum = zero(T)

       # Off-diagonals: x --(k)--> y that remains inside active set
        for (k, ν) in enumerate(model.stoichvecs)
            y = x + ν
            if haskey(state_idx, y)
                α = model.propensities[k](x, rates, t)
                i = state_idx[y]
                if i != j
                    push!(I, i); push!(J, j); push!(V, α)
                    colsum += α
                end
            end
        end

        # Diagonal: negative of kept outflow
        push!(I, j); push!(J, j); push!(V, -colsum)
    end

    return sparse(I, J, V, n, n)
end

"""
    reconstruct_sparse_matrix(builder, A_old, global_old, t; recompute_retained=false)

Reconstruct the restricted generator by **columns** (CSC-friendly):
- Copy retained source columns from `A_old`, filtering rows that are no longer active,
  and recompute the diagonal as `-sum(kept off-diagonals)`.
- Build columns for newly-added sources from scratch.
- If `recompute_retained=true` (e.g. time-dependent propensities with `t` advanced),
  also rebuild the retained columns.

Falls back to a full rebuild if overlap < 30%.
"""
function reconstruct_sparse_matrix(
    builder::MasterEquationBuilder{T},
    A_old::SparseMatrixCSC{T, Int},
    global_old::Vector{Int},
    t::Real;
    recompute_retained::Bool=false
) where {T}

    space = builder.space
    model = builder.model
    rates = builder.rates

    global_new = get_global_indices(space)
    n_new = length(global_new)
    new_compact = Dict(g => i for (i, g) in enumerate(global_new))
    old_compact = Dict(g => i for (i, g) in enumerate(global_old))

    retained   = Set(global_new) ∩ Set(global_old)
    new_cols   = setdiff(Set(global_new), Set(global_old))

    if length(retained) < 0.3 * n_new
        return build_sparse_matrix(builder, t)
    end

    I = Int[]; J = Int[]; V = T[]
    rows = rowvals(A_old); vals = nonzeros(A_old)

    # helper: build a new source column from scratch
    function build_new_column!(src_global::Int, new_j::Int)
        x = space.global_to_state[src_global]
        colsum = zero(T)
        @inbounds for (k, ν) in enumerate(model.stoichvecs)
            y = x + ν
            if haskey(space.state_to_global, y)
                g_tgt = space.state_to_global[y]
                if space.global_active[g_tgt] && haskey(new_compact, g_tgt)
                    i = new_compact[g_tgt]
                    α = model.propensities[k](x, rates, t)
                    if i != new_j
                        push!(I, i); push!(J, new_j); push!(V, α)
                        colsum += α
                    end
                end
            end
        end
        push!(I, new_j); push!(J, new_j); push!(V, -colsum)
        return nothing
    end

    # 1) retained source columns
    @inbounds for src_global in retained
        new_j = new_compact[src_global]

        if recompute_retained
            build_new_column!(src_global, new_j)
            continue
        end

        old_j = old_compact[src_global]
        colsum = zero(T)

        for nz in nzrange(A_old, old_j)
            r_old = rows[nz]
            if r_old == old_j
                continue # skip old diag; we’ll recompute diag to keep exact col-sum zero
            end
            tgt_global = global_old[r_old]
            if haskey(new_compact, tgt_global)
                i_new = new_compact[tgt_global]
                v = vals[nz]
                if i_new != new_j
                    push!(I, i_new); push!(J, new_j); push!(V, v)
                    colsum += v
                end
            end
        end

        push!(I, new_j); push!(J, new_j); push!(V, -colsum)
    end

    # 2) newly added source columns
    @inbounds for src_global in new_cols
        new_j = new_compact[src_global]
        build_new_column!(src_global, new_j)
    end

    return sparse(I, J, V, n_new, n_new)
end
export MasterEquationBuilder, build_sparse_matrix, reconstruct_sparse_matrix
