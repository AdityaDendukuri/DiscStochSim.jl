"""
    build_generator(space, model, rates, t) -> (A, in_flow, out_flow)

Build the CME generator matrix from scratch in CSC-friendly COO format.

For each column j (state x_j), iterate over reactions k with stoich vector ν_k:
  - Off-diagonal: A[i,j] += a_k(x_j) where x_j + ν_k = x_i ∈ active set
  - Diagonal:     A[j,j]  = -sum of off-diagonal column entries

Returns `(A, in_flow, out_flow)` where `A = in_flow + out_flow`,
`out_flow = diag(A)` (negative), `in_flow = A - out_flow` (non-negative).
"""
function build_generator(sp::StateSpace{E,T}, model, rates, t::Real) where {E,T}
    n = length(sp)
    n == 0 && return (spzeros(T, 0, 0), spzeros(T, 0, 0), spzeros(T, 0, 0))

    # Pre-estimate nnz: at most n_reactions off-diagonal entries per column + 1 diagonal
    n_rxn = length(model.stoichvecs)
    sizehint = n * (n_rxn + 1)
    I = Vector{Int}(); J = Vector{Int}(); V = Vector{T}()
    sizehint!(I, sizehint); sizehint!(J, sizehint); sizehint!(V, sizehint)

    @inbounds for j in 1:n
        x = sp.states[j]
        col_sum = zero(T)

        for k in 1:n_rxn
            y = x + model.stoichvecs[k]
            i = get(sp.index, y, 0)
            if i != 0
                α = convert(T, model.propensities[k](x, rates, t))
                if α > 0
                    push!(I, i); push!(J, j); push!(V, α)
                    col_sum += α
                end
            end
        end

        if col_sum > 0
            push!(I, j); push!(J, j); push!(V, -col_sum)
        end
    end

    A = sparse(I, J, V, n, n)
    out_flow = spdiagm(0 => diag(A))
    in_flow = A - out_flow
    (A, in_flow, out_flow)
end

"""
    reconstruct_generator(sp, model, rates, t, A_old, gids_old) -> (A, in_flow, out_flow)

Incremental rebuild: reuse columns from `A_old` for states that were present in the
previous step (matched via global IDs). Build new columns from scratch.

Falls back to `build_generator` if overlap is below 30%.
"""
function reconstruct_generator(sp::StateSpace{E,T}, model, rates, t::Real,
                               A_old::SparseMatrixCSC, gids_old::Vector{Int}) where {E,T}
    n_new = length(sp)
    n_new == 0 && return build_generator(sp, model, rates, t)

    gids_new = get_global_ids(sp)

    # Map old global ID → old column index
    old_gid_to_col = Dict{Int,Int}()
    for (col, gid) in enumerate(gids_old)
        old_gid_to_col[gid] = col
    end

    # Count retained columns
    n_retained = count(gid -> haskey(old_gid_to_col, gid), gids_new)
    if n_retained < 0.3 * n_new
        return build_generator(sp, model, rates, t)
    end

    # Build new-row-index lookup for all active states (by global ID)
    gid_to_new_row = Dict{Int,Int}()
    for (row, gid) in enumerate(gids_new)
        gid_to_new_row[gid] = row
    end

    # Old global ID → old row index
    old_gid_to_row = Dict{Int,Int}()
    for (row, gid) in enumerate(gids_old)
        old_gid_to_row[gid] = row
    end

    n_rxn = length(model.stoichvecs)
    sizehint_n = n_new * (n_rxn + 1)
    I = Vector{Int}(); J = Vector{Int}(); V = Vector{T}()
    sizehint!(I, sizehint_n); sizehint!(J, sizehint_n); sizehint!(V, sizehint_n)

    @inbounds for (j_new, gid) in enumerate(gids_new)
        old_col = get(old_gid_to_col, gid, 0)

        if old_col != 0
            # Retained column: copy old off-diagonal entries + add new neighbors
            x = sp.states[j_new]
            col_sum = zero(T)

            # Collect old off-diagonal targets (by new row index) for dedup
            old_targets = Set{Int}()
            for idx in nzrange(A_old, old_col)
                old_row = A_old.rowval[idx]
                old_row == old_col && continue  # skip diagonal
                val = A_old.nzval[idx]
                old_row_gid = gids_old[old_row]
                new_row = get(gid_to_new_row, old_row_gid, 0)
                if new_row != 0
                    push!(I, new_row); push!(J, j_new); push!(V, val)
                    col_sum += val
                    push!(old_targets, new_row)
                end
            end

            # Check for NEW neighbors not in old state space
            for k in 1:n_rxn
                y = x + model.stoichvecs[k]
                i_new = get(sp.index, y, 0)
                if i_new != 0 && i_new != j_new && !(i_new in old_targets)
                    # This neighbor is new — compute propensity
                    α = convert(T, model.propensities[k](x, rates, t))
                    if α > 0
                        push!(I, i_new); push!(J, j_new); push!(V, α)
                        col_sum += α
                    end
                end
            end

            if col_sum > 0
                push!(I, j_new); push!(J, j_new); push!(V, -col_sum)
            end
        else
            # New column: build from scratch
            x = sp.states[j_new]
            col_sum = zero(T)
            for k in 1:n_rxn
                y = x + model.stoichvecs[k]
                i_new = get(sp.index, y, 0)
                if i_new != 0
                    α = convert(T, model.propensities[k](x, rates, t))
                    if α > 0
                        push!(I, i_new); push!(J, j_new); push!(V, α)
                        col_sum += α
                    end
                end
            end
            if col_sum > 0
                push!(I, j_new); push!(J, j_new); push!(V, -col_sum)
            end
        end
    end

    A = sparse(I, J, V, n_new, n_new)
    out_flow = spdiagm(0 => diag(A))
    in_flow = A - out_flow
    (A, in_flow, out_flow)
end
