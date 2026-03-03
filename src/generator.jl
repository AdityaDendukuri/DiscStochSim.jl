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

    # Check overlap; fall back to full build if too few states are retained
    old_gid_set = Set{Int}(gids_old)
    n_retained = count(gid -> gid in old_gid_set, gids_new)
    if n_retained < 0.3 * n_new
        return build_generator(sp, model, rates, t)
    end

    n_rxn = length(model.stoichvecs)
    sizehint_n = n_new * (n_rxn + 1)
    I = Vector{Int}(); J = Vector{Int}(); V = Vector{T}()
    sizehint!(I, sizehint_n); sizehint!(J, sizehint_n); sizehint!(V, sizehint_n)

    @inbounds for j_new in 1:n_new
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

    A = sparse(I, J, V, n_new, n_new)
    out_flow = spdiagm(0 => diag(A))
    in_flow = A - out_flow
    (A, in_flow, out_flow)
end
