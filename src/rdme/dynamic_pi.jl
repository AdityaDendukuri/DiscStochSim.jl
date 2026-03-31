"""
    compute_dynamic_pi(sp, A; n_max) -> Vector{Vector{Float64}}

Compute dynamic-π prolongation weights from the stationary distribution of the
isolated 2-voxel RDME system with generator `A` on state space `sp`.

Solves  Aᵀ π = 0,  Σπ = 1  via sparse LU (SuiteSparse), then extracts the
conditional distribution for every aggregate count nc:

    pi_table[nc+1][n1+1]  =  p(n1 | n1+n2 = nc)     nc = 0 … 2·n_max

Falls back to Binomial(nc, 1/2) for slices with stationary mass < 1e-16
(states in the deep tail that the system never visits).

Works for **any** 2-voxel RDME (birth-death, Schlögl, …).  The resulting table
can be passed directly to `prolong_dynamic` or `build_schlogl_coarse_system_dynamic`.
"""
function compute_dynamic_pi(sp::StateSpace{CartesianIndex{2}, <:Any},
                             A::SparseMatrixCSC{<:Real};
                             n_max::Int = 200)
    # Stationary distribution satisfies  A π = 0,  Σπ = 1.
    # Replace the last row of A with the normalization constraint [1,…,1]
    # and solve the resulting non-singular linear system.
    n = size(A, 1)
    rows_A, cols_A, vals_A = findnz(A)

    mask  = rows_A .!= n
    rows_B = [rows_A[mask]; fill(n, n)]
    cols_B = [cols_A[mask]; collect(1:n)]
    vals_B = [Float64.(vals_A[mask]); ones(n)]

    B   = sparse(rows_B, cols_B, vals_B, n, n)
    rhs = zeros(n);  rhs[n] = 1.0

    π_stat = B \ rhs
    π_stat = max.(π_stat, 0.0)
    π_stat ./= sum(π_stat)

    # Extract conditional slices  p(n1 | n1+n2 = nc)
    n_coarse_max = 2 * n_max
    pi_table = Vector{Vector{Float64}}(undef, n_coarse_max + 1)

    for nc in 0:n_coarse_max
        slice = zeros(nc + 1)     # index n1+1, for n1 = 0:nc

        for n1 in 0:nc
            n2 = nc - n1
            (n1 <= n_max && n2 <= n_max) || continue
            idx = get(sp.index, CartesianIndex(n1, n2), 0)
            idx != 0 && (slice[n1+1] = π_stat[idx])
        end

        total = sum(slice)
        if total > 1e-16
            pi_table[nc+1] = slice ./ total
        else
            # Fallback: Binomial(nc, 1/2) in log-space
            log_half = -nc * log(2.0)
            b = [exp(log_half +
                     sum(log(i) for i in (nc-n1+1):nc; init=0.0) -
                     sum(log(i) for i in 1:n1;          init=0.0)) for n1 in 0:nc]
            s = sum(b)
            pi_table[nc+1] = s > 0 ? b ./ s : fill(1.0/(nc+1), nc+1)
        end
    end

    return pi_table
end

"""
    prolong_dynamic(sp_coarse, pi_table, n_max) -> StateSpace{CartesianIndex{2·K2}}

Prolongation operator  𝓟_dyn: coarse K2 voxels → fine 2·K2 voxels.

For each coarse state (nc₁, …, nc_{K2}), each pair j is independently
split according to dynamic-π weights:

    p_fine(n1, n2, …)  =  p_coarse(nc) · ∏ⱼ pi_table[ncⱼ+1][n_{2j-1}+1]

where n_{2j-1} + n_{2j} = ncⱼ for every pair j.

Works for K2 = 1 (K=1→2), K2 = 2 (K=2→4), etc.  The Binomial fallback
stored in `pi_table` for negligible-mass nc values ensures numerical stability.
"""
function prolong_dynamic(sp_coarse::StateSpace{CartesianIndex{K2}, T},
                          pi_table::Vector{Vector{Float64}},
                          n_max::Int;
                          weight_tol::Float64 = 1e-14) where {K2, T}
    K = 2 * K2
    sp_fine = StateSpace{CartesianIndex{K}, T}()
    fine_buf = zeros(Int, K)

    for i in eachindex(sp_coarse.states)
        nc_tup = Tuple(sp_coarse.states[i])
        p_nc   = sp_coarse.probs[i]
        p_nc > weight_tol || continue
        _prolong_dyn_recurse!(sp_fine, pi_table, nc_tup, p_nc, fine_buf, 1, K2, n_max,
                               weight_tol)
    end

    sp_fine
end

function _prolong_dyn_recurse!(sp_fine::StateSpace{CartesianIndex{K}, T},
                                pi_table::Vector{Vector{Float64}},
                                nc_tup::NTuple{K2, Int},
                                weight::T,
                                fine_buf::Vector{Int},
                                pair_idx::Int,
                                n_pairs::Int,
                                n_max::Int,
                                weight_tol::Float64) where {K, K2, T}
    if pair_idx > n_pairs
        s   = CartesianIndex(NTuple{K, Int}(fine_buf))
        idx = get(sp_fine.index, s, 0)
        if idx == 0
            add_state!(sp_fine, s, weight)
        else
            sp_fine.probs[idx] += weight
        end
        return
    end

    nc    = nc_tup[pair_idx]
    nc + 1 > length(pi_table) && return
    pi_nc = pi_table[nc + 1]

    for n1 in 0:nc
        n2 = nc - n1
        (n1 <= n_max && n2 <= n_max) || continue
        w = weight * pi_nc[n1 + 1]
        w > weight_tol || continue
        fine_buf[2*pair_idx - 1] = n1
        fine_buf[2*pair_idx]     = n2
        _prolong_dyn_recurse!(sp_fine, pi_table, nc_tup, w,
                               fine_buf, pair_idx + 1, n_pairs, n_max, weight_tol)
    end
end
