# ─── Patch QSD: generalized M-voxel neighborhood prolongation ────────────────
#
# Replaces pair-based Binomial / Multinomial prolongation with a physics-
# informed QSD that accounts for:
#   • The internal adjacency of the fine voxels inside each coarse patch
#   • Bimolecular A+B→∅ reactions that create A-B anti-correlation
#   • Arbitrary patch size M (1D pairs M=2, 2D 2×2 blocks M=4, ...)
#
# Fast-diffusion limit: QSD → Multinomial(nA_c, 1/M,...) ⊗ Multinomial(nB_c, 1/M,...)
# Strong-reaction limit: QSD concentrates on states where A and B are separated

# ─── State enumeration ───────────────────────────────────────────────────────

function _enum_compositions!(out::Vector{Vector{Int}}, buf::Vector{Int},
                              pos::Int, remaining::Int, M::Int)
    if pos == M
        buf[pos] = remaining
        push!(out, copy(buf))
        return
    end
    for n in 0:remaining
        buf[pos] = n
        _enum_compositions!(out, buf, pos + 1, remaining - n, M)
    end
end

function _enum_compositions(N::Int, M::Int)::Vector{Vector{Int}}
    out = Vector{Int}[]
    _enum_compositions!(out, zeros(Int, M), 1, N, M)
    out
end

# ─── PatchQSD struct ─────────────────────────────────────────────────────────

"""
    PatchQSD

Precomputed quasi-stationary distributions for an M-voxel patch, 2 species (A, B).

For each coarse total count (nA_c, nB_c), stores the QSD over all ways to
distribute those molecules among the M fine voxels in the patch, conditioned
on the (nA_c, nB_c) subspace (i.e., treating A+B→∅ as absorbing).

`tables[nA_c+1, nB_c+1]` = sparse list of `(nA_in_patch, nB_in_patch, weight)`,
where `nA_in_patch[m]` and `nB_in_patch[m]` are the molecule counts in
patch-local voxel m (1..M).

Recovered limits:
- k_rxn = 0 : QSD → Multinomial(nA_c, 1/M,...) ⊗ Multinomial(nB_c, 1/M,...)
- d_A, d_B → ∞ relative to k_rxn : same limit
- M = 2, k_rxn = 0 : Binomial(nA_c, 1/2) ⊗ Binomial(nB_c, 1/2)
"""
struct PatchQSD
    M      :: Int
    edges  :: Vector{Tuple{Int,Int}}   # undirected patch-local adjacency (1-indexed)
    nA_max :: Int
    nB_max :: Int
    d_A    :: Float64
    d_B    :: Float64
    k_rxn  :: Float64
    # tables[nA_c+1, nB_c+1] = sparse (nA_in_patch, nB_in_patch, weight)
    tables :: Matrix{Vector{Tuple{Vector{Int}, Vector{Int}, Float64}}}
end

# ─── Local QSD builder ───────────────────────────────────────────────────────

function _build_patch_qsd_2s(
    nA_c::Int, nB_c::Int,
    M::Int, edges::Vector{Tuple{Int,Int}},
    d_A::Float64, d_B::Float64, k_rxn::Float64,
    tol::Float64
)::Vector{Tuple{Vector{Int}, Vector{Int}, Float64}}

    # Trivial: both zero
    if nA_c == 0 && nB_c == 0
        return [(zeros(Int, M), zeros(Int, M), 1.0)]
    end

    A_comps = _enum_compositions(nA_c, M)
    B_comps = _enum_compositions(nB_c, M)
    nS_A    = length(A_comps)
    nS_B    = length(B_comps)
    n_states = nS_A * nS_B

    # States indexed as (iA-1)*nS_B + iB
    # Build reverse lookup: (nA_vec, nB_vec) → linear index
    state_idx = Dict{Tuple{Vector{Int},Vector{Int}}, Int}()
    states    = Vector{Tuple{Vector{Int},Vector{Int}}}(undef, n_states)
    idx = 1
    for nA_vec in A_comps, nB_vec in B_comps
        states[idx] = (nA_vec, nB_vec)
        state_idx[(nA_vec, nB_vec)] = idx
        idx += 1
    end

    if n_states == 1
        s = states[1]
        return [(s[1], s[2], 1.0)]
    end

    # Build local generator Q (column convention: Q[j,i] = rate i→j)
    Q = zeros(n_states, n_states)

    for (i, (nA_vec, nB_vec)) in enumerate(states)
        # Diffusion A (both directions per undirected edge)
        for (ki, kj) in edges
            for (from_k, to_k) in ((ki, kj), (kj, ki))
                r = d_A * nA_vec[from_k]
                iszero(r) && continue
                nA_new      = copy(nA_vec)
                nA_new[from_k] -= 1
                nA_new[to_k]   += 1
                j = get(state_idx, (nA_new, nB_vec), 0)
                j > 0 && (Q[j, i] += r)
                Q[i, i] -= r
            end
        end

        # Diffusion B (both directions per undirected edge)
        for (ki, kj) in edges
            for (from_k, to_k) in ((ki, kj), (kj, ki))
                r = d_B * nB_vec[from_k]
                iszero(r) && continue
                nB_new      = copy(nB_vec)
                nB_new[from_k] -= 1
                nB_new[to_k]   += 1
                j = get(state_idx, (nA_vec, nB_new), 0)
                j > 0 && (Q[j, i] += r)
                Q[i, i] -= r
            end
        end

        # A+B→∅ in each fine voxel k: absorbing (leaves (nA_c,nB_c) subspace)
        if k_rxn > 0
            for k in 1:M
                Q[i, i] -= k_rxn * nA_vec[k] * nB_vec[k]
            end
        end
    end

    # Lead eigenvector (eigenvalue closest to 0) = QSD
    vals, vecs = eigen(Q)
    lead_idx   = argmax(real.(vals))
    π          = real.(vecs[:, lead_idx])
    π          = abs.(π)
    s          = sum(π)
    s > 0 ? (π ./= s) : (π .= 1.0 / n_states)
    π[π .< tol] .= 0.0
    s2 = sum(π)
    s2 > 0 && (π ./= s2)

    result = Tuple{Vector{Int}, Vector{Int}, Float64}[]
    for (i, (nA_vec, nB_vec)) in enumerate(states)
        w = π[i]
        w > tol && push!(result, (nA_vec, nB_vec, w))
    end
    result
end

# ─── Public builder ──────────────────────────────────────────────────────────

"""
    build_patch_qsd(nA_max, nB_max, M, edges, d_A, d_B, k_rxn; tol) → PatchQSD

Precompute the within-patch QSD for all coarse totals (nA_c, nB_c) ∈ [0,nA_max]×[0,nB_max].

Arguments
---------
- `nA_max`, `nB_max` : maximum coarse total count per species to tabulate
- `M`                : patch size (2=1D pair, 4=2D 2×2 block, 8=3D 2×2×2 block)
- `edges`            : undirected patch-local adjacency list (indices 1..M)
- `d_A`, `d_B`       : internal (fine-scale) diffusion rates per species
- `k_rxn`            : A+B→∅ bimolecular rate

Standard edge lists
-------------------
    patch_edges_1d()                       → [(1,2)]              M=2
    patch_edges_2d()                       → [(1,2),(3,4),(1,3),(2,4)]  M=4
    patch_edges_from_cmap(cmap, fine_grid) → derived from CoarseningMap
"""
function build_patch_qsd(
    nA_max::Int, nB_max::Int,
    M::Int,
    edges::Vector{Tuple{Int,Int}},
    d_A::Float64, d_B::Float64, k_rxn::Float64;
    tol::Float64 = 1e-14
)::PatchQSD
    tables = Matrix{Vector{Tuple{Vector{Int},Vector{Int},Float64}}}(undef,
                                                                     nA_max + 1,
                                                                     nB_max + 1)
    for nA_c in 0:nA_max, nB_c in 0:nB_max
        tables[nA_c+1, nB_c+1] = _build_patch_qsd_2s(
            nA_c, nB_c, M, edges, d_A, d_B, k_rxn, tol)
    end
    PatchQSD(M, edges, nA_max, nB_max, d_A, d_B, k_rxn, tables)
end

# ─── Standard patch topologies ───────────────────────────────────────────────

"""Patch-local edges for a 1D pair (M=2): single edge 1↔2."""
patch_edges_1d() = [(1, 2)]

"""Patch-local edges for a 2D 2×2 block (M=4, row-major: TL=1,TR=2,BL=3,BR=4)."""
patch_edges_2d() = [(1,2), (3,4), (1,3), (2,4)]

"""
    patch_edges_from_cmap(cmap, fine_grid) → Vector{Tuple{Int,Int}}

Derive patch-local adjacency from a CoarseningMap and the fine Grid2D.
All patches in a uniform grid share the same topology; computed from the first patch.
"""
function patch_edges_from_cmap(cmap::CoarseningMap, fine_grid::Grid2D)
    patch = cmap.coarse_to_fine[1]
    K_y   = fine_grid.K_y
    edges = Tuple{Int,Int}[]
    M     = length(patch)
    for a in 1:M, b in (a+1):M
        k1, k2 = patch[a], patch[b]
        i1, j1 = (k1 - 1) ÷ K_y + 1, (k1 - 1) % K_y + 1
        i2, j2 = (k2 - 1) ÷ K_y + 1, (k2 - 1) % K_y + 1
        ((i1 == i2 && abs(j1 - j2) == 1) || (j1 == j2 && abs(i1 - i2) == 1)) &&
            push!(edges, (a, b))
    end
    edges
end

# ─── Prolongation with CoarseningMap ─────────────────────────────────────────

"""
    prolong_patch_qsd(sp_coarse, cmap, ::Val{N_F}, qt; weight_tol, prob_tol)
        → StateSpace{CartesianIndex{N_F}}

QSD prolongation for 2-species M-voxel patches using a `PatchQSD` table and
a `CoarseningMap`. Replaces `prolong2s(sp_coarse, cmap, ::Val{N_F})`.

For each coarse state, each patch J is split according to the precomputed QSD
for its (nA_c, nB_c) totals. In the fast-diffusion limit this recovers
Multinomial exactly.
"""
function prolong_patch_qsd(
    sp_coarse::StateSpace{CartesianIndex{N_C}, T},
    cmap::CoarseningMap,
    ::Val{N_F},
    qt::PatchQSD;
    weight_tol::Float64 = 1e-14,
    prob_tol::Float64   = 0.0
) where {N_C, N_F, T}
    N_F == 2 * cmap.n_fine   || error("N_F=$N_F ≠ 2×n_fine=$(2*cmap.n_fine)")
    N_C == 2 * cmap.n_coarse || error("N_C=$N_C ≠ 2×n_coarse=$(2*cmap.n_coarse)")

    sp_fine  = StateSpace{CartesianIndex{N_F}, T}()
    fine_buf = zeros(Int, N_F)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]; prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_pqsd_cmap_recurse!(sp_fine, coarse_tup, prob, fine_buf,
                                     1, cmap, qt, weight_tol, prob_tol)
    end
    sp_fine
end

function _prolong_pqsd_cmap_recurse!(
    sp_fine, coarse_tup::NTuple{N_C,Int}, weight::T,
    fine_buf::Vector{Int}, J::Int,
    cmap, qt, weight_tol, prob_tol
) where {N_C, T}
    if J > cmap.n_coarse
        N_F = length(fine_buf)
        s   = CartesianIndex(NTuple{N_F, Int}(fine_buf))
        idx = get(sp_fine.index, s, 0)
        if idx == 0; add_state!(sp_fine, s, weight)
        else; sp_fine.probs[idx] += weight; end
        return
    end

    nA_c  = coarse_tup[2J - 1]
    nB_c  = coarse_tup[2J]
    patch = cmap.coarse_to_fine[J]

    nA_c <= qt.nA_max && nB_c <= qt.nB_max ||
        error("Count (nA_c=$nA_c, nB_c=$nB_c) exceeds PatchQSD table bounds")

    for (nA_vec, nB_vec, w) in qt.tables[nA_c + 1, nB_c + 1]
        w > prob_tol || continue
        new_wt = weight * w
        new_wt > weight_tol || continue
        # Place fine-voxel counts into fine_buf at global positions
        for (m, k) in enumerate(patch)
            fine_buf[2k - 1] = nA_vec[m]
            fine_buf[2k]     = nB_vec[m]
        end
        _prolong_pqsd_cmap_recurse!(sp_fine, coarse_tup, new_wt, fine_buf,
                                     J + 1, cmap, qt, weight_tol, prob_tol)
    end
end

# ─── 1D pair version (no CoarseningMap) ──────────────────────────────────────

"""
    prolong_patch_qsd(sp_coarse, ::Val{N_F}, qt; weight_tol, prob_tol)
        → StateSpace{CartesianIndex{N_F}}

1D pair version (M=2). Replaces `prolong2s(sp_coarse, ::Val{N_F})`.
State layout: interleaved `(nA_1,nB_1, nA_2,nB_2, …)`.
"""
function prolong_patch_qsd(
    sp_coarse::StateSpace{CartesianIndex{N_C}, T},
    ::Val{N_F},
    qt::PatchQSD;
    weight_tol::Float64 = 1e-14,
    prob_tol::Float64   = 0.0
) where {N_C, N_F, T}
    qt.M == 2 || error("Use the CoarseningMap method for M=$(qt.M)≠2")
    N_F == 2 * N_C || error("N_F=$N_F ≠ 2×N_C=$N_C for M=2 prolongation")
    K2 = N_C ÷ 2   # number of coarse voxels

    sp_fine  = StateSpace{CartesianIndex{N_F}, T}()
    fine_buf = zeros(Int, N_F)

    for i in eachindex(sp_coarse.states)
        prob = sp_coarse.probs[i]; prob > weight_tol || continue
        coarse_tup = Tuple(sp_coarse.states[i])
        _prolong_pqsd_1d_recurse!(sp_fine, coarse_tup, prob, fine_buf,
                                   1, K2, qt, weight_tol, prob_tol)
    end
    sp_fine
end

function _prolong_pqsd_1d_recurse!(
    sp_fine, coarse_tup, weight, fine_buf,
    J, K2, qt, weight_tol, prob_tol
)
    if J > K2
        N_F = length(fine_buf)
        s   = CartesianIndex(NTuple{N_F, Int}(fine_buf))
        idx = get(sp_fine.index, s, 0)
        if idx == 0; add_state!(sp_fine, s, weight)
        else; sp_fine.probs[idx] += weight; end
        return
    end

    nA_c = coarse_tup[2J - 1]
    nB_c = coarse_tup[2J]
    nA_c <= qt.nA_max && nB_c <= qt.nB_max ||
        error("Count (nA_c=$nA_c,nB_c=$nB_c) out of PatchQSD table bounds")

    for (nA_vec, nB_vec, w) in qt.tables[nA_c + 1, nB_c + 1]
        w > prob_tol || continue
        new_wt = weight * w
        new_wt > weight_tol || continue
        # Pair J covers fine voxels (2J-1, 2J); interleaved layout
        fine_buf[4J - 3] = nA_vec[1];  fine_buf[4J - 2] = nB_vec[1]  # voxel 2J-1
        fine_buf[4J - 1] = nA_vec[2];  fine_buf[4J]     = nB_vec[2]  # voxel 2J
        _prolong_pqsd_1d_recurse!(sp_fine, coarse_tup, new_wt, fine_buf,
                                   J + 1, K2, qt, weight_tol, prob_tol)
    end
end

# ─── Diagnostics ─────────────────────────────────────────────────────────────

"""
    patch_qsd_correlation(qt, nA_c, nB_c) → ρ

Compute the Pearson correlation ρ(nA₁, nB₁) under the QSD for patch total
(nA_c, nB_c). Returns 0 for the Multinomial limit (k_rxn=0) and negative
values for strong A+B→∅ (molecules of opposite species repelled from same voxel).
"""
function patch_qsd_correlation(qt::PatchQSD, nA_c::Int, nB_c::Int)
    nA_c <= qt.nA_max && nB_c <= qt.nB_max || return NaN
    tbl = qt.tables[nA_c + 1, nB_c + 1]
    isempty(tbl) && return NaN

    μA = 0.0; μB = 0.0; μAB = 0.0; μA2 = 0.0; μB2 = 0.0
    for (nA_vec, nB_vec, w) in tbl
        a = Float64(nA_vec[1]); b = Float64(nB_vec[1])
        μA  += w * a;    μB  += w * b
        μAB += w * a * b; μA2 += w * a^2; μB2 += w * b^2
    end
    σA = sqrt(max(0.0, μA2 - μA^2));  σB = sqrt(max(0.0, μB2 - μB^2))
    (σA < 1e-12 || σB < 1e-12) ? 0.0 : (μAB - μA*μB) / (σA * σB)
end
