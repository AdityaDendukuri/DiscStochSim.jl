"""
Spatially Adaptive FSP for the 2D birth-death RDME — Galerkin coarse operator.

Three voxel states:
  EMPTY       – P(n=0)=1 implicit; no tracking
  ACTIVE      – full CME distribution; updated either by a local 2x2 V-cycle
                on disjoint fully-active blocks or by the original single-voxel
                expv fallback
  EQUILIBRATED – Galerkin coarse: full distribution, updated with expv at krylov_m_eq
                 (cheap: near-Poisson so small Krylov subspace suffices)

The key upgrades over the mean-field version:
  OLD: equil voxel stores μ_loc = (k_b + d·Σμ_nb)/(k_d + d·|nb|)  [algebraic, frozen at local SS]
  NEW: equil voxel stores p_eq, evolved each step under Q_eff         [Galerkin CME, dynamic]

  OLD: active voxels evolve independently under mean-field-coupled 1-voxel CMEs
  NEW: the full active region can be evolved as one persistent joint FSP state;
       when disabled or too large, a persistent 2x2 block V-cycle fallback remains available

The effective generator for an equilibrated voxel is the Galerkin coarse operator
obtained by marginalising the RDME over all other voxels under the mean-field
(independent-voxel) approximation:

    Q_eff[n→n+1] = k_b + d·Σ_{nb active/equil} μ_nb   (constant gain)
    Q_eff[n→n-1] = (k_d + d·|nb|)·n                    (linear loss)

This is identical to what build_coarse_system / build_rdme_system compute for a
single voxel with effective rates — the Galerkin coarsening of the RDME block.

Benefits over algebraic mean update:
  1. Correct transient: p_eq tracks the actual distribution, not the instantaneous SS.
  2. Accurate coupling: neighbours see voxel_mean(p_eq), which changes smoothly.
  3. Exact total distribution: convolution of actual p_eq, not assumed Poisson(μ_loc).
  4. No staircase artefact in equilibrated means: smooth convergence to global SS.
"""

# ─── Model ───────────────────────────────────────────────────────────────────

abstract type AbstractSpatialRDMEModel end

struct BirthDeathRDME <: AbstractSpatialRDMEModel
    k_b :: Float64
    k_d :: Float64
    D   :: Float64
    dx  :: Float64
end

struct BranchingDeathRDME <: AbstractSpatialRDMEModel
    k_b :: Float64
    k_d :: Float64
    D   :: Float64
    dx  :: Float64
end

jump_rate(m::BirthDeathRDME) = m.D / m.dx^2
ss_mean(m::BirthDeathRDME)   = m.k_b / m.k_d
jump_rate(m::BranchingDeathRDME) = m.D / m.dx^2

has_finite_local_ss(::BirthDeathRDME) = true
has_finite_local_ss(::BranchingDeathRDME) = false

ss_mean(m::BranchingDeathRDME) =
    error("BranchingDeathRDME (X -> 2X, X -> ∅) has no finite nonzero steady mean")

local_birth_rate(m::BirthDeathRDME, n::Int) = m.k_b
local_birth_rate(m::BranchingDeathRDME, n::Int) = m.k_b * n

local_total_birth_rate(m::BirthDeathRDME, n::Int, n_sites::Int=1) = m.k_b * n_sites
local_total_birth_rate(m::BranchingDeathRDME, n::Int, n_sites::Int=1) = m.k_b * n

# ─── State ───────────────────────────────────────────────────────────────────

mutable struct SpatialFSP
    model        :: AbstractSpatialRDMEModel
    K_x          :: Int
    K_y          :: Int
    dists        :: Dict{CartesianIndex{2}, Vector{Float64}}   # ACTIVE: full CME dist
    joint_state  :: Union{Nothing, StateSpace}                  # persistent full joint active-region state
    joint_voxels :: Vector{CartesianIndex{2}}                   # active voxel ordering used by joint_state
    block_states :: Dict{CartesianIndex{2}, StateSpace{CartesianIndex{4}, Float64}}  # persistent 2×2 joint states
    equil_dists  :: Dict{CartesianIndex{2}, Vector{Float64}}   # EQUIL: Galerkin coarse dist
    pending_flux :: Dict{CartesianIndex{2}, Float64}            # cumulative flux → empty
    t            :: Float64
    n_max        :: Int     # truncation for active voxels
    n_max_eq     :: Int     # truncation for equilibrated voxels (can be smaller)
    ε_expand     :: Float64
    ε_equil      :: Float64
    ε_prune      :: Float64
    krylov_m_eq  :: Int     # Krylov dim for equil voxels (small: near-SS distributions)
    use_joint_active :: Bool
    joint_state_tol  :: Float64
    joint_max_states :: Int
    joint_expand_depth :: Int
    joint_expand_threshold :: Float64
    use_block_vcycle :: Bool
    block_state_tol  :: Float64
    block_max_states :: Int
    block_expand_depth :: Int
    block_pre_frac   :: Float64
    block_post_frac  :: Float64
    use_joint_vcycle    :: Bool     # V-cycle on the joint active region
    joint_vcycle_τ_pre  :: Float64  # pre-smooth time (intra-pair diffusion)
    joint_vcycle_τ_post :: Float64  # post-smooth time
end

function SpatialFSP(model::AbstractSpatialRDMEModel, K_x::Int, K_y::Int;
                    n_max::Int        = 40,
                    n_max_eq::Int     = -1,      # default: auto (4·μ_ss + 4)
                    ε_expand::Float64 = 0.3,
                    ε_equil::Float64  = 0.06,
                    ε_prune::Float64  = 1e-14,
                    krylov_m_eq::Int  = 8,        # equil voxels converge fast
                    use_joint_active::Bool = false,
                    joint_state_tol::Float64 = 1e-10,
                    joint_max_states::Int = 250_000,
                    joint_expand_depth::Int = 1,
                    joint_expand_threshold::Float64 = 0.05,
                    use_block_vcycle::Bool = true,
                    block_state_tol::Float64 = 1e-10,
                    block_max_states::Int = 50_000,
                    block_expand_depth::Int = 1,
                    block_pre_frac::Float64 = 0.35,
                    block_post_frac::Float64 = 0.15,
                    use_joint_vcycle::Bool = false,
                    joint_vcycle_τ_pre::Float64 = 0.0,
                    joint_vcycle_τ_post::Float64 = 0.0)
    n_max_eq_val = n_max_eq > 0 ? n_max_eq :
                   has_finite_local_ss(model) ? max(12, round(Int, 4 * ss_mean(model)) + 4) :
                   max(12, n_max)
    SpatialFSP(model, K_x, K_y,
               Dict{CartesianIndex{2}, Vector{Float64}}(),
               nothing,
               CartesianIndex{2}[],
               Dict{CartesianIndex{2}, StateSpace{CartesianIndex{4}, Float64}}(),
               Dict{CartesianIndex{2}, Vector{Float64}}(),
               Dict{CartesianIndex{2}, Float64}(),
               0.0, n_max, n_max_eq_val, ε_expand, ε_equil, ε_prune, krylov_m_eq,
               use_joint_active, joint_state_tol, joint_max_states, joint_expand_depth,
               joint_expand_threshold,
               use_block_vcycle, block_state_tol, block_max_states,
               block_expand_depth, block_pre_frac, block_post_frac,
               use_joint_vcycle, joint_vcycle_τ_pre, joint_vcycle_τ_post)
end

function set_ic!(s::SpatialFSP, idx::CartesianIndex{2}, n0::Int)
    n0 <= s.n_max || error("n0=$n0 > n_max=$(s.n_max)")
    p = zeros(s.n_max + 1);  p[n0 + 1] = 1.0
    s.dists[idx] = p
    s.joint_state = nothing
    empty!(s.joint_voxels)
    empty!(s.block_states)
end

n_active(s::SpatialFSP)       = length(s.dists)
n_equilibrated(s::SpatialFSP) = length(s.equil_dists)
n_empty(s::SpatialFSP)        = s.K_x * s.K_y - n_active(s) - n_equilibrated(s)

# ─── Helpers ─────────────────────────────────────────────────────────────────

voxel_mean(p::AbstractVector) = sum((i - 1) * p[i] for i in eachindex(p))

function grid_neighbors(idx::CartesianIndex{2}, K_x::Int, K_y::Int)
    i, j = Tuple(idx)
    ns = CartesianIndex{2}[]
    i > 1   && push!(ns, CartesianIndex(i - 1, j))
    i < K_x && push!(ns, CartesianIndex(i + 1, j))
    j > 1   && push!(ns, CartesianIndex(i, j - 1))
    j < K_y && push!(ns, CartesianIndex(i, j + 1))
    ns
end

function _single_voxel_generator(model::AbstractSpatialRDMEModel,
                                 in_rate::Float64,
                                 out_coeff::Float64,
                                 n_max::Int;
                                 n_sites::Int = 1)
    n = n_max + 1
    Q = zeros(n, n)
    for col in 1:n
        nm = col - 1
        up = in_rate + local_total_birth_rate(model, nm, n_sites)
        if col < n && up > 0
            Q[col+1, col] += up
            Q[col, col] -= up
        end
        if col > 1
            down = (model.k_d + out_coeff) * nm
            if down > 0
                Q[col-1, col] += down
                Q[col, col] -= down
            end
        end
    end
    Q
end

function _poisson_pmf(λ::Float64, n_max::Int)
    p = zeros(n_max + 1);  p[1] = exp(-λ)
    for n in 1:n_max; p[n+1] = p[n] * λ / n; end
    p
end

poisson_pmf_vec(λ::Float64, n_max::Int) = _poisson_pmf(λ, n_max)

# Normalise and prune in-place
function _clean!(p::Vector{Float64}, ε_prune::Float64)
    p .= max.(0.0, p);  s = sum(p);  s > 0.0 && (p ./= s)
    p[p .< ε_prune] .= 0.0
end

# Truncate/pad a distribution vector to length n_max+1
function _resize_dist(p::Vector{Float64}, n_max::Int)
    n = n_max + 1
    length(p) == n && return copy(p)
    q = zeros(n)
    m = min(length(p), n)
    q[1:m] .= p[1:m]
    s = sum(q);  s > 0.0 && (q ./= s)
    q
end

function _copy_sp_local(sp::StateSpace{E, T}) where {E, T}
    sp2 = StateSpace{E, T}()
    for i in eachindex(sp.states)
        add_state!(sp2, sp.states[i], sp.probs[i])
    end
    sp2
end

_eu_local(i::Int, K::Int) = CartesianIndex(ntuple(j -> j == i ? 1 : 0, K))

_active_voxels(s::SpatialFSP) = sort!(collect(keys(s.dists)); by = Tuple)

function _positive_support(p::Vector{Float64}, tol::Float64)
    supp = findall(p .> tol)
    isempty(supp) ? [argmax(p)] : supp
end

function _product_state_from_dists(dists::Vector{Vector{Float64}},
                                   tol::Float64,
                                   max_states::Int)
    K = length(dists)
    K == 0 && return nothing
    _product_state_from_dists(dists, tol, max_states, Val(K))
end

function _product_state_from_dists(dists::Vector{Vector{Float64}},
                                   tol::Float64,
                                   max_states::Int,
                                   ::Val{K}) where {K}
    supports = [_positive_support(p, tol) for p in dists]
    # Overflow-safe product: use Int128 then compare
    support_prod = Int128(1)
    for supp in supports
        support_prod *= length(supp)
        support_prod > max_states && return nothing
    end

    sp = StateSpace{CartesianIndex{K}, Float64}()
    buf = zeros(Int, K)

    function rec(dim::Int, weight::Float64)
        if dim > K
            add_state!(sp, CartesianIndex(ntuple(i -> buf[i], K)), weight)
            return
        end
        p = dists[dim]
        for idx in supports[dim]
            w = weight * p[idx]
            w > tol || continue
            buf[dim] = idx - 1
            rec(dim + 1, w)
        end
    end

    rec(1, 1.0)
    renormalize!(sp) > 0.0 ? sp : nothing
end

function _marginalize_joint_state(sp_old::StateSpace, keep_positions::Vector{Int})
    K = length(keep_positions)
    K == 0 && return nothing
    _marginalize_joint_state(sp_old, keep_positions, Val(K))
end

function _marginalize_joint_state(sp_old::StateSpace{CartesianIndex{N0}, Float64},
                                  keep_positions::Vector{Int},
                                  ::Val{K}) where {N0, K}
    sp_new = StateSpace{CartesianIndex{K}, Float64}()
    for i in eachindex(sp_old.states)
        p = sp_old.probs[i]
        p > 0.0 || continue
        tup = Tuple(sp_old.states[i])
        s_new = CartesianIndex(ntuple(j -> tup[keep_positions[j]], K))
        _add_block_state!(sp_new, s_new, p)
    end
    renormalize!(sp_new) > 0.0 ? sp_new : nothing
end

function _append_joint_state(sp_base::StateSpace{CartesianIndex{N0}, Float64},
                             add_dists::Vector{Vector{Float64}},
                             tol::Float64,
                             max_states::Int) where {N0}
    n_add = length(add_dists)
    n_add == 0 && return _copy_sp_local(sp_base)
    K = N0 + n_add
    _append_joint_state(sp_base, add_dists, tol, max_states, Val(K))
end

function _append_joint_state(sp_base::StateSpace{CartesianIndex{N0}, Float64},
                             add_dists::Vector{Vector{Float64}},
                             tol::Float64,
                             max_states::Int,
                             ::Val{K}) where {N0, K}
    supports = [_positive_support(p, tol) for p in add_dists]
    support_prod = prod(length(supp) for supp in supports)
    length(sp_base) * support_prod <= max_states || return nothing

    sp_new = StateSpace{CartesianIndex{K}, Float64}()
    buf = zeros(Int, K)
    n_add = length(add_dists)

    function rec_add(dim::Int, weight::Float64)
        if dim > n_add
            add_state!(sp_new, CartesianIndex(ntuple(i -> buf[i], K)), weight)
            return
        end
        p = add_dists[dim]
        for idx in supports[dim]
            w = weight * p[idx]
            w > tol || continue
            buf[N0 + dim] = idx - 1
            rec_add(dim + 1, w)
        end
    end

    for i in eachindex(sp_base.states)
        p = sp_base.probs[i]
        p > 0.0 || continue
        tup = Tuple(sp_base.states[i])
        for j in 1:N0
            buf[j] = tup[j]
        end
        rec_add(1, p)
    end

    renormalize!(sp_new) > 0.0 ? sp_new : nothing
end

function _permute_joint_state(sp_old::StateSpace{CartesianIndex{N}, Float64},
                              perm::Vector{Int}) where {N}
    length(perm) == N || error("Permutation length $(length(perm)) != state dimension $N")
    sp_new = StateSpace{CartesianIndex{N}, Float64}()
    for i in eachindex(sp_old.states)
        p = sp_old.probs[i]
        p > 0.0 || continue
        tup = Tuple(sp_old.states[i])
        s_new = CartesianIndex(ntuple(j -> tup[perm[j]], N))
        _add_block_state!(sp_new, s_new, p)
    end
    renormalize!(sp_new) > 0.0 ? sp_new : nothing
end

function _joint_marginals(sp::StateSpace{CartesianIndex{K}, Float64},
                          n_max::Int) where {K}
    out = [zeros(n_max + 1) for _ in 1:K]
    for i in eachindex(sp.states)
        p = sp.probs[i]
        p > 0.0 || continue
        tup = Tuple(sp.states[i])
        for k in 1:K
            out[k][tup[k] + 1] += p
        end
    end
    for p in out
        _clean!(p, 0.0)
    end
    out
end

function _joint_total_distribution(sp::StateSpace)
    isempty(sp.states) && return [1.0]
    max_total = maximum(sum(Tuple(s)) for s in sp.states)
    out = zeros(max_total + 1)
    for i in eachindex(sp.states)
        out[sum(Tuple(sp.states[i])) + 1] += sp.probs[i]
    end
    out
end

function _block_voxels(anchor::CartesianIndex{2})
    i, j = Tuple(anchor)
    (
        CartesianIndex(i,     j),
        CartesianIndex(i,     j + 1),
        CartesianIndex(i + 1, j),
        CartesianIndex(i + 1, j + 1),
    )
end

function _try_add_block_anchor!(
    anchors::Vector{CartesianIndex{2}},
    covered::Set{CartesianIndex{2}},
    dists::Dict{CartesianIndex{2}, Vector{Float64}},
    anchor::CartesianIndex{2},
    K_x::Int,
    K_y::Int,
)
    i, j = Tuple(anchor)
    1 <= i < K_x || return
    1 <= j < K_y || return

    block = _block_voxels(anchor)
    all(haskey(dists, idx) for idx in block) || return
    any(idx in covered for idx in block) && return

    push!(anchors, anchor)
    foreach(idx -> push!(covered, idx), block)
end

function _select_active_blocks(
    dists::Dict{CartesianIndex{2}, Vector{Float64}},
    K_x::Int,
    K_y::Int,
    preferred::AbstractVector{CartesianIndex{2}} = CartesianIndex{2}[],
)
    anchors = CartesianIndex{2}[]
    covered = Set{CartesianIndex{2}}()

    for anchor in sort!(collect(preferred); by = Tuple)
        _try_add_block_anchor!(anchors, covered, dists, anchor, K_x, K_y)
    end
    for i in 1:(K_x - 1), j in 1:(K_y - 1)
        _try_add_block_anchor!(anchors, covered, dists, CartesianIndex(i, j), K_x, K_y)
    end

    anchors
end

function _build_joint_active_system(model::AbstractSpatialRDMEModel,
                                    voxels::Vector{CartesianIndex{2}},
                                    means::Dict{CartesianIndex{2}, Float64},
                                    K_x::Int,
                                    K_y::Int)
    K = length(voxels)
    K > 0 || error("Active region must be non-empty")
    d = jump_rate(model)
    stoichs      = CartesianIndex{K}[]
    propensities = Function[]
    pos = Dict(voxels[i] => i for i in eachindex(voxels))

    for p in 1:K
        let p = p
            ep = _eu_local(p, K)
            push!(stoichs,  ep); push!(propensities, (x, rv, t) -> local_birth_rate(model, Tuple(x)[p]))
            push!(stoichs, -ep); push!(propensities, (x, rv, t) -> model.k_d * max(0, Tuple(x)[p]))
        end
    end

    for (p, idx) in enumerate(voxels)
        for nb in grid_neighbors(idx, K_x, K_y)
            q = get(pos, nb, 0)
            if q != 0
                p < q || continue
                let p = p, q = q
                    ep = _eu_local(p, K)
                    eq = _eu_local(q, K)
                    push!(stoichs, -ep + eq); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[p]))
                    push!(stoichs,  ep - eq); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[q]))
                end
            else
                μ_nb = get(means, nb, 0.0)
                let p = p, μ_nb = μ_nb
                    ep = _eu_local(p, K)
                    push!(stoichs, -ep)
                    push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[p]))
                    if μ_nb > 0.0
                        push!(stoichs, ep)
                        push!(propensities, (x, rv, t) -> d * μ_nb)
                    end
                end
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

function _joint_active_bc(state, n_max::Int)
    all(0 <= n <= n_max for n in Tuple(state))
end

# ─── Joint active V-cycle ─────────────────────────────────────────────────────

"""
Intra-pair diffusion system: diffusion only within consecutive pairs
(active_voxels[1],active_voxels[2]), (active_voxels[3],active_voxels[4]), ...
Used as the pre/post-smooth operator. Pairs that are not spatially adjacent on
the 2D grid are silently skipped (no smoothing needed if fast mixing doesn't apply).
"""
function _build_joint_intra_pair_system(
    model::AbstractSpatialRDMEModel,
    active_voxels::Vector{CartesianIndex{2}},
    K_x::Int, K_y::Int,
    ::Val{K_FINE}
) where K_FINE
    d = jump_rate(model)
    stoichs      = CartesianIndex{K_FINE}[]
    propensities = Function[]
    for j in 1:(K_FINE ÷ 2)
        v1 = active_voxels[2j-1]; v2 = active_voxels[2j]
        any(==(v2), grid_neighbors(v1, K_x, K_y)) || continue
        let p = 2j-1, q = 2j
            ep = _eu_local(p, K_FINE); eq = _eu_local(q, K_FINE)
            push!(stoichs, -ep + eq); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[p]))
            push!(stoichs,  ep - eq); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[q]))
        end
    end
    DiscreteStochasticSystem{CartesianIndex{K_FINE}}(stoichs, propensities)
end

"""
Coarse system for K_COARSE super-voxels, each covering 2 consecutive fine voxels.
Rates are derived from the fine-grid 2D geometry assuming fast intra-pair mixing
(each molecule is at each of the 2 fine positions with probability 1/2).

  Birth:  2·local_birth_rate per super-voxel  (2 fine voxels)
  Death:  k_d per molecule
  Diffusion j→k:  d·(edge_count_jk / 2)·N_j  (averaged over uniform pair distribution)
  Outflow to inactive:  d·(ext_edges_j / 2)·N_j
  Inflow from inactive: d·Σ μ_nb  (summed over all inactive neighbours of either fine voxel)
"""
function _build_joint_pair_coarse_system(
    model::AbstractSpatialRDMEModel,
    active_voxels::Vector{CartesianIndex{2}},
    means::Dict{CartesianIndex{2}, Float64},
    K_x::Int, K_y::Int,
    ::Val{K_COARSE}
) where K_COARSE
    K_FINE = 2 * K_COARSE
    length(active_voxels) == K_FINE || error("Need $(K_FINE) fine voxels, got $(length(active_voxels))")
    d = jump_rate(model)
    stoichs      = CartesianIndex{K_COARSE}[]
    propensities = Function[]

    active_set      = Set(active_voxels)
    fine_to_coarse  = Dict(active_voxels[k] => div(k-1, 2) + 1 for k in 1:K_FINE)

    # ── Reactions ──────────────────────────────────────────────────────────────
    for j in 1:K_COARSE
        let j = j
            ej = _eu_local(j, K_COARSE)
            push!(stoichs,  ej)
            push!(propensities, (x, rv, t) -> 2 * local_birth_rate(model, Tuple(x)[j]))
            push!(stoichs, -ej)
            push!(propensities, (x, rv, t) -> model.k_d * max(0, Tuple(x)[j]))
        end
    end

    # ── Inter-pair diffusion (count boundary edges between pairs) ─────────────
    pair_edges = Dict{Tuple{Int,Int}, Int}()
    for k_fine in 1:K_FINE
        v = active_voxels[k_fine]; j = div(k_fine - 1, 2) + 1
        for nb in grid_neighbors(v, K_x, K_y)
            nb ∈ active_set || continue
            l = fine_to_coarse[nb]
            l == j && continue
            key = (min(j, l), max(j, l))
            pair_edges[key] = get(pair_edges, key, 0) + 1
        end
    end
    for ((j, k), n_edges) in pair_edges
        let j = j, k = k, r = d * n_edges / 2   # 1/pair_size = 1/2
            ej = _eu_local(j, K_COARSE); ek = _eu_local(k, K_COARSE)
            push!(stoichs, -ej + ek); push!(propensities, (x, rv, t) -> r * max(0, Tuple(x)[j]))
            push!(stoichs,  ej - ek); push!(propensities, (x, rv, t) -> r * max(0, Tuple(x)[k]))
        end
    end

    # ── Boundary with inactive neighbours ─────────────────────────────────────
    out_rate = zeros(K_COARSE)
    in_rate  = zeros(K_COARSE)
    for k_fine in 1:K_FINE
        v = active_voxels[k_fine]; j = div(k_fine - 1, 2) + 1
        for nb in grid_neighbors(v, K_x, K_y)
            nb ∈ active_set && continue
            out_rate[j] += d / 2                        # per-molecule outflow
            in_rate[j]  += d * get(means, nb, 0.0)     # inflow (not per-molecule)
        end
    end
    for j in 1:K_COARSE
        let j = j, r_out = out_rate[j], r_in = in_rate[j]
            ej = _eu_local(j, K_COARSE)
            r_out > 0 && (push!(stoichs, -ej); push!(propensities, (x, rv, t) -> r_out * max(0, Tuple(x)[j])))
            r_in  > 0 && (push!(stoichs,  ej); push!(propensities, (x, rv, t) -> r_in))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_COARSE}}(stoichs, propensities)
end

"""
    _binomial_split_marginal(p_N, n_max_out, tol) -> Vector{Float64}

Given the coarse total-count distribution P(N=n), compute the marginal
P(n_1 = k) assuming n_1 ~ Binomial(N, 1/2). By symmetry P(n_2) is identical.
Valid when fast intra-pair diffusion has equilibrated the within-pair positions.
"""
function _binomial_split_marginal(p_N::Vector{Float64}, n_max_out::Int, tol::Float64)
    p_out = zeros(n_max_out + 1)
    for N in 0:length(p_N)-1
        p_N[N+1] > tol || continue
        log_half_N = -N * log(2.0)
        for k in 0:min(N, n_max_out)
            log_bc = (sum(log(Float64(i)) for i in (N-k+1):N; init=0.0) -
                      sum(log(Float64(i)) for i in 1:k;      init=0.0))
            w = p_N[N+1] * exp(log_half_N + log_bc)
            w > tol && (p_out[k+1] += w)
        end
    end
    s = sum(p_out); s > 0.0 && (p_out ./= s)
    p_out
end

"""
    _evolve_joint_active_vcycle!(s, means, dt; krylov_m)

Coarse-first marginal V-cycle — no fine CartesianIndex state is ever formed.

  1. Coarse marginals: convolve consecutive pairs → K_act/2 1-D distributions
  2. Coarse product state: CartesianIndex{K_act/2}, ≈(2·n_max)^(K_act/2) states
  3. Coarse solve: expand! + expv under the pair coarse generator
  4. Extract coarse 1-D marginals from coarse joint state
  5. Binomial-split each coarse marginal → two per-voxel marginals P(n_1), P(n_2)
     [Optional: 2-voxel post-smooth via intra-pair CME when τ_post > 0]
  6. Write per-voxel marginals back to s.dists

State space at coarse level: (2·n_max+1)^(K_act/2) — tractable up to K_act ≈ 16
with n_max = 10. Falls back to _evolve_joint_active! for odd K_act or < 2.
"""
function _evolve_joint_active_vcycle!(s::SpatialFSP,
                                       means::Dict{CartesianIndex{2}, Float64},
                                       dt::Float64;
                                       krylov_m::Int)
    active_voxels = _active_voxels(s)
    K_act = length(active_voxels)

    # Fallback: odd or too-small → direct joint expv (cheap for K_act ≤ 3)
    (K_act >= 2 && iseven(K_act)) || return _evolve_joint_active!(s, means, dt; krylov_m)

    K_coarse = K_act ÷ 2
    n_max_c  = 2 * s.n_max

    # ── 1. Coarse marginals: convolve pairs ───────────────────────────────────
    coarse_marginals = Vector{Float64}[]
    for j in 1:K_coarse
        p1 = s.dists[active_voxels[2j-1]]
        p2 = s.dists[active_voxels[2j]]
        p_pair = zeros(length(p1) + length(p2) - 1)
        for i1 in eachindex(p1), i2 in eachindex(p2)
            p_pair[i1 + i2 - 1] += p1[i1] * p2[i2]
        end
        push!(coarse_marginals, p_pair)
    end

    # ── 2. Coarse product state ───────────────────────────────────────────────
    sp_coarse = _product_state_from_dists(coarse_marginals, s.joint_state_tol, s.joint_max_states)
    sp_coarse === nothing && return _evolve_joint_active!(s, means, dt; krylov_m)

    # ── 3. Coarse solve ───────────────────────────────────────────────────────
    sys_coarse = _build_joint_pair_coarse_system(s.model, active_voxels, means,
                                                  s.K_x, s.K_y, Val(K_coarse))
    expand!(sp_coarse, sys_coarse, x -> _joint_active_bc(x, n_max_c); depth = 1)
    length(sp_coarse) <= s.joint_max_states || return _evolve_joint_active!(s, means, dt; krylov_m)
    A_c, = build_generator(sp_coarse, sys_coarse, nothing, s.t)
    sp_coarse.probs .= expv(Float64(dt), A_c, sp_coarse.probs;
                             m = min(krylov_m, size(A_c, 1)))
    prune_threshold!(sp_coarse, s.joint_state_tol)
    renormalize!(sp_coarse) > 0.0 || return false

    # ── 4. Coarse 1-D marginals from joint coarse state ───────────────────────
    coarse_marg_post = _joint_marginals(sp_coarse, n_max_c)

    # ── 5. Split each coarse marginal → two per-voxel distributions ──────────
    τ_post = s.joint_vcycle_τ_post
    d = jump_rate(s.model)

    for j in 1:K_coarse
        v1 = active_voxels[2j-1]; v2 = active_voxels[2j]
        p_N = coarse_marg_post[j]

        if τ_post > 0.0 && any(==(v2), grid_neighbors(v1, s.K_x, s.K_y))
            # 2-voxel post-smooth: build CartesianIndex{2} joint, run intra-pair CME
            p1_init = _binomial_split_marginal(p_N, s.n_max, s.joint_state_tol)
            p2_init = copy(p1_init)
            # Build 2-voxel joint state from product p1 ⊗ p2
            sp2 = StateSpace{CartesianIndex{2}, Float64}()
            for i1 in eachindex(p1_init), i2 in eachindex(p2_init)
                w = p1_init[i1] * p2_init[i2]
                w > s.joint_state_tol || continue
                add_state!(sp2, CartesianIndex(i1-1, i2-1), w)
            end
            if !isempty(sp2.states)
                e1 = CartesianIndex(1,0); e2 = CartesianIndex(0,1)
                stoichs2 = [e2 - e1, e1 - e2]
                props2   = Function[
                    (x, rv, t) -> d * max(0, Tuple(x)[1]),
                    (x, rv, t) -> d * max(0, Tuple(x)[2]),
                ]
                intra2 = DiscreteStochasticSystem{CartesianIndex{2}}(stoichs2, props2)
                A2, = build_generator(sp2, intra2, nothing, s.t)
                sp2.probs .= expv(Float64(τ_post), A2, sp2.probs; m = min(krylov_m, size(A2,1)))
                _clean!(sp2.probs, s.ε_prune)
                p1 = zeros(s.n_max + 1); p2 = zeros(s.n_max + 1)
                for idx in eachindex(sp2.states)
                    tup = Tuple(sp2.states[idx]); pr = sp2.probs[idx]
                    tup[1] <= s.n_max && (p1[tup[1]+1] += pr)
                    tup[2] <= s.n_max && (p2[tup[2]+1] += pr)
                end
                _clean!(p1, 0.0); _clean!(p2, 0.0)
                s.dists[v1] = p1; s.dists[v2] = p2
                continue
            end
        end

        # Default: symmetric Binomial split
        p_split = _binomial_split_marginal(p_N, s.n_max, s.joint_state_tol)
        s.dists[v1] = copy(p_split)
        s.dists[v2] = copy(p_split)
    end

    # Joint state not maintained (only marginals): clear to avoid stale cache
    s.joint_state  = nothing
    empty!(s.joint_voxels)
    true
end

function _sync_joint_state!(s::SpatialFSP, active_voxels::Vector{CartesianIndex{2}})
    isempty(active_voxels) && begin
        s.joint_state = nothing
        empty!(s.joint_voxels)
        return false
    end

    # Fast cap: if K_act is so large that even the smallest product state would exceed
    # max_states, bail immediately before computing supports.
    length(active_voxels) > 60 && begin   # 2^60 >> any max_states
        s.joint_state = nothing; empty!(s.joint_voxels); return false
    end

    if s.joint_state === nothing || isempty(s.joint_voxels)
        sp = _product_state_from_dists([s.dists[idx] for idx in active_voxels],
                                       s.joint_state_tol, s.joint_max_states)
        sp === nothing && return false
        s.joint_state = sp
        s.joint_voxels = copy(active_voxels)
        return true
    end

    s.joint_voxels == active_voxels && return true

    prev_voxels = copy(s.joint_voxels)
    prev_pos = Dict(prev_voxels[i] => i for i in eachindex(prev_voxels))
    keep_positions = [i for i in eachindex(prev_voxels) if prev_voxels[i] in active_voxels]
    work_voxels = prev_voxels[keep_positions]

    sp_work = isempty(keep_positions) ? nothing : _marginalize_joint_state(s.joint_state, keep_positions)
    if sp_work === nothing
        add_dists = [s.dists[idx] for idx in active_voxels]
        sp_new = _product_state_from_dists(add_dists, s.joint_state_tol, s.joint_max_states)
        sp_new === nothing && return false
        s.joint_state = sp_new
        s.joint_voxels = copy(active_voxels)
        return true
    end

    added_voxels = [idx for idx in active_voxels if !haskey(prev_pos, idx)]
    if !isempty(added_voxels)
        add_dists = [s.dists[idx] for idx in added_voxels]
        sp_work = sp_work === nothing ?
            _product_state_from_dists(add_dists, s.joint_state_tol, s.joint_max_states) :
            _append_joint_state(sp_work, add_dists, s.joint_state_tol, s.joint_max_states)
        sp_work === nothing && return false
        append!(work_voxels, added_voxels)
    end

    if work_voxels != active_voxels
        work_pos = Dict(work_voxels[i] => i for i in eachindex(work_voxels))
        perm = [work_pos[idx] for idx in active_voxels]
        sp_work = _permute_joint_state(sp_work, perm)
        sp_work === nothing && return false
    end

    s.joint_state = sp_work
    s.joint_voxels = copy(active_voxels)
    true
end

function _evolve_joint_active!(s::SpatialFSP,
                               means::Dict{CartesianIndex{2}, Float64},
                               dt::Float64;
                               krylov_m::Int)
    active_voxels = _active_voxels(s)
    _sync_joint_state!(s, active_voxels) || begin
        s.joint_state = nothing
        empty!(s.joint_voxels)
        return false
    end

    sp = _copy_sp_local(s.joint_state)
    sys = _build_joint_active_system(s.model, active_voxels, means, s.K_x, s.K_y)

    if s.joint_expand_depth > 0
        if s.joint_expand_threshold > 0.0
            expand_flux!(sp, sys, nothing, s.t, x -> _joint_active_bc(x, s.n_max);
                         depth = s.joint_expand_depth, threshold = s.joint_expand_threshold)
        else
            expand!(sp, sys, x -> _joint_active_bc(x, s.n_max); depth = s.joint_expand_depth)
        end
        length(sp) <= s.joint_max_states || begin
            s.joint_state = nothing
            empty!(s.joint_voxels)
            return false
        end
    end

    A, = build_generator(sp, sys, nothing, s.t)
    sp.probs .= expv(Float64(dt), A, sp.probs; m = min(krylov_m, size(A, 1)))
    prune_threshold!(sp, s.joint_state_tol)
    renormalize!(sp) > 0.0 || begin
        s.joint_state = nothing
        empty!(s.joint_voxels)
        return false
    end

    s.joint_state = sp
    s.joint_voxels = copy(active_voxels)

    marginals = _joint_marginals(sp, s.n_max)
    empty!(s.dists)
    for (idx, p) in zip(active_voxels, marginals)
        s.dists[idx] = p
    end
    true
end

function _build_block_system(model::AbstractSpatialRDMEModel,
                             block::NTuple{4, CartesianIndex{2}},
                             means::Dict{CartesianIndex{2}, Float64},
                             K_x::Int,
                             K_y::Int)
    d = jump_rate(model)
    stoichs      = CartesianIndex{4}[]
    propensities = Function[]
    block_set    = Set(block)

    for p in 1:4
        let p = p
            ep = _eu_local(p, 4)
            push!(stoichs,  ep); push!(propensities, (x, rv, t) -> local_birth_rate(model, Tuple(x)[p]))
            push!(stoichs, -ep); push!(propensities, (x, rv, t) -> model.k_d * max(0, Tuple(x)[p]))
        end
    end

    for (a, b) in ((1, 2), (1, 3), (2, 4), (3, 4))
        let a = a, b = b
            ea = _eu_local(a, 4)
            eb = _eu_local(b, 4)
            push!(stoichs, -ea + eb); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[a]))
            push!(stoichs,  ea - eb); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[b]))
        end
    end

    ext_mean_sum = 0.0
    n_boundary_edges = 0
    for (p, idx) in enumerate(block)
        for nb in grid_neighbors(idx, K_x, K_y)
            nb in block_set && continue
            n_boundary_edges += 1
            μ_nb = get(means, nb, 0.0)
            ext_mean_sum += μ_nb
            let p = p, μ_nb = μ_nb
                ep = _eu_local(p, 4)
                push!(stoichs, -ep)
                push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[p]))
                μ_nb > 0.0 || continue
                push!(stoichs, ep)
                push!(propensities, (x, rv, t) -> d * μ_nb)
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{4}}(stoichs, propensities),
    ext_mean_sum,
    n_boundary_edges
end

function _build_block_intra_system(model::AbstractSpatialRDMEModel)
    d = jump_rate(model)
    stoichs      = CartesianIndex{4}[]
    propensities = Function[]

    for (a, b) in ((1, 2), (1, 3), (2, 4), (3, 4))
        let a = a, b = b
            ea = _eu_local(a, 4)
            eb = _eu_local(b, 4)
            push!(stoichs, -ea + eb); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[a]))
            push!(stoichs,  ea - eb); push!(propensities, (x, rv, t) -> d * max(0, Tuple(x)[b]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{4}}(stoichs, propensities)
end

function _block_bc(state::CartesianIndex{4}, n_max::Int)
    all(0 <= n <= n_max for n in Tuple(state))
end

function _build_block_product_state(p_vox::NTuple{4, Vector{Float64}},
                                    weight_tol::Float64)
    supports = ntuple(k -> begin
        supp = findall(p_vox[k] .> weight_tol)
        isempty(supp) ? [argmax(p_vox[k])] : supp
    end, 4)

    sp = StateSpace{CartesianIndex{4}, Float64}()
    for i1 in supports[1], i2 in supports[2], i3 in supports[3], i4 in supports[4]
        w = p_vox[1][i1] * p_vox[2][i2] * p_vox[3][i3] * p_vox[4][i4]
        w > weight_tol || continue
        add_state!(sp, CartesianIndex(i1 - 1, i2 - 1, i3 - 1, i4 - 1), w)
    end
    renormalize!(sp)
    sp
end

function _init_block_state(s::SpatialFSP, block::NTuple{4, CartesianIndex{2}})
    p_vox = (s.dists[block[1]], s.dists[block[2]], s.dists[block[3]], s.dists[block[4]])
    support_sizes = ntuple(k -> count(>(s.block_state_tol), p_vox[k]), 4)
    support_prod  = prod(max(1, n) for n in support_sizes)
    support_prod <= s.block_max_states || return nothing

    sp = _build_block_product_state(p_vox, s.block_state_tol)
    length(sp) > 0 || return nothing
    sp
end

function _sync_block_cache!(s::SpatialFSP, anchors::Vector{CartesianIndex{2}})
    allowed = Set(anchors)
    for anchor in collect(keys(s.block_states))
        anchor in allowed || delete!(s.block_states, anchor)
    end

    for anchor in anchors
        haskey(s.block_states, anchor) && continue
        block = _block_voxels(anchor)
        sp = _init_block_state(s, block)
        sp === nothing || (s.block_states[anchor] = sp)
    end
end

function _restrict_block_total(sp::StateSpace{CartesianIndex{4}, Float64}, total_n_max::Int)
    p_total = zeros(total_n_max + 1)
    for i in eachindex(sp.states)
        total = sum(Tuple(sp.states[i]))
        p_total[total + 1] += sp.probs[i]
    end
    p_total
end

function _add_block_state!(sp::StateSpace{E, T},
                           state::E,
                           weight::Float64) where {E, T<:Real}
    weight > 0.0 || return
    idx = get(sp.index, state, 0)
    if idx == 0
        add_state!(sp, state, weight)
    else
        sp.probs[idx] += weight
    end
end

function _logfactorials(n::Int)
    out = zeros(n + 1)
    for k in 1:n
        out[k + 1] = out[k] + log(Float64(k))
    end
    out
end

function _add_multinomial_block_states!(sp::StateSpace{CartesianIndex{4}, Float64},
                                        N::Int,
                                        total_prob::Float64,
                                        n_max::Int,
                                        logfacts::Vector{Float64},
                                        prob_tol::Float64)
    log_quarter = log(4.0)
    for n1 in 0:min(n_max, N)
        rem1 = N - n1
        for n2 in 0:min(n_max, rem1)
            rem2 = rem1 - n2
            for n3 in 0:min(n_max, rem2)
                n4 = rem2 - n3
                n4 <= n_max || continue
                logw = logfacts[N + 1] - logfacts[n1 + 1] - logfacts[n2 + 1] -
                       logfacts[n3 + 1] - logfacts[n4 + 1] - N * log_quarter
                w = total_prob * exp(logw)
                w > prob_tol || continue
                _add_block_state!(sp, CartesianIndex(n1, n2, n3, n4), w)
            end
        end
    end
end

function _block_marginals(sp::StateSpace{CartesianIndex{4}, Float64}, n_max::Int)
    out = ntuple(_ -> zeros(n_max + 1), 4)
    for i in eachindex(sp.states)
        p = sp.probs[i]
        p > 0.0 || continue
        t = Tuple(sp.states[i])
        for k in 1:4
            out[k][t[k] + 1] += p
        end
    end
    for k in 1:4
        _clean!(out[k], 0.0)
    end
    out
end

function _evolve_block_vcycle(s::SpatialFSP,
                              block::NTuple{4, CartesianIndex{2}},
                              sp_seed::StateSpace{CartesianIndex{4}, Float64},
                              means::Dict{CartesianIndex{2}, Float64},
                              dt::Float64;
                              krylov_m::Int)
    length(sp_seed) > 0 || return nothing
    length(sp_seed) <= s.block_max_states || return nothing

    sp_fine = _copy_sp_local(sp_seed)
    renormalize!(sp_fine) > 0.0 || return nothing

    sys_full, ext_mean_sum, n_boundary_edges =
        _build_block_system(s.model, block, means, s.K_x, s.K_y)

    if s.block_expand_depth > 0
        expand!(sp_fine, sys_full, x -> _block_bc(x, s.n_max); depth = s.block_expand_depth)
        length(sp_fine) <= s.block_max_states || return nothing
    end

    sp_smoothed = _copy_sp_local(sp_fine)
    τ_pre = s.block_pre_frac > 0.0 ? s.block_pre_frac * dt : 0.0
    if τ_pre > 0.0
        intra_sys = _build_block_intra_system(s.model)
        A_intra, = build_generator(sp_smoothed, intra_sys, nothing, s.t)
        sp_smoothed.probs .= expv(Float64(τ_pre), A_intra, sp_smoothed.probs;
                                  m = min(krylov_m, size(A_intra, 1)))
        _clean!(sp_smoothed.probs, s.ε_prune)
    end

    total_n_max = 4 * s.n_max
    p_total_pre = _restrict_block_total(sp_smoothed, total_n_max)

    # Coarse correction acts on the total count in the 2x2 block. Fast internal
    # diffusion is handled by smoothing; births/deaths and exterior exchange are
    # collapsed to a single coarse birth-death process.
    in_rate = jump_rate(s.model) * ext_mean_sum
    out_coeff = jump_rate(s.model) * (n_boundary_edges / 4)
    Q_c = _single_voxel_generator(s.model, in_rate, out_coeff, total_n_max; n_sites = 4)
    p_total_post = expv(Float64(dt), Q_c, p_total_pre; m = min(krylov_m, size(Q_c, 1)))
    _clean!(p_total_post, s.ε_prune)

    sp_new = StateSpace{CartesianIndex{4}, Float64}()
    for i in eachindex(sp_smoothed.states)
        p = sp_smoothed.probs[i]
        p > 0.0 || continue
        total = sum(Tuple(sp_smoothed.states[i]))
        p_pre = p_total_pre[total + 1]
        p_post = p_total_post[total + 1]
        p_new = p_pre > 1e-300 ? p * (p_post / p_pre) : 0.0
        p_new > 0.0 && add_state!(sp_new, sp_smoothed.states[i], p_new)
    end

    logfacts = _logfactorials(total_n_max)
    for total in 0:total_n_max
        p_pre = p_total_pre[total + 1]
        p_post = p_total_post[total + 1]
        p_pre <= s.block_state_tol && p_post > s.block_state_tol || continue
        _add_multinomial_block_states!(sp_new, total, p_post, s.n_max, logfacts, s.block_state_tol)
    end

    if s.block_expand_depth > 0
        intra_sys = _build_block_intra_system(s.model)
        expand!(sp_new, intra_sys, x -> _block_bc(x, s.n_max); depth = 1)
        length(sp_new) <= s.block_max_states || return nothing
    end

    τ_post = s.block_post_frac > 0.0 ? s.block_post_frac * dt : 0.0
    if τ_post > 0.0
        intra_sys = _build_block_intra_system(s.model)
        A_intra, = build_generator(sp_new, intra_sys, nothing, s.t)
        sp_new.probs .= expv(Float64(τ_post), A_intra, sp_new.probs;
                             m = min(krylov_m, size(A_intra, 1)))
        _clean!(sp_new.probs, s.ε_prune)
    end

    prune_threshold!(sp_new, s.block_state_tol)
    renormalize!(sp_new) > 0.0 || return nothing
    sp_new
end

# ─── Main time step ───────────────────────────────────────────────────────────

"""
    step!(s, dt; krylov_m=30)

Galerkin-coarse step:
  1. Means from BOTH active (voxel_mean(p)) and equilibrated (voxel_mean(p_eq)).
  2. Evolve active voxels under either:
       - a persistent fully joint active-region FSP, or
       - a persistent local 2x2 block V-cycle, or
       - the original 1-voxel Galerkin CME fallback.
  3. Evolve equilibrated voxels under Q_eff_eq (Galerkin coarse CME, krylov_m_eq).
     This replaces the old algebraic μ_loc update.
  4. Accumulate flux into empty voxels from active + equilibrated.
  5. Activate (prolong) empty voxels crossing threshold.
  6. Deactivate (restrict) active voxels near local Poisson SS;
     store their CURRENT distribution as the initial equil_dist.
"""
function step!(s::SpatialFSP, dt::Float64; krylov_m::Int = 30)
    d = jump_rate(s.model)

    # ── 1. Effective means: active → voxel_mean(p),  equil → voxel_mean(p_eq) ─
    means = Dict{CartesianIndex{2}, Float64}()
    for (idx, p)  in s.dists;       means[idx] = voxel_mean(p);    end
    for (idx, p)  in s.equil_dists; means[idx] = voxel_mean(p);    end

    # ── 2. Evolve active voxels (joint active-region FSP, then block V-cycle, then 1-voxel fallback) ─
    new_dists = Dict{CartesianIndex{2}, Vector{Float64}}()
    handled = Set{CartesianIndex{2}}()

    if s.use_joint_active && !isempty(s.dists)
        succeeded = s.use_joint_vcycle ?
                    _evolve_joint_active_vcycle!(s, means, dt; krylov_m) :
                    _evolve_joint_active!(s, means, dt; krylov_m)
        if succeeded
            new_dists = copy(s.dists)
            union!(handled, keys(new_dists))
            empty!(s.block_states)
        end
    end

    if isempty(handled) && s.use_block_vcycle
        active_blocks = _select_active_blocks(s.dists, s.K_x, s.K_y, collect(keys(s.block_states)))
        _sync_block_cache!(s, active_blocks)

        for anchor in active_blocks
            haskey(s.block_states, anchor) || continue
            block = _block_voxels(anchor)
            sp_new_block = _evolve_block_vcycle(s, block, s.block_states[anchor], means, dt; krylov_m)
            if sp_new_block === nothing
                delete!(s.block_states, anchor)
                continue
            end

            s.block_states[anchor] = sp_new_block
            p_new_block = _block_marginals(sp_new_block, s.n_max)
            for k in 1:4
                new_dists[block[k]] = p_new_block[k]
                push!(handled, block[k])
            end
        end
    end

    for (idx, p_old) in s.dists
        idx in handled && continue
        ns       = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in     = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        in_rate  = d * μ_in
        out_coeff = d * length(ns)
        Q        = _single_voxel_generator(s.model, in_rate, out_coeff, s.n_max)
        p_new    = expv(Float64(dt), Q, p_old; m = min(krylov_m, size(Q, 1)))
        _clean!(p_new, s.ε_prune)
        new_dists[idx] = p_new
    end
    s.dists = new_dists

    # ── 3. Update equilibrated voxels ────────────────────────────────────────
    #
    # Two sub-cases based on neighbourhood:
    #
    #  BOUNDARY equil (≥1 active neighbour):
    #    Effective rates are rapidly changing (active neighbour mean is moving).
    #    Use Galerkin expv — cheap because near-Poisson → tiny Krylov dim.
    #
    #  INTERIOR equil (all neighbours equil or empty):
    #    Effective rates change only as slow background equil means drift.
    #    Galerkin expv is wasteful: distribution barely moves in one step.
    #    Use algebraic mean → Poisson update instead (O(n_max_eq), no matrix exp).
    #
    # Cost breakdown per step:
    #   Boundary equil: O(wavefront circumference) ≈ O(√K)  × expv(n_max_eq, krylov_m_eq)
    #   Interior equil: O(K_equil)                           × O(n_max_eq)  algebraic
    #
    if has_finite_local_ss(s.model)
        new_equil = Dict{CartesianIndex{2}, Vector{Float64}}()
        for (idx, p_eq) in s.equil_dists
            ns      = grid_neighbors(idx, s.K_x, s.K_y)
            μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
            in_rate = d * μ_in
            out_coeff = d * length(ns)

            if any(nb -> haskey(s.dists, nb), ns)
                # ── BOUNDARY: Galerkin expv ───────────────────────────────────────
                Q_eq     = _single_voxel_generator(s.model, in_rate, out_coeff, s.n_max_eq)
                p_old_eq = _resize_dist(p_eq, s.n_max_eq)
                p_new_eq = expv(Float64(dt), Q_eq, p_old_eq; m = min(s.krylov_m_eq, size(Q_eq, 1)))
                _clean!(p_new_eq, s.ε_prune)
                new_equil[idx] = p_new_eq
            else
                # ── INTERIOR: algebraic mean → Poisson (no matrix exp) ───────────
                μ_loc = (s.model.k_b + in_rate) / (s.model.k_d + out_coeff)
                new_equil[idx] = _poisson_pmf(μ_loc, s.n_max_eq)
            end
        end
        s.equil_dists = new_equil
    else
        empty!(s.equil_dists)
    end

    # ── 4. Accumulate flux into empty voxels ─────────────────────────────────
    for (idx, μ) in means
        μ > 0.0 || continue
        for nb in grid_neighbors(idx, s.K_x, s.K_y)
            haskey(s.dists, nb)       && continue
            haskey(s.equil_dists, nb) && continue
            s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + d * μ * dt
        end
    end

    # ── 5. Activate (prolong) empty voxels ────────────────────────────────────
    for (nb, fl) in s.pending_flux
        if fl >= s.ε_expand && !haskey(s.dists, nb) && !haskey(s.equil_dists, nb)
            p0 = zeros(s.n_max + 1);  p0[1] = 1.0
            s.dists[nb] = p0
            delete!(s.pending_flux, nb)
        end
    end

    # ── 6. Restrict: deactivate voxels near their local Poisson SS ───────────
    # Store the CURRENT active distribution as the initial equil_dist.
    # No information is lost: the restriction is exact given the local-SS criterion.
    if has_finite_local_ss(s.model)
        to_deactivate = CartesianIndex{2}[]
        for (idx, p) in s.dists
            ns      = grid_neighbors(idx, s.K_x, s.K_y)
            μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
            μ_loc   = (s.model.k_b + d * μ_in) / (s.model.k_d + d * length(ns))

            μ = voxel_mean(p)
            abs(μ - μ_loc) / (μ_loc + 1e-10) < s.ε_equil || continue
            p_loc = _poisson_pmf(μ_loc, s.n_max)
            tv = 0.5 * sum(abs(p[i] - p_loc[i]) for i in eachindex(p_loc))
            tv < 0.08 || continue

            push!(to_deactivate, idx)
            s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
        end
        for idx in to_deactivate; delete!(s.dists, idx); end
    end

    if s.use_joint_active
        _sync_joint_state!(s, _active_voxels(s)) || begin
            s.joint_state = nothing
            empty!(s.joint_voxels)
        end
        empty!(s.block_states)
    elseif s.use_block_vcycle
        next_blocks = _select_active_blocks(s.dists, s.K_x, s.K_y, collect(keys(s.block_states)))
        _sync_block_cache!(s, next_blocks)
    else
        s.joint_state = nothing
        empty!(s.joint_voxels)
        empty!(s.block_states)
    end

    s.t += dt
    return s
end

# ─── Grid views ──────────────────────────────────────────────────────────────

function status_grid(s::SpatialFSP)
    g = zeros(Int, s.K_x, s.K_y)
    for idx in keys(s.equil_dists); g[idx] = 2; end
    for idx in keys(s.dists);       g[idx] = 1; end
    g
end

function mean_grid(s::SpatialFSP)
    g = zeros(s.K_x, s.K_y)
    for (idx, p) in s.equil_dists; g[idx] = voxel_mean(p); end
    for (idx, p) in s.dists;       g[idx] = voxel_mean(p); end
    g
end

"""
total_distribution:
- if a joint active-region state is present, use its exact total-count law;
- otherwise convolve active voxel marginals;
- then convolve equilibrated voxel distributions.
"""
function total_distribution(s::SpatialFSP)
    function _conv!(a, b)
        n1, n2 = length(a), length(b)
        c = zeros(n1 + n2 - 1)
        for i in 1:n1, j in 1:n2; c[i+j-1] += a[i]*b[j]; end
        resize!(a, length(c));  a .= c
    end

    p_tot = if s.use_joint_active && s.joint_state !== nothing && !isempty(s.joint_voxels)
        _joint_total_distribution(s.joint_state)
    else
        [1.0]
    end

    if !(s.use_joint_active && s.joint_state !== nothing && !isempty(s.joint_voxels))
        for (_, p) in s.dists; _conv!(p_tot, p); end
    end
    for (_, p) in s.equil_dists; _conv!(p_tot, p); end
    p_tot
end

function equilibration_fraction(s::SpatialFSP; rtol::Float64 = 0.15)
    has_finite_local_ss(s.model) || return 0.0
    μ_ss  = ss_mean(s.model)
    n_tot = s.K_x * s.K_y
    n_done = 0
    for (_, p) in s.equil_dists
        abs(voxel_mean(p) - μ_ss)/μ_ss < rtol && (n_done += 1)
    end
    for (_, p) in s.dists
        abs(voxel_mean(p) - μ_ss)/μ_ss < rtol && (n_done += 1)
    end
    n_done / n_tot
end
