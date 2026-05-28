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
    use_multigrid       :: Bool     # true multigrid: generalized spatial V-cycle
    multigrid_levels    :: Int      # number of spatial levels (1=pairs, 2=quads, 3=8-blocks…)
    use_block_exact         :: Bool     # exact Galerkin pairs: full 2-body joint CME, ADI H/V passes
    use_block_vcycle_4body  :: Bool     # 2×2 block V-cycle: compact rstep restrict + correction prolong
    block_vcycle_4body_max  :: Int      # max rstep-expanded states per block (skip dense interior)
    use_amg_block           :: Bool     # no mean-field: 2×2 active + frontier joint, rstep state space
    amg_block_max_states    :: Int      # max rstep states per extended block
    amg_activate_tol        :: Float64  # activate frontier voxel if P(n≥1) > this
    # Organic activation: bypass pending_flux accumulation entirely.
    # Activate empty neighbor nb when any adjacent tracked voxel k has P(n_k>0) > ε_organic.
    # Eliminates activation delay (ε_expand) — the main source of wave speed error.
    organic_activation :: Bool
    ε_organic          :: Float64
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
                    joint_vcycle_τ_post::Float64 = 0.0,
                    use_multigrid::Bool = false,
                    multigrid_levels::Int = 3,
                    use_block_exact::Bool = false,
                    use_block_vcycle_4body::Bool = false,
                    block_vcycle_4body_max::Int = 3000,
                    use_amg_block::Bool = false,
                    amg_block_max_states::Int = 30_000,
                    amg_activate_tol::Float64 = 1e-4,
                    organic_activation::Bool = false,
                    ε_organic::Float64 = 1e-4)
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
               use_joint_vcycle, joint_vcycle_τ_pre, joint_vcycle_τ_post,
               use_multigrid, multigrid_levels,
               use_block_exact, use_block_vcycle_4body, block_vcycle_4body_max,
               use_amg_block, amg_block_max_states, amg_activate_tol,
               organic_activation, ε_organic)
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

"""
    _single_voxel_generator_range(model, in_rate, out_coeff, n_lo, n_hi)

Build a single-voxel CME generator restricted to molecule counts [n_lo, n_hi].
Output is (n_hi-n_lo+1) × (n_hi-n_lo+1) — smaller than the full [0, n_max] matrix.
Used by _expv_voxel to restrict computation to the active support of p.
"""
function _single_voxel_generator_range(model::AbstractSpatialRDMEModel,
                                         in_rate::Float64, out_coeff::Float64,
                                         n_lo::Int, n_hi::Int;
                                         n_sites::Int = 1)
    nv = n_hi - n_lo + 1
    Q  = zeros(Float64, nv, nv)
    for col in 1:nv
        n = col - 1 + n_lo   # actual molecule count
        up = in_rate + local_total_birth_rate(model, n, n_sites)
        if col < nv && up > 0
            Q[col+1, col] += up;  Q[col, col] -= up
        end
        down = _voxel_down_rate(model, n, out_coeff)
        if col > 1 && down > 0
            Q[col-1, col] += down;  Q[col, col] -= down
        end
    end
    Q
end

# Per-model downward rate (excludes diffusion-out which is added as out_coeff×n)
_voxel_down_rate(m::BirthDeathRDME,    n, out_coeff) = (m.k_d + out_coeff) * n
_voxel_down_rate(m::BranchingDeathRDME, n, out_coeff) = (m.k_d + out_coeff) * n
# SchloglModel1D dispatch defined in schlogl_spatial_fsp.jl (included after rdme_model.jl)

# ── Sparse Tridiagonal overrides for BD models ────────────────────────────────
# BirthDeath and BranchingDeath are also birth-death chains (tridiagonal).
# Returning Tridiagonal instead of a dense matrix reduces expv mv cost
# from O(n²) to O(n), giving the same speedup as for Schlögl.
function _single_voxel_generator_range(model::BirthDeathRDME,
                                         in_rate::Float64, out_coeff::Float64,
                                         n_lo::Int, n_hi::Int;
                                         n_sites::Int = 1)
    nv = n_hi - n_lo + 1
    dl = zeros(Float64, nv - 1)
    d  = zeros(Float64, nv)
    du = zeros(Float64, nv - 1)
    for col in 1:nv
        n    = col - 1 + n_lo
        up   = in_rate + local_total_birth_rate(model, n, n_sites)
        down = _voxel_down_rate(model, n, out_coeff)
        if col < nv && up > 0;   dl[col]    = up;   d[col]   -= up;   end
        if col > 1 && down > 0;  du[col-1]  = down; d[col]   -= down; end
    end
    LinearAlgebra.Tridiagonal(dl, d, du)
end

function _single_voxel_generator_range(model::BranchingDeathRDME,
                                         in_rate::Float64, out_coeff::Float64,
                                         n_lo::Int, n_hi::Int;
                                         n_sites::Int = 1)
    nv = n_hi - n_lo + 1
    dl = zeros(Float64, nv - 1)
    d  = zeros(Float64, nv)
    du = zeros(Float64, nv - 1)
    for col in 1:nv
        n    = col - 1 + n_lo
        up   = in_rate + local_total_birth_rate(model, n, n_sites)
        down = _voxel_down_rate(model, n, out_coeff)
        if col < nv && up > 0;   dl[col]    = up;   d[col]   -= up;   end
        if col > 1 && down > 0;  du[col-1]  = down; d[col]   -= down; end
    end
    LinearAlgebra.Tridiagonal(dl, d, du)
end

"""
    _expv_voxel(model, in_rate, out_coeff, p, dt, krylov_m; margin, prune_tol) → p_new

Evolve a per-voxel distribution using expv on the ACTIVE SUPPORT of p only.
Equivalent to rstep for dense 1D distributions:
  - Finds [lo, hi] = range where p > prune_tol
  - Builds generator only for [lo-margin, hi+margin] ← small matrix
  - Runs expv on that compact system
  - Returns result in original-length vector

For Schlögl at n_hi=162:  ~40 states instead of 221  (5× speedup)
For BirthDeath at n≈3:    ~15 states instead of  20  (modest)
The margin ensures boundary losses are negligible.
"""
function _expv_voxel(model::AbstractSpatialRDMEModel,
                      in_rate::Float64, out_coeff::Float64,
                      p::Vector{Float64}, dt::Float64,
                      krylov_m::Int;
                      n_sites::Int   = 1,
                      margin::Int    = 8,
                      n_include::Int = -1,     # force this molecule count into support
                      lo_override::Int = -1,   # override lower bound of support
                      hi_override::Int = -1,   # override upper bound of support
                      prune_tol::Float64 = 1e-10)
    n_max  = length(p) - 1
    lo_idx = findfirst(>(prune_tol), p)
    hi_idx = findlast(>(prune_tol), p)
    isnothing(lo_idx) && return copy(p)   # zero distribution

    lo_mol = max(0,     lo_idx - 1 - margin)
    hi_mol = min(n_max, hi_idx - 1 + margin)

    if n_include >= 0
        lo_mol = min(lo_mol, max(0,     n_include - margin))
        hi_mol = max(hi_mol, min(n_max, n_include + margin))
    end
    # Hard overrides: for bistable models with large drifts, the support must span
    # the full state space or probability escapes the compact support in one step.
    if lo_override >= 0; lo_mol = lo_override; end
    if hi_override >= 0; hi_mol = hi_override; end

    Q_sub  = _single_voxel_generator_range(model, in_rate, out_coeff, lo_mol, hi_mol;
                                            n_sites)
    p_sub  = p[lo_mol+1 : hi_mol+1]
    p_sub_new = expv(Float64(dt), Q_sub, p_sub; m = min(krylov_m, size(Q_sub, 1)))

    p_new = zeros(Float64, n_max + 1)
    p_new[lo_mol+1 : hi_mol+1] = max.(0.0, p_sub_new)
    p_new
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

function _build_pair_system_full(model :: AbstractSpatialRDMEModel,
                                v1    :: CartesianIndex{2},
                                v2    :: CartesianIndex{2},
                                means :: Dict{CartesianIndex{2}, Float64},
                                K_x   :: Int,
                                K_y   :: Int)
    d = jump_rate(model)
    stoichs      = CartesianIndex{2}[]
    propensities = Function[]
    pair_set     = (v1, v2)

    # Local reactions for each voxel
    for p in 1:2
        let p = p
            ep = _eu_local(p, 2)
            push!(stoichs,  ep); push!(propensities, (x,r,t) -> local_birth_rate(model, Tuple(x)[p]))
            push!(stoichs, -ep); push!(propensities, (x,r,t) -> model.k_d * max(0, Tuple(x)[p]))
        end
    end

    # Intra-pair diffusion
    push!(stoichs, CartesianIndex(-1, 1))
    push!(propensities, (x,r,t) -> d * max(0, Tuple(x)[1]))
    push!(stoichs, CartesianIndex( 1,-1))
    push!(propensities, (x,r,t) -> d * max(0, Tuple(x)[2]))

    # External coupling: each voxel's non-pair neighbors
    for (p, vox) in ((1, v1), (2, v2))
        for nb in grid_neighbors(vox, K_x, K_y)
            (nb === v1 || nb === v2) && continue
            μ_nb = get(means, nb, 0.0)
            let p = p, μ_nb = μ_nb
                ep = _eu_local(p, 2)
                push!(stoichs, -ep)
                push!(propensities, (x,r,t) -> d * max(0, Tuple(x)[p]))
                μ_nb > 0.0 || continue
                push!(stoichs,  ep)
                push!(propensities, (x,r,t) -> d * μ_nb)
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{2}}(stoichs, propensities)
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

# ─── Exact-Galerkin joint evolution (no total-count restriction) ─────────────
#
# Pair (2-body, ADI): Q_pair = L₁ + L₂ + D₁₂ — 1D coarsening, alternating H/V
#
# Block V-cycle (4-body): proper 2-level multigrid
#   Restriction:  compact rstep seed (top-k states/voxel) + depth-d expansion
#   Coarse solve: full 4-body joint CME with all 4 intra-block diffusion directions
#   Prolongation: multiplicative correction  p_new = p_pre * (marg_final/marg_init)
#                 preserves pre-smooth information, corrects only what the joint adds

"""
    _evolve_block_vcycle_4body(model, block, p_vox, means, K_x, K_y, dt, n_max; ...)

2×2 block V-cycle for the RDME:
  1. Pre-smooth: per-voxel expv for τ = dt * tau_pre_frac
  2. Restrict:   top seed_k states per voxel → product → rstep expand (depth rstep_depth)
  3. Project:    per-voxel probability onto rstep subspace (independence approx)
  4. Coarse solve: evolve joint 4-body CME for τ_coarse = dt * (1 - 2*tau_pre_frac)
  5. Prolongate:  p_new[n] = p_pre[n] * (marg_final[n] / marg_init[n])
                  (states outside subspace keep p_pre unchanged)

Returns NTuple{4, Vector{Float64}} of updated per-voxel distributions, or nothing on fail.
"""
function _evolve_block_vcycle_4body(model        :: AbstractSpatialRDMEModel,
                                     block        :: NTuple{4, CartesianIndex{2}},
                                     p_vox        :: NTuple{4, Vector{Float64}},
                                     means        :: Dict{CartesianIndex{2}, Float64},
                                     K_x          :: Int, K_y :: Int,
                                     dt           :: Float64,
                                     n_max        :: Int;
                                     krylov_m     :: Int     = 30,
                                     state_tol    :: Float64 = 1e-10,
                                     max_states   :: Int     = 3000,
                                     tau_pre_frac :: Float64 = 0.3,
                                     seed_k       :: Int     = 3,
                                     rstep_depth  :: Int     = 2)

    d    = jump_rate(model)
    τ_pre    = dt * tau_pre_frac
    τ_coarse = dt - τ_pre            # joint CME covers the rest; no post-smooth
    τ_coarse ≤ 0.0 && return nothing

    block_set = Set(block)

    # ── 1. Pre-smooth (external coupling only) ───────────────────────────────
    # Intra-block diffusion is handled exclusively by the joint CME.
    # Including it here would double-count it and homogenize the block.
    p_pre = ntuple(k -> begin
        idx      = block[k]
        ext_ns   = filter(nb -> !(nb in block_set), grid_neighbors(idx, K_x, K_y))
        μ_in     = sum(get(means, nb, 0.0) for nb in ext_ns; init = 0.0)
        in_rate  = d * μ_in
        out_coef = d * length(ext_ns)
        isempty(ext_ns) ?
            copy(p_vox[k]) :
            _expv_voxel(model, in_rate, out_coef, p_vox[k], Float64(τ_pre), krylov_m)
    end, 4)

    # ── 2. Compact rstep restriction: top-seed_k states per voxel ────────────
    top_idx = ntuple(k -> begin
        p   = p_pre[k]
        nz  = count(>(state_tol), p)
        k_s = max(1, min(seed_k, nz))
        partialsortperm(p, 1:k_s; rev = true)
    end, 4)

    sp = StateSpace{CartesianIndex{4}, Float64}()
    for i1 in top_idx[1], i2 in top_idx[2], i3 in top_idx[3], i4 in top_idx[4]
        add_state!(sp, CartesianIndex(i1-1, i2-1, i3-1, i4-1), 1.0)
    end
    isempty(sp.states) && return nothing

    sys_full, _, _ = _build_block_system(model, block, means, K_x, K_y)
    expand!(sp, sys_full, x -> _block_bc(x, n_max); depth = rstep_depth)
    length(sp) > max_states && return nothing

    # ── 3. Project pre-smoothed probability onto rstep subspace ─────────────
    for (i, ci) in enumerate(sp.states)
        t = Tuple(ci)
        sp.probs[i] = prod(k -> p_pre[k][t[k]+1], 1:4)
    end
    s_init = sum(sp.probs)
    s_init ≤ 0.0 && return nothing
    sp.probs ./= s_init

    marg_init = _block_marginals(sp, n_max)   # marginals before coarse solve

    # ── 4. Coarse solve ──────────────────────────────────────────────────────
    A, = build_generator(sp, sys_full, nothing, 0.0)
    sp.probs .= expv(Float64(τ_coarse), A, sp.probs; m = min(krylov_m, size(A, 1)))
    _clean!(sp.probs, state_tol)
    renormalize!(sp) > 0.0 || return nothing

    marg_final = _block_marginals(sp, n_max)  # marginals after coarse solve

    # ── 5. Multiplicative correction prolongation ────────────────────────────
    # p_new[n] = p_pre[n] * (marg_final[n] / marg_init[n])
    # States outside subspace (marg_init[n]=0): keep p_pre unchanged.
    ntuple(k -> begin
        p_out = copy(p_pre[k])
        mi = marg_init[k];  mf = marg_final[k]
        for n in 0:n_max
            mi[n+1] > state_tol || continue
            p_out[n+1] = p_pre[k][n+1] * (mf[n+1] / mi[n+1])
        end
        _clean!(p_out, state_tol)
        p_out
    end, 4)
end

function _step_block_vcycle_4body!(s        :: SpatialFSP,
                                    means    :: Dict{CartesianIndex{2}, Float64},
                                    dt       :: Float64;
                                    krylov_m :: Int = 30,
                                    step_idx :: Int = 0)
    # 4-step offset cycle so every row-pair and col-pair gets V-cycle coupled:
    #   pass 0: anchors at (odd_row,  odd_col)   → (1,1),(1,3),(3,1),(3,3)...
    #   pass 1: anchors at (even_row, even_col)  → (2,2),(2,4),(4,2),(4,4)...
    #   pass 2: anchors at (odd_row,  even_col)  → (1,2),(1,4),(3,2),(3,4)...
    #   pass 3: anchors at (even_row, odd_col)   → (2,1),(2,3),(4,1),(4,3)...
    pass    = step_idx % 4
    i_off   = pass ∈ (1, 3) ? 1 : 0   # row offset: 0=odd start, 1=even start
    j_off   = pass ∈ (1, 2) ? 1 : 0   # col offset
    K_x, K_y = s.K_x, s.K_y
    preferred = [CartesianIndex(i, j)
                 for i in (1+i_off):2:(K_x-1) for j in (1+j_off):2:(K_y-1)]
    anchors   = _select_active_blocks(s.dists, K_x, K_y, preferred)
    handled   = Set{CartesianIndex{2}}()

    for anchor in anchors
        block = _block_voxels(anchor)
        p_vox = ntuple(k -> s.dists[block[k]], 4)

        result = _evolve_block_vcycle_4body(s.model, block, p_vox, means,
                                            s.K_x, s.K_y, dt, s.n_max;
                                            krylov_m,
                                            state_tol    = s.block_state_tol,
                                            max_states   = s.block_vcycle_4body_max,
                                            tau_pre_frac = s.block_pre_frac)
        result === nothing && continue

        for (k, idx) in enumerate(block)
            m = result[k];  sm = sum(m)
            sm > 0 || continue   # degenerate → fall through to per-voxel
            s.dists[idx] = m ./ sm
            push!(handled, idx)
        end
    end
    handled
end

function _evolve_pair_exact(model    :: AbstractSpatialRDMEModel,
                             v1      :: CartesianIndex{2},
                             v2      :: CartesianIndex{2},
                             p1      :: Vector{Float64},
                             p2      :: Vector{Float64},
                             means   :: Dict{CartesianIndex{2}, Float64},
                             K_x     :: Int, K_y :: Int,
                             dt      :: Float64,
                             n_max   :: Int;
                             krylov_m  :: Int   = 30,
                             state_tol :: Float64 = 1e-10,
                             max_states :: Int  = 2000)

    # Build sparse product state
    supp1 = findall(>(state_tol), p1); isempty(supp1) && (supp1 = [argmax(p1)])
    supp2 = findall(>(state_tol), p2); isempty(supp2) && (supp2 = [argmax(p2)])
    length(supp1) * length(supp2) > max_states && return nothing

    sp = StateSpace{CartesianIndex{2}, Float64}()
    for i1 in supp1, i2 in supp2
        w = p1[i1] * p2[i2]; w > state_tol || continue
        add_state!(sp, CartesianIndex(i1 - 1, i2 - 1), w)
    end
    isempty(sp.states) && return nothing
    renormalize!(sp)

    # Full pair system: reactions + intra-diffusion + external coupling
    sys = _build_pair_system_full(model, v1, v2, means, K_x, K_y)
    bc  = x -> all(c -> 0 ≤ c ≤ n_max, Tuple(x))

    expand!(sp, sys, bc; depth = 1)
    length(sp) > max_states && return nothing

    # Evolve full joint CME — exact Galerkin, no restriction
    A, = build_generator(sp, sys, nothing, 0.0)
    sp.probs .= expv(Float64(dt), A, sp.probs; m = min(krylov_m, size(A, 1)))
    _clean!(sp.probs, state_tol)

    # Extract marginals
    p1_new = zeros(n_max + 1); p2_new = zeros(n_max + 1)
    for (ci, prob) in zip(sp.states, sp.probs)
        t = Tuple(ci)
        t[1] ≤ n_max && (p1_new[t[1] + 1] += prob)
        t[2] ≤ n_max && (p2_new[t[2] + 1] += prob)
    end
    s1 = sum(p1_new); s2 = sum(p2_new)
    s1 > 0 && (p1_new ./= s1); s2 > 0 && (p2_new ./= s2)
    (p1_new, p2_new)
end

function _step_pair_exact!(s        :: SpatialFSP,
                            means    :: Dict{CartesianIndex{2}, Float64},
                            dt       :: Float64;
                            krylov_m :: Int  = 30,
                            row_pass :: Bool = true)
    n_max     = s.n_max
    state_tol = s.block_state_tol
    max_s     = s.block_max_states
    Δ = row_pass ? CartesianIndex(0, 1) : CartesianIndex(1, 0)

    handled = Set{CartesianIndex{2}}()
    for v1 in keys(s.dists)
        v1 ∈ handled && continue
        v2 = v1 + Δ
        haskey(s.dists, v2) || continue
        v2 ∈ handled && continue

        result = _evolve_pair_exact(s.model, v1, v2,
                                    s.dists[v1], s.dists[v2],
                                    means, s.K_x, s.K_y, dt, n_max;
                                    krylov_m, state_tol, max_states = max_s)
        result === nothing && continue
        s.dists[v1], s.dists[v2] = result
        push!(handled, v1); push!(handled, v2)
    end
    handled
end

# ─── Multi-level stochastic prolongation for pair V-cycle ────────────────────
#
# Three-level hierarchy along the split axis n₁|N:
#   Level 0 (fine):     p(n₁|N)        — exact conditional from pre-smooth
#   Level 1 (mid):      p(⌊n₁/2⌋|N)   — bin-pair coarsening of split
#   Level 2 (coarsest): p(N) only      — total count (from 1D CME coarse solve)
#
# Prolongation (level 2 → 0):
#   Existing N (p_pre[N]>0): multiplicative correction on level-0 fiber — exact
#   New N:  1. propagate level-1 fiber from nearest known N via coarse L/R
#           2. refine level-1 → level-0 using within-bin weights from nearest
#              known level-0 fiber (borrowed across N — changes slowly)

# Build level-0 and level-1 fibers from the pre-smooth joint sp.
# fiber0[N][n₁]   = p(n₁|N)        — exact conditional
# fiber1[N][⌊n₁/2⌋] = p(⌊n₁/2⌋|N) — coarsened conditional
function _build_pair_fibers(sp::StateSpace{CartesianIndex{2},Float64},
                             nmax2::Int, state_tol::Float64)
    fibers0 = Dict{Int, Dict{Int,Float64}}()
    fibers1 = Dict{Int, Dict{Int,Float64}}()
    for (ci, prob) in zip(sp.states, sp.probs)
        prob < state_tol && continue
        n1, n2 = Tuple(ci)
        N = n1 + n2
        N > nmax2 && continue
        d0 = get!(fibers0, N, Dict{Int,Float64}())
        d0[n1] = get(d0, n1, 0.0) + prob
        d1 = get!(fibers1, N, Dict{Int,Float64}())
        m  = n1 ÷ 2
        d1[m] = get(d1, m, 0.0) + prob
    end
    for d in values(fibers0); s = sum(values(d)); s > 1e-300 && map!(v -> v/s, values(d)); end
    for d in values(fibers1); s = sum(values(d)); s > 1e-300 && map!(v -> v/s, values(d)); end
    fibers0, fibers1
end

# Lift level-1 fiber from N to N+1 (birth event, coarse L/R operator).
# Uses uniform within-bin weights — only called for new N values where
# fine structure is unknown. Accurate when propensities change slowly.
function _lift_fiber1(fiber1::Dict{Int,Float64}, N::Int,
                       in1::Float64, in2::Float64,
                       model::AbstractSpatialRDMEModel, n_max::Int)
    new_f = Dict{Int,Float64}()
    for (m, prob) in fiber1
        prob < 1e-300 && continue
        for n1 in (2m, 2m+1)           # both fine states in coarse bin m
            (n1 > min(N, n_max) || N - n1 < 0) && continue
            n2  = N - n1
            b1  = local_total_birth_rate(model, n1, 1) + in1
            b2  = local_total_birth_rate(model, n2, 1) + in2
            tot = b1 + b2; tot < 1e-300 && continue
            w   = 0.5 * prob            # uniform within-bin weight
            if n1 + 1 <= n_max          # born at v1: n₁ → n₁+1
                m1 = (n1 + 1) ÷ 2
                new_f[m1] = get(new_f, m1, 0.0) + w * b1/tot
            end
            new_f[m] = get(new_f, m, 0.0) + w * b2/tot   # born at v2: n₁ unchanged
        end
    end
    s = sum(values(new_f)); s > 1e-300 && map!(v -> v/s, values(new_f))
    new_f
end

# Lower level-1 fiber from N to N-1 (death event, coarse L/R operator).
function _lower_fiber1(fiber1::Dict{Int,Float64}, N::Int,
                        out1::Float64, out2::Float64,
                        model::AbstractSpatialRDMEModel)
    N <= 0 && return Dict{Int,Float64}()
    new_f = Dict{Int,Float64}()
    for (m, prob) in fiber1
        prob < 1e-300 && continue
        for n1 in (2m, 2m+1)
            n1 > N && continue
            n2  = N - n1
            d1  = _voxel_down_rate(model, n1, out1)
            d2  = _voxel_down_rate(model, n2, out2)
            tot = d1 + d2; tot < 1e-300 && continue
            w   = 0.5 * prob
            if n1 > 0                   # death at v1: n₁ → n₁-1
                m1 = (n1 - 1) ÷ 2
                new_f[m1] = get(new_f, m1, 0.0) + w * d1/tot
            end
            new_f[m] = get(new_f, m, 0.0) + w * d2/tot   # death at v2: n₁ unchanged
        end
    end
    s = sum(values(new_f)); s > 1e-300 && map!(v -> v/s, values(new_f))
    new_f
end

# Refine: level-1 → level-0 by splitting each coarse bin using within-bin
# weights borrowed from a reference level-0 fiber at the nearest known N.
function _refine_fiber1(fiber1::Dict{Int,Float64},
                         ref0  ::Dict{Int,Float64},
                         n_max ::Int)
    fiber0 = Dict{Int,Float64}()
    for (m, cp) in fiber1
        cp < 1e-300 && continue
        n1_lo, n1_hi = 2m, 2m+1
        p_lo = get(ref0, n1_lo, 0.0)
        p_hi = get(ref0, n1_hi, 0.0)
        tot  = p_lo + p_hi
        w_lo, w_hi = tot > 1e-300 ? (p_lo/tot, p_hi/tot) : (0.5, 0.5)
        n1_lo ≤ n_max && (fiber0[n1_lo] = get(fiber0, n1_lo, 0.0) + cp * w_lo)
        n1_hi ≤ n_max && (fiber0[n1_hi] = get(fiber0, n1_hi, 0.0) + cp * w_hi)
    end
    s = sum(values(fiber0)); s > 1e-300 && map!(v -> v/s, values(fiber0))
    fiber0
end

# Multi-level prolongation: assemble sp_new from fibers and p_post weights.
# For existing N: exact multiplicative correction (fiber0 × p_post[N]).
# For new N:      coarse L/R propagation (fiber1) → level-0 refinement.
function _prolong_multilevel(sp        ::StateSpace{CartesianIndex{2},Float64},
                              p_pre     ::Vector{Float64},
                              p_post    ::Vector{Float64},
                              in1       ::Float64, in2  ::Float64,
                              out1      ::Float64, out2 ::Float64,
                              model     ::AbstractSpatialRDMEModel,
                              n_max     ::Int,
                              state_tol ::Float64,
                              max_states::Int)
    nmax2 = 2 * n_max
    fibers0, fibers1 = _build_pair_fibers(sp, nmax2, state_tol)
    isempty(fibers0) && return nothing

    known_Ns   = sort(collect(keys(fibers0)))
    derived0   = Dict{Int, Dict{Int,Float64}}()   # cache for new-N fiber0
    derived1   = Dict{Int, Dict{Int,Float64}}()   # cache for new-N fiber1

    # Derive fiber0 at a new N via: coarse L/R on fiber1 → refine → cache.
    function derive_fiber0(N::Int)
        haskey(fibers0,  N) && return fibers0[N]
        haskey(derived0, N) && return derived0[N]
        isempty(known_Ns) && return nothing

        # Nearest known anchor for both fibers
        idx      = searchsortedfirst(known_Ns, N)
        N_anchor = idx > length(known_Ns) ? known_Ns[end] :
                   idx == 1               ? known_Ns[1]   :
                   (N - known_Ns[idx-1] ≤ known_Ns[idx] - N) ? known_Ns[idx-1] : known_Ns[idx]

        # Propagate level-1 fiber step by step from anchor → N (coarse, cheap)
        cur1 = get(fibers1, N_anchor, Dict{Int,Float64}())
        cur_N = N_anchor
        step  = cur_N < N ? 1 : -1
        while cur_N != N
            next1 = step > 0 ? _lift_fiber1(cur1,  cur_N, in1, in2, model, n_max) :
                                _lower_fiber1(cur1, cur_N, out1, out2, model)
            cur_N += step
            isempty(next1) && return nothing
            if !haskey(fibers1, cur_N) && !haskey(derived1, cur_N)
                derived1[cur_N] = next1
            end
            cur1 = get(fibers1, cur_N, get(derived1, cur_N, next1))
        end

        # Refine level-1 → level-0 using within-bin weights from anchor fiber0
        ref0  = get(fibers0, N_anchor, Dict{Int,Float64}())
        f0new = _refine_fiber1(cur1, ref0, n_max)
        derived0[N] = f0new
        f0new
    end

    # Build sp_new
    sp_new = StateSpace{CartesianIndex{2}, Float64}()
    for N in 0:nmax2
        p_post[N+1] > state_tol || continue
        fiber0 = (haskey(fibers0, N) && p_pre[N+1] > state_tol) ?
                 fibers0[N] : derive_fiber0(N)
        fiber0 === nothing && continue
        isempty(fiber0)   && continue
        for (n1, cp) in fiber0
            n2 = N - n1
            (0 ≤ n1 ≤ n_max && 0 ≤ n2 ≤ n_max) || continue
            w = cp * p_post[N+1]
            w > state_tol || continue
            ci  = CartesianIndex(n1, n2)
            idx = get(sp_new.index, ci, 0)
            if idx == 0; add_state!(sp_new, ci, w)
            else;        sp_new.probs[idx] += w; end
        end
        length(sp_new) > max_states && return nothing
    end
    sp_new
end

# ─── ADI pair V-cycle (K=2, correct 2D, O(N) new states vs O(N³) for K=4) ──

"""
    _evolve_pair_vcycle(model, v1, v2, p1, p2, means, K_x, K_y, dt, n_max; ...)
        → (p1_new, p2_new) or nothing

K=2 V-cycle for a single adjacent pair of active voxels.

Structure (same as K=4 block but 50× cheaper for new states):
  1. Product state IC from p1⊗p2 (sparse, state_tol filtered)
  2. Pre-smooth: pair diffusion only (τ_pre fraction of dt)
  3. Restrict: p_total[N] = Σ_{n1+n2=N} p(n1,n2)  →  1D vector (2n_max+1 states)
  4. Coarse solve: 1D CME on total N, includes external boundary flux
  5. Prolong: multi-level stochastic interpolation
               — existing N: exact multiplicative correction (fiber0 × p_post/p_pre)
               — new N:      coarse L/R on fiber1 → refine to fiber0
  6. Post-smooth: pair diffusion (τ_post fraction of dt)
  7. Extract marginals → p1_new, p2_new

ADI (alternating row/column sweeps) captures all 4 nearest-neighbor correlations
in 2D by alternating between horizontal and vertical pair sweeps each step.
"""
function _evolve_pair_vcycle(
        model   :: AbstractSpatialRDMEModel,
        v1      :: CartesianIndex{2},
        v2      :: CartesianIndex{2},
        p1      :: Vector{Float64},
        p2      :: Vector{Float64},
        means   :: Dict{CartesianIndex{2}, Float64},
        K_x     :: Int, K_y :: Int,
        dt      :: Float64,
        n_max   :: Int;
        krylov_m :: Int    = 30,
        state_tol :: Float64 = 1e-8,
        pre_frac  :: Float64 = 0.1,
        post_frac :: Float64 = 0.1,
        max_states :: Int   = 2000)

    d = jump_rate(model)

    # ── 1. Build sparse product state p1⊗p2 ──────────────────────────────────
    supp1 = findall(>(state_tol), p1); isempty(supp1) && (supp1 = [argmax(p1)])
    supp2 = findall(>(state_tol), p2); isempty(supp2) && (supp2 = [argmax(p2)])
    length(supp1) * length(supp2) > max_states && return nothing

    sp = StateSpace{CartesianIndex{2}, Float64}()
    for i1 in supp1, i2 in supp2
        w = p1[i1] * p2[i2]
        w > state_tol || continue
        add_state!(sp, CartesianIndex(i1-1, i2-1), w)
    end
    isempty(sp.states) && return nothing
    renormalize!(sp)

    # ── 1b. Organic rstep expansion of the fine-level state space ─────────────
    # Expand from the product state using the pair RDME system.
    # This is the organic bound: only states reachable by one reaction from
    # the current support are added (no fixed n_max pre-allocation).
    pair_intra = DiscreteStochasticSystem{CartesianIndex{2}}(
        [CartesianIndex(-1,1), CartesianIndex(1,-1)],
        [(x,r,t) -> d*max(0,Tuple(x)[1]),
         (x,r,t) -> d*max(0,Tuple(x)[2])])
    bc2 = x -> all(c -> 0 ≤ c ≤ n_max, Tuple(x))
    expand!(sp, pair_intra, bc2; depth=1)   # rstep: organic, not n_max-bounded
    length(sp) > max_states && return nothing

    # ── 2. Pre-smooth: intra-pair diffusion only ──────────────────────────────
    τ_pre = pre_frac * dt
    if τ_pre > 0.0 && !isempty(sp.states)
        A, = build_generator(sp, pair_intra, nothing, 0.0)
        sp.probs .= expv(τ_pre, A, sp.probs; m=min(krylov_m,size(A,1)))
        _clean!(sp.probs, state_tol)
    end

    # ── 3. Restrict: 1D total count distribution ──────────────────────────────
    nmax2 = 2 * n_max
    p_pre = zeros(nmax2 + 1)
    for (ci, prob) in zip(sp.states, sp.probs)
        N = sum(Tuple(ci)); N <= nmax2 && (p_pre[N+1] += prob)
    end
    p_pre ./= max(sum(p_pre), 1e-100)

    # ── 4. Coarse solve: 1D CME on total N ───────────────────────────────────
    # Boundary contributions: external equil/empty neighbours of v1 and v2
    in1 = sum(d * get(means, nb, 0.0)
              for nb in grid_neighbors(v1, K_x, K_y)
              if nb != v2 && !haskey(means, nb); init=0.0)
    in1 += sum(d * get(means, nb, 0.0)
               for nb in grid_neighbors(v1, K_x, K_y)
               if nb != v2 && haskey(means, nb); init=0.0)
    in2  = sum(d * get(means, nb, 0.0)
               for nb in grid_neighbors(v2, K_x, K_y)
               if nb != v1; init=0.0)
    in_rate   = in1 + in2

    out1 = sum(d for nb in grid_neighbors(v1, K_x, K_y) if nb != v2 && !haskey(means, nb); init=0.0)
    out1 += sum(d for nb in grid_neighbors(v1, K_x, K_y) if nb != v2 && haskey(means, nb); init=0.0)
    out2  = sum(d for nb in grid_neighbors(v2, K_x, K_y) if nb != v1; init=0.0)
    out_coeff = (out1 + out2) / 2   # average per-molecule out rate

    Q_c = _single_voxel_generator(model, in_rate, out_coeff, nmax2; n_sites=2)
    p_post = expv(Float64(dt), Q_c, p_pre; m=min(krylov_m, size(Q_c,1)))
    _clean!(p_post, state_tol)

    # ── 5. Prolong: multi-level stochastic interpolation ─────────────────────
    #
    # Three-level hierarchy (fine → coarse split → total N):
    #   Existing N: exact multiplicative correction using pre-smooth fiber0
    #   New N:      coarse L/R propagation on fiber1 → refine to fiber0
    # See _prolong_multilevel and helpers above for full derivation.

    sp_new = _prolong_multilevel(sp, p_pre, p_post, in1, in2, out1, out2,
                                  model, n_max, state_tol, max_states)
    sp_new === nothing && return nothing
    length(sp_new) > max_states && return nothing

    # ── 6. Post-smooth: intra-pair diffusion ──────────────────────────────────
    τ_post = post_frac * dt
    if τ_post > 0.0 && !isempty(sp_new.states)
        intra = DiscreteStochasticSystem{CartesianIndex{2}}(
            [CartesianIndex(-1,1), CartesianIndex(1,-1)],
            [(x,r,t) -> d*max(0,Tuple(x)[1]),
             (x,r,t) -> d*max(0,Tuple(x)[2])])
        expand!(sp_new, intra, x -> all(c->0≤c≤n_max,Tuple(x)); depth=1)
        A, = build_generator(sp_new, intra, nothing, 0.0)
        sp_new.probs .= expv(τ_post, A, sp_new.probs; m=min(krylov_m,size(A,1)))
        _clean!(sp_new.probs, state_tol)
    end

    renormalize!(sp_new) > 0.0 || return nothing

    # ── 7. Extract updated per-voxel marginals ────────────────────────────────
    p1_new = zeros(n_max + 1); p2_new = zeros(n_max + 1)
    for (ci, prob) in zip(sp_new.states, sp_new.probs)
        t = Tuple(ci)
        t[1] ≤ n_max && (p1_new[t[1]+1] += prob)
        t[2] ≤ n_max && (p2_new[t[2]+1] += prob)
    end
    s1 = sum(p1_new); s1 > 0 && (p1_new ./= s1)
    s2 = sum(p2_new); s2 > 0 && (p2_new ./= s2)

    (p1_new, p2_new)
end

# ─── Three-level multigrid V-cycle for RDME active voxels ───────────────────
#
# Hierarchy:
#   L0  individual voxels    — w states each
#   L1  adjacent pairs       — 2w total-count states   (restrict = convolution)
#   L2  2×2 quads of pairs   — 4w total-count states   (coarsest, CME solve here)
#
# V-cycle (down then up):
#   L0 pre-smooth → restrict L0→L1 → restrict L1→L2
#   → coarse CME solve at L2
#   → multiplicative correction L2→L1 → multiplicative correction L1→L0
#   → L0 post-smooth
#
# The multiplicative correction at each level:
#   p_fine_new(n_a) ∝ p_fine(n_a) × E[C(N_coarse) | n_a]
#                   = p_fine(n_a) × Σ_{n_b} C(n_a+n_b) × p_other(n_b)
# where C(N) = p_coarse_new[N] / p_coarse_old[N].  This is exact for existing N
# values and falls back to identity for new N values created by the coarse solve.

# Convolve two 1D distributions → total-count distribution (restriction operator)
function _mg_convolve(p1::Vector{Float64}, p2::Vector{Float64}, stol::Float64)
    n1 = length(p1)-1;  n2 = length(p2)-1
    out = zeros(n1+n2+1)
    for i in 1:length(p1)
        p1[i] < stol && continue
        for j in 1:length(p2)
            p2[j] < stol && continue
            out[i+j-1] += p1[i]*p2[j]
        end
    end
    s = sum(out); s > 1e-300 && (out ./= s)
    out
end

# Apply multiplicative coarse correction (prolongation operator).
# Given p_this (this group's distribution), p_other (sibling group's distribution),
# p_coarse_old and p_coarse_new (before/after coarse solve of the merged group):
#   p_this_new(N_a) ∝ p_this(N_a) × Σ_{N_b} C(N_a+N_b) × p_other(N_b)
# where C(N) = p_coarse_new[N]/p_coarse_old[N]   (1.0 when old≈0, no info)
function _mg_correct(p_this::Vector{Float64}, p_other::Vector{Float64},
                     p_old::Vector{Float64},  p_new::Vector{Float64},
                     stol::Float64)
    na = length(p_this)-1;  nb = length(p_other)-1
    out = zeros(na+1)
    for Na in 0:na
        p_this[Na+1] < stol && continue
        c = 0.0
        for Nb in 0:nb
            p_other[Nb+1] < stol && continue
            Nc = Na+Nb
            Nc < length(p_old) || continue
            C  = p_old[Nc+1] > 1e-300 ? p_new[Nc+1]/p_old[Nc+1] : 1.0
            c += C * p_other[Nb+1]
        end
        out[Na+1] = p_this[Na+1] * c
    end
    s = sum(out)
    s > 1e-300 ? out./s : copy(p_this)   # fallback: identity if correction collapses
end

# Galerkin prolongation: reconstruct per-voxel marginals from the updated
# pair total-count distribution p_new(N) using the conditional
#   φ_N(n1) ∝ p1_ref(n1) * p2_ref(N-n1)
# This is the Galerkin-optimal prolongation: the coarse variable is N=n1+n2,
# the intra-pair transitions at fixed N are symmetric diffusion hops whose
# stationary is Binomial(N,1/2), and the marginal product p1_ref*p2_ref
# encodes the current spatial distribution within each N-fiber (from the
# pre-smooth step).  Falls back to the uniform (binomial) conditional when
# the marginal product has no support at a given N.
function _prolong_galerkin(p_new  :: Vector{Float64},
                            p1_ref :: Vector{Float64},
                            p2_ref :: Vector{Float64},
                            n_max  :: Int,
                            stol   :: Float64)
    p1_out = zeros(n_max + 1)
    p2_out = zeros(n_max + 1)
    for N in 0:2*n_max
        N + 1 > length(p_new) && continue
        p_new[N+1] < stol && continue
        n1_lo = max(0, N - n_max)
        n1_hi = min(N, n_max)
        # Conditional φ_N(n1) ∝ p1_ref(n1) * p2_ref(N-n1)
        s = 0.0
        for n1 in n1_lo:n1_hi
            s += p1_ref[n1+1] * p2_ref[N-n1+1]
        end
        if s < stol
            # Fallback: uniform over valid N-fiber (≡ Binomial(N,1/2) weight)
            cnt = n1_hi - n1_lo + 1
            cnt <= 0 && continue
            w = p_new[N+1] / cnt
            for n1 in n1_lo:n1_hi
                p1_out[n1+1] += w
                p2_out[N-n1+1] += w
            end
        else
            w = p_new[N+1] / s
            for n1 in n1_lo:n1_hi
                c = w * p1_ref[n1+1] * p2_ref[N-n1+1]
                p1_out[n1+1] += c
                p2_out[N-n1+1] += c
            end
        end
    end
    s1 = sum(p1_out); s2 = sum(p2_out)
    s1 > 1e-300 && (p1_out ./= s1)
    s2 > 1e-300 && (p2_out ./= s2)
    p1_out, p2_out
end

# Build the joint product state p1⊗p2, apply intra-pair diffusion for τ,
# then return (p1_smoothed, p2_smoothed, p_total_count).  Used at ℓ=1 in
# _step_multigrid! so the L1 restriction carries intra-pair correlations.
# Returns nothing if the pair state space is too large.
function _presmooth_pair(
        d         :: Float64,
        n_max     :: Int,
        p1        :: Vector{Float64},
        p2        :: Vector{Float64},
        τ         :: Float64;
        stol      :: Float64 = 1e-8,
        krylov_m  :: Int     = 30,
        max_states :: Int    = 2000)

    supp1 = findall(>(stol), p1); isempty(supp1) && (supp1 = [argmax(p1)])
    supp2 = findall(>(stol), p2); isempty(supp2) && (supp2 = [argmax(p2)])
    length(supp1) * length(supp2) > max_states && return nothing

    sp = StateSpace{CartesianIndex{2}, Float64}()
    for i1 in supp1, i2 in supp2
        w = p1[i1] * p2[i2]
        w > stol || continue
        add_state!(sp, CartesianIndex(i1-1, i2-1), w)
    end
    isempty(sp.states) && return nothing
    renormalize!(sp)

    pair_intra = DiscreteStochasticSystem{CartesianIndex{2}}(
        [CartesianIndex(-1,1), CartesianIndex(1,-1)],
        [(x,r,t) -> d*max(0,Tuple(x)[1]),
         (x,r,t) -> d*max(0,Tuple(x)[2])])
    bc2 = x -> all(c -> 0 ≤ c ≤ n_max, Tuple(x))
    expand!(sp, pair_intra, bc2; depth=1)
    length(sp) > max_states && return nothing

    if τ > 0.0 && !isempty(sp.states)
        A, = build_generator(sp, pair_intra, nothing, 0.0)
        sp.probs .= expv(τ, A, sp.probs; m=min(krylov_m, size(A,1)))
        _clean!(sp.probs, stol)
    end

    p1_new = zeros(n_max + 1); p2_new = zeros(n_max + 1)
    for (ci, prob) in zip(sp.states, sp.probs)
        t = Tuple(ci)
        t[1] ≤ n_max && (p1_new[t[1]+1] += prob)
        t[2] ≤ n_max && (p2_new[t[2]+1] += prob)
    end
    s1 = sum(p1_new); s1 > 0 && (p1_new ./= s1)
    s2 = sum(p2_new); s2 > 0 && (p2_new ./= s2)

    nmax2 = 2 * n_max
    p_total = zeros(nmax2 + 1)
    for (ci, prob) in zip(sp.states, sp.probs)
        N = sum(Tuple(ci)); N <= nmax2 && (p_total[N+1] += prob)
    end
    p_total ./= max(sum(p_total), 1e-100)

    (p1_new, p2_new, p_total)
end

"""
    _step_multigrid!(s, means, dt; krylov_m, row_pass) → handled::Set

Three-level multigrid V-cycle for the RDME active wavefront.

Levels: L0 (per-voxel) → L1 (pairs, total N₂) → L2 (quads, total N₄).
Each level coarsens the spatial scale by 2×; the coarse CME solve at the
deepest level captures collective molecular dynamics and feeds multiplicative
corrections back down to individual voxels.

Pairing direction alternates each level (ADI-style):
  odd  levels → Δ_par  (row-pass: horizontal, col-pass: vertical)
  even levels → Δ_perp (perpendicular)

Between-representative offset at level ℓ: Δ(ℓ+1) × 2^(ℓ÷2)
This gives: L1=pairs, L2=2×2 quads, L3=2×4 blocks, L4=4×4 blocks, …
"""
function _step_multigrid!(s        :: SpatialFSP,
                           means    :: Dict{CartesianIndex{2},Float64},
                           dt       :: Float64;
                           krylov_m :: Int    = 30,
                           row_pass :: Bool   = true,
                           n_levels :: Int    = 3,
                           α_pre    :: Float64 = 0.08,
                           α_post   :: Float64 = 0.08)

    d      = jump_rate(s.model)
    n_max  = s.n_max;  stol = s.block_state_tol
    _n_un  = s.model isa SchloglModel1D ? schlogl_fixed_points(s.model)[2] : -1
    τ_pre  = α_pre  * dt
    τ_mid  = (1.0 - α_pre - α_post) * dt
    τ_post = α_post * dt
    n_lev  = clamp(n_levels, 1, 8)

    active_keys = collect(keys(s.dists))
    isempty(active_keys) && return Set{CartesianIndex{2}}()
    active_set  = Set(active_keys)
    equil_set   = Set(keys(s.equil_dists))

    # Classify active voxels:
    #   interior — all neighbours active or equil of the SAME basin (safe to coarsen)
    #   edge     — ≥1 empty neighbour, OR (Schlögl) adjacent to equil of opposite basin
    #
    # For bistable Schlögl the wavefront voxels sit between equil-hi and equil-lo.
    # Total-count restriction loses the bimodal structure for these voxels, so they
    # are kept at full per-voxel resolution (edge) rather than coarsened (interior).
    n_un_mg = _n_un > 0 ? Float64(_n_un) : -1.0
    interior = CartesianIndex{2}[]
    edge      = CartesianIndex{2}[]
    for idx in active_keys
        ns = grid_neighbors(idx, s.K_x, s.K_y)
        is_edge = if n_un_mg > 0.0
            mean_idx = get(means, idx, 0.0)
            any(ns) do nb
                (!haskey(s.dists, nb) && !haskey(s.equil_dists, nb)) ||
                (haskey(s.equil_dists, nb) &&
                 xor(mean_idx >= n_un_mg, get(means, nb, 0.0) >= n_un_mg))
            end
        else
            !all(nb -> nb in active_set || nb in equil_set, ns)
        end
        push!(is_edge ? edge : interior, idx)
    end
    interior_set = Set(interior)

    _n_hi_mg = s.model isa SchloglModel1D ? schlogl_fixed_points(s.model)[3] : -1
    # Per-voxel expv helper — support must span n_hi for Schlögl (same reason as step!)
    function smooth!(idx, τ)
        τ ≤ 0.0 && return
        ns  = grid_neighbors(idx, s.K_x, s.K_y)
        μin = sum(get(means,nb,0.0) for nb in ns; init=0.0)
        p   = _expv_voxel(s.model, d*μin, d*length(ns), s.dists[idx],
                           τ, krylov_m; n_include = _n_hi_mg > 0 ? _n_hi_mg : _n_un)
        _clean!(p, stol); s.dists[idx] = p
    end

    # ── Pre-smooth (all active) ───────────────────────────────────────────────
    for idx in active_keys; smooth!(idx, τ_pre); end

    # ── ADI direction vectors ─────────────────────────────────────────────────
    Δ_par  = row_pass ? CartesianIndex(0,1) : CartesianIndex(1,0)
    Δ_perp = row_pass ? CartesianIndex(1,0) : CartesianIndex(0,1)

    # Edge voxels are handled by pre/post smooth; multigrid covers interior only
    frontier_handled = Set{CartesianIndex{2}}()
    for idx in edge; smooth!(idx, τ_mid); end

    interior_mg = interior

    # ── Build hierarchy over interior voxels only ─────────────────────────────
    # Alternating pairing directions; at level ℓ the expected rep offset is
    # Δ(ℓ+1) × 2^(ℓ÷2) where Δ(ℓ) = Δ_par for odd ℓ, Δ_perp for even ℓ.
    Δdir(ℓ) = isodd(ℓ) ? Δ_par : Δ_perp
    function pair_delta(ℓ)   # offset between reps when pairing level-ℓ groups
        stride = 2^((ℓ-1)÷2)
        CartesianIndex(Tuple(Δdir(ℓ+1)) .* stride)
    end

    # level_reps[ℓ]   — representative voxel (min CI) for each group
    # level_voxels[ℓ] — all voxels in each group (for external-rate computation)
    # level_dists[ℓ]  — total-count distribution for each group
    # level_pairs[ℓ]  — (i,j) child-index pairs that form level ℓ+1

    level_reps   = [interior_mg]                                 # level_reps[ℓ][g]   = rep of group g at level ℓ
    level_voxels = [[[v] for v in interior_mg]]                  # level_voxels[ℓ][g] = voxels in group g at level ℓ
    level_dists  = [[copy(s.dists[v]) for v in interior_mg]]     # level_dists[ℓ][g]  = total-count dist of group g
    level_pairs  = Vector{Tuple{Int,Int}}[]                 # level_pairs[ℓ]     = (i,j) forming level ℓ+1

    for ℓ in 1:n_lev
        prev_reps = level_reps[end]
        isempty(prev_reps) && break
        rep_idx   = Dict(r=>i for (i,r) in enumerate(prev_reps))
        δ         = pair_delta(ℓ)

        used = Set{Int}()
        pairs_ℓ   = Tuple{Int,Int}[]
        new_reps   = CartesianIndex{2}[]
        new_vox    = Vector{Vector{CartesianIndex{2}}}()
        new_dists  = Vector{Float64}[]

        for (i,rep) in enumerate(prev_reps)
            i in used && continue
            rep2 = rep + δ
            j    = get(rep_idx, rep2, 0)
            (j==0 || j in used) && continue
            push!(pairs_ℓ,   (i,j))
            push!(used, i); push!(used, j)
            push!(new_reps,  min(rep, rep2))
            push!(new_vox,   vcat(level_voxels[end][i], level_voxels[end][j]))
            # At every level: presmooth with effective inter-group coupling.
            # k_ℓ = 1/2^((ℓ-1)÷2): fraction of boundary voxels per child group,
            # so each particle in a child group exits to the sibling at rate d*k_ℓ.
            # At ℓ=1: k=1 (single voxels, direct edge), identical to old code.
            # At ℓ=2: k=1 (2-voxel pairs, full column boundary with sibling pair).
            # At ℓ=3+: k<1 as interior voxels accumulate in each child group.
            k_ℓ     = 1.0 / 2.0^((ℓ-1)÷2)
            d_eff   = d * k_ℓ
            n_max_ℓ = length(level_dists[end][i]) - 1
            ps = _presmooth_pair(d_eff, n_max_ℓ,
                                 level_dists[end][i], level_dists[end][j],
                                 τ_pre; stol, krylov_m)
            if ps !== nothing
                p_i_s, p_j_s, p_total = ps
                level_dists[end][i] = p_i_s   # update reference for Galerkin prolongation
                level_dists[end][j] = p_j_s
                push!(new_dists, p_total)
            else
                push!(new_dists, _mg_convolve(level_dists[end][i],
                                               level_dists[end][j], stol))
            end
        end

        isempty(pairs_ℓ) && break
        push!(level_pairs,  pairs_ℓ)
        push!(level_reps,   new_reps)
        push!(level_voxels, new_vox)
        push!(level_dists,  new_dists)
    end

    n_built = length(level_pairs)   # actual levels built (≤ n_lev)

    # ── Coarse solve at deepest level ─────────────────────────────────────────
    if n_built > 0
        old_coarse = copy.(level_dists[end])
        new_coarse = similar.(old_coarse)
        n_vox_coarse = length(level_voxels[end][1])

        for (g, vxs) in enumerate(level_voxels[end])
            nv  = length(vxs)
            nmg = nv * n_max
            in_g = 0.0; n_ext = 0
            for vi in vxs, nb in grid_neighbors(vi, s.K_x, s.K_y)
                nb in vxs && continue
                in_g += d * get(means,nb,0.0); n_ext += 1
            end
            Q_g = _single_voxel_generator(s.model, in_g, d*n_ext/nv, nmg; n_sites=nv)
            p_g = expv(Float64(τ_mid), Q_g, level_dists[end][g];
                        m=min(krylov_m, size(Q_g,1)))
            _clean!(p_g, stol); new_coarse[g] = p_g
        end

        # ── Correction phase (deep → shallow) ─────────────────────────────────
        # current_corr[ℓ] = corrected distributions at level ℓ (1-indexed, 1=voxels)
        current_corr = [copy.(d) for d in level_dists]
        current_corr[end] = new_coarse   # deepest level = coarse solve output

        # Full Galerkin cascade: deepest level → voxels.
        # Before each Galerkin step at level ℓ, groups at level ℓ+1 that are
        # NOT covered by any deeper-level solve (orphaned groups) receive a
        # direct coarse solve at their own level so they also get τ_mid dynamics.
        for ℓ in n_built:-1:1
            pairs_ℓ = level_pairs[ℓ]
            old_ℓ   = level_dists[ℓ]
            new_c   = current_corr[ℓ+1]
            n_max_ℓ = length(old_ℓ[1]) - 1

            # Orphan coarse solve: groups at level ℓ+1 not in any level-(ℓ+2) group
            if ℓ < n_built
                covered_above = Set{Int}()
                for (i,j) in level_pairs[ℓ+1]
                    push!(covered_above, i); push!(covered_above, j)
                end
                for k in 1:length(new_c)
                    k in covered_above && continue
                    vxs = level_voxels[ℓ+1][k]
                    nv  = length(vxs);  nmg = nv * n_max
                    ig  = 0.0;          ne  = 0
                    for vi in vxs, nb in grid_neighbors(vi, s.K_x, s.K_y)
                        nb in vxs && continue
                        ig += d * get(means,nb,0.0);  ne += 1
                    end
                    Q_k = _single_voxel_generator(s.model, ig, d*ne/nv, nmg; n_sites=nv)
                    p_k = expv(Float64(τ_mid), Q_k, level_dists[ℓ+1][k];
                                m=min(krylov_m, size(Q_k,1)))
                    _clean!(p_k, stol);  current_corr[ℓ+1][k] = p_k
                end
            end

            corrected = copy(old_ℓ)
            for (g,(i,j)) in enumerate(pairs_ℓ)
                corrected[i], corrected[j] =
                    _prolong_galerkin(new_c[g], old_ℓ[i], old_ℓ[j], n_max_ℓ, stol)
            end
            current_corr[ℓ] = corrected
        end

        # Apply corrected voxel distributions (current_corr[1]) to s.dists
        handled = union(frontier_handled, Set{CartesianIndex{2}}(edge))
        for (g,(i,j)) in enumerate(level_pairs[1])
            vi = interior_mg[i]; vj = interior_mg[j]
            p_i = current_corr[1][i]; p_j = current_corr[1][j]
            _clean!(p_i, stol); _clean!(p_j, stol)
            si = sum(p_i); sj = sum(p_j)
            si > 1e-300 && (s.dists[vi] = p_i./si; push!(handled, vi))
            sj > 1e-300 && (s.dists[vj] = p_j./sj; push!(handled, vj))
        end

        # Interior voxels not in any pair: per-voxel expv for τ_mid
        for idx in interior
            idx in handled && continue
            ns  = grid_neighbors(idx, s.K_x, s.K_y)
            μin = sum(get(means,nb,0.0) for nb in ns; init=0.0)
            p   = _expv_voxel(s.model, d*μin, d*length(ns), s.dists[idx],
                               τ_mid, krylov_m; n_include=_n_un)
            _clean!(p, stol); s.dists[idx] = p; push!(handled, idx)
        end
    else
        # No interior pairs formed: per-voxel for remaining interior voxels
        handled = union(frontier_handled, Set{CartesianIndex{2}}(edge))
        for idx in interior_mg
            ns  = grid_neighbors(idx, s.K_x, s.K_y)
            μin = sum(get(means,nb,0.0) for nb in ns; init=0.0)
            p   = _expv_voxel(s.model, d*μin, d*length(ns), s.dists[idx],
                               τ_mid, krylov_m; n_include=_n_un)
            _clean!(p, stol); s.dists[idx] = p; push!(handled, idx)
        end
    end

    # ── Post-smooth (all active) ──────────────────────────────────────────────
    for idx in active_keys; smooth!(idx, τ_post); end

    handled
end

"""
    _step_adi_pair!(s, means, dt; krylov_m, row_pass)

ADI pair V-cycle: sweep all horizontal pairs (row_pass=true) or all vertical
pairs (row_pass=false) of active voxels using `_evolve_pair_vcycle`.

Alternating between row_pass=true and row_pass=false each step captures all
4 nearest-neighbor 2D correlations at 2× the cost of one direction.
"""
function _step_adi_pair!(s::SpatialFSP,
                          means     :: Dict{CartesianIndex{2}, Float64},
                          dt        :: Float64;
                          krylov_m  :: Int  = 30,
                          row_pass  :: Bool = true)
    n_max     = s.n_max
    state_tol = s.block_state_tol
    max_s     = s.block_max_states

    # Direction: row-adjacent (Δi=0, Δj=1) or column-adjacent (Δi=1, Δj=0)
    Δ = row_pass ? CartesianIndex(0,1) : CartesianIndex(1,0)

    handled = Set{CartesianIndex{2}}()
    for v1 in keys(s.dists)
        v1 ∈ handled && continue
        v2 = v1 + Δ
        haskey(s.dists, v2) || continue
        v2 ∈ handled && continue

        result = _evolve_pair_vcycle(
            s.model, v1, v2,
            s.dists[v1], s.dists[v2],
            means, s.K_x, s.K_y, dt, n_max;
            krylov_m, state_tol, max_states = max_s,
            pre_frac  = s.block_pre_frac,
            post_frac = s.block_post_frac)

        if result !== nothing
            s.dists[v1], s.dists[v2] = result
            push!(handled, v1); push!(handled, v2)
        end
    end
    handled   # return set of voxels processed; unmatched fall through to per-voxel
end

# ─── AMG block: no mean-field within 2×2 active + frontier ──────────────────
#
# All other methods approximate the coupling from a neighbour voxel as d·μ_nb
# (its mean only).  Here, the 2×2 active block plus every empty/equil neighbour
# of that block are all placed into a single joint CME.  Coupling between any
# two voxels that are BOTH in this extended set uses the exact rate d·n_i.
# Mean-field is applied only at the outer boundary (outside the extended set).
#
# State space is built by rstep from the active product seed; the frontier
# voxels start at n=0 in the seed and their n>0 states are generated organically
# by the rstep expansion through seeding reactions.  No (n_max+1)^N enumeration.
#
# After the solve, any frontier voxel whose marginal has P(n≥1) > amg_activate_tol
# is immediately promoted to active (bypassing the pending_flux threshold).
# This is the "prolongation that preserves pre-restricted state": the voxel
# starts its active career with the physically-correct solved distribution.

function _evolve_amg_block(
    model,
    all_voxels  :: Vector{CartesianIndex{2}},
    all_dists   :: Vector{Vector{Float64}},
    outer_means :: Dict{CartesianIndex{2}, Float64},
    K_x :: Int, K_y :: Int,
    dt  :: Float64,
    n_max :: Int;
    seed_k      :: Int     = 2,
    rstep_depth :: Int     = 2,
    max_states  :: Int     = 30_000,
    krylov_m    :: Int     = 30,
    state_tol   :: Float64 = 1e-10
)
    N        = length(all_voxels)
    N_active = 4

    # Exact joint CME for all_voxels; mean-field only at outer boundary
    sys = _build_joint_active_system(model, all_voxels, outer_means, K_x, K_y)
    bc  = x -> all(0 ≤ c ≤ n_max for c in Tuple(x))

    # Seed: product of active top-seed_k states; frontier fixed at n=0
    active_supps = [let p = all_dists[k]
                        nz = count(>(state_tol), p)
                        partialsortperm(p, 1:max(1, min(seed_k, nz)); rev=true)
                    end for k in 1:N_active]

    sp = StateSpace{CartesianIndex{N}, Float64}()
    for combo in Iterators.product(active_supps...)
        n_a = Tuple(c - 1 for c in combo)
        n_f = ntuple(_ -> 0, N - N_active)
        ci  = CartesianIndex(n_a..., n_f...)
        w   = prod(k -> all_dists[k][combo[k]], 1:N_active)
        w > state_tol && add_state!(sp, ci, w)
    end
    isempty(sp.states) && return nothing

    # rstep expansion: frontier n>0 states generated by seeding reactions
    expand!(sp, sys, bc; depth = rstep_depth)
    length(sp) > max_states && return nothing

    # Initialise from product distribution.
    # Frontier states at n>0 get p=0 (correct: frontier starts empty).
    # expv will populate them via seeding transitions.
    for i in eachindex(sp.states)
        t = Tuple(sp.states[i])
        sp.probs[i] = prod(k -> begin
            nk = t[k];  nk < length(all_dists[k]) ? all_dists[k][nk + 1] : 0.0
        end, 1:N)
    end
    s = sum(sp.probs);  s ≤ 0.0 && return nothing
    sp.probs ./= s

    A, = build_generator(sp, sys, nothing, 0.0)
    sp.probs .= expv(Float64(dt), A, sp.probs; m = min(krylov_m, size(A, 1)))
    _clean!(sp.probs, state_tol)
    renormalize!(sp) > 0.0 || return nothing

    # Extract per-voxel marginals
    marginals = [zeros(n_max + 1) for _ in 1:N]
    for i in eachindex(sp.states)
        t = Tuple(sp.states[i]);  pr = sp.probs[i]
        for k in 1:N
            nk = t[k];  nk ≤ n_max && (marginals[k][nk + 1] += pr)
        end
    end
    for m in marginals
        s = sum(m);  s > 0.0 && (m ./= s)
    end
    marginals
end

function _step_amg_block!(
    s        :: SpatialFSP,
    means    :: Dict{CartesianIndex{2}, Float64},
    dt       :: Float64;
    krylov_m :: Int = 30,
    step_idx :: Int = 0
)
    pass  = step_idx % 4
    i_off = pass ∈ (1, 3) ? 1 : 0
    j_off = pass ∈ (1, 2) ? 1 : 0
    preferred  = [CartesianIndex(i, j)
                  for i in (1+i_off):2:(s.K_x-1) for j in (1+j_off):2:(s.K_y-1)]
    anchors    = _select_active_blocks(s.dists, s.K_x, s.K_y, preferred)
    active_set = Set(keys(s.dists))

    handled = Set{CartesianIndex{2}}()

    for anchor in anchors
        block = _block_voxels(anchor)
        all(b -> b ∈ active_set, block)  || continue
        any(b -> b ∈ handled, block)     && continue

        # Frontier: EMPTY voxels adjacent to block only (not active, not equil).
        # Equil voxels are handled via mean-field in outer_means (step 3 evolves
        # them separately; including them here would re-activate saturated interior).
        seen     = Set{CartesianIndex{2}}(block)
        frontier = CartesianIndex{2}[]
        for v in block
            for nb in grid_neighbors(v, s.K_x, s.K_y)
                nb ∈ seen                  && continue
                nb ∈ active_set            && continue   # active → mean-field
                haskey(s.equil_dists, nb)  && continue   # equil  → mean-field
                1 ≤ nb[1] ≤ s.K_x && 1 ≤ nb[2] ≤ s.K_y || continue
                push!(seen, nb);  push!(frontier, nb)
            end
        end

        all_voxels = [block..., frontier...]

        p_block    = [s.dists[b] for b in block]
        p_frontier = [vcat(1.0, zeros(s.n_max)) for _ in frontier]  # all empty: δ(n=0)
        all_dists  = [p_block..., p_frontier...]

        # outer_means: mean-field for everything outside all_voxels
        # (other active, equil, and far-empty voxels)
        all_set     = Set(all_voxels)
        outer_means = Dict{CartesianIndex{2}, Float64}()
        for v in all_voxels
            for nb in grid_neighbors(v, s.K_x, s.K_y)
                nb ∈ all_set && continue
                outer_means[nb] = get(means, nb, 0.0)
            end
        end

        result = _evolve_amg_block(
            s.model, all_voxels, all_dists, outer_means,
            s.K_x, s.K_y, dt, s.n_max;
            seed_k      = 2,
            rstep_depth = 2,
            max_states  = s.amg_block_max_states,
            krylov_m,
            state_tol   = s.block_state_tol)
        result === nothing && continue

        # Update active block voxels
        for (k, v) in enumerate(block)
            m = result[k];  sm = sum(m)
            sm > 0 || continue
            s.dists[v] = m ./ sm
            push!(handled, v)
        end

        # Activate empty frontier voxels where P(n≥1) > amg_activate_tol.
        # Frontier is EMPTY only (equil excluded above), so no equil_dists cleanup.
        for (k, v) in enumerate(frontier)
            m  = result[4 + k]
            sm = sum(m);  sm > 0 || continue
            m ./= sm
            1.0 - m[1] > s.amg_activate_tol || continue
            haskey(s.dists, v) && continue   # edge case: activated by another block
            s.dists[v] = m
            delete!(s.pending_flux, v)
        end
    end

    handled
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

    row_pass = isodd(round(Int, s.t / dt))

    if isempty(handled) && s.use_multigrid
        # Generalized multigrid V-cycle: n_levels spatial hierarchy.
        mg_handled = _step_multigrid!(s, means, dt; krylov_m, row_pass,
                                       n_levels = s.multigrid_levels)
        for idx in mg_handled; new_dists[idx] = s.dists[idx]; end
        union!(handled, mg_handled)
        empty!(s.block_states)
    end

    if isempty(handled) && s.use_block_vcycle_4body
        # 2×2 block V-cycle: compact rstep restriction + 4-body joint CME + correction prolongation.
        step_idx = round(Int, s.t / dt)
        vcb_handled = _step_block_vcycle_4body!(s, means, dt; krylov_m, step_idx)
        for idx in vcb_handled; new_dists[idx] = s.dists[idx]; end
        union!(handled, vcb_handled)
        empty!(s.block_states)
    end

    if isempty(handled) && s.use_amg_block
        # AMG block: no mean-field within 2×2 active + frontier joint state.
        # Frontier (empty/equil neighbours) included in joint CME via rstep.
        # Frontier voxels with P(n≥1) > amg_activate_tol promoted to active.
        step_idx    = round(Int, s.t / dt)
        amg_handled = _step_amg_block!(s, means, dt; krylov_m, step_idx)
        for idx in amg_handled; new_dists[idx] = s.dists[idx]; end
        union!(handled, amg_handled)
        empty!(s.block_states)
    end

    if isempty(handled) && s.use_block_exact
        # Exact-Galerkin pair joint CME: evolve full 2-body generator,
        # no total-count restriction.  ADI alternating row/column passes.
        exact_handled = _step_pair_exact!(s, means, dt; krylov_m, row_pass)
        for idx in exact_handled; new_dists[idx] = s.dists[idx]; end
        union!(handled, exact_handled)
        empty!(s.block_states)
    end

    if isempty(handled) && s.use_block_vcycle
        # ADI pair V-cycle: K=2 pairs, alternating row/column each step.
        pair_handled = _step_adi_pair!(s, means, dt; krylov_m, row_pass)
        for idx in pair_handled
            new_dists[idx] = s.dists[idx]
        end
        union!(handled, pair_handled)
        empty!(s.block_states)
    end

    # Schlögl: support must span [n_lo, n_hi] so strong diffusion inflow from an
    # equil-hi neighbour can push probability all the way to n_hi in one step.
    # Additionally, use time-splitting when λ·dt > krylov_m (the Krylov dimension
    # needed for expv scales with λ·dt; exceeding krylov_m causes mass loss).
    # Each sub-step keeps λ·dt_sub ≤ krylov_m so m=30 is sufficient.
    _n_un = s.model isa SchloglModel1D ? schlogl_fixed_points(s.model)[2] : -1
    _n_hi_sp = s.model isa SchloglModel1D ? schlogl_fixed_points(s.model)[3] : -1
    for (idx, p_old) in s.dists
        idx in handled && continue
        ns        = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in      = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        in_rate   = d * μ_in
        out_coeff = d * length(ns)
        n_inc     = _n_hi_sp > 0 ? _n_hi_sp : _n_un
        # Time-splitting: split into sub-steps so λ·dt_sub ≤ krylov_m
        max_rate  = in_rate + (_n_hi_sp > 0 ? out_coeff * _n_hi_sp : 0.0)
        n_sub     = max_rate > 0.0 ? max(1, ceil(Int, max_rate * Float64(dt) / krylov_m)) : 1
        dt_sub    = Float64(dt) / n_sub
        p_cur = p_old
        for _ in 1:n_sub
            p_cur = _expv_voxel(s.model, in_rate, out_coeff, p_cur, dt_sub, krylov_m;
                                n_include = n_inc)
            _clean!(p_cur, s.ε_prune)
        end
        new_dists[idx] = p_cur
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
                p_old_eq = _resize_dist(p_eq, s.n_max_eq)
                p_new_eq = _expv_voxel(s.model, in_rate, out_coeff, p_old_eq,
                                        Float64(dt), s.krylov_m_eq)
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
        # BranchingDeathRDME: boundary-only expv.
        # SchloglModel1D: equil voxels are FROZEN — they stay at their current basin
        #   (delta at n_lo or n_hi) until cross-basin flux from an active wave front
        #   re-activates them in step 5. This is the exact bistable analogue of the
        #   BD "empty voxel": both are background states that only become active when
        #   a wave front crosses them. No ad-hoc step needed — the pending_flux
        #   mechanism (steps 4–5) handles both cases uniformly.
        if !(s.model isa SchloglModel1D)
            new_equil = Dict{CartesianIndex{2}, Vector{Float64}}()
            for (idx, p_eq) in s.equil_dists
                ns = grid_neighbors(idx, s.K_x, s.K_y)
                if any(nb -> haskey(s.dists, nb), ns)
                    μ_in      = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
                    in_rate   = d * μ_in
                    out_coeff = d * length(ns)
                    p_old_eq  = _resize_dist(p_eq, s.n_max_eq)
                    p_new_eq  = _expv_voxel(s.model, in_rate, out_coeff, p_old_eq,
                                            Float64(dt), s.krylov_m_eq)
                    _clean!(p_new_eq, s.ε_prune)
                    new_equil[idx] = p_new_eq
                else
                    new_equil[idx] = p_eq
                end
            end
            s.equil_dists = new_equil
        end
        # SchloglModel1D: s.equil_dists unchanged (frozen background)
    end

    # ── 3a. Schlögl: early deactivation (BEFORE cross-basin flux) ────────────
    # A voxel deactivates and snaps to its nearest stable basin when either:
    #   (a) p(n_un) < ε_equil  [distribution has committed to one basin], or
    #   (b) |mean - n_un| > snap_margin  [mean-field fixed-point criterion]
    #
    # Criterion (b) resolves a mean-field artefact: leading-edge voxels with
    # predominantly equil-lo neighbours develop a spurious stable fixed point
    # near n_un (e.g. mean ≈ 68 for 1 equil-hi + 3 equil-lo neighbours with
    # D=2). The true joint CME is bimodal, but the per-voxel mean-field
    # approximation creates a single attractor just above n_un, so p(n_un)
    # never drops below ε_equil and the wavefront stalls.  Snapping when the
    # mean has moved >snap_margin from n_un resolves the trapping: after one
    # expv step from δ(n_lo) the mean already overshoots n_un (e.g. ≈69.5),
    # so the voxel immediately commits to equil-hi and drives the next ring.
    #
    # The snap stores δ(n_hi) or δ(n_lo) — consistent with the frozen-
    # background assumption — rather than the raw evolved distribution.
    if s.model isa SchloglModel1D
        n_lo_e, n_un_e, n_hi_e = schlogl_fixed_points(s.model)
        ε_out      = s.ε_equil
        snap_margin = 4.0   # commit when mean is >4 away from n_un
        to_deact_early = CartesianIndex{2}[]
        for (idx, p) in s.dists
            nv = length(p)
            p_at_nun = n_un_e + 1 <= nv ? p[n_un_e + 1] : 0.0
            mean_p   = voxel_mean(p)
            deactivate = p_at_nun < ε_out || abs(mean_p - Float64(n_un_e)) > snap_margin
            deactivate || continue
            push!(to_deact_early, idx)
            basin_n = mean_p >= Float64(n_un_e) ? n_hi_e : n_lo_e
            p_snap  = zeros(s.n_max_eq + 1);  p_snap[basin_n + 1] = 1.0
            s.equil_dists[idx] = p_snap
        end
        for idx in to_deact_early; delete!(s.dists, idx); end
    end

    # ── 4. Accumulate flux into target voxels ────────────────────────────────
    # BD:      total flux  D·μ·dt  from any source → empty targets (μ=0 baseline).
    # Schlögl: cross-basin flux from EQUIL (settled) sources only.
    #   Only equil voxels that have fully committed to one basin generate flux —
    #   active (transitioning) voxels don't. This is the bistable analog of the
    #   BD rule "only settled voxels radiate flux into empty neighbors."
    #   Once an equil-hi voxel is deactivated, it radiates cross-basin flux to
    #   adjacent equil-lo voxels, driving the next layer of wave propagation.
    if s.model isa SchloglModel1D
        _, n_un_schlogl, _ = schlogl_fixed_points(s.model)
        n_un_f = Float64(n_un_schlogl)
        # Sources: equil voxels only (settled, committed to one basin)
        for (idx, p_eq) in s.equil_dists
            μ = get(means, idx, voxel_mean(p_eq))
            μ > 0.0 || continue
            for nb in grid_neighbors(idx, s.K_x, s.K_y)
                haskey(s.dists, nb)       && continue  # active targets excluded
                haskey(s.equil_dists, nb) || continue  # empty (Schlögl has none)
                μ_nb  = get(means, nb, -1.0)
                # Only equil-hi drives equil-lo: hi wave propagates outward.
                # Equil-lo must NOT generate flux back to equil-hi — that would
                # re-activate stable hi voxels and cause the hi region to collapse.
                cross = (μ_nb >= 0.0 && μ_nb < n_un_f && μ > n_un_f) ?
                            d * (μ - n_un_f) * dt : 0.0
                cross > 0.0 && (s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + cross)
            end
        end
    else
        for (idx, μ) in means
            μ > 0.0 || continue
            for nb in grid_neighbors(idx, s.K_x, s.K_y)
                haskey(s.dists, nb)       && continue
                haskey(s.equil_dists, nb) && continue
                s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + d * μ * dt
            end
        end
    end

    # ── 5. Activate (prolong) target voxels ──────────────────────────────────
    # Empty target:  start from delta(0).
    # Equil target (Schlögl bistable): start from the stored basin distribution.
    for (nb, fl) in s.pending_flux
        if fl >= s.ε_expand && !haskey(s.dists, nb)
            if haskey(s.equil_dists, nb)
                s.dists[nb] = _resize_dist(s.equil_dists[nb], s.n_max)
                delete!(s.equil_dists, nb)
            else
                p0 = zeros(s.n_max + 1);  p0[1] = 1.0
                s.dists[nb] = p0
            end
            delete!(s.pending_flux, nb)
        end
    end

    # ── 5b. Organic activation (optional) ────────────────────────────────────
    # Activate empty neighbor nb immediately when any adjacent tracked voxel k
    # has P(n_k > 0) = 1 - p_k[1] > ε_organic.  Bypasses pending_flux entirely,
    # eliminating the activation delay that causes wave speed error.
    if s.organic_activation
        for (idx, p) in s.dists
            (1.0 - p[1]) > s.ε_organic || continue
            for nb in grid_neighbors(idx, s.K_x, s.K_y)
                haskey(s.dists, nb) && continue
                if haskey(s.equil_dists, nb)
                    s.dists[nb] = _resize_dist(s.equil_dists[nb], s.n_max)
                    delete!(s.equil_dists, nb)
                else
                    p0 = zeros(s.n_max + 1);  p0[1] = 1.0
                    s.dists[nb] = p0
                end
                delete!(s.pending_flux, nb)
            end
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
    else
        to_deactivate = CartesianIndex{2}[]
        if s.model isa SchloglModel1D
            # Schlögl: deactivation was handled in step 3a (early deactivation).
            # Nothing to do here — all converged voxels already moved to equil_dists.
        else
            # BranchingDeathRDME: saturation-based equil (mean > 50% of n_max)
            sat_thresh = 0.50 * Float64(s.n_max)
            for (idx, p) in s.dists
                voxel_mean(p) >= sat_thresh || continue
                ns = grid_neighbors(idx, s.K_x, s.K_y)
                any(nb -> !haskey(s.dists,nb) && !haskey(s.equil_dists,nb), ns) && continue
                push!(to_deactivate, idx)
                s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
            end
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

"""
    step!(s, dt, fg; krylov_m=30, flux_max_depth=6)

Flux-vcycle variant of step!: replaces the active-voxel solve with
`evolve_flux_vcycle!` (adjacency-graph adaptive coarsening), then runs the
same equilibration/activation/restriction lifecycle as the standard `step!`.

Falls back to per-voxel CME for any active voxel not handled by the vcycle.
"""
function step!(s::SpatialFSP, dt::Float64, fg::Any;
               krylov_m::Int = 30, flux_max_depth::Int = 6)
    d = jump_rate(s.model)

    # ── 1. Means ─────────────────────────────────────────────────────────────
    means = Dict{CartesianIndex{2}, Float64}()
    for (idx, p) in s.dists;       means[idx] = voxel_mean(p); end
    for (idx, p) in s.equil_dists; means[idx] = voxel_mean(p); end

    # ── 2. Evolve active voxels via flux V-cycle ──────────────────────────────
    handled = Set{CartesianIndex{2}}()
    if !isempty(s.dists)
        ok = evolve_flux_vcycle!(s, fg, means, dt;
                                 krylov_m = krylov_m, max_depth = flux_max_depth)
        ok && union!(handled, keys(s.dists))
    end

    # fallback for any voxel the vcycle didn't cover
    new_dists = copy(s.dists)
    for (idx, p_old) in s.dists
        idx in handled && continue
        ns        = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in      = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        in_rate   = d * μ_in
        out_coeff = d * length(ns)
        _n_un2 = s.model isa SchloglModel1D ? schlogl_fixed_points(s.model)[2] : -1
        p_new     = _expv_voxel(s.model, in_rate, out_coeff, p_old, Float64(dt), krylov_m;
                               n_include = _n_un2)
        _clean!(p_new, s.ε_prune)
        new_dists[idx] = p_new
    end
    s.dists = new_dists
    s.joint_state = nothing; empty!(s.joint_voxels); empty!(s.block_states)

    # ── 3. Update equilibrated voxels (identical to standard step!) ───────────
    if has_finite_local_ss(s.model)
        new_equil = Dict{CartesianIndex{2}, Vector{Float64}}()
        for (idx, p_eq) in s.equil_dists
            ns        = grid_neighbors(idx, s.K_x, s.K_y)
            μ_in      = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
            in_rate   = d * μ_in
            out_coeff = d * length(ns)
            if any(nb -> haskey(s.dists, nb), ns)
                Q_eq     = _single_voxel_generator(s.model, in_rate, out_coeff, s.n_max_eq)
                p_old_eq = _resize_dist(p_eq, s.n_max_eq)
                p_new_eq = expv(Float64(dt), Q_eq, p_old_eq;
                                m = min(s.krylov_m_eq, size(Q_eq, 1)))
                _clean!(p_new_eq, s.ε_prune)
                new_equil[idx] = p_new_eq
            else
                μ_loc = (s.model.k_b + in_rate) / (s.model.k_d + out_coeff)
                new_equil[idx] = _poisson_pmf(μ_loc, s.n_max_eq)
            end
        end
        s.equil_dists = new_equil
    else
        # BranchingDeathRDME: boundary-only expv.
        # SchloglModel1D: frozen (see main step! for rationale).
        if !(s.model isa SchloglModel1D)
            new_equil = Dict{CartesianIndex{2}, Vector{Float64}}()
            for (idx, p_eq) in s.equil_dists
                ns = grid_neighbors(idx, s.K_x, s.K_y)
                if any(nb -> haskey(s.dists, nb), ns)
                    μ_in      = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
                    in_rate   = d * μ_in
                    out_coeff = d * length(ns)
                    p_old_eq  = _resize_dist(p_eq, s.n_max_eq)
                    p_new_eq  = _expv_voxel(s.model, in_rate, out_coeff, p_old_eq,
                                            Float64(dt), s.krylov_m_eq)
                    _clean!(p_new_eq, s.ε_prune)
                    new_equil[idx] = p_new_eq
                else
                    new_equil[idx] = p_eq
                end
            end
            s.equil_dists = new_equil
        end
    end

    # ── 3a. Schlögl: early deactivation (BEFORE cross-basin flux) ────────────
    if s.model isa SchloglModel1D
        _, n_un_e2, _ = schlogl_fixed_points(s.model)
        ε_out2 = s.ε_equil
        to_deact_early2 = CartesianIndex{2}[]
        for (idx, p) in s.dists
            nv = length(p)
            p_at_nun = n_un_e2 + 1 <= nv ? p[n_un_e2 + 1] : 0.0
            p_at_nun < ε_out2 || continue
            push!(to_deact_early2, idx)
            s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
        end
        for idx in to_deact_early2; delete!(s.dists, idx); end
    end

    # ── 4. Accumulate flux into target voxels ────────────────────────────────
    if s.model isa SchloglModel1D
        _, n_un_schlogl, _ = schlogl_fixed_points(s.model)
        n_un_f = Float64(n_un_schlogl)
        for (idx, p_eq) in s.equil_dists
            μ = get(means, idx, voxel_mean(p_eq))
            μ > 0.0 || continue
            for nb in grid_neighbors(idx, s.K_x, s.K_y)
                haskey(s.dists, nb)       && continue
                haskey(s.equil_dists, nb) || continue
                μ_nb  = get(means, nb, -1.0)
                cross = if μ_nb >= 0.0 && μ_nb < n_un_f && μ > n_un_f
                    d * (μ - n_un_f) * dt
                elseif μ_nb >= n_un_f && μ < n_un_f
                    d * (n_un_f - μ) * dt
                else
                    0.0
                end
                cross > 0.0 && (s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + cross)
            end
        end
    else
        for (idx, μ) in means
            μ > 0.0 || continue
            for nb in grid_neighbors(idx, s.K_x, s.K_y)
                haskey(s.dists, nb)       && continue
                haskey(s.equil_dists, nb) && continue
                s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + d * μ * dt
            end
        end
    end

    # ── 5. Activate target voxels ─────────────────────────────────────────────
    for (nb, fl) in s.pending_flux
        if fl >= s.ε_expand && !haskey(s.dists, nb)
            if haskey(s.equil_dists, nb)
                s.dists[nb] = _resize_dist(s.equil_dists[nb], s.n_max)
                delete!(s.equil_dists, nb)
            else
                p0 = zeros(s.n_max + 1); p0[1] = 1.0
                s.dists[nb] = p0
            end
            delete!(s.pending_flux, nb)
        end
    end

    # ── 6. Restrict active voxels near local Poisson SS ───────────────────────
    if has_finite_local_ss(s.model)
        to_deact = CartesianIndex{2}[]
        for (idx, p) in s.dists
            ns      = grid_neighbors(idx, s.K_x, s.K_y)
            μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
            μ_loc   = (s.model.k_b + d * μ_in) / (s.model.k_d + d * length(ns))
            μ       = voxel_mean(p)
            abs(μ - μ_loc) / (μ_loc + 1e-10) < s.ε_equil || continue
            p_loc   = _poisson_pmf(μ_loc, s.n_max)
            tv      = 0.5 * sum(abs(p[i] - p_loc[i]) for i in eachindex(p_loc))
            tv < 0.08 || continue
            push!(to_deact, idx)
            s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
        end
        for idx in to_deact; delete!(s.dists, idx); end
    else
        to_deact = CartesianIndex{2}[]
        if s.model isa SchloglModel1D
            nothing  # deactivated in step 3a above
        else
            sat_thresh = 0.50 * Float64(s.n_max)
            for (idx, p) in s.dists
                voxel_mean(p) >= sat_thresh || continue
                ns = grid_neighbors(idx, s.K_x, s.K_y)
                any(nb -> !haskey(s.dists,nb) && !haskey(s.equil_dists,nb), ns) && continue
                push!(to_deact, idx); s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
            end
        end
        for idx in to_deact; delete!(s.dists, idx); end
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
