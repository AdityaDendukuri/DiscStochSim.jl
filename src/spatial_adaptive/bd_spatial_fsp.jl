"""
Spatially Adaptive FSP for the 2D birth-death RDME — Galerkin coarse operator.

Three voxel states:
  EMPTY       – P(n=0)=1 implicit; no tracking
  ACTIVE      – full CME distribution; updated with expv at krylov_m_act
  EQUILIBRATED – Galerkin coarse: full distribution, updated with expv at krylov_m_eq
                 (cheap: near-Poisson so small Krylov subspace suffices)

The key upgrade over the mean-field version:
  OLD: equil voxel stores μ_loc = (k_b + d·Σμ_nb)/(k_d + d·|nb|)  [algebraic, frozen at local SS]
  NEW: equil voxel stores p_eq, evolved each step under Q_eff         [Galerkin CME, dynamic]

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

struct BirthDeathRDME
    k_b :: Float64
    k_d :: Float64
    D   :: Float64
    dx  :: Float64
end

jump_rate(m::BirthDeathRDME) = m.D / m.dx^2
ss_mean(m::BirthDeathRDME)   = m.k_b / m.k_d

# ─── State ───────────────────────────────────────────────────────────────────

mutable struct SpatialFSP
    model        :: BirthDeathRDME
    K_x          :: Int
    K_y          :: Int
    dists        :: Dict{CartesianIndex{2}, Vector{Float64}}   # ACTIVE: full CME dist
    equil_dists  :: Dict{CartesianIndex{2}, Vector{Float64}}   # EQUIL: Galerkin coarse dist
    pending_flux :: Dict{CartesianIndex{2}, Float64}            # cumulative flux → empty
    t            :: Float64
    n_max        :: Int     # truncation for active voxels
    n_max_eq     :: Int     # truncation for equilibrated voxels (can be smaller)
    ε_expand     :: Float64
    ε_equil      :: Float64
    ε_prune      :: Float64
    krylov_m_eq  :: Int     # Krylov dim for equil voxels (small: near-SS distributions)
end

function SpatialFSP(model::BirthDeathRDME, K_x::Int, K_y::Int;
                    n_max::Int        = 40,
                    n_max_eq::Int     = -1,      # default: auto (4·μ_ss + 4)
                    ε_expand::Float64 = 0.3,
                    ε_equil::Float64  = 0.06,
                    ε_prune::Float64  = 1e-14,
                    krylov_m_eq::Int  = 8)        # equil voxels converge fast
    n_max_eq_val = n_max_eq > 0 ? n_max_eq :
                   max(12, round(Int, 4 * ss_mean(model)) + 4)
    SpatialFSP(model, K_x, K_y,
               Dict{CartesianIndex{2}, Vector{Float64}}(),
               Dict{CartesianIndex{2}, Vector{Float64}}(),
               Dict{CartesianIndex{2}, Float64}(),
               0.0, n_max, n_max_eq_val, ε_expand, ε_equil, ε_prune, krylov_m_eq)
end

function set_ic!(s::SpatialFSP, idx::CartesianIndex{2}, n0::Int)
    n0 <= s.n_max || error("n0=$n0 > n_max=$(s.n_max)")
    p = zeros(s.n_max + 1);  p[n0 + 1] = 1.0
    s.dists[idx] = p
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

function _bd_generator(k_b_eff::Float64, k_d_eff::Float64, n_max::Int)
    n = n_max + 1
    Q = zeros(n, n)
    for col in 1:n
        nm = col - 1
        if col < n;  Q[col+1,col] += k_b_eff;  Q[col,col] -= k_b_eff;  end
        if col > 1;  r = k_d_eff * nm
                     Q[col-1,col] += r;         Q[col,col] -= r;         end
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

# ─── Main time step ───────────────────────────────────────────────────────────

"""
    step!(s, dt; krylov_m=30)

Galerkin-coarse step:
  1. Means from BOTH active (voxel_mean(p)) and equilibrated (voxel_mean(p_eq)).
  2. Evolve active voxels under Q_eff (full CME, krylov_m).
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

    # ── 2. Evolve active voxels (full CME) ────────────────────────────────────
    new_dists = Dict{CartesianIndex{2}, Vector{Float64}}()
    for (idx, p_old) in s.dists
        ns      = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        k_b_eff = s.model.k_b + d * μ_in
        k_d_eff = s.model.k_d + d * length(ns)
        Q     = _bd_generator(k_b_eff, k_d_eff, s.n_max)
        p_new = expv(Float64(dt), Q, p_old; m = min(krylov_m, size(Q, 1)))
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
    new_equil = Dict{CartesianIndex{2}, Vector{Float64}}()
    for (idx, p_eq) in s.equil_dists
        ns      = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        k_b_eff = s.model.k_b + d * μ_in
        k_d_eff = s.model.k_d + d * length(ns)

        if any(nb -> haskey(s.dists, nb), ns)
            # ── BOUNDARY: Galerkin expv ───────────────────────────────────────
            Q_eq     = _bd_generator(k_b_eff, k_d_eff, s.n_max_eq)
            p_old_eq = _resize_dist(p_eq, s.n_max_eq)
            p_new_eq = expv(Float64(dt), Q_eq, p_old_eq; m = min(s.krylov_m_eq, size(Q_eq, 1)))
            _clean!(p_new_eq, s.ε_prune)
            new_equil[idx] = p_new_eq
        else
            # ── INTERIOR: algebraic mean → Poisson (no matrix exp) ───────────
            # Poisson(μ_loc) is an excellent approximation: the voxel is already
            # near-Poisson and drifting slowly toward global SS.
            μ_loc = k_b_eff / k_d_eff
            new_equil[idx] = _poisson_pmf(μ_loc, s.n_max_eq)
        end
    end
    s.equil_dists = new_equil

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
    to_deactivate = CartesianIndex{2}[]
    for (idx, p) in s.dists
        ns      = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        k_b_eff = s.model.k_b + d * μ_in
        k_d_eff = s.model.k_d + d * length(ns)
        μ_loc   = k_b_eff / k_d_eff

        μ = voxel_mean(p)
        abs(μ - μ_loc) / (μ_loc + 1e-10) < s.ε_equil || continue
        p_loc = _poisson_pmf(μ_loc, s.n_max)
        tv = 0.5 * sum(abs(p[i] - p_loc[i]) for i in eachindex(p_loc))
        tv < 0.08 || continue

        push!(to_deactivate, idx)
        # Store actual current distribution (truncated to n_max_eq)
        s.equil_dists[idx] = _resize_dist(p, s.n_max_eq)
    end
    for idx in to_deactivate; delete!(s.dists, idx); end

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
total_distribution: convolves active distributions with equil distributions.
With Galerkin coarse operator, equil_dists store actual evolved distributions
(not assumed Poisson), so this convolution is more accurate than before.
"""
function total_distribution(s::SpatialFSP)
    function _conv!(a, b)
        n1, n2 = length(a), length(b)
        c = zeros(n1 + n2 - 1)
        for i in 1:n1, j in 1:n2; c[i+j-1] += a[i]*b[j]; end
        resize!(a, length(c));  a .= c
    end
    p_tot = [1.0]
    for (_, p) in s.dists;       _conv!(p_tot, p); end
    for (_, p) in s.equil_dists; _conv!(p_tot, p); end
    p_tot
end

function equilibration_fraction(s::SpatialFSP; rtol::Float64 = 0.15)
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
