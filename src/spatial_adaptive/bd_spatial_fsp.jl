"""
Spatially Adaptive FSP for the 2D birth-death RDME.

Three voxel states:
  EMPTY       – not yet reached; P(n=0)=1 implicit
  ACTIVE      – tracked by FSP; stores full distribution p
  EQUILIBRATED – deactivated; stored mean μ_eq = LOCAL SS given current neighbours

Deactivation uses the LOCAL steady-state criterion: a voxel deactivates when
its distribution is near Poisson(μ_ss_local), where

    μ_ss_local = (k_b + d·Σ_nb μ_nb) / (k_d + d·|neighbours|)

This fires EARLY: the center reaches its local SS (≈ k_b/k_d_eff ≈ 0.22 for
four empty neighbours) within a fraction of one equilibration timescale.

After deactivation, μ_eq is updated each step as neighbours change:

    μ_eq_new = (k_b + d·Σ_nb μ_nb) / (k_d + d·|neighbours|)

Equilibrated voxels still drive the pending_flux for their EMPTY neighbours,
so the wavefront keeps propagating even when all voxels behind it are merged.

At long time, all μ_eq → k_b/k_d (global SS) and the total count distribution
collapses to Poisson(K²·k_b/k_d).  The active count stays small throughout —
only the wavefront ring is ever tracked simultaneously.
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
    dists        :: Dict{CartesianIndex{2}, Vector{Float64}}
    equil_means  :: Dict{CartesianIndex{2}, Float64}   # local SS mean, evolves each step
    pending_flux :: Dict{CartesianIndex{2}, Float64}   # cumulative flux into EMPTY voxels
    t            :: Float64
    n_max        :: Int
    ε_expand     :: Float64   # cumulative molecules to activate an EMPTY voxel
    ε_equil      :: Float64   # relative tolerance for LOCAL SS deactivation
    ε_prune      :: Float64
end

function SpatialFSP(model::BirthDeathRDME, K_x::Int, K_y::Int;
                    n_max::Int        = 40,
                    ε_expand::Float64 = 0.3,
                    ε_equil::Float64  = 0.06,
                    ε_prune::Float64  = 1e-14)
    SpatialFSP(model, K_x, K_y,
               Dict{CartesianIndex{2}, Vector{Float64}}(),
               Dict{CartesianIndex{2}, Float64}(),
               Dict{CartesianIndex{2}, Float64}(),
               0.0, n_max, ε_expand, ε_equil, ε_prune)
end

function set_ic!(s::SpatialFSP, idx::CartesianIndex{2}, n0::Int)
    n0 <= s.n_max || error("n0=$n0 > n_max=$(s.n_max)")
    p = zeros(s.n_max + 1);  p[n0 + 1] = 1.0
    s.dists[idx] = p
end

n_active(s::SpatialFSP)       = length(s.dists)
n_equilibrated(s::SpatialFSP) = length(s.equil_means)
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

# Effective mean for any voxel: active → from dist, equilibrated → stored, empty → 0
function _eff_mean(s::SpatialFSP, idx::CartesianIndex{2},
                   means::Dict{CartesianIndex{2}, Float64})
    get(means, idx, get(s.equil_means, idx, 0.0))
end

function _bd_generator(k_b_eff::Float64, k_d_eff::Float64, n_max::Int)
    n = n_max + 1
    Q = zeros(n, n)
    for col in 1:n
        nm = col - 1
        if col < n;  Q[col+1,col] += k_b_eff;  Q[col,col] -= k_b_eff;  end
        if col > 1;  r = k_d_eff * nm
                     Q[col-1,col] += r;       Q[col,col] -= r;           end
    end
    Q
end

function _poisson_pmf(λ::Float64, n_max::Int)
    p = zeros(n_max + 1);  p[1] = exp(-λ)
    for n in 1:n_max; p[n+1] = p[n] * λ / n; end
    p
end

poisson_pmf_vec(λ::Float64, n_max::Int) = _poisson_pmf(λ, n_max)

# ─── Main step ───────────────────────────────────────────────────────────────

"""
    step!(s, dt; krylov_m=30)

1. Build effective means: active dists + stored equil_means.
2. Evolve active voxels with mean-field coupling from ALL known neighbours.
3. Update equil_means: each equilibrated voxel's local SS tracks its neighbours.
4. Accumulate pending_flux from BOTH active AND equilibrated voxels into empty neighbours.
5. Activate empty voxels that have accumulated ≥ ε_expand expected molecules.
6. Deactivate active voxels near their LOCAL Poisson SS.
"""
function step!(s::SpatialFSP, dt::Float64; krylov_m::Int = 30)
    d = jump_rate(s.model)

    # ── 1. Effective means ────────────────────────────────────────────────────
    means = Dict{CartesianIndex{2}, Float64}()
    for (idx, p) in s.dists;         means[idx] = voxel_mean(p);    end
    for (idx, μ) in s.equil_means;   means[idx] = μ;                end

    # ── 2. Evolve active voxels ───────────────────────────────────────────────
    new_dists = Dict{CartesianIndex{2}, Vector{Float64}}()
    for (idx, p_old) in s.dists
        ns      = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        k_b_eff = s.model.k_b + d * μ_in
        k_d_eff = s.model.k_d + d * length(ns)
        Q     = _bd_generator(k_b_eff, k_d_eff, s.n_max)
        p_new = expv(Float64(dt), Q, p_old; m = min(krylov_m, size(Q, 1)))
        p_new .= max.(0.0, p_new);  p_new ./= sum(p_new)
        p_new[p_new .< s.ε_prune] .= 0.0
        new_dists[idx] = p_new
    end
    s.dists = new_dists

    # ── 3. Update equilibrated voxels' stored means (local SS relaxation) ─────
    for (idx, _) in s.equil_means
        ns      = grid_neighbors(idx, s.K_x, s.K_y)
        μ_in    = sum(get(means, nb, 0.0) for nb in ns; init = 0.0)
        k_b_eff = s.model.k_b + d * μ_in
        k_d_eff = s.model.k_d + d * length(ns)
        s.equil_means[idx] = k_b_eff / k_d_eff
    end

    # ── 4. Accumulate flux into empty voxels from ACTIVE + EQUILIBRATED ───────
    for (idx, μ) in means          # covers active + equil
        μ > 0.0 || continue
        for nb in grid_neighbors(idx, s.K_x, s.K_y)
            haskey(s.dists, nb)       && continue
            haskey(s.equil_means, nb) && continue
            s.pending_flux[nb] = get(s.pending_flux, nb, 0.0) + d * μ * dt
        end
    end

    # ── 5. Activate empty voxels ──────────────────────────────────────────────
    for (nb, fl) in s.pending_flux
        if fl >= s.ε_expand && !haskey(s.dists, nb) && !haskey(s.equil_means, nb)
            p0 = zeros(s.n_max + 1);  p0[1] = 1.0
            s.dists[nb] = p0
            delete!(s.pending_flux, nb)
        end
    end

    # ── 6. Deactivate active voxels near their LOCAL Poisson SS ──────────────
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
        s.equil_means[idx] = μ_loc
    end
    for idx in to_deactivate; delete!(s.dists, idx); end

    s.t += dt
    return s
end

# ─── Grid views ──────────────────────────────────────────────────────────────

"""
status_grid: 0=empty, 1=active, 2=equilibrated
"""
function status_grid(s::SpatialFSP)
    g = zeros(Int, s.K_x, s.K_y)
    for idx in keys(s.equil_means); g[idx] = 2; end
    for idx in keys(s.dists);       g[idx] = 1; end
    g
end

"""
mean_grid: E[n] per voxel. Active→from dist. Equilibrated→stored mean. Empty→0.
"""
function mean_grid(s::SpatialFSP)
    g = zeros(s.K_x, s.K_y)
    for (idx, μ) in s.equil_means; g[idx] = μ; end
    for (idx, p) in s.dists;       g[idx] = voxel_mean(p); end
    g
end

"""
total_distribution: convolves active distributions with Poisson(equil_means[k]).
Correct when voxels are independent (exact at RDME SS for birth-death).
"""
function total_distribution(s::SpatialFSP)
    function _conv!(a, b)
        n1, n2 = length(a), length(b)
        c = zeros(n1 + n2 - 1)
        for i in 1:n1, j in 1:n2; c[i+j-1] += a[i]*b[j]; end
        resize!(a, length(c));  a .= c
    end
    p_tot = [1.0]
    for (_, p)  in s.dists;       _conv!(p_tot, p); end
    for (_, μ)  in s.equil_means; _conv!(p_tot, _poisson_pmf(μ, s.n_max)); end
    p_tot
end

function equilibration_fraction(s::SpatialFSP; rtol::Float64 = 0.15)
    μ_ss   = ss_mean(s.model)
    n_tot  = s.K_x * s.K_y
    n_done = count(μ -> abs(μ - μ_ss)/μ_ss < rtol, values(s.equil_means))
    for (_, p) in s.dists
        abs(voxel_mean(p) - μ_ss)/μ_ss < rtol && (n_done += 1)
    end
    n_done / n_tot
end
