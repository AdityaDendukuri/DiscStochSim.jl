"""
joint_rdme_fsp.jl

Backbone RDME FSP: single joint CME over all tracked voxels.

State: Vector{Int16} of length (K_f + 1), where
  - dims 1..K_f  : particle counts for frontier (fine) voxels
  - dim  K_f+1   : N_int = total count in coarse interior region

Algorithm per step:
  1. diffuse_expand!   — add untracked spatial neighbors as new fine dimensions
  2. kemeny_restrict!  — demote frontier voxels (all neighbors tracked) to coarse interior
  3. rstep_expand!     — add stoichiometric neighbors of current support
  4. evolve!           — expv on joint CME generator
  5. flux_prune!       — remove states with p < ε_flux

Frontier voxels: those with at least one untracked spatial neighbor.
Interior voxels: all tracked voxels with all spatial neighbors tracked.
Interior is represented as a single coarse dimension N_int = total particle count.
This is exact under fast intra-interior mixing (Kemeny-Snell lumpability in the
limit of high interior diffusion relative to birth-death timescale).
"""

# ─── Struct ───────────────────────────────────────────────────────────────────

"""
    JointRDMEFSP

Joint CME over frontier voxels plus a coarse interior dimension.

State vector: (n_{v1}, ..., n_{vK_f}, N_int)
  - First K_f components: fine frontier voxel particle counts
  - Last component (index K_f+1): total particles in coarse interior

`fine_voxels[k]` is the CartesianIndex of the k-th fine dimension.
`interior_voxels` is the set of interior voxels aggregated into N_int.
"""
mutable struct JointRDMEFSP
    model            :: AbstractSpatialRDMEModel
    K_x              :: Int
    K_y              :: Int
    fine_voxels      :: Vector{CartesianIndex{2}}   # active fine voxels (dims 1..K_f)
    interior_voxels  :: Set{CartesianIndex{2}}       # coarse equilibrated interior
    n_interior       :: Int                          # number of interior sub-voxels (for coupling rate)
    ss               :: StateSpace{Vector{Int16}, Float64}
    t                :: Float64
    n_max            :: Int
    ε_flux           :: Float64
    ε_restrict       :: Float64   # coarsen voxel k when P(n_k=0) < ε_restrict (almost certainly occupied)
    krylov_m         :: Int
end

function JointRDMEFSP(model::AbstractSpatialRDMEModel, K_x::Int, K_y::Int;
                       n_max::Int       = 20,
                       ε_flux::Float64  = 1e-6,
                       ε_restrict::Float64 = 0.01,
                       krylov_m::Int    = 30)
    ss = StateSpace{Vector{Int16}, Float64}()
    JointRDMEFSP(model, K_x, K_y,
                 CartesianIndex{2}[],
                 Set{CartesianIndex{2}}(),
                 0,
                 ss, 0.0, n_max, ε_flux, ε_restrict, krylov_m)
end

# ─── IC ───────────────────────────────────────────────────────────────────────

"""
    set_ic!(jf::JointRDMEFSP, idx, n0)

Add a voxel to the initial condition with n0 particles.
All IC voxels start as fine (frontier classification happens after all ICs are set).
The new voxel becomes dimension K_f+1 of the joint state; all existing states are
lifted by appending n0 to each.
"""
function set_ic!(jf::JointRDMEFSP, idx::CartesianIndex{2}, n0::Int)
    idx ∈ jf.fine_voxels && error("Voxel $idx already in fine set")
    push!(jf.fine_voxels, idx)

    if length(jf.ss) == 0
        # First voxel — initial state: (n0, 0) where last dim is N_int=0
        add_state!(jf.ss, Int16[n0, 0], 1.0)
    else
        # Lift: insert new dimension before the last (N_int) dimension
        K_f_old = length(jf.fine_voxels) - 1  # before push
        new_ss = StateSpace{Vector{Int16}, Float64}()
        for i in eachindex(jf.ss.states)
            s = jf.ss.states[i]
            # Insert n0 at position K_f_old+1 (before N_int)
            s_new = vcat(s[1:K_f_old], Int16(n0), s[K_f_old+1:end])
            add_state!(new_ss, s_new, jf.ss.probs[i])
        end
        jf.ss = new_ss
    end
end

# State dimension helpers
_dim_Kf(jf::JointRDMEFSP)    = length(jf.fine_voxels)
_dim_total(jf::JointRDMEFSP)  = length(jf.fine_voxels) + 1  # K_f + N_int

# ─── Spatial helpers ──────────────────────────────────────────────────────────

"""
    _is_equilibrated(k, jf)

Returns true iff fine voxel `k` is eligible for coarsening into N_int:
  1. All its spatial neighbors are already tracked (fine or interior).
  2. P(n_k = 0) < ε_restrict — the voxel is almost certainly occupied.

Condition 2 is the Kemeny-Snell criterion for BranchingDeath: deep interior
voxels have exponentially growing mean counts so P(empty) → 0.  Wavefront
voxels and state-space-frontier voxels (added by rstep but not yet reached
by the wave) have P(empty) ≫ ε_restrict and must remain fine.
"""
function _is_equilibrated(k::Int, jf::JointRDMEFSP)
    idx = jf.fine_voxels[k]
    # Condition 1: all neighbors tracked
    for nb in grid_neighbors(idx, jf.K_x, jf.K_y)
        (nb ∈ jf.fine_voxels || nb ∈ jf.interior_voxels) || return false
    end
    # Condition 2: P(n_k = 0) < ε_restrict
    p_empty = 0.0
    for i in eachindex(jf.ss.states)
        jf.ss.states[i][k] == 0 && (p_empty += jf.ss.probs[i])
    end
    p_empty < jf.ε_restrict
end

# ─── Spatial expansion ────────────────────────────────────────────────────────

"""
    _diffuse_expand!(jf)

Add untracked spatial neighbors as new fine dimensions whenever any state
with non-negligible probability has n_k > 0 for an adjacent fine voxel k.
This is the organic spatial expansion: add dimension u when particles at k
can actually diffuse there (conservative: no threshold delay).
"""
function _diffuse_expand!(jf::JointRDMEFSP)
    K_f     = _dim_Kf(jf)
    tracked = union(Set(jf.fine_voxels), jf.interior_voxels)

    # Compute E[n_k] for each fine voxel; add untracked neighbors when E[n_k] > ε_flux.
    # Using mean flux D·E[n_k] > ε_flux keeps the state space frontier at the physical
    # wavefront, preventing premature coarsening of still-active boundary voxels.
    mean_n = zeros(K_f)
    for i in eachindex(jf.ss.states)
        p = jf.ss.probs[i]; p <= 0.0 && continue
        s = jf.ss.states[i]
        for k in 1:K_f
            mean_n[k] += p * s[k]
        end
    end

    to_add = CartesianIndex{2}[]
    for k in 1:K_f
        mean_n[k] * jump_rate(jf.model) < jf.ε_flux && continue
        for nb in grid_neighbors(jf.fine_voxels[k], jf.K_x, jf.K_y)
            (nb ∈ tracked || nb ∈ to_add) && continue
            push!(to_add, nb)
        end
    end

    for nb in to_add
        _add_fine_dim!(jf, nb)
    end
end

"""
    _add_fine_dim!(jf, idx)

Add untracked voxel `idx` as a new fine dimension (inserted before N_int).
All existing states are lifted: the new dimension starts at n=0.
"""
function _add_fine_dim!(jf::JointRDMEFSP, idx::CartesianIndex{2})
    K_f = _dim_Kf(jf)
    push!(jf.fine_voxels, idx)

    # Lift all states: insert 0 at position K_f+1 (before N_int which is at K_f+1)
    new_ss = StateSpace{Vector{Int16}, Float64}()
    for i in eachindex(jf.ss.states)
        s = jf.ss.states[i]
        s_new = vcat(s[1:K_f], Int16(0), s[K_f+1:end])
        add_state!(new_ss, s_new, jf.ss.probs[i])
    end
    jf.ss = new_ss
end

# ─── Kemeny-Snell restriction ─────────────────────────────────────────────────

"""
    _kemeny_restrict!(jf)

Demote fine voxels whose ALL spatial neighbors are now tracked into the coarse
interior. This is exact lumpability: total count N_int absorbs the demoted voxel.

For each demoted voxel k: marginalise dim k into N_int by summing
  p_new[(..., N_int')] += p_old[(..., n_k, ..., N_int)]  for all n_k (N_int' = N_int + n_k)
"""
function _kemeny_restrict!(jf::JointRDMEFSP)
    changed = true
    while changed
        changed = false
        for k in reverse(eachindex(jf.fine_voxels))   # reverse so indices stay valid
            _is_equilibrated(k, jf) || continue
            _demote_fine_to_interior!(jf, k)
            changed = true
        end
    end
end

function _demote_fine_to_interior!(jf::JointRDMEFSP, k::Int)
    K_f     = _dim_Kf(jf)
    dim_int = K_f + 1   # N_int dimension index

    new_ss = StateSpace{Vector{Int16}, Float64}()
    for i in eachindex(jf.ss.states)
        s = jf.ss.states[i]; p_old = jf.ss.probs[i]
        n_k   = s[k]
        N_int = s[dim_int]
        # New state: remove dim k, add n_k to N_int
        s_new = Vector{Int16}(undef, K_f)   # length K_f = (K_f+1) - 1
        pos = 0
        for j in 1:K_f
            j == k && continue
            pos += 1
            s_new[pos] = s[j]
        end
        s_new[pos+1] = N_int + n_k   # last entry is new N_int

        idx_new = get(new_ss.index, s_new, 0)
        if idx_new == 0
            add_state!(new_ss, s_new, p_old)
        else
            new_ss.probs[idx_new] += p_old
        end
    end

    push!(jf.interior_voxels, jf.fine_voxels[k])
    jf.n_interior += 1
    deleteat!(jf.fine_voxels, k)
    jf.ss = new_ss
end

# ─── Rstep expand ────────────────────────────────────────────────────────────

"""
    _rstep_expand!(jf)

Add all stoichiometric neighbors of the current support (birth, death,
intra-fine diffusion, fine↔coarse exchange) that satisfy the boundary condition.
"""
function _rstep_expand!(jf::JointRDMEFSP)
    K_f     = _dim_Kf(jf)
    n_max   = jf.n_max
    d       = jump_rate(jf.model)
    pos     = Dict(jf.fine_voxels[k] => k for k in 1:K_f)
    dim_int = K_f + 1
    n_int_max = jf.n_max * max(1, jf.n_interior)  # generous upper bound for N_int

    bc = s -> (length(s) == dim_int &&
               all(Int(s[k]) >= 0 && Int(s[k]) <= n_max for k in 1:K_f) &&
               Int(s[dim_int]) >= 0 && Int(s[dim_int]) <= n_int_max)

    current = copy(jf.ss.states)
    for s in current
        for k in 1:K_f
            # Birth at fine voxel k
            if Int(s[k]) < n_max
                s_new = copy(s); s_new[k] += Int16(1)
                bc(s_new) && !haskey(jf.ss.index, s_new) && add_state!(jf.ss, s_new, 0.0)
            end
            # Death at fine voxel k
            if Int(s[k]) > 0
                s_new = copy(s); s_new[k] -= Int16(1)
                bc(s_new) && !haskey(jf.ss.index, s_new) && add_state!(jf.ss, s_new, 0.0)
            end
            # Intra-fine diffusion k → q
            Int(s[k]) > 0 || continue
            for nb in grid_neighbors(jf.fine_voxels[k], jf.K_x, jf.K_y)
                q = get(pos, nb, 0); q == 0 && continue
                Int(s[q]) < n_max || continue
                s_new = copy(s); s_new[k] -= Int16(1); s_new[q] += Int16(1)
                bc(s_new) && !haskey(jf.ss.index, s_new) && add_state!(jf.ss, s_new, 0.0)
            end
            # Fine → interior (if fine voxel k is adjacent to an interior voxel)
            _has_interior_nb(jf.fine_voxels[k], jf) || continue
            Int(s[dim_int]) < n_int_max || continue
            s_new = copy(s); s_new[k] -= Int16(1); s_new[dim_int] += Int16(1)
            bc(s_new) && !haskey(jf.ss.index, s_new) && add_state!(jf.ss, s_new, 0.0)
        end
        # Interior → fine (N_int > 0 and interior touches this fine voxel)
        Int(s[dim_int]) > 0 || continue
        jf.n_interior > 0  || continue
        for k in 1:K_f
            _has_interior_nb(jf.fine_voxels[k], jf) || continue
            Int(s[k]) < n_max || continue
            s_new = copy(s); s_new[k] += Int16(1); s_new[dim_int] -= Int16(1)
            bc(s_new) && !haskey(jf.ss.index, s_new) && add_state!(jf.ss, s_new, 0.0)
        end
    end
end

function _has_interior_nb(idx::CartesianIndex{2}, jf::JointRDMEFSP)
    for nb in grid_neighbors(idx, jf.K_x, jf.K_y)
        nb ∈ jf.interior_voxels && return true
    end
    false
end

# ─── Generator ───────────────────────────────────────────────────────────────

"""
    _build_joint_rdme_system(jf) → DiscreteStochasticSystem{Vector{Int16}}

Build the RDME stochastic system for the current fine voxels + coarse interior.

Reactions:
  - Birth at each fine voxel k (rate: local_birth_rate(n_k))
  - Death at each fine voxel k (rate: k_d × n_k)
  - Birth in coarse interior (rate: local_birth_rate(N_int, n_interior))
  - Death in coarse interior (rate: k_d × N_int)
  - Intra-fine diffusion k↔q for adjacent fine voxels (rate: D × n_k or D × n_q)
  - Fine k → interior (rate: D × n_k × n_interior_neighbors_k)
  - Interior → fine k (rate: D × N_int / n_interior × n_interior_neighbors_k)
    (uniform mixing assumption: each interior voxel equally likely to hold a particle)
  - Fine k → untracked (rate: D × n_k × n_untracked_neighbors_k)  [absorbed: not added]
"""
function _build_joint_rdme_system(jf::JointRDMEFSP)
    K_f     = _dim_Kf(jf)
    d       = jump_rate(jf.model)
    pos     = Dict(jf.fine_voxels[k] => k for k in 1:K_f)
    dim_int = K_f + 1
    n_int   = max(1, jf.n_interior)  # avoid divide-by-zero

    stoichs      = Vector{Vector{Int16}}()
    propensities = Function[]

    function e_k(k)   # unit vector: +1 at dim k, 0 elsewhere
        v = zeros(Int16, dim_int); v[k] = Int16(1); v
    end

    # Birth / death at each fine voxel
    for k in 1:K_f
        let k = k
            push!(stoichs,  e_k(k)); push!(propensities, (x, rv, t) -> local_birth_rate(jf.model, Int(x[k])))
            push!(stoichs, -e_k(k)); push!(propensities, (x, rv, t) -> jf.model.k_d * max(0, Int(x[k])))
        end
    end

    # Birth / death in coarse interior
    if jf.n_interior > 0
        push!(stoichs,  e_k(dim_int))
        push!(propensities, (x, rv, t) -> local_total_birth_rate(jf.model, Int(x[dim_int]), n_int))
        push!(stoichs, -e_k(dim_int))
        push!(propensities, (x, rv, t) -> jf.model.k_d * max(0, Int(x[dim_int])))
    end

    # Intra-fine diffusion k ↔ q (add each pair once)
    for k in 1:K_f
        for nb in grid_neighbors(jf.fine_voxels[k], jf.K_x, jf.K_y)
            q = get(pos, nb, 0); (q == 0 || q <= k) && continue
            let k = k, q = q
                v_kq = zeros(Int16, dim_int); v_kq[k] = Int16(-1); v_kq[q] = Int16(1)
                v_qk = zeros(Int16, dim_int); v_qk[q] = Int16(-1); v_qk[k] = Int16(1)
                push!(stoichs, v_kq); push!(propensities, (x, rv, t) -> d * max(0, Int(x[k])))
                push!(stoichs, v_qk); push!(propensities, (x, rv, t) -> d * max(0, Int(x[q])))
            end
        end
    end

    # Fine k ↔ interior (one reaction per interior neighbor of k)
    for k in 1:K_f
        n_int_nb = count(nb -> nb ∈ jf.interior_voxels,
                         grid_neighbors(jf.fine_voxels[k], jf.K_x, jf.K_y))
        n_int_nb == 0 && continue
        let k = k, n_int_nb = n_int_nb, n_int = n_int
            # Fine k → interior: rate D × n_k × n_int_nb
            v_to_int = zeros(Int16, dim_int); v_to_int[k] = Int16(-1); v_to_int[dim_int] = Int16(1)
            push!(stoichs, v_to_int)
            push!(propensities, (x, rv, t) -> d * n_int_nb * max(0, Int(x[k])))
            # Interior → fine k: rate D × N_int / n_int × n_int_nb
            # (uniform mixing: each interior sub-voxel has prob 1/n_int of holding a particle)
            v_from_int = zeros(Int16, dim_int); v_from_int[k] = Int16(1); v_from_int[dim_int] = Int16(-1)
            push!(stoichs, v_from_int)
            push!(propensities, (x, rv, t) -> d * n_int_nb / n_int * max(0, Int(x[dim_int])))
        end
    end

    # Note: fine k → untracked exterior reactions are intentionally omitted.
    # With a conservative generator (absorbing=false), probability reaching the
    # boundary of the state space is not drained.  _diffuse_expand! proactively
    # adds the next ring of voxels so particles are never blocked at the frontier.

    DiscreteStochasticSystem{Vector{Int16}}(stoichs, propensities)
end

# ─── Step ────────────────────────────────────────────────────────────────────

function step!(jf::JointRDMEFSP, dt::Float64)
    isempty(jf.fine_voxels) && return

    # 1. Organic spatial expansion
    _diffuse_expand!(jf)

    # 2. Kemeny-Snell restriction (interior voxels → coarse)
    _kemeny_restrict!(jf)

    # 3. Rstep expand (birth/death/diffusion neighbors of current support)
    _rstep_expand!(jf)

    # 4. Build generator and evolve
    K_f      = _dim_Kf(jf)
    dim_int  = K_f + 1
    n_int    = max(1, jf.n_interior)
    n_int_max = jf.n_max * n_int

    sys = _build_joint_rdme_system(jf)
    A,  = build_generator(jf.ss, sys, nothing, jf.t; absorbing=false)

    km = min(jf.krylov_m, length(jf.ss))
    km = max(km, 4)
    jf.ss.probs .= expv(dt, A, jf.ss.probs; m = km)
    jf.ss.probs  = max.(0.0, jf.ss.probs)

    # 5. Probability prune: remove states below threshold
    prune_threshold!(jf.ss, jf.ε_flux)

    jf.t += dt
end

# ─── Diagnostics ─────────────────────────────────────────────────────────────

n_fine(jf::JointRDMEFSP)     = length(jf.fine_voxels)
n_interior(jf::JointRDMEFSP) = length(jf.interior_voxels)
n_tracked(jf::JointRDMEFSP)  = n_fine(jf) + n_interior(jf)
n_joint_states(jf::JointRDMEFSP) = length(jf.ss)

"""
    wavefront_radius(jf, ci, cj; p_thresh=0.01) → Float64

Maximum distance from (ci, cj) of any fine voxel where P(n_k > 0) ≥ p_thresh.

With p_thresh ≈ 1/N_traj this matches the SSA definition: max radius that at
least one of N_traj trajectories visits.  Default 0.01 ≈ 1/100 trajectories.
"""
function wavefront_radius(jf::JointRDMEFSP, ci::Int, cj::Int; p_thresh::Float64=0.01)
    isempty(jf.fine_voxels) && return 0.0
    K_f  = _dim_Kf(jf)
    max_r = 0.0
    for k in 1:K_f
        # P(n_k > 0) = sum of p[s] where s[k] > 0
        p_pos = 0.0
        for i in eachindex(jf.ss.states)
            jf.ss.states[i][k] > 0 && (p_pos += jf.ss.probs[i])
        end
        p_pos < p_thresh && continue
        v = jf.fine_voxels[k]
        r = sqrt(Float64(v[1]-ci)^2 + Float64(v[2]-cj)^2)
        max_r = max(max_r, r)
    end
    max_r
end

"""
    wavefront_radius_mean(jf, ci, cj) → Float64

E[max_r(t)] = Σ_s p[s] × max_{k: s[k]>0} r_k  over fine voxels.

This is the correct SSA-comparable wavefront metric: it matches the mean of
max occupied radius across trajectories, rather than the frontier of the state
space (which grows at rstep speed and is not physically meaningful).
"""
function wavefront_radius_mean(jf::JointRDMEFSP, ci::Int, cj::Int)
    isempty(jf.fine_voxels) && return 0.0
    K_f = _dim_Kf(jf)

    # Precompute radius of each fine voxel
    radii = [sqrt(Float64(jf.fine_voxels[k][1]-ci)^2 + Float64(jf.fine_voxels[k][2]-cj)^2)
             for k in 1:K_f]

    E_max_r = 0.0
    for i in eachindex(jf.ss.states)
        p = jf.ss.probs[i]; p > 0 || continue
        s = jf.ss.states[i]
        max_r = 0.0
        for k in 1:K_f
            s[k] > 0 && (max_r = max(max_r, radii[k]))
        end
        E_max_r += p * max_r
    end
    E_max_r
end

"""
    voxel_mean(jf, k) → Float64

Expected particle count in fine voxel k.
"""
function voxel_mean(jf::JointRDMEFSP, k::Int)
    μ = 0.0
    for i in eachindex(jf.ss.states)
        μ += jf.ss.states[i][k] * jf.ss.probs[i]
    end
    μ
end

"""
    interior_mean(jf) → Float64

Expected total particle count in the coarse interior.
"""
function interior_mean(jf::JointRDMEFSP)
    dim_int = _dim_Kf(jf) + 1
    μ = 0.0
    for i in eachindex(jf.ss.states)
        μ += jf.ss.states[i][dim_int] * jf.ss.probs[i]
    end
    μ
end
