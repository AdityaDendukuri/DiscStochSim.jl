# ── Flux-based adaptive V-cycle for 2D RDME front voxels ─────────────────────
#
# Key idea: maintain a weighted graph on active voxels where
#   w_{ij} = EMA of  d · min(μ_i, μ_j) · dt
#
# High w_{ij}  → both voxels have similar non-zero means → fast intra-pair
#                mixing → Multinomial prolongation is valid → MERGE
# Low w_{ij}   → large gradient across edge (front boundary) → KEEP FINE
#
# At each V-cycle level:
#   1. Greedy max-weight matching → groups (pairs + singletons)
#   2. Convolve group marginals   → coarse marginals
#   3. Recurse at coarser level
#   4. Multinomial-split coarse marginals → per-voxel marginals (weighted by μ)

# ─── FluxGraph ────────────────────────────────────────────────────────────────

mutable struct FluxGraph
    voxels  :: Vector{CartesianIndex{2}}
    edges   :: Vector{Tuple{Int,Int}}    # (i,j) i<j in local voxel index
    weights :: Vector{Float64}           # EMA flux weights per edge
    α       :: Float64                   # EMA decay  (new = α·old + (1-α)·current)
end

"""
    build_flux_graph(active_voxels, K_x, K_y; α=0.7) -> FluxGraph

Build a flux graph on the active voxels with zero initial weights.
"""
function build_flux_graph(voxels::Vector{CartesianIndex{2}}, K_x::Int, K_y::Int;
                           α::Float64 = 0.7)
    pos = Dict(v => k for (k, v) in enumerate(voxels))
    edges = Tuple{Int,Int}[]
    for (k, v) in enumerate(voxels)
        for nb in grid_neighbors(v, K_x, K_y)
            l = get(pos, nb, 0)
            l > 0 && k < l && push!(edges, (k, l))
        end
    end
    FluxGraph(voxels, edges, zeros(length(edges)), α)
end

"""
    update_flux!(fg, dists, d, dt)

Update edge weights using EMA:  w = α·w_old + (1-α)·d·min(μ_i,μ_j)·dt
"""
function update_flux!(fg::FluxGraph,
                      dists::Dict{CartesianIndex{2}, Vector{Float64}},
                      d::Float64, dt::Float64)
    for (e, (i, j)) in enumerate(fg.edges)
        μi = haskey(dists, fg.voxels[i]) ? voxel_mean(dists[fg.voxels[i]]) : 0.0
        μj = haskey(dists, fg.voxels[j]) ? voxel_mean(dists[fg.voxels[j]]) : 0.0
        w_cur = d * min(μi, μj) * dt
        fg.weights[e] = fg.α * fg.weights[e] + (1.0 - fg.α) * w_cur
    end
end

# ─── Greedy max-weight matching ───────────────────────────────────────────────

"""
    greedy_flux_matching(fg) -> (pairs, singletons)

Greedy maximum-weight matching: sort edges by weight (high = good merge),
greedily assign non-overlapping pairs.  Unmatched voxels become singletons.
"""
function greedy_flux_matching(fg::FluxGraph)
    n      = length(fg.voxels)
    used   = fill(false, n)
    order  = sortperm(fg.weights; rev = true)
    pairs  = Tuple{Int,Int}[]
    for e in order
        i, j = fg.edges[e]
        used[i] || used[j] || begin
            push!(pairs, (i, j)); used[i] = used[j] = true
        end
    end
    pairs, findall(.!used)
end

# ─── Coarse system for arbitrary groups ──────────────────────────────────────

"""
    _build_group_coarse_sys(model, fine_voxels, groups, means, K_x, K_y)

Build the coarse RDME for K_G super-voxels whose membership is described by
`groups[j]` = vector of local indices into `fine_voxels`.

Birth:           λ · |group_j| per super-voxel j
Death:           k_d · N_j per super-voxel j
Diffusion j→k:   d · (# fine-grid edges between groups j,k) / |group_j|  per molecule
Outflow to bg:   d · (# ext edges from group_j) / |group_j|              per molecule
Inflow from bg:  d · Σ μ_nb  (sum over inactive neighbours of any voxel in j)
"""
function _build_group_coarse_sys(
    model       :: AbstractSpatialRDMEModel,
    fine_voxels :: Vector{CartesianIndex{2}},
    groups      :: Vector{Vector{Int}},
    means       :: Dict{CartesianIndex{2}, Float64},
    K_x         :: Int,
    K_y         :: Int
)
    K_G      = length(groups)
    d        = jump_rate(model)
    stoichs  = CartesianIndex{K_G}[]
    props    = Function[]

    active_set  = Set(fine_voxels)
    fine_to_grp = Dict{CartesianIndex{2}, Int}()
    for (j, g) in enumerate(groups)
        for k in g; fine_to_grp[fine_voxels[k]] = j; end
    end
    gsz = [length(g) for g in groups]

    eu(j) = CartesianIndex(ntuple(k -> k == j ? 1 : 0, K_G))

    # Birth + death
    for j in 1:K_G
        let j=j, sz=Float64(gsz[j]), kd=model.k_d
            ej = eu(j)
            push!(stoichs,  ej); push!(props, (x,rv,t) -> sz * local_birth_rate(model, Tuple(x)[j]))
            push!(stoichs, -ej); push!(props, (x,rv,t) -> kd * max(0, Tuple(x)[j]))
        end
    end

    # Inter-group diffusion
    inter = Dict{Tuple{Int,Int}, Int}()
    for (j, g) in enumerate(groups), k_fine in g
        for nb in grid_neighbors(fine_voxels[k_fine], K_x, K_y)
            nb ∈ active_set || continue
            l = fine_to_grp[nb]
            l == j && continue
            key = minmax(j, l)
            inter[key] = get(inter, key, 0) + 1
        end
    end
    for ((j, l), ne) in inter
        let j=j, l=l, rjl=d*ne/gsz[j], rlj=d*ne/gsz[l]
            ej = eu(j); el = eu(l)
            push!(stoichs, -ej+el); push!(props, (x,rv,t) -> rjl * max(0, Tuple(x)[j]))
            push!(stoichs,  ej-el); push!(props, (x,rv,t) -> rlj * max(0, Tuple(x)[l]))
        end
    end

    # Boundary with inactive
    r_out = zeros(K_G); r_in = zeros(K_G)
    for (j, g) in enumerate(groups), k_fine in g
        for nb in grid_neighbors(fine_voxels[k_fine], K_x, K_y)
            nb ∈ active_set && continue
            r_out[j] += d / gsz[j]
            r_in[j]  += d * get(means, nb, 0.0)
        end
    end
    for j in 1:K_G
        let j=j, ro=r_out[j], ri=r_in[j]
            ej = eu(j)
            ro > 0 && (push!(stoichs,-ej); push!(props,(x,rv,t)->ro*max(0,Tuple(x)[j])))
            ri > 0 && (push!(stoichs, ej); push!(props,(x,rv,t)->ri))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_G}}(stoichs, props)
end

# ─── Helpers ──────────────────────────────────────────────────────────────────

# Convolve a list of 1-D probability vectors (total-count distribution)
function _convolve_group(ps::Vector{Vector{Float64}})
    out = copy(ps[1])
    for i in 2:length(ps)
        p2  = ps[i]
        n1, n2 = length(out), length(p2)
        c = zeros(n1 + n2 - 1)
        for a in 1:n1, b in 1:n2; c[a+b-1] += out[a] * p2[b]; end
        out = c
    end
    out
end

# Extract 1-D marginal of dimension k from a joint StateSpace
function _marginal_k(sp::StateSpace{CartesianIndex{K}, Float64}, k::Int, n_max::Int) where K
    p = zeros(n_max + 1)
    for i in eachindex(sp.states)
        sp.probs[i] > 0.0 || continue
        n = Tuple(sp.states[i])[k]
        n <= n_max && (p[n+1] += sp.probs[i])
    end
    s = sum(p); s > 0 && (p ./= s)
    p
end

"""
    _multinomial_split(p_N, vox_means, n_max_out, tol) -> Vector{Vector}

Split total-count distribution P(N=n) into per-voxel marginals using
Multinomial(N, w_1/W,...,w_M/W) where w_k ∝ vox_means[k].
Reduces to Binomial(N,½) when all means are equal.
"""
function _multinomial_split(p_N::Vector{Float64}, vox_means::Vector{Float64},
                             n_max_out::Int, tol::Float64)
    M = length(vox_means)
    W = max(sum(vox_means), 1e-10)
    probs = [max(w, 1e-10) / W for w in vox_means]
    out   = [zeros(n_max_out + 1) for _ in 1:M]
    buf   = zeros(Int, M)
    for N in 0:length(p_N)-1
        p_N[N+1] > tol || continue
        _ms_rec!(out, buf, N, p_N[N+1], probs, n_max_out, tol, 1)
    end
    for p in out; s = sum(p); s > 0 && (p ./= s); end
    out
end

function _ms_rec!(out, buf, n_rem, w, probs, n_max, tol, m)
    M = length(out)
    if m == M
        n_rem <= n_max || return
        buf[m] = n_rem
        for k in 1:M; out[k][buf[k]+1] += w; end
        return
    end
    p_m   = probs[m]
    p_rest = sum(probs[m+1:end])
    p_tot  = p_m + p_rest
    p_cond = p_tot > 0 ? p_m / p_tot : 0.0
    lnf = sum(log(Float64(i)) for i in 1:n_rem; init=0.0)
    for k in 0:min(n_rem, n_max)
        bw = if p_cond <= 0.0
                 k == 0 ? 1.0 : 0.0
             elseif p_cond >= 1.0
                 k == n_rem ? 1.0 : 0.0
             else
                 lk  = sum(log(Float64(i)) for i in 1:k;       init=0.0)
                 lnk = sum(log(Float64(i)) for i in 1:n_rem-k; init=0.0)
                 exp(lnf - lk - lnk + k*log(p_cond) + (n_rem-k)*log(1-p_cond))
             end
        wk = w * bw; wk > tol || continue
        buf[m] = k
        _ms_rec!(out, buf, n_rem-k, wk, probs, n_max, tol, m+1)
    end
end

# Build a coarser flux graph from groups (representative = first voxel of each group)
function _coarsen_flux_graph(fg::FluxGraph, groups::Vector{Vector{Int}},
                              rep_voxels::Vector{CartesianIndex{2}}, K_x::Int, K_y::Int)
    # Build coarser version: only keep edges between representatives of different groups
    new_fg = build_flux_graph(rep_voxels, K_x, K_y; α = fg.α)
    # Transfer weights: edge (j_rep, l_rep) weight = max weight of fine edges between groups j and l
    grp_of = Dict{Int,Int}()
    for (j, g) in enumerate(groups), k in g; grp_of[k] = j; end
    for (e, (i, j)) in enumerate(fg.edges)
        gi = grp_of[i]; gj = grp_of[j]
        gi == gj && continue
        # Find the corresponding edge in new_fg
        for (e2, (ri, rj)) in enumerate(new_fg.edges)
            (ri == gi && rj == gj) || (ri == gj && rj == gi) || continue
            new_fg.weights[e2] = max(new_fg.weights[e2], fg.weights[e])
            break
        end
    end
    new_fg
end

# ─── Recursive V-cycle ────────────────────────────────────────────────────────

"""
    _flux_vcycle!(dists, fine_voxels, fg, means, model, K_x, K_y, dt,
                  n_max, tol, ε_prune, krylov_m, max_states, depth, t)

Recursive flux-based V-cycle on `length(dists)` active voxels.

dists   – Vector of per-voxel marginals P(n_k=·) (modified in-place)
fg      – FluxGraph for this level
depth   – remaining recursion levels (0 = base: direct expv)

Returns true on success, false if state space would explode (caller should fall back).
"""
function _flux_vcycle!(
    dists        :: Vector{Vector{Float64}},
    fine_voxels  :: Vector{CartesianIndex{2}},
    fg           :: FluxGraph,
    means        :: Dict{CartesianIndex{2}, Float64},
    model        :: AbstractSpatialRDMEModel,
    K_x          :: Int, K_y :: Int,
    dt           :: Float64,
    n_max        :: Int,
    tol          :: Float64,
    ε_prune      :: Float64,
    krylov_m     :: Int,
    max_states   :: Int,
    depth        :: Int,
    t            :: Float64
)
    K = length(dists)
    K == length(fine_voxels) || error("dists/voxels length mismatch")

    # ── Base case: direct joint expv ─────────────────────────────────────────
    if depth == 0 || K <= 2
        sp = _product_state_from_dists(dists, tol, max_states)
        sp === nothing && return false
        groups_id = [[k] for k in 1:K]
        sys = _build_group_coarse_sys(model, fine_voxels, groups_id, means, K_x, K_y)
        bc  = x -> _joint_active_bc(x, n_max)
        expand!(sp, sys, bc; depth=1)
        length(sp) > max_states && return false
        A, = build_generator(sp, sys, nothing, t)
        sp.probs .= expv(Float64(dt), A, sp.probs; m=min(krylov_m, size(A,1)))
        prune_threshold!(sp, tol); renormalize!(sp)
        margs = _joint_marginals(sp, n_max)
        for k in 1:K; dists[k] = margs[k]; end
        return true
    end

    # ── Greedy matching ───────────────────────────────────────────────────────
    pairs, singletons = greedy_flux_matching(fg)
    groups = vcat([[p[1], p[2]] for p in pairs], [[k] for k in singletons])
    K_G    = length(groups)

    # ── Restrict: save local means for prolongation, convolve marginals ───────
    grp_means  = [[voxel_mean(dists[k]) for k in g] for g in groups]
    coarse_d   = [_convolve_group([dists[k] for k in g]) for g in groups]

    # Coarse n_max: use actual support of coarse distributions, not 2*n_max.
    # This prevents exponential blowup at deeper recursion levels.
    n_max_c = maximum(
        findlast(>(tol), p) - 1
        for p in coarse_d
        if any(>(tol), p);
        init = 2 * n_max
    )
    n_max_c = max(n_max_c, n_max)   # never smaller than fine n_max

    # ── Coarsen flux graph ────────────────────────────────────────────────────
    rep_voxels = [fine_voxels[g[1]] for g in groups]
    fg_c       = _coarsen_flux_graph(fg, groups, rep_voxels, K_x, K_y)

    # ── Recurse ───────────────────────────────────────────────────────────────
    sys_c = _build_group_coarse_sys(model, fine_voxels, groups, means, K_x, K_y)
    bc_c  = x -> _joint_active_bc(x, n_max_c)

    sp_c = _product_state_from_dists(coarse_d, tol, max_states)
    sp_c === nothing && return false
    expand!(sp_c, sys_c, bc_c; depth=1)
    length(sp_c) > max_states && return false

    success = _flux_vcycle!(coarse_d, rep_voxels, fg_c, means, model, K_x, K_y,
                             dt, n_max_c, tol, ε_prune, krylov_m, max_states, depth-1, t)

    if !success
        # Fallback: direct solve at this coarse level
        A_c, = build_generator(sp_c, sys_c, nothing, t)
        sp_c.probs .= expv(Float64(dt), A_c, sp_c.probs; m=min(krylov_m, size(A_c,1)))
        prune_threshold!(sp_c, tol); renormalize!(sp_c)
        coarse_d = _joint_marginals(sp_c, n_max_c)
    end
    # coarse_d now holds the post-solve K_G coarse marginals

    # ── Prolong: Multinomial split weighted by pre-solve means ────────────────
    for (j, g) in enumerate(groups)
        fine_margs = _multinomial_split(coarse_d[j], grp_means[j], n_max, tol)
        for (i, k) in enumerate(g); dists[k] = fine_margs[i]; end
    end

    true
end

# ─── Main entry point ─────────────────────────────────────────────────────────

"""
    evolve_flux_vcycle!(s, fg, means, dt; krylov_m, max_depth)

Apply the flux-based adaptive V-cycle to all active voxels in `s`.

1. Updates FluxGraph weights from current per-voxel means.
2. Runs recursive V-cycle:
   - Greedy max-weight matching → groups
   - Convolve group marginals → coarse marginals
   - Recurse (depth levels) → coarser solve
   - Multinomial-split → per-voxel marginals
3. Writes updated marginals back to s.dists.

The FluxGraph is maintained externally so it accumulates history across steps.
"""
function evolve_flux_vcycle!(
    s          :: SpatialFSP,
    fg         :: FluxGraph,
    means      :: Dict{CartesianIndex{2}, Float64},
    dt         :: Float64;
    krylov_m   :: Int   = 30,
    max_depth  :: Int   = 6
)
    active_voxels = _active_voxels(s)
    isempty(active_voxels) && return false

    # Rebuild flux graph if active set changed
    if fg.voxels != active_voxels
        new_fg = build_flux_graph(active_voxels, s.K_x, s.K_y; α = fg.α)
        # Transfer weights for edges that still exist
        old_pos = Dict(v => k for (k, v) in enumerate(fg.voxels))
        old_emap = Dict((min(i,j), max(i,j)) => e for (e,(i,j)) in enumerate(fg.edges))
        for (e2, (i2, j2)) in enumerate(new_fg.edges)
            vi = new_fg.voxels[i2]; vj = new_fg.voxels[j2]
            oi = get(old_pos, vi, 0); oj = get(old_pos, vj, 0)
            oi == 0 || oj == 0 && continue
            key = (min(oi,oj), max(oi,oj))
            e_old = get(old_emap, key, 0)
            e_old > 0 && (new_fg.weights[e2] = fg.weights[e_old])
        end
        copy!(fg.voxels, new_fg.voxels)
        copy!(fg.edges,   new_fg.edges)
        resize!(fg.weights, length(new_fg.weights))
        copy!(fg.weights, new_fg.weights)
    end

    update_flux!(fg, s.dists, jump_rate(s.model), dt)

    dists = [s.dists[v] for v in active_voxels]

    success = _flux_vcycle!(dists, active_voxels, fg, means, s.model,
                             s.K_x, s.K_y, dt, s.n_max, s.joint_state_tol,
                             s.ε_prune, krylov_m, s.joint_max_states, max_depth, s.t)

    success || return false

    for (k, v) in enumerate(active_voxels)
        s.dists[v] = dists[k]
    end
    s.joint_state = nothing; empty!(s.joint_voxels)
    true
end
