"""
Graph-based RDME multigrid FSP
==============================

Generalises the 1-D Schlögl multigrid to an *arbitrary* voxel graph.
The spatial structure is described by a `VoxelGraph` (nodes + edges);
the same V-cycle machinery then works unchanged for 1-D chains, 2-D
grids, 3-D lattices, or any hand-crafted topology.

Public API
----------
  VoxelGraph(n, edges)            – arbitrary graph
  chain(K)                        – 1-D chain (K voxels)
  grid_2d(Kx, Ky)                 – 2-D rectangular grid
  grid_3d(Kx, Ky, Kz)             – 3-D rectangular grid

  build_schlogl_graph_system(model, graph) → DiscreteStochasticSystem
  solve_rdme_graph(model, graph, u0, tspan, alg) → RDMESolution

The V-cycle pairing is determined at runtime by a greedy maximum-weight
matching on diffusion-flux edge weights, identical to flux_vcycle.jl.
"""

# ── VoxelGraph ────────────────────────────────────────────────────────────────

"""
    VoxelGraph

Describes the spatial connectivity of a voxel RDME.

Fields
------
  n_voxels  – total number of voxels
  edges     – list of (i,j) adjacency pairs, i < j, 1-indexed
  positions – optional K×D matrix of voxel coordinates (for visualisation)
"""
struct VoxelGraph
    n_voxels  :: Int
    edges     :: Vector{Tuple{Int,Int}}
    positions :: Union{Nothing, Matrix{Float64}}
end

VoxelGraph(n::Int, edges::Vector{Tuple{Int,Int}}) = VoxelGraph(n, edges, nothing)

"""1-D chain of K voxels."""
function chain(K::Int)
    edges = [(k, k+1) for k in 1:K-1]
    pos   = reshape(Float64.(1:K), K, 1)
    VoxelGraph(K, edges, pos)
end

"""2-D rectangular K_x × K_y grid."""
function grid_2d(K_x::Int, K_y::Int)
    N   = K_x * K_y
    idx = (i,j) -> (i-1)*K_y + j
    edges = Tuple{Int,Int}[]
    for i in 1:K_x, j in 1:K_y
        j < K_y && push!(edges, (idx(i,j), idx(i,j+1)))
        i < K_x && push!(edges, (idx(i,j), idx(i+1,j)))
    end
    pos = Matrix{Float64}(undef, N, 2)
    for i in 1:K_x, j in 1:K_y
        pos[idx(i,j), :] = [Float64(i), Float64(j)]
    end
    VoxelGraph(N, edges, pos)
end

"""3-D rectangular K_x × K_y × K_z grid."""
function grid_3d(K_x::Int, K_y::Int, K_z::Int)
    N   = K_x * K_y * K_z
    idx = (i,j,k) -> (i-1)*K_y*K_z + (j-1)*K_z + k
    edges = Tuple{Int,Int}[]
    for i in 1:K_x, j in 1:K_y, k in 1:K_z
        k < K_z && push!(edges, (idx(i,j,k), idx(i,j,k+1)))
        j < K_y && push!(edges, (idx(i,j,k), idx(i,j+1,k)))
        i < K_x && push!(edges, (idx(i,j,k), idx(i+1,j,k)))
    end
    pos = Matrix{Float64}(undef, N, 3)
    for i in 1:K_x, j in 1:K_y, k in 1:K_z
        pos[idx(i,j,k), :] = [Float64(i), Float64(j), Float64(k)]
    end
    VoxelGraph(N, edges, pos)
end

"""Degree (number of neighbours) of each voxel."""
function degrees(g::VoxelGraph)
    deg = zeros(Int, g.n_voxels)
    for (i,j) in g.edges
        deg[i] += 1; deg[j] += 1
    end
    deg
end

"""Adjacency list: neighbours[k] = list of voxel indices adjacent to k."""
function adjacency_list(g::VoxelGraph)
    nb = [Int[] for _ in 1:g.n_voxels]
    for (i,j) in g.edges
        push!(nb[i], j); push!(nb[j], i)
    end
    nb
end

# ── Generator ─────────────────────────────────────────────────────────────────

"""
    build_schlogl_graph_system(model, graph) → DiscreteStochasticSystem

Build the Schlögl RDME generator for `model` on `graph`.

State space: CartesianIndex{K} where K = graph.n_voxels.
Component k of the index is the molecule count in voxel k.

Reactions per voxel k (identical Schlögl kinetics, independent of geometry):
  ∅  →  X_k          at rate  c3
  X_k →  ∅           at rate  c4 · n_k
  2X_k → 3X_k        at rate  c1 · n_k(n_k-1)/2
  3X_k → 2X_k        at rate  c2 · n_k(n_k-1)(n_k-2)/6

Diffusion for each edge (i,j) in graph.edges:
  X_i → X_j          at rate  D · n_i
  X_j → X_i          at rate  D · n_j

where D = model.D (the diffusion coefficient; dx=1 assumed throughout,
so the jump rate d = D/dx² = D).
"""
function build_schlogl_graph_system(model::SchloglModel1D, graph::VoxelGraph;
                                    Dedge::Union{Nothing,AbstractVector} = nothing)
    K  = graph.n_voxels
    D  = model.D          # jump rate d = D (dx=1); overridden per-edge by Dedge
    c1 = model.c1; c2 = model.c2; c3 = model.c3; c4 = model.c4
    Dedge === nothing || length(Dedge) == length(graph.edges) ||
        error("Dedge must have one entry per graph edge ($(length(graph.edges)))")

    stoichs      = CartesianIndex{K}[]
    propensities = Function[]

    eu(k) = CartesianIndex(ntuple(i -> i == k ? 1 : 0, K))

    # Per-voxel Schlögl reactions
    for k in 1:K
        let k=k, ek=eu(k)
            push!(stoichs,  ek); push!(propensities, (x,rv,t) -> c3)
            push!(stoichs, -ek); push!(propensities, (x,rv,t) -> c4 * max(0, Tuple(x)[k]))
            push!(stoichs,  ek); push!(propensities,
                (x,rv,t) -> begin nk=Tuple(x)[k]; c1*max(0,nk)*max(0,nk-1)/2; end)
            push!(stoichs, -ek); push!(propensities,
                (x,rv,t) -> begin nk=Tuple(x)[k]; c2*max(0,nk)*max(0,nk-1)*max(0,nk-2)/6; end)
        end
    end

    # Diffusion along each edge (bidirectional); per-edge rate if Dedge given
    for (e, (i,j)) in enumerate(graph.edges)
        let i=i, j=j, ei=eu(i), ej=eu(j), De = Dedge === nothing ? D : float(Dedge[e])
            push!(stoichs, -ei+ej); push!(propensities,
                (x,rv,t) -> De * max(0, Tuple(x)[i]))
            push!(stoichs,  ei-ej); push!(propensities,
                (x,rv,t) -> De * max(0, Tuple(x)[j]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

# ── Graph-aware V-cycle pairing ───────────────────────────────────────────────

"""
    graph_flux_pairs(means, graph, dt, α) → (pairs, singletons)

Build flux-weighted adjacency weights  w_{ij} = D · min(μ_i, μ_j) · dt
and greedily select a maximum-weight non-overlapping matching.

Returns (pairs, singletons) as lists of voxel-index tuples/integers.
"""
function graph_flux_pairs(means::Vector{Float64}, graph::VoxelGraph,
                           D::Float64, dt::Float64)
    K      = graph.n_voxels
    n_e    = length(graph.edges)
    used   = fill(false, K)
    weights = [D * min(means[i], means[j]) * dt for (i,j) in graph.edges]
    order   = sortperm(weights; rev=true)

    pairs     = Tuple{Int,Int}[]
    for e in order
        i, j = graph.edges[e]
        used[i] || used[j] || begin
            push!(pairs, (i,j)); used[i] = used[j] = true
        end
    end
    pairs, findall(.!used)
end

# ── Coarse graph ──────────────────────────────────────────────────────────────

"""
    coarsen_graph(graph, pairs, singletons) → (coarse_graph, groups, node_map)

Build the coarse graph by contracting each pair to a super-voxel.

groups[j]   = vector of fine voxel indices in super-voxel j  (length 2 or 1)
node_map[k] = super-voxel index of fine voxel k
"""
function coarsen_graph(graph::VoxelGraph,
                        pairs::Vector{Tuple{Int,Int}},
                        singletons::Vector{Int})
    groups   = vcat([[p[1],p[2]] for p in pairs], [[s] for s in singletons])
    K_c      = length(groups)
    node_map = zeros(Int, graph.n_voxels)
    for (j,g) in enumerate(groups), k in g; node_map[k] = j; end

    # Coarse edges: (j,l) if any fine edge connects group j to group l
    coarse_edge_set = Set{Tuple{Int,Int}}()
    for (i,j_fine) in graph.edges
        gi = node_map[i]; gj = node_map[j_fine]
        gi != gj && push!(coarse_edge_set, minmax(gi,gj))
    end

    coarse_edges = collect(coarse_edge_set)

    # Coarse positions: centroid of constituent fine voxels
    coarse_pos = if graph.positions !== nothing
        D_pos = size(graph.positions, 2)
        pos_c = zeros(K_c, D_pos)
        for (j,g) in enumerate(groups)
            pos_c[j,:] = mean(graph.positions[k,:] for k in g)
        end
        pos_c
    else
        nothing
    end

    VoxelGraph(K_c, coarse_edges, coarse_pos), groups, node_map
end

# ── Graph RDME solve ──────────────────────────────────────────────────────────

"""
    GraphRDMESolution

Lightweight result type returned by `solve_rdme_graph`.
"""
struct GraphRDMESolution
    t                :: Vector{Float64}
    snapshots        :: Vector{Tuple{Vector,Vector{Float64}}}
    state_space_sizes:: Vector{Int}
    graph            :: VoxelGraph
end

"""
    mean_voxel_counts(sol::GraphRDMESolution) → Matrix  (n_t × K)
"""
function mean_voxel_counts(sol::GraphRDMESolution)
    K   = sol.graph.n_voxels
    out = zeros(length(sol.t), K)
    for (ti, (states, probs)) in enumerate(sol.snapshots)
        for (s, pr) in zip(states, probs)
            t = Tuple(s)
            for k in 1:K; out[ti,k] += t[k] * pr; end
        end
    end
    out
end

"""
    config_probs(sol::GraphRDMESolution, n_thresh, ti) → Vector  (length K+1)

P(exactly m voxels with n > n_thresh) at snapshot index ti.
"""
function config_probs(sol::GraphRDMESolution, n_thresh::Int, ti::Int)
    K = sol.graph.n_voxels
    p = zeros(K+1)
    states, probs = sol.snapshots[ti]
    for (s, pr) in zip(states, probs)
        m = count(n -> n > n_thresh, Tuple(s))
        p[m+1] += pr
    end
    p
end

"""
    solve_rdme_graph(model, graph, u0, tspan, alg) → GraphRDMESolution

Run the multigrid FSP on an arbitrary VoxelGraph.

Arguments
---------
  model  – SchloglModel1D (kinetics + diffusion coefficient)
  graph  – VoxelGraph describing spatial connectivity
  u0     – CartesianIndex{K} initial condition (molecule counts per voxel)
  tspan  – (t0, t_end)
  alg    – RDMEMultigridFSP parameter struct (dt, n_max, prune_tol, etc.)

The V-cycle pairing is recomputed at each step from current marginal means,
using greedy maximum-weight matching on diffusion-flux edge weights.

At each time step:
  1. Compute per-voxel marginal means from the joint FSP state.
  2. Build flux weights: w_{ij} = D · min(μ_i, μ_j) · dt.
  3. Greedy matching → pairs + singletons.
  4. Restrict: convolve within-pair marginals → coarse 1-D dists.
  5. Build coarse graph system on K_c super-voxels.
  6. Solve coarse FSP with expv.
  7. Prolong: multinomial split → per-voxel marginals.
  8. Reconstruct joint state via product approximation.
"""
function solve_rdme_graph(model::SchloglModel1D,
                           graph::VoxelGraph,
                           u0::CartesianIndex,
                           tspan::Tuple{Float64,Float64},
                           alg::RDMEMultigridFSP)
    K      = graph.n_voxels
    length(Tuple(u0)) == K || error("u0 dimension $(length(Tuple(u0))) ≠ graph.n_voxels $K")
    deg    = degrees(graph)
    nb     = adjacency_list(graph)
    D      = model.D
    dt     = Float64(alg.dt)
    n_max  = alg.n_max
    t0, tf = tspan

    sys_fine = build_schlogl_graph_system(model, graph)

    # Boundary condition: all voxel counts in [0, n_max]
    bc = x -> all(0 <= n <= n_max for n in Tuple(x))

    # Initialise state
    sp = StateSpace{CartesianIndex{K}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, sys_fine, bc; depth = alg.coarse_d_depth)
    prune_threshold!(sp, alg.prune_tol)

    t    = t0
    step = 0
    save_every = alg.save_every

    snapshots  = Tuple{Vector,Vector{Float64}}[]
    times      = Float64[]
    ss_sizes   = Int[]

    function save_snap!()
        push!(times, t)
        push!(snapshots, (copy(sp.states), copy(sp.probs)))
        push!(ss_sizes, length(sp))
    end

    save_snap!()

    while t < tf - dt/2
        dt_s = min(dt, tf - t)

        # ── Per-voxel marginal means ───────────────────────────────────────
        means = zeros(K)
        for (s, pr) in zip(sp.states, sp.probs)
            t_  = Tuple(s)
            for k in 1:K; means[k] += t_[k] * pr; end
        end

        # ── Graph flux pairing ─────────────────────────────────────────────
        pairs, singletons = graph_flux_pairs(means, graph, D, dt_s)

        # ── V-cycle step ───────────────────────────────────────────────────
        success = _graph_vcycle_step!(sp, sys_fine, pairs, singletons, graph,
                                       model, means, bc, dt_s, n_max,
                                       alg.max_states, alg.prune_tol,
                                       alg.weight_tol, alg.τ_pre)

        if !success
            # Fallback: direct expv on fine joint state
            isempty(sp.states) && add_state!(sp, u0, 1.0)   # safety: re-seed if cleared
            length(sp) < alg.max_states && expand!(sp, sys_fine, bc; depth=1)
            m_kry = min(30, length(sp))
            if !isempty(sp.states) && length(sp) <= alg.max_states && m_kry > 0
                A, = build_generator(sp, sys_fine, nothing, t)
                sp.probs .= expv(dt_s, A, sp.probs; m = m_kry)
            end
        end

        prune_threshold!(sp, alg.prune_tol)
        renormalize!(sp)
        t += dt_s; step += 1

        step % save_every == 0 && save_snap!()
    end

    t ≈ times[end] || save_snap!()

    GraphRDMESolution(times, snapshots, ss_sizes, graph)
end

# ── Internal: graph V-cycle step ─────────────────────────────────────────────

function _graph_vcycle_step!(sp::StateSpace{CartesianIndex{K}, Float64},
                              sys_fine,
                              pairs::Vector{Tuple{Int,Int}},
                              singletons::Vector{Int},
                              graph::VoxelGraph,
                              model::SchloglModel1D,
                              means::Vector{Float64},
                              bc,
                              dt::Float64,
                              n_max::Int,
                              max_states::Int,
                              prune_tol::Float64,
                              weight_tol::Float64,
                              τ_pre::Float64) where K

    isempty(pairs) && return false   # no pairs → fall through to direct solve

    groups   = vcat([[p[1],p[2]] for p in pairs], [[s] for s in singletons])
    K_c      = length(groups)
    node_map = zeros(Int, K)
    for (j,g) in enumerate(groups), k in g; node_map[k] = j; end

    gsz = [length(g) for g in groups]   # super-voxel sizes (1 or 2)

    # ── 1. Extract per-voxel marginals from joint state ───────────────────
    fine_margs = [zeros(n_max+1) for _ in 1:K]
    for (s, pr) in zip(sp.states, sp.probs)
        t_ = Tuple(s)
        for k in 1:K
            n = t_[k]; n <= n_max && (fine_margs[k][n+1] += pr)
        end
    end
    for p in fine_margs; s=sum(p); s>0 && (p./=s); end

    # ── 2. Restrict: convolve marginals within pairs → coarse 1-D dists ──
    coarse_d = Vector{Float64}[]
    pre_means = Vector{Float64}[]   # per-fine-voxel means before step (for prolong)
    for g in groups
        ps = [fine_margs[k] for k in g]
        μs = [means[k] for k in g]
        push!(pre_means, μs)
        if length(g) == 2
            p1, p2 = ps[1], ps[2]
            p_pair = zeros(length(p1)+length(p2)-1)
            for i1 in eachindex(p1), i2 in eachindex(p2)
                p_pair[i1+i2-1] += p1[i1]*p2[i2]
            end
            push!(coarse_d, p_pair)
        else
            push!(coarse_d, copy(ps[1]))
        end
    end

    # ── 3. Build coarse graph and system ──────────────────────────────────
    coarse_graph, _, _ = coarsen_graph(graph, pairs, singletons)
    coarse_model = SchloglModel1D(model.D/2, model.c1, model.c2,
                                   model.c3, model.c4)  # D/2 per fine voxel in pair

    # Effective coarse Schlögl model: reactions scale with group size
    # Build coarse system manually using group-aware generator
    n_max_c = minimum(length(p)-1 for p in coarse_d)   # actual coarse support
    n_max_c = max(n_max_c, n_max)

    # Build coarse product state
    coarse_tol = max(prune_tol, weight_tol)
    supp = [findall(>(coarse_tol), p) for p in coarse_d]
    prod_size = prod(Int128(length(s)) for s in supp; init=Int128(1))
    prod_size > max_states && return false

    sp_c = StateSpace{CartesianIndex{K_c}, Float64}()
    _build_coarse_product_state!(sp_c, coarse_d, supp, coarse_tol)
    isempty(sp_c.states) && return false

    # Coarse generator using group-aware coarse system
    sys_c = _build_graph_coarse_system(model, graph, groups, node_map, gsz, K_c)
    bc_c  = x -> all(0 <= n <= n_max_c for n in Tuple(x))
    expand!(sp_c, sys_c, bc_c; depth=1)
    length(sp_c) > max_states && return false

    A_c, = build_generator(sp_c, sys_c, nothing, 0.0)
    sp_c.probs .= expv(dt, A_c, sp_c.probs; m=min(30, size(A_c,1)))
    prune_threshold!(sp_c, coarse_tol)
    renormalize!(sp_c) > 0 || return false

    # ── 4. Extract coarse marginals ────────────────────────────────────────
    coarse_post = [zeros(n_max_c+1) for _ in 1:K_c]
    for (s,pr) in zip(sp_c.states, sp_c.probs)
        t_ = Tuple(s)
        for j in 1:K_c; n=t_[j]; n<=n_max_c && (coarse_post[j][n+1]+=pr); end
    end
    for p in coarse_post; s=sum(p); s>0 && (p./=s); end

    # ── 5. Prolong: multinomial split → per-voxel marginals ───────────────
    fine_post = copy(fine_margs)
    for (j,g) in enumerate(groups)
        p_N = coarse_post[j]
        μs  = pre_means[j]
        if length(g) == 2
            # Weighted binomial split proportional to pre-solve means
            W = sum(μs); W < 1e-10 && (W = 1.0)
            probs_split = [max(μ,1e-10)/W for μ in μs]
            p1 = zeros(n_max+1); p2 = zeros(n_max+1)
            for N in 0:length(p_N)-1
                p_N[N+1] > coarse_tol || continue
                for k1 in 0:min(N, n_max)
                    k2 = N - k1; k2 <= n_max || continue
                    binom_w = begin
                        p = probs_split[1]
                        if p <= 0.0; (k1==0 ? 1.0 : 0.0)
                        elseif p >= 1.0; (k1==N ? 1.0 : 0.0)
                        else
                            lnf = sum(log(Float64(i)) for i in 1:N; init=0.0)
                            lk1 = sum(log(Float64(i)) for i in 1:k1; init=0.0)
                            lk2 = sum(log(Float64(i)) for i in 1:k2; init=0.0)
                            exp(lnf - lk1 - lk2 + k1*log(p) + k2*log(1-p))
                        end
                    end
                    w = p_N[N+1] * binom_w
                    w > coarse_tol || continue
                    p1[k1+1] += w; p2[k2+1] += w
                end
            end
            for p in (p1,p2); s=sum(p); s>0 && (p./=s); end
            fine_post[g[1]] = p1; fine_post[g[2]] = p2
        else
            # Singleton: trim to n_max
            p = p_N[1:min(length(p_N), n_max+1)]
            if length(p) < n_max+1; append!(p, zeros(n_max+1-length(p))); end
            s=sum(p); s>0 && (p./=s)
            fine_post[g[1]] = p
        end
    end

    # ── 6. Reconstruct joint state via product approximation ──────────────
    new_sp = StateSpace{CartesianIndex{K}, Float64}()
    _product_joint_from_marginals!(new_sp, fine_post, n_max, prune_tol, max_states)
    isempty(new_sp.states) && return false

    expand!(new_sp, sys_fine, bc; depth=1)
    length(new_sp) > max_states && return false

    # Update in-place
    empty!(sp.states); empty!(sp.probs); empty!(sp.index)
    for (s,p) in zip(new_sp.states, new_sp.probs)
        add_state!(sp, s, p)
    end
    true
end

# ── Helpers ───────────────────────────────────────────────────────────────────

function _build_coarse_product_state!(sp_c::StateSpace{CartesianIndex{Kc}, Float64},
                                       dists, supp, tol) where Kc
    buf = zeros(Int, Kc)
    function rec(dim, w)
        dim > Kc && (add_state!(sp_c, CartesianIndex(ntuple(i->buf[i],Kc)), w); return)
        for idx in supp[dim]
            ww = w * dists[dim][idx]; ww > tol || continue
            buf[dim] = idx-1; rec(dim+1, ww)
        end
    end
    rec(1, 1.0)
    renormalize!(sp_c)
end

function _product_joint_from_marginals!(sp::StateSpace{CartesianIndex{K}, Float64},
                                         margs, n_max, tol, max_states) where K
    supp = [findall(>(tol), p) for p in margs]
    prod_size = prod(Int128(length(s)) for s in supp; init=Int128(1))
    prod_size > max_states && return
    buf = zeros(Int, K)
    function rec(dim, w)
        dim > K && (add_state!(sp, CartesianIndex(ntuple(i->buf[i],K)), w); return)
        for idx in supp[dim]
            ww = w * margs[dim][idx]; ww > tol || continue
            buf[dim] = idx-1; rec(dim+1, ww)
        end
    end
    rec(1, 1.0)
    renormalize!(sp)
end

function _build_graph_coarse_system(model::SchloglModel1D,
                                     fine_graph::VoxelGraph,
                                     groups::Vector{Vector{Int}},
                                     node_map::Vector{Int},
                                     gsz::Vector{Int},
                                     K_c::Int)
    D  = model.D
    c1 = model.c1; c2 = model.c2; c3 = model.c3; c4 = model.c4

    stoichs = CartesianIndex{K_c}[]
    props   = Function[]
    eu(j)   = CartesianIndex(ntuple(i -> i==j ? 1 : 0, K_c))

    # Per-super-voxel reactions (scale birth by group size)
    for j in 1:K_c
        let j=j, ej=eu(j), sz=Float64(gsz[j])
            push!(stoichs,  ej); push!(props, (x,rv,t) -> sz*c3)
            push!(stoichs, -ej); push!(props, (x,rv,t) -> c4*max(0,Tuple(x)[j]))
            push!(stoichs,  ej); push!(props,
                (x,rv,t) -> begin nj=Tuple(x)[j]; c1*max(0,nj)*max(0,nj-1)/2; end)
            push!(stoichs, -ej); push!(props,
                (x,rv,t) -> begin nj=Tuple(x)[j];
                    c2*max(0,nj)*max(0,nj-1)*max(0,nj-2)/6; end)
        end
    end

    # Inter-super-voxel diffusion: count fine edges between groups
    inter_edges = Dict{Tuple{Int,Int}, Int}()
    for (i,j) in fine_graph.edges
        gi = node_map[i]; gj = node_map[j]; gi == gj && continue
        key = minmax(gi, gj)
        inter_edges[key] = get(inter_edges, key, 0) + 1
    end
    for ((gi,gj), n_e) in inter_edges
        let gi=gi, gj=gj, ri=D*n_e/gsz[gi], rj=D*n_e/gsz[gj]
            ei=eu(gi); ej=eu(gj)
            push!(stoichs,-ei+ej); push!(props,(x,rv,t)->ri*max(0,Tuple(x)[gi]))
            push!(stoichs, ei-ej); push!(props,(x,rv,t)->rj*max(0,Tuple(x)[gj]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_c}}(stoichs, props)
end

# ── mean helper ───────────────────────────────────────────────────────────────
_mean_vec(v) = length(v) == 0 ? 0.0 : sum(v) / length(v)
