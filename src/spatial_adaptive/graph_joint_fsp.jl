"""
graph_joint_fsp.jl

Persisted joint CME FSP over an arbitrary `VoxelGraph` (shape-agnostic).

A single joint `StateSpace` over `CartesianIndex{K}` (K = graph.n_voxels) is grown
by rstep (`expand!`) and trimmed to its realised support (`prune_threshold!`).
There is NO per-voxel `n_max` closure and NO product/marginal reconstruction: the
genuine joint distribution — and therefore every voxel–voxel correlation — is
carried across time steps. Geometry enters only through `graph.edges`, so the same
solver runs on a ring, a 2-D grid, a 3-D lattice, or any hand-built topology.

This is the fine-level kernel of the graph-multigrid hybrid. Phase 2 layers
graph-aggregation coarsening + reservoir lumping on top of this core to scale to
node counts whose full joint is intractable.
"""

# ─── Struct ───────────────────────────────────────────────────────────────────

mutable struct GraphJointFSP{K}
    graph    :: VoxelGraph
    sys      :: DiscreteStochasticSystem{CartesianIndex{K}}
    sp       :: StateSpace{CartesianIndex{K}, Float64}
    bc       :: Function          # boundary condition for rstep (default: nonnegativity only)
    t        :: Float64
    ε_prune  :: Float64
    krylov_m :: Int
end

"""
    GraphJointFSP(graph, sys; bc=nonneg, ε_prune=1e-8, krylov_m=30)

Build a persisted-joint FSP for the channel system `sys` (e.g. from
`build_schlogl_graph_system`) on `graph`. `bc` is the rstep boundary condition;
the default admits any nonnegative state (support is bounded by pruning, not by a
fixed n_max). Pass a conservation predicate (e.g. `x -> sum(Tuple(x)) <= N`) for
mass-conserved systems — that is conservation, not an n_max cap.
"""
function GraphJointFSP(graph::VoxelGraph,
                       sys::DiscreteStochasticSystem{CartesianIndex{K}};
                       bc::Function = x -> all(>=(0), Tuple(x)),
                       ε_prune::Float64 = 1e-8,
                       krylov_m::Int = 30) where {K}
    graph.n_voxels == K || error("graph.n_voxels=$(graph.n_voxels) ≠ system dim K=$K")
    sp = StateSpace{CartesianIndex{K}, Float64}()
    GraphJointFSP{K}(graph, sys, sp, bc, 0.0, ε_prune, krylov_m)
end

# ─── IC ─────────────────────────────────────────────────────────────────────

"""
    set_ic!(f, u0::CartesianIndex{K}, p=1.0)

Seed the joint state `u0` with probability `p`.
"""
function set_ic!(f::GraphJointFSP{K}, u0::CartesianIndex{K}, p::Float64 = 1.0) where {K}
    add_state!(f.sp, u0, p)
    f
end

# ─── Step ─────────────────────────────────────────────────────────────────────

"""
    step!(f::GraphJointFSP, dt)

One persisted-joint step: rstep-expand the support, build the joint generator,
evolve with `expv`, clip negatives, prune to support, renormalize. No product
reconstruction — correlations persist by construction.
"""
function step!(f::GraphJointFSP{K}, dt::Float64) where {K}
    isempty(f.sp.states) && return f

    expand!(f.sp, f.sys, f.bc; depth = 1)

    A, = build_generator(f.sp, f.sys, nothing, f.t; absorbing = false)
    m  = clamp(f.krylov_m, 4, length(f.sp))
    f.sp.probs .= expv(dt, A, f.sp.probs; m = m)
    f.sp.probs .= max.(0.0, f.sp.probs)

    prune_threshold!(f.sp, f.ε_prune)
    renormalize!(f.sp)
    f.t += dt
    f
end

# ─── Observables (read straight from the joint; means/cov come OUT, never in) ──

n_states(f::GraphJointFSP) = length(f.sp)

"""
    graph_voxel_means(f) -> Vector{Float64}  (length K)

Exact per-voxel mean ⟨n_v⟩ = Σ n·p from the realised joint.
"""
function graph_voxel_means(f::GraphJointFSP{K}) where {K}
    μ = zeros(K)
    for (s, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(s)
        @inbounds for v in 1:K
            μ[v] += t[v] * p
        end
    end
    μ
end

"""
    graph_voxel_cov(f) -> Matrix{Float64}  (K × K)

Exact voxel–voxel covariance `Cov(n_i,n_j)` from the joint. A product/marginal
solver returns this as identically zero off-diagonal; a non-zero off-diagonal
here is direct evidence the joint (the multigrid target) is being represented.
"""
function graph_voxel_cov(f::GraphJointFSP{K}) where {K}
    μ  = zeros(K)
    M2 = zeros(K, K)
    for (s, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(s)
        @inbounds for i in 1:K
            ni = t[i]
            μ[i] += ni * p
            for j in i:K
                M2[i, j] += ni * t[j] * p
            end
        end
    end
    C = zeros(K, K)
    @inbounds for i in 1:K, j in i:K
        c = M2[i, j] - μ[i] * μ[j]
        C[i, j] = c; C[j, i] = c
    end
    C
end

"""
    graph_voxel_marginal(f, v; n_max=nothing) -> Vector{Float64}  (length n_max+1)

Exact per-voxel marginal `P_v(n) = Σ_{other voxels} P(joint)`, read straight from the
realised joint and indexed `n = 0 … n_max` (defaults to the largest count voxel `v`
attains in the current support). This is the bimodality diagnostic: a voxel holding
mass in BOTH basins (lo and hi) returns a two-peaked marginal — the quantity the
restriction story is tested against, and the per-voxel view the click-a-voxel widget
draws. The forward map joint→marginal is exact but lossy (drops correlation); it does
NOT invert (rebuilding a joint from marginals is the product closure, the q(Y) failure).
"""
function graph_voxel_marginal(f::GraphJointFSP{K}, v::Int;
                              n_max::Union{Int, Nothing} = nothing) where {K}
    isempty(f.sp.states) && return Float64[]
    nm = n_max === nothing ? maximum(s[v] for s in f.sp.states) : n_max
    m  = zeros(nm + 1)
    for (s, p) in zip(f.sp.states, f.sp.probs)
        n = s[v]
        n <= nm && (m[n + 1] += p)
    end
    m
end

"""
    graph_voxel_marginals(f; n_max=nothing) -> Vector{Vector{Float64}}  (length K)

All `K` per-voxel marginals in a single pass over the joint, on a SHARED `n_max` so
they are directly comparable / stackable for the widget. `n_max` defaults to the
global max count over all voxels in the current support.
"""
function graph_voxel_marginals(f::GraphJointFSP{K};
                               n_max::Union{Int, Nothing} = nothing) where {K}
    isempty(f.sp.states) && return [Float64[] for _ in 1:K]
    nm = n_max === nothing ? maximum(maximum(Tuple(s)) for s in f.sp.states) : n_max
    M  = [zeros(nm + 1) for _ in 1:K]
    for (s, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(s)
        @inbounds for v in 1:K
            n = t[v]
            n <= nm && (M[v][n + 1] += p)
        end
    end
    M
end

"""
    graph_config_prob(f, n_thresh) -> Vector{Float64}  (length K+1)

P(exactly m voxels have n > n_thresh). The bistable/front configuration
indicator; only available from the joint.
"""
function graph_config_prob(f::GraphJointFSP{K}, n_thresh::Int) where {K}
    p = zeros(K + 1)
    for (s, pr) in zip(f.sp.states, f.sp.probs)
        m = count(>(n_thresh), Tuple(s))
        p[m + 1] += pr
    end
    p
end
