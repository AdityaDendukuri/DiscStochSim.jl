"""
Voxel graph
===========

Spatial connectivity for graph-general RDME: a `VoxelGraph` (nodes + edges)
plus the standard chain / 2-D / 3-D constructors. The joint RDME generator on
any such graph is assembled by `build_rdme_joint` (unified_rdme.jl); the SA-AMG
stepper (`UnifiedFSP`) then drives it unchanged.

Public API
----------
  VoxelGraph(n, edges)            – arbitrary graph
  chain(K)                        – 1-D chain (K voxels)
  grid_2d(Kx, Ky)                 – 2-D rectangular grid
  grid_3d(Kx, Ky, Kz)             – 3-D rectangular grid
  degrees, adjacency_list         – graph queries
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
