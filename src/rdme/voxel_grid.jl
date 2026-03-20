"""
    VoxelGrid

1D uniform spatial grid for the RDME.  The voxel side length is `dx` and
there are `n_voxels` voxels at refinement `level` (0 = finest).
"""
struct VoxelGrid
    n_voxels::Int
    dx::Float64
    level::Int
end

"""
    coarsen(g) -> VoxelGrid

Merge adjacent pairs: returns a grid with half as many voxels and double `dx`.
`n_voxels` must be even.
"""
function coarsen(g::VoxelGrid)
    iseven(g.n_voxels) || error("n_voxels = $(g.n_voxels) must be even to coarsen")
    VoxelGrid(g.n_voxels ÷ 2, g.dx * 2, g.level + 1)
end

"""
    build_hierarchy(g, n_levels) -> Vector{VoxelGrid}

Return a vector of `n_levels + 1` grids: `levels[1]` is the finest (= `g`),
`levels[end]` is the coarsest.
"""
function build_hierarchy(g::VoxelGrid, n_levels::Int)
    levels = Vector{VoxelGrid}(undef, n_levels + 1)
    levels[1] = g
    for ℓ in 2:n_levels+1
        levels[ℓ] = coarsen(levels[ℓ-1])
    end
    levels
end

"""
    diffusion_rate(D, g) -> Float64

Discrete diffusion jump rate `D / dx²` for diffusion coefficient `D` on grid `g`.
"""
diffusion_rate(D::Real, g::VoxelGrid) = D / g.dx^2
