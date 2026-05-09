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

# ─── 2D grid ──────────────────────────────────────────────────────────────────

"""
    Grid2D

2D uniform spatial grid for the RDME.
Voxels indexed (i,j) with i in 1..K_x, j in 1..K_y.
Linear index: k = (i-1)*K_y + j.
"""
struct Grid2D
    K_x::Int
    K_y::Int
    dx::Float64   # voxel side length (square voxels)
    level::Int
end

n_voxels(g::Grid2D) = g.K_x * g.K_y

function coarsen2d(g::Grid2D)
    iseven(g.K_x) && iseven(g.K_y) || error("K_x=$(g.K_x), K_y=$(g.K_y) must both be even")
    Grid2D(g.K_x ÷ 2, g.K_y ÷ 2, g.dx * 2, g.level + 1)
end

diffusion_rate(D::Real, g::Grid2D) = D / g.dx^2

function build_hierarchy2d(g::Grid2D, n_levels::Int)
    levels = Vector{Grid2D}(undef, n_levels + 1)
    levels[1] = g
    for ℓ in 2:n_levels+1
        levels[ℓ] = coarsen2d(levels[ℓ-1])
    end
    levels
end

# ─── Coarsening map ───────────────────────────────────────────────────────────

"""
    CoarseningMap

Encodes which fine voxels form each coarse voxel (partition-based coarsening).

Fields:
- `n_fine`          : total number of fine voxels
- `n_coarse`        : total number of coarse voxels
- `patch_size`      : fine voxels per coarse voxel (2 for 1D pairs, 4 for 2D 2×2 blocks)
- `fine_to_coarse`  : length n_fine; fine_to_coarse[k] = coarse voxel index for fine voxel k
- `coarse_to_fine`  : length n_coarse; coarse_to_fine[J] = sorted list of fine voxel indices in patch J
"""
struct CoarseningMap
    n_fine::Int
    n_coarse::Int
    patch_size::Int
    fine_to_coarse::Vector{Int}
    coarse_to_fine::Vector{Vector{Int}}
end

"""
    CoarseningMap1D(K) -> CoarseningMap

Pairs adjacent 1D voxels: (1,2)→1, (3,4)→2, ...
"""
function CoarseningMap1D(K::Int)
    iseven(K) || error("K=$K must be even")
    K2 = K ÷ 2
    ftc = [ceil(Int, k / 2) for k in 1:K]
    ctf = [[2j-1, 2j] for j in 1:K2]
    CoarseningMap(K, K2, 2, ftc, ctf)
end

"""
    CoarseningMap2D(K_x, K_y) -> CoarseningMap

Groups 2×2 blocks of fine voxels into coarse voxels.
Fine voxel (i,j) → linear k = (i-1)*K_y + j.
Coarse voxel (I,J) → linear Jc = (I-1)*(K_y÷2) + J.
"""
function CoarseningMap2D(K_x::Int, K_y::Int)
    iseven(K_x) && iseven(K_y) || error("K_x=$K_x, K_y=$K_y must both be even")
    K2_x, K2_y = K_x ÷ 2, K_y ÷ 2
    n_fine = K_x * K_y
    n_coarse = K2_x * K2_y
    ftc = zeros(Int, n_fine)
    ctf = [Int[] for _ in 1:n_coarse]

    for i in 1:K_x, j in 1:K_y
        k  = (i - 1) * K_y + j
        I  = (i + 1) ÷ 2
        J  = (j + 1) ÷ 2
        Jc = (I - 1) * K2_y + J
        ftc[k] = Jc
        push!(ctf[Jc], k)
    end

    CoarseningMap(n_fine, n_coarse, 4, ftc, ctf)
end
