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
    CoarseningMapFull(K) -> CoarseningMap

Merges all K voxels into a single super-voxel.
Useful for adaptive coarsening when all voxels are in the same phase
(fast-diffusion limit: all voxels correlated).
"""
function CoarseningMapFull(K::Int)
    CoarseningMap(K, 1, K, fill(1, K), [collect(1:K)])
end

"""
    build_adaptive_cmap_1d(means, n_lo, n_hi; θ_lo, θ_hi, min_front) -> CoarseningMap

Build an adaptive 1D CoarseningMap from per-voxel means.

Voxels with mean near n_lo (within θ_lo) or near n_hi (within θ_hi) are "stable"
and merged into contiguous super-voxels.  Voxels in the transition region
(n_lo + θ_lo ≤ μ_k ≤ n_hi - θ_hi) are kept at fine resolution as the front.
`min_front` ensures at least this many voxels around the steepest gradient are kept fine.
"""
function build_adaptive_cmap_1d(means::Vector{Float64}, n_lo::Int, n_hi::Int;
                                  θ_lo::Float64 = 20.0, θ_hi::Float64 = 20.0,
                                  min_front::Int = 2,
                                  max_front::Int = typemax(Int))
    K = length(means)
    labels = similar(means, Symbol)   # :lo, :hi, :front

    for k in 1:K
        μ = means[k]
        if μ ≤ n_lo + θ_lo
            labels[k] = :lo
        elseif μ ≥ n_hi - θ_hi
            labels[k] = :hi
        else
            labels[k] = :front
        end
    end

    # Ensure at least min_front voxels around each lo→hi or hi→lo transition are :front
    for k in 1:K-1
        if labels[k] != labels[k+1] && labels[k] != :front && labels[k+1] != :front
            for δ in max(1,k-min_front+1):min(K,k+min_front)
                labels[δ] = :front
            end
        end
    end

    # Assign region indices
    assignments = zeros(Int, K)
    region = 0
    prev = :none

    # Count front voxels to compute group size for max_front cap
    n_front = count(==(:(front)), labels)
    front_group_size = max_front < typemax(Int) && n_front > 0 ?
                       ceil(Int, n_front / min(max_front, n_front)) : 1
    front_counter = 0

    for k in 1:K
        if labels[k] == :front
            front_counter += 1
            # Start a new region at the beginning of each group
            if front_counter == 1 || (front_counter - 1) % front_group_size == 0
                region += 1
            end
            assignments[k] = region
            prev = :front
        elseif labels[k] != prev
            region += 1
            assignments[k] = region
            prev = labels[k]
        else
            assignments[k] = region
        end
    end

    CoarseningMapRegions(assignments)
end

"""
    CoarseningMapRegions(assignments) -> CoarseningMap

Build a CoarseningMap from a region-assignment vector.
`assignments[k]` = coarse voxel index (1-based) that fine voxel k belongs to.
Regions must be contiguous integers starting at 1.
"""
function CoarseningMapRegions(assignments::Vector{Int})
    K = length(assignments)
    K_c = maximum(assignments)
    ctf = [Int[] for _ in 1:K_c]
    for (k, j) in enumerate(assignments)
        push!(ctf[j], k)
    end
    patch_size = maximum(length(r) for r in ctf)
    CoarseningMap(K, K_c, patch_size, assignments, ctf)
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

"""
    build_adaptive_cmap_2d(means, K_x, K_y; θ_lo, θ_hi) -> CoarseningMap

Build an adaptive 2D CoarseningMap from per-voxel means.

Voxels are classified as :lo (mean ≤ θ_lo), :hi (mean ≥ θ_hi), or :front.
Connected components of same-label voxels (4-connectivity) are merged into
coarse super-voxels.  Fine-resolution is kept only at the active wavefront
(the :front voxels and their immediate neighbourhood).

Linear index: k = (i-1)*K_y + j  for voxel (i,j).
"""
function build_adaptive_cmap_2d(means::Matrix{Float64}, K_x::Int, K_y::Int;
                                  θ_lo::Float64 = 1.0, θ_hi::Float64 = 8.5)
    # ── classify ──────────────────────────────────────────────────────────────
    labels = similar(means, Symbol)
    for i in 1:K_x, j in 1:K_y
        μ = means[i,j]
        labels[i,j] = μ ≤ θ_lo ? :lo : μ ≥ θ_hi ? :hi : :front
    end

    # ── connected-component BFS (4-connectivity, same label) ──────────────────
    # All three labels are merged by connected component:
    # :lo/:hi  → merged (low net activity: empty or at SS)
    # :front   → also merged by connectivity (the wavefront is one connected region)
    # For narrow-wavefront models (point source, bistable), :front forms a thin band
    # → small connected component → few extra regions.
    # For broad-wavefront models (global birth), :front covers most of grid
    # → one large connected component → K_eff stays small.
    component = zeros(Int, K_x, K_y)
    region = 0
    for i0 in 1:K_x, j0 in 1:K_y
        component[i0,j0] != 0 && continue
        region += 1
        lab = labels[i0,j0]
        queue = Tuple{Int,Int}[(i0, j0)]
        component[i0,j0] = region
        while !isempty(queue)
            (i, j) = popfirst!(queue)
            for (ni,nj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
                1 ≤ ni ≤ K_x && 1 ≤ nj ≤ K_y || continue
                component[ni,nj] != 0 && continue
                labels[ni,nj] == lab   || continue
                component[ni,nj] = region
                push!(queue, (ni,nj))
            end
        end
    end

    # ── build CoarseningMap ───────────────────────────────────────────────────
    K = K_x * K_y
    ftc = zeros(Int, K)
    for i in 1:K_x, j in 1:K_y
        ftc[(i-1)*K_y + j] = component[i,j]
    end
    K_c = region
    ctf = [Int[] for _ in 1:K_c]
    for k in 1:K; push!(ctf[ftc[k]], k); end
    CoarseningMap(K, K_c, maximum(length, ctf), ftc, ctf)
end
