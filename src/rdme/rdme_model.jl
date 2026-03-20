"""
    RDMEModel1D

1D single-species RDME with birth-death reactions and diffusion.

State type: `CartesianIndex{K}` where component `k` holds the molecule count
in voxel `k`.  This slots directly into the existing `StateSpace{CartesianIndex{K},T}`
and `build_generator` machinery.

Fields:
- `D`          : diffusion coefficient [length²/time]
- `birth_rate` : zero-order production rate per voxel [molecules/time]
- `death_rate` : first-order degradation rate per molecule [1/time]
"""
struct RDMEModel1D
    D::Float64
    birth_rate::Float64
    death_rate::Float64
end

# ─── unit stoich vector helpers ───────────────────────────────────────────────

"""Return CartesianIndex{K} with +1 in position `k`, 0 elsewhere."""
@inline _e(k::Int, ::Val{K}) where {K} =
    CartesianIndex(ntuple(j -> j == k ? 1 : 0, Val(K)))

"""Boundary condition: all counts in [0, n_max]."""
function rdme_bc(state::CartesianIndex{K}, n_max::Int) where {K}
    all(c -> 0 ≤ c ≤ n_max, Tuple(state))
end

# ─── full RDME system (reactions + all diffusion) ─────────────────────────────

"""
    build_rdme_system(model, grid) -> DiscreteStochasticSystem{CartesianIndex{K}}

Build the `DiscreteStochasticSystem` representing the full RDME:
  - Birth  ∅ → A  in each voxel at rate `birth_rate`
  - Death  A → ∅  in each voxel at rate `death_rate · nₖ`
  - Diffusion A_k ⇌ A_{k±1}  at rate `D/dx²`

The returned system can be passed directly to `build_generator`.
"""
function build_rdme_system(model::RDMEModel1D, grid::VoxelGrid)
    K   = grid.n_voxels
    d   = diffusion_rate(model.D, grid)
    k_b = model.birth_rate
    k_d = model.death_rate

    stoichs     = CartesianIndex{K}[]
    propensities = Function[]

    # ── per-voxel reactions ──────────────────────────────────────────────────
    for k in 1:K
        let k = k
            ek = _e(k, Val(K))

            # Birth: +eₖ
            push!(stoichs, ek)
            push!(propensities, (x, rates, t) -> k_b)

            # Death: −eₖ
            push!(stoichs, -ek)
            push!(propensities, (x, rates, t) -> k_d * max(0, Tuple(x)[k]))
        end
    end

    # ── diffusion between adjacent fine voxels ───────────────────────────────
    for k in 1:K-1
        let k = k
            ek  = _e(k,   Val(K))
            ek1 = _e(k+1, Val(K))

            # k → k+1
            push!(stoichs, -ek + ek1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[k]))

            # k+1 → k
            push!(stoichs, ek - ek1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[k+1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

# ─── intra-coarse-pair diffusion system (smoother) ────────────────────────────

"""
    build_intra_system(model, fine_grid, coarse_grid)
        -> DiscreteStochasticSystem{CartesianIndex{K}}

Build a `DiscreteStochasticSystem` containing only the diffusion jumps
*within* each coarse voxel pair  {2j-1, 2j}.  Used as the pre-/post-smoother
in the multigrid V-cycle: evolving this operator for time τ equilibrates the
within-pair molecule distribution towards Binomial(n̄ⱼ, 1/2).
"""
function build_intra_system(model::RDMEModel1D, fine_grid::VoxelGrid,
                             coarse_grid::VoxelGrid)
    K  = fine_grid.n_voxels
    K2 = coarse_grid.n_voxels
    @assert K == 2 * K2
    d = diffusion_rate(model.D, fine_grid)  # fine-grid rate D/h²

    stoichs      = CartesianIndex{K}[]
    propensities = Function[]

    for j in 1:K2
        let j = j
            kL = 2j - 1   # left  fine voxel in pair j
            kR = 2j       # right fine voxel in pair j
            eL = _e(kL, Val(K))
            eR = _e(kR, Val(K))

            # kL → kR
            push!(stoichs, -eL + eR)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[kL]))

            # kR → kL
            push!(stoichs, eL - eR)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[kR]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

# ─── coarse-level RDME system ─────────────────────────────────────────────────

"""
    build_coarse_system(model, coarse_grid)
        -> DiscreteStochasticSystem{CartesianIndex{K2}}

Build the coarse-level RDME system on the coarse grid.  The coarse diffusion
rate is `D / (2h)²` (automatically correct since `coarse_grid.dx = 2h`).
Birth-death rates are unchanged (linear reactions coarsen exactly,
Proposition 4.1 of the paper).
"""
function build_coarse_system(model::RDMEModel1D, coarse_grid::VoxelGrid)
    # Structurally identical to the fine system, just on the coarser grid.
    # `build_rdme_system` uses `grid.dx` for the diffusion rate, so
    # coarse_grid (dx = 2h) automatically gives the correct rate D/(2h)².
    build_rdme_system(model, coarse_grid)
end
