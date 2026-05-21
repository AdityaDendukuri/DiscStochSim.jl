"""
    SchloglModel1D

1D single-species Schlögl RDME: bistable autocatalytic reactions + diffusion.

Per-voxel reactions (X = molecule count in voxel k):
  ∅  →^{c3}          X        (zeroth-order production)
  X  →^{c4}          ∅        (first-order death)
  2X →^{c1·n(n-1)/2} 3X       (bimolecular autocatalysis)
  3X →^{c2·n(n-1)(n-2)/6} 2X  (trimolecular suppression)

Fields:
- `D`  : diffusion coefficient [length²/time]
- `c1` : autocatalysis rate constant   (propensity = c1·n(n-1)/2)
- `c2` : reverse rate constant         (propensity = c2·n(n-1)(n-2)/6)
- `c3` : zeroth-order production rate per voxel
- `c4` : first-order death rate per molecule

Default parameters give bistability with single-voxel stable states at
n_low ≈ 31 and n_high ≈ 148 (unstable fixed point at n ≈ 78).
"""
struct SchloglModel1D <: AbstractSpatialRDMEModel
    D::Float64
    c1::Float64
    c2::Float64
    c3::Float64
    c4::Float64
end

SchloglModel1D(D) = SchloglModel1D(D, 0.028, 3.2e-4, 19.5, 1.0)

"""
    schlogl_rate_eq(model, n) -> Float64

Deterministic rate equation f(n) = net production rate at molecule count n.
Roots are the stable (f'<0) and unstable (f'>0) fixed points.
"""
function schlogl_rate_eq(m::SchloglModel1D, n::Real)
    m.c1*n*(n-1)/2 - m.c2*n*(n-1)*(n-2)/6 + m.c3 - m.c4*n
end

"""
    schlogl_fixed_points(model; n_max=250) -> (n_low, n_unstable, n_high)

Find the three integer fixed points of the Schlögl rate equation by scanning
and locating sign changes.  Returns `nothing` for any missing root.
"""
function schlogl_fixed_points(m::SchloglModel1D; n_max::Int=250)
    f = n -> schlogl_rate_eq(m, n)
    roots = Int[]
    for n in 1:n_max-1
        if f(n) * f(n+1) < 0
            # Bisect to nearest integer
            push!(roots, abs(f(n)) < abs(f(n+1)) ? n : n+1)
        end
    end
    length(roots) == 3 ? Tuple(roots) : (get(roots,1,nothing), get(roots,2,nothing), get(roots,3,nothing))
end

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

"""
    coarsen_model(model, r) -> model

Scale model parameters for a grid coarsened by ratio `r = coarse_dx / fine_dx`.
"""
coarsen_model(m::RDMEModel1D, r::Real) = RDMEModel1D(m.D, m.birth_rate * r, m.death_rate)
coarsen_model(m::SchloglModel1D, r::Real) = SchloglModel1D(m.D, m.c1 / r, m.c2 / r^2, m.c3 * r, m.c4)

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
    build_coarse_system(model, coarse_grid, fine_grid)
        -> DiscreteStochasticSystem{CartesianIndex{K2}}

Build the coarse-level RDME system on the coarse grid.

Coarsening rules (Proposition 4.1 of the paper):
- **Diffusion**: rate D/(2h)² is automatic since `coarse_grid.dx = 2h`.
- **Zero-order birth** (∅ → A, rate k_b per voxel): each coarse voxel
  represents `r = coarse_dx / fine_dx` fine voxels, so the coarse birth
  rate is `r * k_b`.  For one coarsening level r = 2.
- **First-order death** (A → ∅, rate k_d per molecule): propensity is
  k_d * n̄ where n̄ is the coarse molecule count — unchanged.
"""
function build_coarse_system(model::RDMEModel1D,
                              coarse_grid::VoxelGrid,
                              fine_grid::VoxelGrid)
    r = coarse_grid.dx / fine_grid.dx   # coarsening ratio (= 2 for one level)
    scaled_model = RDMEModel1D(model.D, model.birth_rate * r, model.death_rate)
    build_rdme_system(scaled_model, coarse_grid)
end

# ─── Schlögl RDME system ───────────────────────────────────────────────────────

"""
    build_schlogl_rdme_system(model, grid) -> DiscreteStochasticSystem{CartesianIndex{K}}

Build the full fine-grid RDME for the Schlögl model:
  - Per voxel: zeroth-order birth, first-order death, autocatalysis, reverse
  - Between adjacent voxels: diffusion at rate D/dx²
"""
function build_schlogl_rdme_system(model::SchloglModel1D, grid::VoxelGrid)
    K  = grid.n_voxels
    d  = diffusion_rate(model.D, grid)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4

    stoichs      = CartesianIndex{K}[]
    propensities = Function[]

    for k in 1:K
        let k = k
            ek = _e(k, Val(K))

            # Production: ∅ → X
            push!(stoichs, ek)
            push!(propensities, (x, rates, t) -> c3)

            # Death: X → ∅
            push!(stoichs, -ek)
            push!(propensities, (x, rates, t) -> c4 * max(0, Tuple(x)[k]))

            # Autocatalysis: 2X → 3X  (net +1)
            push!(stoichs, ek)
            push!(propensities, (x, rates, t) -> begin
                n = Tuple(x)[k]; c1 * max(0, n*(n-1)) / 2
            end)

            # Reverse: 3X → 2X  (net −1)
            push!(stoichs, -ek)
            push!(propensities, (x, rates, t) -> begin
                n = Tuple(x)[k]; c2 * max(0, n*(n-1)*(n-2)) / 6
            end)
        end
    end

    for k in 1:K-1
        let k = k
            ek  = _e(k,   Val(K))
            ek1 = _e(k+1, Val(K))
            push!(stoichs, -ek + ek1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[k]))
            push!(stoichs,  ek - ek1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[k+1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

"""
    build_schlogl_adaptive_system(model, fine_grid, cmap, ::Val{K_eff};
                                   region_birth=nothing)
        -> DiscreteStochasticSystem{CartesianIndex{K_eff}}

Build a K_eff-voxel Schlögl system for an adaptive CoarseningMap with
heterogeneous patch sizes M_j = length(cmap.coarse_to_fine[j]).

Reactions within region j use Binomial-corrected rates:
  production:    M_j · c3   (or region_birth[j] if supplied)
  death:         c4 · n̄_j
  autocatalysis: (c1/M_j) · n̄_j(n̄_j-1)/2
  reverse:       (c2/M_j²) · n̄_j(n̄_j-1)(n̄_j-2)/6

Diffusion between adjacent regions j and j+1 uses the boundary-molecule
approximation (fast within-patch equilibrium):
  j → j+1:  (d/M_j)   · n̄_j    [right boundary of patch j]
  j+1 → j:  (d/M_{j+1}) · n̄_{j+1} [left boundary of patch j+1]

`region_birth`: optional length-K_eff vector of per-region production rates.
  Use to implement localized sources (set rate=0 for non-source regions).
  When nothing, defaults to M_j·c3 for each region.

For M_j=1 (a fine voxel) all rates reduce to the standard Schlögl RDME.
"""
function build_schlogl_adaptive_system(
    model::SchloglModel1D,
    fine_grid::VoxelGrid,
    cmap::CoarseningMap,
    ::Val{K_eff};
    region_birth::Union{Nothing, Vector{Float64}} = nothing
) where {K_eff}
    K_eff == cmap.n_coarse || error("K_eff=$K_eff ≠ cmap.n_coarse=$(cmap.n_coarse)")
    d            = diffusion_rate(model.D, fine_grid)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
    patch_sizes  = [length(cmap.coarse_to_fine[j]) for j in 1:K_eff]
    birth_rates  = region_birth === nothing ?
                   [Float64(patch_sizes[j]) * c3 for j in 1:K_eff] :
                   region_birth

    stoichs      = CartesianIndex{K_eff}[]
    propensities = Function[]

    for j in 1:K_eff
        let j = j, Mj = patch_sizes[j], bj = birth_rates[j]
            ej = _e(j, Val(K_eff))

            # Production: ∅ →(bj)→ n̄_j + 1
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> bj)

            # Death: n̄_j →(c4·n̄_j)→ n̄_j - 1
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> c4 * max(0, Tuple(x)[j]))

            # Autocatalysis: (c1/Mj)·n̄_j(n̄_j-1)/2
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> begin
                n = Tuple(x)[j]; (c1/Mj) * max(0, n*(n-1)) / 2
            end)

            # Reverse: (c2/Mj²)·n̄_j(n̄_j-1)(n̄_j-2)/6
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> begin
                n = Tuple(x)[j]; (c2/Mj^2) * max(0, n*(n-1)*(n-2)) / 6
            end)
        end
    end

    for j in 1:K_eff-1
        let j = j, Mj = patch_sizes[j], Mj1 = patch_sizes[j+1]
            ej  = _e(j,   Val(K_eff))
            ej1 = _e(j+1, Val(K_eff))

            # j → j+1: boundary-molecule rate (d/Mj)·n̄_j
            push!(stoichs, ej1 - ej)
            push!(propensities, (x, rates, t) -> (d/Mj) * max(0, Tuple(x)[j]))

            # j+1 → j: boundary-molecule rate (d/Mj1)·n̄_{j+1}
            push!(stoichs, ej - ej1)
            push!(propensities, (x, rates, t) -> (d/Mj1) * max(0, Tuple(x)[j+1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_eff}}(stoichs, propensities)
end

"""
    build_rdme_adaptive_2d(model, grid2d, cmap, ::Val{K_eff}; region_birth=nothing)
        -> DiscreteStochasticSystem{CartesianIndex{K_eff}}

Build a K_eff-region birth-death RDME system for an adaptive 2D CoarseningMap.

Reactions within region J use patch-size-corrected rates:
  production:    birth_J  (defaults to M_J · c3, or region_birth[J] if supplied)
  death:         c4 · n̄_J

Diffusion between adjacent regions uses the boundary-molecule approximation:
  J → J':  (d · n_adj(J,J') / M_J)  · n̄_J
where n_adj(J,J') is the number of fine voxel pairs crossing the J–J' boundary.
This generalises `build_schlogl_adaptive_system` to arbitrary 2D patch shapes.
"""
function build_rdme_adaptive_2d(
    model::SchloglModel1D,
    grid::Grid2D,
    cmap::CoarseningMap,
    ::Val{K_eff};
    region_birth::Union{Nothing, Vector{Float64}} = nothing
) where {K_eff}
    K_eff == cmap.n_coarse || error("K_eff=$K_eff ≠ cmap.n_coarse=$(cmap.n_coarse)")
    d   = diffusion_rate(model.D, grid)
    c3, c4 = model.c3, model.c4
    K_y = grid.K_y
    patch_sizes = [length(cmap.coarse_to_fine[j]) for j in 1:K_eff]
    birth_rates = region_birth === nothing ?
                  [Float64(patch_sizes[j]) * c3 for j in 1:K_eff] :
                  region_birth

    stoichs      = CartesianIndex{K_eff}[]
    propensities = Function[]

    # ── per-region reactions ──────────────────────────────────────────────────
    for j in 1:K_eff
        let j=j, bj=birth_rates[j]
            ej = _e(j, Val(K_eff))
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> bj)
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> c4 * max(0, Tuple(x)[j]))
        end
    end

    # ── inter-region diffusion via boundary-molecule approximation ────────────
    # Count shared boundary fine-voxel pairs for every ordered pair (J, J').
    # Two fine voxels k=(i-1)*K_y+j and k'=(i'-1)*K_y+j' are adjacent iff
    # they differ by exactly 1 in one coordinate.
    K_y_grid = grid.K_y
    n_adj = zeros(Int, K_eff, K_eff)
    for k in 1:cmap.n_fine
        i = (k - 1) ÷ K_y_grid + 1
        j = (k - 1) % K_y_grid + 1
        Jk = cmap.fine_to_coarse[k]
        # right neighbour
        if j < K_y_grid
            k2 = k + 1;  Jk2 = cmap.fine_to_coarse[k2]
            if Jk != Jk2; n_adj[Jk, Jk2] += 1; n_adj[Jk2, Jk] += 1; end
        end
        # down neighbour
        if i < grid.K_x
            k2 = k + K_y_grid;  Jk2 = cmap.fine_to_coarse[k2]
            if Jk != Jk2; n_adj[Jk, Jk2] += 1; n_adj[Jk2, Jk] += 1; end
        end
    end

    for ja in 1:K_eff, jb in (ja+1):K_eff
        n_adj[ja,jb] == 0 && continue
        let ja=ja, jb=jb, na=Float64(n_adj[ja,jb]), Ma=Float64(patch_sizes[ja]), Mb=Float64(patch_sizes[jb])
            eja = _e(ja, Val(K_eff)); ejb = _e(jb, Val(K_eff))
            push!(stoichs, ejb - eja)
            push!(propensities, (x, rates, t) -> (d * na / Ma) * max(0, Tuple(x)[ja]))
            push!(stoichs, eja - ejb)
            push!(propensities, (x, rates, t) -> (d * na / Mb) * max(0, Tuple(x)[jb]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K_eff}}(stoichs, propensities)
end

"""
    build_schlogl_coarse_system(model, coarse_grid, fine_grid)
        -> DiscreteStochasticSystem{CartesianIndex{K2}}

Galerkin coarse-grid Schlögl system: A_c = R A P where P uses Binomial(n_c, 1/2)
prolongation weights (valid when diffusion ≫ reactions within each pair).

Coarse propensities for r = coarse_dx/fine_dx = 2 (one coarsening level):

  a1_coarse(nc) = c1 · nc(nc-1)/4       [E_{Binom}(n_L(n_L-1)/2 + n_R(n_R-1)/2)]
  a2_coarse(nc) = c2 · nc(nc-1)(nc-2)/24 [E_{Binom}(n_L(n_L-1)(n_L-2)/6 + …)]
  a3_coarse     = r · c3                 [r fine voxels each producing at rate c3]
  a4_coarse(nc) = c4 · nc               [first-order, unchanged per molecule]

The derivation uses falling factorial moments of Binom(nc, 1/2):
  E[n_L^{(m)}] = nc^{(m)} / 2^m,  so  E[n_L(n_L-1)] = nc(nc-1)/4, etc.

These coarse rates differ from the fine rates at the Schlögl bistable front,
where the actual conditional distribution is bimodal (not Binomial).
"""
function build_schlogl_coarse_system(model::SchloglModel1D,
                                      coarse_grid::VoxelGrid,
                                      fine_grid::VoxelGrid)
    r  = coarse_grid.dx / fine_grid.dx    # = 2 for one coarsening level
    K2 = coarse_grid.n_voxels
    d  = diffusion_rate(model.D, coarse_grid)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4

    stoichs      = CartesianIndex{K2}[]
    propensities = Function[]

    for j in 1:K2
        let j = j
            ej = _e(j, Val(K2))

            # Production (r fine voxels per coarse voxel)
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> r * c3)

            # Death (unchanged per molecule)
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> c4 * max(0, Tuple(x)[j]))

            # Autocatalysis: Binomial-corrected  → c1·nc(nc-1)/4
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> begin
                nc = Tuple(x)[j]; c1 * max(0, nc*(nc-1)) / 4
            end)

            # Reverse trimolecular: Binomial-corrected → c2·nc(nc-1)(nc-2)/24
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> begin
                nc = Tuple(x)[j]; c2 * max(0, nc*(nc-1)*(nc-2)) / 24
            end)
        end
    end

    for j in 1:K2-1
        let j = j
            ej  = _e(j,   Val(K2))
            ej1 = _e(j+1, Val(K2))
            push!(stoichs, -ej + ej1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[j]))
            push!(stoichs,  ej - ej1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[j+1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K2}}(stoichs, propensities)
end

# ─── adaptive mixed-grid Schlögl system ──────────────────────────────────────

"""
    build_schlogl_mixed_system(model, fine_grid, levels)
        -> DiscreteStochasticSystem{CartesianIndex{K_mixed}}

Build the Schlögl RDME on a multi-level mixed-resolution grid.

For each pair j:
  - levels[j] = 0 → fine: two voxels, fine reactions + intra-pair diffusion
  - levels[j] = 1 → pair-coarse: one voxel nc_j, Binomial-corrected rates (r=2)
  - levels[j] = 2 → quad-coarse: levels[j] and levels[j+1] must both be 2;
                    one voxel ñ representing 4 fine voxels, Binomial-corrected
                    rates (r=4)

**Inter-boundary diffusion** uses first-moment closure:
  fine voxel:      propensity = d_fine × n
  level-1 voxel:   propensity = d_fine × nc / 2
  level-2 voxel:   propensity = d_fine × ñ  / 4

**Level-2 coarse propensities** (Binomial(ñ, 1/4) per fine voxel, r = 4):
  Production:     4 c3
  Death:          c4 ñ
  Autocatalysis:  c1 ñ(ñ-1) / 8
  Reverse:        c2 ñ(ñ-1)(ñ-2) / 96
"""
function build_schlogl_mixed_system(model::SchloglModel1D,
                                     fine_grid::VoxelGrid,
                                     levels::Vector{Int})
    K = fine_grid.n_voxels
    iseven(K) || error("K=$K must be even")
    n_pairs = K ÷ 2
    length(levels) == n_pairs || error("levels length $(length(levels)) ≠ K/2=$n_pairs")

    d_fine = diffusion_rate(model.D, fine_grid)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4

    KM            = _mixed_dim(levels)
    mixed_offsets = _mixed_offsets(levels)

    stoichs      = CartesianIndex{KM}[]
    propensities = Function[]

    # ── per-voxel reactions ───────────────────────────────────────────────────
    j = 1
    while j <= n_pairs
        mp = mixed_offsets[j]
        l  = levels[j]

        if l == 0
            # Fine pair: two independent fine voxels
            for p in (mp, mp+1)
                let p = p
                    ep = _e(p, Val(KM))
                    push!(stoichs,  ep); push!(propensities, (x,rv,t) -> c3)
                    push!(stoichs, -ep); push!(propensities, (x,rv,t) -> c4*max(0,Tuple(x)[p]))
                    push!(stoichs,  ep); push!(propensities, (x,rv,t) -> begin
                        n=Tuple(x)[p]; c1*max(0,n*(n-1))/2 end)
                    push!(stoichs, -ep); push!(propensities, (x,rv,t) -> begin
                        n=Tuple(x)[p]; c2*max(0,n*(n-1)*(n-2))/6 end)
                end
            end
            let p=mp, p1=mp+1
                ep=_e(p,Val(KM)); ep1=_e(p1,Val(KM))
                push!(stoichs,-ep+ep1); push!(propensities,(x,rv,t)->d_fine*max(0,Tuple(x)[p]))
                push!(stoichs, ep-ep1); push!(propensities,(x,rv,t)->d_fine*max(0,Tuple(x)[p1]))
            end
            j += 1

        elseif l == 1
            # Level-1 pair: one coarse voxel nc, Binomial-corrected (r=2)
            let p = mp
                ep = _e(p, Val(KM))
                push!(stoichs,  ep); push!(propensities, (x,rv,t) -> 2*c3)
                push!(stoichs, -ep); push!(propensities, (x,rv,t) -> c4*max(0,Tuple(x)[p]))
                push!(stoichs,  ep); push!(propensities, (x,rv,t) -> begin
                    nc=Tuple(x)[p]; c1*max(0,nc*(nc-1))/4 end)
                push!(stoichs, -ep); push!(propensities, (x,rv,t) -> begin
                    nc=Tuple(x)[p]; c2*max(0,nc*(nc-1)*(nc-2))/24 end)
            end
            j += 1

        else
            # Level-2 quad: one coarse voxel ñ representing 4 fine voxels (r=4)
            let p = mp
                ep = _e(p, Val(KM))
                push!(stoichs,  ep); push!(propensities, (x,rv,t) -> 4*c3)
                push!(stoichs, -ep); push!(propensities, (x,rv,t) -> c4*max(0,Tuple(x)[p]))
                push!(stoichs,  ep); push!(propensities, (x,rv,t) -> begin
                    nq=Tuple(x)[p]; c1*max(0,nq*(nq-1))/8 end)
                push!(stoichs, -ep); push!(propensities, (x,rv,t) -> begin
                    nq=Tuple(x)[p]; c2*max(0,nq*(nq-1)*(nq-2))/96 end)
            end
            j += 2  # skip the paired level-2 entry
        end
    end

    # ── inter-pair diffusion ──────────────────────────────────────────────────
    # Iterate over boundaries between adjacent groups.
    # For each boundary j | j+1, identify the rightmost coord of group j and
    # leftmost coord of group j+1, and the effective voxel divisor for closure.
    # Divisor: level-0 → 1, level-1 → 2, level-2 → 4.
    j = 1
    while j <= n_pairs
        l = levels[j]
        j_next = l == 2 ? j + 2 : j + 1   # index of next group
        j_next > n_pairs && break

        # Rightmost position and divisor for group j
        if l == 0
            p_r   = mixed_offsets[j] + 1   # second fine voxel
            div_r = 1
        elseif l == 1
            p_r   = mixed_offsets[j]
            div_r = 2
        else  # level 2: the shared quad coord
            p_r   = mixed_offsets[j]
            div_r = 4
        end

        # Leftmost position and divisor for group j_next
        l_next = levels[j_next]
        p_l    = mixed_offsets[j_next]
        div_l  = l_next == 0 ? 1 : l_next == 1 ? 2 : 4

        let p_r=p_r, p_l=p_l, div_r=div_r, div_l=div_l
            er = _e(p_r, Val(KM)); el = _e(p_l, Val(KM))
            push!(stoichs, -er+el)
            push!(propensities, (x,rv,t) -> d_fine*max(0,Tuple(x)[p_r])/div_r)
            push!(stoichs,  er-el)
            push!(propensities, (x,rv,t) -> d_fine*max(0,Tuple(x)[p_l])/div_l)
        end

        j = l == 2 ? j + 2 : j + 1
    end

    DiscreteStochasticSystem{CartesianIndex{KM}}(stoichs, propensities)
end

# ─── dynamic-π coarse Schlögl system ──────────────────────────────────────────

"""
    build_schlogl_coarse_system_dynamic(model, coarse_grid, fine_grid, pi_table)
        -> DiscreteStochasticSystem{CartesianIndex{K2}}

Coarse-grid Schlögl system using dynamic-π prolongation weights (from
`compute_dynamic_pi`) instead of the Binomial(nc, 1/2) approximation.

Coarse propensities are computed as expectations under the precomputed
conditional distribution π(n1 | nc):

  a_prod(nc)  = r · c3                            [zeroth-order, exact]
  a_death(nc) = c4 · nc                           [first-order, exact]
  a_auto(nc)  = Σ π(n1|nc) · [c1·n1(n1-1)/2 + c1·n2(n2-1)/2]
  a_rev(nc)   = Σ π(n1|nc) · [c2·n1(n1-1)(n1-2)/6 + c2·n2(n2-1)(n2-2)/6]

The expectations are **precomputed once** from pi_table (not re-evaluated at
every generator build), giving O(1) propensity evaluation.

Diffusion conserves nc for K2=1, so no intra-coarse-voxel diffusion term
appears.  Inter-coarse-voxel diffusion (K2 > 1) uses D/(2h)².
"""
function build_schlogl_coarse_system_dynamic(model::SchloglModel1D,
                                              coarse_grid::VoxelGrid,
                                              fine_grid::VoxelGrid,
                                              pi_table::Vector{Vector{Float64}})
    r  = coarse_grid.dx / fine_grid.dx
    K2 = coarse_grid.n_voxels
    d  = diffusion_rate(model.D, coarse_grid)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
    n_coarse_max = length(pi_table) - 1   # = 2 · n_max_fine

    # Precompute expectation tables indexed by nc = 0:n_coarse_max
    auto_rate = zeros(n_coarse_max + 1)
    rev_rate  = zeros(n_coarse_max + 1)
    for nc in 0:n_coarse_max
        pi_nc = pi_table[nc+1]
        for n1 in 0:nc
            n2 = nc - n1
            w = pi_nc[n1+1]
            auto_rate[nc+1] += w * c1 * (n1*(n1-1)/2  + n2*(n2-1)/2)
            rev_rate[nc+1]  += w * c2 * (n1*(n1-1)*(n1-2)/6 + n2*(n2-1)*(n2-2)/6)
        end
    end

    stoichs      = CartesianIndex{K2}[]
    propensities = Function[]

    for j in 1:K2
        let j = j, auto_rate = auto_rate, rev_rate = rev_rate
            ej = _e(j, Val(K2))

            # Production (r fine voxels per coarse voxel)
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> r * c3)

            # Death (first-order: exact regardless of π)
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> c4 * max(0, Tuple(x)[j]))

            # Autocatalysis: table lookup
            push!(stoichs, ej)
            push!(propensities, (x, rates, t) -> begin
                nc = Tuple(x)[j]
                0 <= nc <= n_coarse_max ? auto_rate[nc+1] : c1*nc*(nc-1)/4
            end)

            # Reverse trimolecular: table lookup
            push!(stoichs, -ej)
            push!(propensities, (x, rates, t) -> begin
                nc = Tuple(x)[j]
                0 <= nc <= n_coarse_max ? rev_rate[nc+1] : c2*nc*(nc-1)*(nc-2)/24
            end)
        end
    end

    # Inter-coarse-voxel diffusion (no-op for K2=1 since there are no pairs)
    for j in 1:K2-1
        let j = j
            ej  = _e(j,   Val(K2))
            ej1 = _e(j+1, Val(K2))
            push!(stoichs, -ej + ej1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[j]))
            push!(stoichs,  ej - ej1)
            push!(propensities, (x, rates, t) -> d * max(0, Tuple(x)[j+1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K2}}(stoichs, propensities)
end

# ─── Two-species bottleneck RDME ──────────────────────────────────────────────

"""
    BottleneckModel1D

1D two-species RDME:  ∅ →(k_prod) A →(k_AB) B →(k_deg) ∅  with diffusion.

State type: `CartesianIndex{2K}` with interleaved layout
`(nA_1, nB_1, nA_2, nB_2, …, nA_K, nB_K)` where `nA_k` / `nB_k` are the
A / B molecule counts in voxel `k`.

Analytical steady state (linear kinetics):
  μ_A = k_prod / k_AB,  μ_B = k_prod / k_deg
  P(state) = ∏_k Pois(μ_A)(nA_k) · Pois(μ_B)(nB_k)
"""
struct BottleneckModel1D
    D_A::Float64
    D_B::Float64
    k_prod::Float64
    k_AB::Float64
    k_deg::Float64
end

coarsen_model(m::BottleneckModel1D, r::Real) =
    BottleneckModel1D(m.D_A, m.D_B, m.k_prod * r, m.k_AB, m.k_deg)

"""Boundary condition for 2-species RDME: all counts in [0, n_max]."""
function bottleneck_bc(state::CartesianIndex{N}, n_max::Int) where {N}
    all(c -> 0 ≤ c ≤ n_max, Tuple(state))
end

"""
    build_bottleneck_system(model, grid) -> DiscreteStochasticSystem{CartesianIndex{2K}}

Full RDME for the bottleneck model. State layout: interleaved `(nA_1, nB_1, nA_2, nB_2, …)`.
"""
function build_bottleneck_system(model::BottleneckModel1D, grid::VoxelGrid)
    K  = grid.n_voxels
    N  = 2 * K
    dA = diffusion_rate(model.D_A, grid)
    dB = diffusion_rate(model.D_B, grid)
    k_prod, k_AB, k_deg = model.k_prod, model.k_AB, model.k_deg

    stoichs      = CartesianIndex{N}[]
    propensities = Function[]

    for k in 1:K
        let k = k
            iA = 2k - 1
            iB = 2k
            eA = _e(iA, Val(N))
            eB = _e(iB, Val(N))

            push!(stoichs, eA);        push!(propensities, (x, rv, t) -> k_prod)
            push!(stoichs, -eA + eB); push!(propensities, (x, rv, t) -> k_AB * max(0, Tuple(x)[iA]))
            push!(stoichs, -eB);      push!(propensities, (x, rv, t) -> k_deg * max(0, Tuple(x)[iB]))
        end
    end

    for k in 1:(K-1)
        let k = k
            iAk  = 2k - 1;  iAk1 = 2k + 1
            iBk  = 2k;      iBk1 = 2k + 2
            eAk  = _e(iAk,  Val(N));  eAk1 = _e(iAk1, Val(N))
            eBk  = _e(iBk,  Val(N));  eBk1 = _e(iBk1, Val(N))

            push!(stoichs, -eAk + eAk1); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAk]))
            push!(stoichs,  eAk - eAk1); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAk1]))
            push!(stoichs, -eBk + eBk1); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBk]))
            push!(stoichs,  eBk - eBk1); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBk1]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, propensities)
end

"""
    build_intra_system_2s(model, fine_grid, coarse_grid) -> DiscreteStochasticSystem{CartesianIndex{2K}}

Within-pair diffusion smoother for 2-species RDME.
Intra-pair A and B diffusion only; excludes reactions and inter-pair diffusion.
"""
function build_intra_system_2s(model::BottleneckModel1D,
                                fine_grid::VoxelGrid,
                                coarse_grid::VoxelGrid)
    K  = fine_grid.n_voxels
    K2 = coarse_grid.n_voxels
    @assert K == 2 * K2
    N  = 2 * K
    dA = diffusion_rate(model.D_A, fine_grid)
    dB = diffusion_rate(model.D_B, fine_grid)

    stoichs      = CartesianIndex{N}[]
    propensities = Function[]

    for j in 1:K2
        let j = j
            kL = 2j - 1;  kR = 2j
            iAL = 2kL - 1;  iAR = 2kR - 1
            iBL = 2kL;      iBR = 2kR
            eAL = _e(iAL, Val(N));  eAR = _e(iAR, Val(N))
            eBL = _e(iBL, Val(N));  eBR = _e(iBR, Val(N))

            push!(stoichs, -eAL + eAR); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAL]))
            push!(stoichs,  eAL - eAR); push!(propensities, (x, rv, t) -> dA * max(0, Tuple(x)[iAR]))
            push!(stoichs, -eBL + eBR); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBL]))
            push!(stoichs,  eBL - eBR); push!(propensities, (x, rv, t) -> dB * max(0, Tuple(x)[iBR]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, propensities)
end

# ─── 2D two-species bottleneck RDME ──────────────────────────────────────────

"""
    build_bottleneck_system_2d(model, grid; source_voxels=nothing)
        -> DiscreteStochasticSystem{CartesianIndex{2*K_x*K_y}}

2D RDME for the bottleneck model on a Grid2D.
State: CartesianIndex{2K} interleaved (nA_1,nB_1,...,nA_K,nB_K) with K=K_x*K_y,
       linear index k=(i-1)*K_y+j for voxel (i,j).

source_voxels: Set or Vector of linear voxel indices where ∅→A production occurs.
               If nothing, production occurs in all voxels.
"""
function build_bottleneck_system_2d(model::BottleneckModel1D, grid::Grid2D;
                                     source_voxels=nothing)
    K_x, K_y = grid.K_x, grid.K_y
    K  = K_x * K_y
    N  = 2 * K
    dA = diffusion_rate(model.D_A, grid)
    dB = diffusion_rate(model.D_B, grid)
    k_prod, k_AB, k_deg = model.k_prod, model.k_AB, model.k_deg

    src = source_voxels === nothing ? (1:K) : source_voxels
    src_set = Set(src)

    stoichs      = CartesianIndex{N}[]
    propensities = Function[]

    # ── per-voxel reactions ──────────────────────────────────────────────────
    for k in 1:K
        let k = k
            iA = 2k - 1; iB = 2k
            eA = _e(iA, Val(N)); eB = _e(iB, Val(N))

            if k in src_set
                push!(stoichs, eA)
                push!(propensities, (x, rv, t) -> k_prod)
            end
            push!(stoichs, -eA + eB)
            push!(propensities, (x, rv, t) -> k_AB * max(0, Tuple(x)[iA]))
            push!(stoichs, -eB)
            push!(propensities, (x, rv, t) -> k_deg * max(0, Tuple(x)[iB]))
        end
    end

    # ── 4-connected diffusion ────────────────────────────────────────────────
    for i in 1:K_x, j in 1:K_y
        let i=i, j=j
            k = (i-1)*K_y + j
            iAk = 2k-1; iBk = 2k
            eAk = _e(iAk, Val(N)); eBk = _e(iBk, Val(N))

            # Right neighbor (j+1)
            if j < K_y
                k2 = k + 1
                iAk2=2k2-1; iBk2=2k2
                eAk2=_e(iAk2, Val(N)); eBk2=_e(iBk2, Val(N))
                push!(stoichs, -eAk+eAk2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iAk]))
                push!(stoichs,  eAk-eAk2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iAk2]))
                push!(stoichs, -eBk+eBk2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iBk]))
                push!(stoichs,  eBk-eBk2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iBk2]))
            end
            # Down neighbor (i+1)
            if i < K_x
                k2 = k + K_y
                iAk2=2k2-1; iBk2=2k2
                eAk2=_e(iAk2, Val(N)); eBk2=_e(iBk2, Val(N))
                push!(stoichs, -eAk+eAk2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iAk]))
                push!(stoichs,  eAk-eAk2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iAk2]))
                push!(stoichs, -eBk+eBk2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iBk]))
                push!(stoichs,  eBk-eBk2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iBk2]))
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, propensities)
end

"""
    build_intra_system_2d(model, fine_grid, cmap)
        -> DiscreteStochasticSystem{CartesianIndex{2*K}}

Within-patch diffusion smoother for 2D 2-species RDME.
Only includes diffusion between voxels within the same coarse cell.
"""
function build_intra_system_2d(model::BottleneckModel1D,
                                fine_grid::Grid2D,
                                cmap::CoarseningMap)
    K  = n_voxels(fine_grid)
    N  = 2K
    dA = diffusion_rate(model.D_A, fine_grid)
    dB = diffusion_rate(model.D_B, fine_grid)
    K_y = fine_grid.K_y

    # Build adjacency set: (k1,k2) are neighbors if |k1-k2|==1 or |k1-k2|==K_y AND same row check
    stoichs      = CartesianIndex{N}[]
    propensities = Function[]

    for J in 1:cmap.n_coarse
        patch = cmap.coarse_to_fine[J]
        # Emit diffusion for each adjacent pair within the patch
        for a in 1:length(patch), b in (a+1):length(patch)
            k1, k2 = patch[a], patch[b]
            # Check if k1 and k2 are spatially adjacent (horizontal or vertical)
            i1, j1 = (k1-1) ÷ K_y + 1, (k1-1) % K_y + 1
            i2, j2 = (k2-1) ÷ K_y + 1, (k2-1) % K_y + 1
            adjacent = (i1==i2 && abs(j1-j2)==1) || (j1==j2 && abs(i1-i2)==1)
            adjacent || continue
            let k1=k1, k2=k2
                iA1=2k1-1; iB1=2k1; iA2=2k2-1; iB2=2k2
                eA1=_e(iA1,Val(N)); eA2=_e(iA2,Val(N))
                eB1=_e(iB1,Val(N)); eB2=_e(iB2,Val(N))
                push!(stoichs, -eA1+eA2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iA1]))
                push!(stoichs,  eA1-eA2); push!(propensities, (x,rv,t)->dA*max(0,Tuple(x)[iA2]))
                push!(stoichs, -eB1+eB2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iB1]))
                push!(stoichs,  eB1-eB2); push!(propensities, (x,rv,t)->dB*max(0,Tuple(x)[iB2]))
            end
        end
    end

    DiscreteStochasticSystem{CartesianIndex{N}}(stoichs, propensities)
end
