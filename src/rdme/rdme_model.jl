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
struct SchloglModel1D
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
