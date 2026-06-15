"""
Graded-grid variable-width RDME operator
========================================

Shape-agnostic spatial coarsening for the RDME.  A settled block of fine voxels
is *collapsed into one super-voxel of larger physical width*.  Nothing about the
geometry of the active region is assumed — the only thing that changes when a
block is merged is **diffusion**, because the discrete diffusion rate is set by
cell width through a finite-volume face flux.

Finite-volume face flux between two cells `a`, `b`:

    g_ab = D · ℓ_ab / L_ab                       (conductance)
    k_{a→b} = g_ab / V_a ,   k_{b→a} = g_ab / V_b (per-molecule hop rates)

with
    ℓ_ab = length of the shared face (physical),
    L_ab = distance between the two cell centres (physical),
    V_a  = area (volume) of cell a.

For two equal fine cells (`ℓ = L = h`, `V = h²`) this collapses to the usual
`D / h²` jump rate.  For a wide super-voxel `V_a ∝ W²`, the per-molecule rate out
of it is `∝ 1/W` — **wider cell ⇒ slower diffusion**, exactly the merge rule.
The flux at uniform concentration is `g_ab·c` in both directions ⇒ no net flux,
so the operator is mass-conserving and leaves a flat field invariant.

Reactions are rescaled with the cell volume `Ω` (number of fine sub-voxels) by
the standard system-size scaling, so a settled super-voxel sits at the same
basin *concentration* as the fine voxels it replaced:
    ∅→X  : Ω·c3        X→∅  : c4·n
    2X→3X: c1·n(n-1)/2 / Ω      3X→2X: c2·n(n-1)(n-2)/6 / Ω²

Public API
----------
  GradedCell, GradedGrid
  build_graded_grid(Kx, Ky, active; h, max_width, grade) → GradedGrid
  active_from_means(means, n_lo, n_hi; tol_lo, tol_hi)   → BitMatrix
  graded_diffusion_rates(grid, D)        → (edges, k_fwd, k_rev)
  build_graded_schlogl_system(grid, model) → DiscreteStochasticSystem
  to_voxel_graph(grid)                   → VoxelGraph (super-voxel connectivity)
"""

# ── cell / grid types ──────────────────────────────────────────────────────────

"""
    GradedCell

One super-voxel: a `w × w` square block of fine voxels (`w = 1` is a fine voxel).
`i0:i1, j0:j1` is the (inclusive) fine-index bounding box; `cx,cy` the centre in
fine-index units; `area` = number of fine voxels = w².
"""
struct GradedCell
    i0::Int; i1::Int; j0::Int; j1::Int
    w::Int
    cx::Float64; cy::Float64
    area::Float64
end

"""
    GradedGrid

A partition of a `Kx × Ky` fine grid into square super-voxels of varying width,
with precomputed adjacency, shared-face lengths and centre distances (all in
fine-index units; multiply by `h` for physical lengths).
"""
struct GradedGrid
    Kx::Int; Ky::Int
    h::Float64
    cells::Vector{GradedCell}
    cell_of::Matrix{Int}            # Kx×Ky → cell index
    edges::Vector{Tuple{Int,Int}}   # adjacent cell pairs (a<b)
    face_len::Vector{Float64}       # shared-face length (fine units) per edge
    center_dist::Vector{Float64}    # centre-to-centre distance (fine units) per edge
end

n_cells(g::GradedGrid) = length(g.cells)

# ── activity (settledness) mask ─────────────────────────────────────────────────

"""
    active_from_means(means, n_lo, n_hi; tol_lo, tol_hi) → BitMatrix

A voxel is *active* (unsettled) when either
  (i)  its mean count sits between the two basins (within neither `tol_lo` of
       `n_lo` nor `tol_hi` of `n_hi`), or
  (ii) it is settled in one basin but borders a voxel settled in the *other*
       basin — i.e. it lies on a basin interface, where diffusion does work.

Criterion (ii) is what keeps a sharp front fine: with a hard lo/hi field every
voxel is at a basin, so (i) alone flags nothing.  Both tests are settledness
tests on the basin label — they make no assumption about where the front is or
what shape it has; the active set (the front) emerges from the field.
"""
function active_from_means(means::AbstractMatrix{<:Real}, n_lo::Real, n_hi::Real;
                           tol_lo::Real = 0.15 * (n_hi - n_lo),
                           tol_hi::Real = 0.15 * (n_hi - n_lo))
    Kx, Ky = size(means)
    # basin label: -1 = lo, +1 = hi, 0 = between (intrinsically active)
    label = zeros(Int, Kx, Ky)
    for i in 1:Kx, j in 1:Ky
        μ = means[i, j]
        label[i, j] = μ ≤ n_lo + tol_lo ?  -1 :
                      μ ≥ n_hi - tol_hi ?  +1 : 0
    end
    act = falses(Kx, Ky)
    for i in 1:Kx, j in 1:Ky
        if label[i, j] == 0
            act[i, j] = true                      # (i) between basins
            continue
        end
        for (ni, nj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))   # (ii) basin interface
            1 ≤ ni ≤ Kx && 1 ≤ nj ≤ Ky || continue
            if label[ni, nj] != 0 && label[ni, nj] != label[i, j]
                act[i, j] = true; break
            end
        end
    end
    act
end

# ── graded partition (shape-agnostic, quadtree-style block sizes) ───────────────

# largest power of two ≤ x (and ≤ cap)
_pow2_floor(x::Int, cap::Int) = (s = 1; while 2s ≤ x && 2s ≤ cap; s *= 2; end; s)

"""
    build_graded_grid(Kx, Ky, active; h=1.0, max_width=8, grade=1) → GradedGrid

Build a graded square-block partition of the `Kx × Ky` fine grid.

The target width of a voxel grows with its (Chebyshev) distance to the nearest
*active* voxel: `target = 1 + dist ÷ grade`, snapped down to a power of two and
capped at `max_width`.  Active voxels stay at width 1.  Blocks are then placed
greedily; a block never coarsens a voxel below its own target, so the resolution
degrades smoothly away from the active set.  Purely distance-driven ⇒ shape-blind.
"""
function build_graded_grid(Kx::Int, Ky::Int, active::AbstractMatrix{Bool};
                           h::Float64 = 1.0, max_width::Int = 8, grade::Int = 1)
    size(active) == (Kx, Ky) || error("active must be $Kx×$Ky")

    # ── Chebyshev distance to nearest active voxel (two-pass) ──────────────────
    INF = Kx + Ky + 1
    dist = fill(INF, Kx, Ky)
    for i in 1:Kx, j in 1:Ky
        active[i, j] && (dist[i, j] = 0)
    end
    for i in 1:Kx, j in 1:Ky
        d = dist[i, j]
        i > 1 && (d = min(d, dist[i-1, j] + 1))
        j > 1 && (d = min(d, dist[i, j-1] + 1))
        dist[i, j] = d
    end
    for i in Kx:-1:1, j in Ky:-1:1
        d = dist[i, j]
        i < Kx && (d = min(d, dist[i+1, j] + 1))
        j < Ky && (d = min(d, dist[i, j+1] + 1))
        dist[i, j] = d
    end

    # ── per-voxel target block width (power of two) ────────────────────────────
    target = Matrix{Int}(undef, Kx, Ky)
    for i in 1:Kx, j in 1:Ky
        t = 1 + dist[i, j] ÷ max(grade, 1)
        target[i, j] = _pow2_floor(t, max_width)
    end

    # ── greedy square-block placement ──────────────────────────────────────────
    cell_of = zeros(Int, Kx, Ky)
    cells = GradedCell[]
    for i0 in 1:Kx, j0 in 1:Ky
        cell_of[i0, j0] == 0 || continue
        # largest square that fits the bounds, is unassigned, and does not
        # over-coarsen any member (s ≤ min target in the block)
        s = min(target[i0, j0], Kx - i0 + 1, Ky - j0 + 1)
        while s > 1
            ok = true
            for i in i0:i0+s-1, j in j0:j0+s-1
                if cell_of[i, j] != 0 || target[i, j] < s
                    ok = false; break
                end
            end
            ok && break
            s ÷= 2
        end
        idx = length(cells) + 1
        for i in i0:i0+s-1, j in j0:j0+s-1
            cell_of[i, j] = idx
        end
        push!(cells, GradedCell(i0, i0+s-1, j0, j0+s-1, s,
                                (i0 + i0+s-1) / 2, (j0 + j0+s-1) / 2, float(s*s)))
    end

    # ── adjacency, shared-face length, centre distance ─────────────────────────
    edge_set = Dict{Tuple{Int,Int}, Float64}()   # (a<b) → shared face length
    for i in 1:Kx, j in 1:Ky
        a = cell_of[i, j]
        # only look right and down across fine-voxel faces; accumulate face length
        if i < Kx
            b = cell_of[i+1, j]
            a != b && (key = minmax(a, b); edge_set[key] = get(edge_set, key, 0.0) + 1.0)
        end
        if j < Ky
            b = cell_of[i, j+1]
            a != b && (key = minmax(a, b); edge_set[key] = get(edge_set, key, 0.0) + 1.0)
        end
    end

    edges = Tuple{Int,Int}[]
    face_len = Float64[]
    center_dist = Float64[]
    for ((a, b), ℓ) in edge_set
        ca = cells[a]; cb = cells[b]
        L = hypot(ca.cx - cb.cx, ca.cy - cb.cy)
        push!(edges, (a, b)); push!(face_len, ℓ); push!(center_dist, L)
    end

    GradedGrid(Kx, Ky, h, cells, cell_of, edges, face_len, center_dist)
end

# ── diffusion: face-flux directional per-molecule rates ─────────────────────────

"""
    graded_diffusion_rates(grid, D) → (edges, k_fwd, k_rev)

For each cell edge `(a,b)` return the directional per-molecule hop rates
`k_fwd = k_{a→b} = g_ab / V_a` and `k_rev = k_{b→a} = g_ab / V_b`,
where `g_ab = D · ℓ_ab / L_ab` is the finite-volume conductance and the lengths
are made physical with the fine width `grid.h`.
"""
function graded_diffusion_rates(grid::GradedGrid, D::Real)
    h = grid.h
    ne = length(grid.edges)
    k_fwd = Vector{Float64}(undef, ne)
    k_rev = Vector{Float64}(undef, ne)
    for e in 1:ne
        a, b = grid.edges[e]
        ℓ = grid.face_len[e] * h
        L = grid.center_dist[e] * h
        g = D * ℓ / L
        Va = grid.cells[a].area * h^2
        Vb = grid.cells[b].area * h^2
        k_fwd[e] = g / Va
        k_rev[e] = g / Vb
    end
    grid.edges, k_fwd, k_rev
end

# ── super-voxel connectivity as a VoxelGraph (for visualisation/reuse) ──────────

"""
    to_voxel_graph(grid) → VoxelGraph

Super-voxel connectivity with cell centres as positions.
"""
function to_voxel_graph(grid::GradedGrid)
    pos = Matrix{Float64}(undef, n_cells(grid), 2)
    for (c, cell) in enumerate(grid.cells)
        pos[c, :] = [cell.cx, cell.cy]
    end
    VoxelGraph(n_cells(grid), copy(grid.edges), pos)
end

# ── generator: volume-scaled Schlögl reactions + graded diffusion ───────────────

"""
    build_graded_schlogl_system(grid, model) → DiscreteStochasticSystem

Schlögl RDME generator on the graded super-voxel grid.  Reactions are scaled by
each cell's volume `Ω` (number of fine sub-voxels); diffusion uses the
finite-volume face-flux directional rates from [`graded_diffusion_rates`](@ref).
State index component `c` is the *total* molecule count in cell `c`.
"""
function build_graded_schlogl_system(grid::GradedGrid, model::SchloglModel1D)
    Nc = n_cells(grid)
    c1 = model.c1; c2 = model.c2; c3 = model.c3; c4 = model.c4
    eu(c) = CartesianIndex(ntuple(i -> i == c ? 1 : 0, Nc))

    stoichs      = CartesianIndex{Nc}[]
    propensities = Function[]

    # volume-scaled reactions per cell
    for c in 1:Nc
        Ω = grid.cells[c].area
        let c=c, ec=eu(c), Ω=Ω, invΩ=1.0/Ω, invΩ2=1.0/Ω^2
            push!(stoichs,  ec); push!(propensities, (x,rv,t) -> Ω * c3)
            push!(stoichs, -ec); push!(propensities, (x,rv,t) -> c4 * max(0, Tuple(x)[c]))
            push!(stoichs,  ec); push!(propensities,
                (x,rv,t) -> begin n=Tuple(x)[c]; c1*max(0,n)*max(0,n-1)/2 * invΩ; end)
            push!(stoichs, -ec); push!(propensities,
                (x,rv,t) -> begin n=Tuple(x)[c]; c2*max(0,n)*max(0,n-1)*max(0,n-2)/6 * invΩ2; end)
        end
    end

    # graded diffusion (directional, face-flux)
    edges, k_fwd, k_rev = graded_diffusion_rates(grid, model.D)
    for e in eachindex(edges)
        a, b = edges[e]
        let a=a, b=b, ea=eu(a), eb=eu(b), kab=k_fwd[e], kba=k_rev[e]
            push!(stoichs, -ea+eb); push!(propensities, (x,rv,t) -> kab * max(0, Tuple(x)[a]))
            push!(stoichs,  ea-eb); push!(propensities, (x,rv,t) -> kba * max(0, Tuple(x)[b]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{Nc}}(stoichs, propensities)
end
