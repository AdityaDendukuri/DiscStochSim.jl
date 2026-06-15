# ─── Block-decomposition FSP for the weakly-correlated regime ────────────────
#
# Some spatial systems (e.g. 2-D birth–death with independent filling) have a
# joint support that grows ≈ ∏(per-voxel supports): the full joint blows up
# geometrically even though every marginal is small.  For these the global-joint
# UnifiedFSP is intractable.  This module solves them by DECOMPOSITION.
#
# The voxel set is partitioned into blocks of uniform size `Kb`.  Each block is a
# small *exact* joint CME (solved by the same SA-AMG backward-Euler step as
# UnifiedFSP).  Diffusion edges *within* a block are kept exactly.  Diffusion
# edges *across* a block boundary are split into:
#
#     outflux  a → ext :  D · n_a            (depends on the block's own state — exact)
#     influx   ext → a :  D · ⟨n_g⟩          (neighbour's true marginal mean — closure)
#
# Because diffusion is LINEAR in the count, coupling through the marginal mean
# D·⟨n_g⟩ is the EXACT first-moment consequence of treating blocks as
# independent — it is NOT a Poisson assumption and NOT a mean-field ODE.  The
# marginal means ⟨n_g⟩ = Σ n_g · p are read straight from the neighbour blocks'
# realised distributions (see [[feedback_no_meanfield]]).  The only approximation
# is neglecting cross-block correlation; for weakly-correlated systems that is
# negligible, and it vanishes as block size grows.
#
# The frozen-mean operator splitting introduces a time-integration error that is
# O(dt); it is controlled by `nsub` boundary-coupling sub-cycles per step
# (recompute ⟨n_g⟩ between sub-steps), which converges the block means to the
# full-joint means as nsub → ∞ / dt → 0.

"""
    block_subsystem(model::GraphRDMEModel, ::Val{Kb}, gvox, internal_edges, boundary_edges, μ_ext)

Build a single block's local `DiscreteStochasticSystem` over `CartesianIndex{Kb}`
configurations (one component per voxel in the block).

- `gvox[a]` : global voxel id of local voxel `a` (selects per-voxel rates)
- per-voxel birth `birth_rates[g]` and death `death_rates[g] · n`
- `internal_edges :: (a, b, D)` — bidirectional diffusion between local voxels `a,b`
- `boundary_edges :: (a, g, D)` — local voxel `a` ↔ external global voxel `g`:
  outflux `D·n_a` (exact) and influx `D·μ_ext[g]` (mean closure)
"""
function block_subsystem(model::GraphRDMEModel, ::Val{Kb}, gvox,
                         internal_edges, boundary_edges, μ_ext) where Kb
    stoich = CartesianIndex{Kb}[]; props = Function[]
    ek(k) = CartesianIndex(ntuple(j -> j == k ? 1 : 0, Kb))
    for a in 1:Kb
        let ea = ek(a), kb = model.birth_rates[gvox[a]], kd = model.death_rates[gvox[a]], a = a
            push!(stoich,  ea); push!(props, (x, r, t) -> kb)
            push!(stoich, -ea); push!(props, (x, r, t) -> kd * max(0, Tuple(x)[a]))
        end
    end
    for (a, b, D) in internal_edges
        let ea = ek(a), eb = ek(b), D = D
            push!(stoich, -ea + eb); push!(props, (x, r, t) -> D * max(0, Tuple(x)[a]))
            push!(stoich,  ea - eb); push!(props, (x, r, t) -> D * max(0, Tuple(x)[b]))
        end
    end
    for (a, g, D) in boundary_edges
        let ea = ek(a), rin = D * μ_ext[g], D = D, a = a
            push!(stoich, -ea); push!(props, (x, r, t) -> D * max(0, Tuple(x)[a]))  # outflux (exact)
            push!(stoich,  ea); push!(props, (x, r, t) -> rin)                       # influx (closure)
        end
    end
    DiscreteStochasticSystem{CartesianIndex{Kb}}(stoich, props)
end

"""
    partition_mesh_edges(mesh, voxel_block, nblocks)
        -> (block_voxels, local_idx, internal, boundary)

Split a `VoxelMesh`'s diffusion edges by block membership.

- `block_voxels[b]` : global voxel ids in block `b`
- `local_idx[v]`    : local component index of global voxel `v` within its block
- `internal[b]`     : `(local_a, local_b, D)` edges with both endpoints in block `b`
- `boundary[b]`     : `(local_a, global_g, D)` edges leaving block `b`
"""
function partition_mesh_edges(mesh, voxel_block::AbstractVector{Int}, nblocks::Int)
    block_voxels = [Int[] for _ in 1:nblocks]
    for v in 1:mesh.n_voxels; push!(block_voxels[voxel_block[v]], v); end
    local_idx = zeros(Int, mesh.n_voxels)
    for b in 1:nblocks, (li, v) in enumerate(block_voxels[b]); local_idx[v] = li; end
    internal = [Tuple{Int,Int,Float64}[] for _ in 1:nblocks]
    boundary = [Tuple{Int,Int,Float64}[] for _ in 1:nblocks]
    for (e, (u, v)) in enumerate(mesh.edges)
        D = mesh.D_rates[e]; bu = voxel_block[u]; bv = voxel_block[v]
        if bu == bv
            push!(internal[bu], (local_idx[u], local_idx[v], D))
        else
            push!(boundary[bu], (local_idx[u], v, D))
            push!(boundary[bv], (local_idx[v], u, D))
        end
    end
    block_voxels, local_idx, internal, boundary
end

"""
    BlockFSP{Kb}

Block-decomposition FSP state over blocks of uniform size `Kb`.  Each block is a
`StateSpace{CartesianIndex{Kb},Float64}` evolved as an exact local joint; blocks
are coupled through marginal means at their shared boundaries.
"""
mutable struct BlockFSP{Kb}
    model        :: GraphRDMEModel
    blocks       :: Vector{StateSpace{CartesianIndex{Kb}, Float64}}
    block_voxels :: Vector{Vector{Int}}
    local_idx    :: Vector{Int}
    internal     :: Vector{Vector{Tuple{Int,Int,Float64}}}
    boundary     :: Vector{Vector{Tuple{Int,Int,Float64}}}
    n_voxels     :: Int
    nmax         :: Int
    dt           :: Float64
    nsub         :: Int
    expand_depth :: Int
    prune_tol    :: Float64
    n_cycles     :: Int
    rtol         :: Float64
    t            :: Float64
end

"""
    BlockFSP(model, mesh, voxel_block, ::Val{Kb}; dt, nmax, nsub=1,
             expand_depth=1, prune_tol=1e-10, n_cycles=30, rtol=1e-10) -> BlockFSP{Kb}

Construct a block solver.  `voxel_block[v]` assigns each global voxel to a block;
every block must contain exactly `Kb` voxels.  Each block is initialised to the
all-zero configuration with probability 1.  `nsub` boundary-coupling sub-cycles
per `dt` control the operator-splitting error.
"""
function BlockFSP(model::GraphRDMEModel, mesh, voxel_block::AbstractVector{Int},
                  ::Val{Kb}; dt::Float64, nmax::Int, nsub::Int=1,
                  expand_depth::Int=1, prune_tol::Float64=1e-10,
                  n_cycles::Int=30, rtol::Float64=1e-10) where Kb
    nblocks = maximum(voxel_block)
    bvs, lidx, intl, bnd = partition_mesh_edges(mesh, voxel_block, nblocks)
    for b in 1:nblocks
        length(bvs[b]) == Kb || error("block $b has $(length(bvs[b])) voxels, expected Kb=$Kb")
    end
    blocks = Vector{StateSpace{CartesianIndex{Kb}, Float64}}(undef, nblocks)
    for b in 1:nblocks
        sp = StateSpace{CartesianIndex{Kb}, Float64}()
        add_state!(sp, CartesianIndex(ntuple(_ -> 0, Val(Kb))), 1.0)
        blocks[b] = sp
    end
    BlockFSP{Kb}(model, blocks, bvs, lidx, intl, bnd, mesh.n_voxels,
                 nmax, dt, nsub, expand_depth, prune_tol, n_cycles, rtol, 0.0)
end

"""
    block_means(f::BlockFSP) -> Vector{Float64}

Per-voxel marginal mean ⟨n_v⟩ = Σ n_v · p, read from the realised block
distributions (exact, no mean-field).
"""
function block_means(f::BlockFSP{Kb}) where Kb
    μ = zeros(f.n_voxels)
    for b in 1:length(f.blocks)
        bv = f.block_voxels[b]
        for (ci, p) in zip(f.blocks[b].states, f.blocks[b].probs)
            t = Tuple(ci)
            @inbounds for (li, gv) in enumerate(bv); μ[gv] += t[li] * p; end
        end
    end
    μ
end

n_active(f::BlockFSP) = sum(length, f.blocks)
max_block_size(f::BlockFSP) = maximum(length, f.blocks)

"""
    step!(f::BlockFSP) -> f

Advance every block by `dt` using `nsub` frozen-mean sub-cycles: recompute the
boundary marginal means, rebuild each block's local generator with the frozen
closure, and take a backward-Euler SA-AMG sub-step.
"""
function step!(f::BlockFSP{Kb}) where Kb
    bc = x -> all(c -> 0 <= c <= f.nmax, Tuple(x))
    h  = f.dt / f.nsub
    for _ in 1:f.nsub
        μ = block_means(f)
        for b in 1:length(f.blocks)
            sys = block_subsystem(f.model, Val(Kb), f.block_voxels[b], f.internal[b], f.boundary[b], μ)
            sp  = f.blocks[b]
            expand!(sp, sys, bc; depth=f.expand_depth)
            A, = build_generator(sp, sys, Float64[], f.t)
            sp.probs .= amg_be_step(A, sp.probs, h; n_cycles=f.n_cycles, rtol=f.rtol)
            sp.probs .= max.(0.0, sp.probs)
            renormalize!(sp); prune_threshold!(sp, f.prune_tol)
        end
    end
    f.t += f.dt
    f
end
