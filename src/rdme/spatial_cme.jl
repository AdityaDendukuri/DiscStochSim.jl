"""
Spatial CME with multigrid coarsening/prolongation on the voxel axis
====================================================================

The full joint CME for a K-voxel RDME (single species) lives on the state space
    StateSpace{CartesianIndex{K}, Float64}
where component k of the CartesianIndex holds n_k, the molecule count in voxel k.

Multigrid coarsening merges adjacent voxel pairs: K → K' < K.

Restriction of the state space: CartesianIndex{K} → CartesianIndex{K'} by
summing molecule counts of merged voxels (exact: no information loss in probability).

Restriction of the model: recompute birth rates for merged super-voxels
(birth rate of super-voxel = sum of constituent fine-voxel birth rates);
death rate is per-molecule so it is unchanged.

Prolongation: coarse CartesianIndex{K'} → fine CartesianIndex{K} by splitting
merged molecule counts between fine voxels.

All existing CME machinery (build_generator, expv, expand!) runs unchanged —
coarsening just reduces the dimension K.

This matches the approach used by build_coarse_system / build_intra_system
in rdme_model.jl, generalised to arbitrary graph topologies.

Key types
---------
  VoxelMesh         — spatial graph (voxels + edges + per-edge diffusion rates)
  GraphRDMEModel    — per-voxel birth rates + shared death rate
  CoarseningPlan    — greedy pairing of voxels, mapping fine→coarse indices
  SpatialCMEState   — full joint distribution + model + mesh at one level

Key functions
-------------
  build_graph_rdme_system   — construct DiscreteStochasticSystem{CI{K}}
  coarsen_model             — scale birth rates for merged super-voxels
  coarsen_mesh              — merge edges for the coarse graph
  make_coarsening_plan      — flux-based greedy voxel pairing
  restrict                  — map StateSpace{CI{K}} → StateSpace{CI{K'}}
  prolong                   — map StateSpace{CI{K'}} → StateSpace{CI{K}}
  spatial_vcycle!           — one V-cycle step on a SpatialCMEState
  mean_counts               — exact FSP mean per voxel (no mean-field)
"""

# ── VoxelMesh ──────────────────────────────────────────────────────────────────

"""
    VoxelMesh

Spatial graph for a graph-based RDME.  Voxels are nodes; undirected edges
carry per-edge diffusion rates D/dx².

Fields
------
- `n_voxels`   : total number of voxels K
- `edges`      : undirected edges as (u, v) with u < v
- `D_rates`    : diffusion rate D/dx² for edge e (scalar per edge)
- `positions`  : optional K×ndim matrix of voxel centre coordinates
"""
struct VoxelMesh
    n_voxels  :: Int
    edges     :: Vector{Tuple{Int,Int}}
    D_rates   :: Vector{Float64}
    positions :: Union{Nothing, Matrix{Float64}}
end

VoxelMesh(n, edges, D_rates) = VoxelMesh(n, edges, D_rates, nothing)

"""Build a VoxelMesh from a VoxelGraph (graph_rdme.jl) with uniform D/dx²."""
function VoxelMesh(g::VoxelGraph, D::Float64, dx::Float64)
    rate = D / dx^2
    VoxelMesh(g.n_voxels, g.edges, fill(rate, length(g.edges)), g.positions)
end

"""Build a Kx×Ky rectangular 4-connected VoxelMesh."""
function rect_voxel_mesh(Kx::Int, Ky::Int, D::Float64; dx::Float64=1.0)
    rate = D / dx^2
    edges = Tuple{Int,Int}[]
    idx(i,j) = (i-1)*Ky + j
    for i in 1:Kx, j in 1:Ky
        j < Ky && push!(edges, (idx(i,j), idx(i,j+1)))
        i < Kx && push!(edges, (idx(i,j), idx(i+1,j)))
    end
    pos = Matrix{Float64}(undef, Kx*Ky, 2)
    for i in 1:Kx, j in 1:Ky
        pos[idx(i,j), :] = [Float64(i), Float64(j)]
    end
    VoxelMesh(Kx*Ky, edges, fill(rate, length(edges)), pos)
end

# ── GraphRDMEModel ─────────────────────────────────────────────────────────────

"""
    GraphRDMEModel

Model parameters for a single-species birth-death RDME on a graph.
Diffusion rates live in the `VoxelMesh`; this struct holds reaction rates only.

Fields
------
- `birth_rates` : zeroth-order birth rate per voxel (length K vector)
- `death_rate`  : first-order death rate per molecule (scalar)
"""
struct GraphRDMEModel
    birth_rates :: Vector{Float64}   # per voxel, length K
    death_rates :: Vector{Float64}   # per voxel, length K (typically all k_d, but boundary voxels may differ)
end

"""Uniform birth and death rates across all K voxels."""
GraphRDMEModel(K::Int, k_b::Float64, k_d::Float64) =
    GraphRDMEModel(fill(k_b, K), fill(k_d, K))

"""Uniform death rate from a scalar."""
GraphRDMEModel(birth_rates::Vector{Float64}, k_d::Float64) =
    GraphRDMEModel(birth_rates, fill(k_d, length(birth_rates)))

# ── Build RDME system ──────────────────────────────────────────────────────────

"""
    build_graph_rdme_system(model, mesh) → DiscreteStochasticSystem{CartesianIndex{K}}

Build the full RDME `DiscreteStochasticSystem` for a graph topology:
  - Birth  ∅ → A  in voxel k at rate `birth_rates[k]`
  - Death  A → ∅  in voxel k at rate `death_rate × n_k`
  - Diffusion A_u ⇌ A_v along each edge at rate `D_rates[e] × n_u` (or n_v)

The returned system passes directly to `build_generator`.
Propensities capture model parameters in closures; pass `Float64[]` as rates.
"""
function build_graph_rdme_system(model::GraphRDMEModel, mesh::VoxelMesh)
    K = mesh.n_voxels
    @assert length(model.birth_rates) == K

    stoichs      = CartesianIndex{K}[]
    propensities = Function[]

    # Per-voxel birth and death
    for k in 1:K
        let k=k, kb=model.birth_rates[k], kd=model.death_rates[k]
            ek = _e(k, Val(K))
            push!(stoichs, ek);  push!(propensities, (x, r, t) -> kb)
            push!(stoichs, -ek); push!(propensities, (x, r, t) -> kd * max(0, Tuple(x)[k]))
        end
    end

    # Diffusion along each edge (both directions)
    for (e, (u, v)) in enumerate(mesh.edges)
        let u=u, v=v, d=mesh.D_rates[e]
            eu = _e(u, Val(K)); ev = _e(v, Val(K))
            push!(stoichs, -eu + ev); push!(propensities, (x, r, t) -> d * max(0, Tuple(x)[u]))
            push!(stoichs,  eu - ev); push!(propensities, (x, r, t) -> d * max(0, Tuple(x)[v]))
        end
    end

    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, propensities)
end

# ── CoarseningPlan ─────────────────────────────────────────────────────────────

"""
    CoarseningPlan

Greedy pairing of voxels, collapsing K fine voxels to K' super-voxels.

Fields
------
- `voxel_map` : fine voxel index → coarse super-voxel index (length K)
- `K_coarse`  : number of super-voxels K' ≤ K
- `pairs`     : merged (fine_a, fine_b) pairs
- `singles`   : fine voxels that stay as singleton super-voxels
"""
struct CoarseningPlan
    voxel_map :: Vector{Int}
    K_coarse  :: Int
    pairs     :: Vector{Tuple{Int,Int}}
    singles   :: Vector{Int}
end

"""
    make_coarsening_plan(mesh, mean_counts; threshold=0.0)

Flux-based greedy coarsening: sort edges by weight `D_rate × min(μ_u, μ_v)`,
greedily merge voxel pairs with the highest flux.

Voxels with low/zero flux (empty voxels) become singletons.
Returns a `CoarseningPlan`.
"""
function make_coarsening_plan(mesh::VoxelMesh,
                              mean_counts::AbstractVector{Float64};
                              threshold::Float64 = 0.0)
    K = mesh.n_voxels
    weights = [mesh.D_rates[e] * min(mean_counts[u], mean_counts[v])
               for (e, (u,v)) in enumerate(mesh.edges)]

    used   = fill(false, K)
    pairs  = Tuple{Int,Int}[]
    singles= Int[]

    for e in sortperm(weights; rev=true)
        weights[e] < threshold && break
        u, v = mesh.edges[e]
        (used[u] || used[v]) && continue
        push!(pairs, (u, v))
        used[u] = used[v] = true
    end
    for v in 1:K; used[v] || push!(singles, v); end

    vmap = zeros(Int, K)
    g = 0
    for (a,b) in pairs; g += 1; vmap[a] = g; vmap[b] = g; end
    for s in singles;   g += 1; vmap[s] = g; end

    CoarseningPlan(vmap, g, pairs, singles)
end

# ── Model and mesh coarsening ─────────────────────────────────────────────────

"""
    coarsen_model(model, plan) → GraphRDMEModel

Coarse birth rates = sum of constituent fine-voxel birth rates.
Death rate is unchanged (per-molecule, no coarsening needed).
"""
function coarsen_model(model::GraphRDMEModel, plan::CoarseningPlan)
    coarse_births = zeros(plan.K_coarse)
    coarse_deaths = zeros(plan.K_coarse)
    cnts          = zeros(Int, plan.K_coarse)
    for (v, kb) in enumerate(model.birth_rates)
        g = plan.voxel_map[v]
        coarse_births[g] += kb
        coarse_deaths[g] += model.death_rates[v]
        cnts[g] += 1
    end
    # Death rate: average over merged voxels (per-molecule rate is intensive)
    coarse_deaths ./= max.(cnts, 1)
    GraphRDMEModel(coarse_births, coarse_deaths)
end

"""
    coarsen_mesh(mesh, plan) → VoxelMesh

Drop intra-super-voxel edges; keep inter-super-voxel edges (deduplicated).
Average voxel positions for each super-voxel.
"""
function coarsen_mesh(mesh::VoxelMesh, plan::CoarseningPlan)
    K_c  = plan.K_coarse
    seen = Set{Tuple{Int,Int}}()
    edges  = Tuple{Int,Int}[]
    drates = Float64[]

    for (e, (u,v)) in enumerate(mesh.edges)
        gu, gv = plan.voxel_map[u], plan.voxel_map[v]
        gu == gv && continue
        key = gu < gv ? (gu,gv) : (gv,gu)
        key ∈ seen && continue
        push!(seen, key); push!(edges, key); push!(drates, mesh.D_rates[e])
    end

    positions = if isnothing(mesh.positions)
        nothing
    else
        pos  = zeros(K_c, size(mesh.positions, 2))
        cnts = zeros(Int, K_c)
        for v in 1:mesh.n_voxels
            g = plan.voxel_map[v]
            pos[g, :] .+= mesh.positions[v, :]
            cnts[g] += 1
        end
        pos ./ max.(cnts, 1)
    end

    VoxelMesh(K_c, edges, drates, positions)
end

# ── State-space restriction ────────────────────────────────────────────────────

"""
    restrict_state(state, plan) → CartesianIndex{K'}

Map fine CartesianIndex{K} → coarse CartesianIndex{K'} by summing molecule
counts of voxels within each super-voxel.
"""
@inline function restrict_state(state::CartesianIndex{K},
                                 plan::CoarseningPlan) where K
    t      = Tuple(state)
    coarse = zeros(Int, plan.K_coarse)
    @inbounds for v in 1:K
        coarse[plan.voxel_map[v]] += t[v]
    end
    CartesianIndex(Tuple(coarse))
end

"""
    restrict(sp, plan) → StateSpace{CartesianIndex{K'}, T}

Exact restriction: map fine state space to coarse by summing probabilities
of fine states that share the same coarse representative.
"""
function restrict(sp::StateSpace{CartesianIndex{K},T},
                  plan::CoarseningPlan) where {K,T}
    K_c    = plan.K_coarse
    coarse = StateSpace{CartesianIndex{K_c},T}()
    for (st, prob) in zip(sp.states, sp.probs)
        cs  = restrict_state(st, plan)
        idx = get(coarse.index, cs, 0)
        if idx == 0
            add_state!(coarse, cs, prob)
        else
            coarse.probs[idx] += prob
        end
    end
    coarse
end

# ── State-space prolongation ──────────────────────────────────────────────────

"""
    prolong(coarse_sp, plan, fine_K; split=:equal) → StateSpace{CartesianIndex{K}}

Prolong coarse state space back to fine level.

Each merged-pair super-voxel with total count N is split between two fine voxels:
- `:equal` (default): deterministic midpoint N÷2, N-N÷2
- `:all`: enumerate all Binomial(N,½) splits (exact but expensive for large N)
"""
function prolong(coarse_sp::StateSpace{CartesianIndex{K_c},T},
                 plan::CoarseningPlan,
                 fine_K::Int;
                 split::Symbol = :equal) where {K_c, T}

    # Inverse map: coarse super-voxel → list of fine voxels
    sv_contents = [Int[] for _ in 1:plan.K_coarse]
    for v in 1:fine_K
        push!(sv_contents[plan.voxel_map[v]], v)
    end

    fine_sp = StateSpace{CartesianIndex{fine_K},T}()

    for (cs, cp) in zip(coarse_sp.states, coarse_sp.probs)
        ct = Tuple(cs)

        if split === :equal
            fine_t = zeros(Int, fine_K)
            for g in 1:plan.K_coarse
                fvs = sv_contents[g]
                if length(fvs) == 1
                    fine_t[fvs[1]] = ct[g]
                else                         # merged pair: midpoint split
                    N = ct[g]
                    fine_t[fvs[1]] = N ÷ 2
                    fine_t[fvs[2]] = N - N ÷ 2
                end
            end
            fs  = CartesianIndex(Tuple(fine_t))
            idx = get(fine_sp.index, fs, 0)
            if idx == 0; add_state!(fine_sp, fs, T(cp))
            else;        fine_sp.probs[idx] += T(cp); end

        elseif split === :all
            states_weights = Tuple{Vector{Int},Float64}[(zeros(Int, fine_K), 1.0)]
            for g in 1:plan.K_coarse
                fvs = sv_contents[g]
                if length(fvs) == 1
                    for (ft,_) in states_weights; ft[fvs[1]] = ct[g]; end
                else
                    N = ct[g]; coeff = 2.0^(-N)
                    new_sw = Tuple{Vector{Int},Float64}[]
                    for (ft, w) in states_weights
                        for na in 0:N
                            ft2 = copy(ft); ft2[fvs[1]] = na; ft2[fvs[2]] = N - na
                            push!(new_sw, (ft2, w * coeff * binomial(N, na)))
                        end
                    end
                    states_weights = new_sw
                end
            end
            for (ft, w) in states_weights
                fs  = CartesianIndex(Tuple(ft))
                idx = get(fine_sp.index, fs, 0)
                if idx == 0; add_state!(fine_sp, fs, T(cp * w))
                else;        fine_sp.probs[idx] += T(cp * w); end
            end
        else
            error("Unknown split mode: $split")
        end
    end

    fine_sp
end

# ── SpatialCMEState ────────────────────────────────────────────────────────────

"""
    SpatialCMEState{K}

Full joint CME distribution over K voxels + associated model + spatial graph.

Fields
------
- `sp`    : StateSpace{CartesianIndex{K}, Float64}
- `model` : per-voxel birth rates + death rate
- `mesh`  : spatial graph at this level
- `level` : multigrid level (0 = finest)
"""
struct SpatialCMEState{K}
    sp    :: StateSpace{CartesianIndex{K}, Float64}
    model :: GraphRDMEModel
    mesh  :: VoxelMesh
    level :: Int
end

"""
    SpatialCMEState(model, mesh; n_max)

Initialise with all-zero state and expand once using the RDME reactions.
"""
function SpatialCMEState(model::GraphRDMEModel, mesh::VoxelMesh; n_max::Int=40)
    K   = mesh.n_voxels
    @assert length(model.birth_rates) == K
    sp  = StateSpace{CartesianIndex{K}, Float64}()
    sys = build_graph_rdme_system(model, mesh)
    bc  = x -> all(c -> 0 ≤ c ≤ n_max, Tuple(x))
    add_state!(sp, CartesianIndex(ntuple(_ -> 0, Val(K))), 1.0)
    expand!(sp, sys, bc; depth=1)
    SpatialCMEState{K}(sp, model, mesh, 0)
end

"""Set IC: all probability at molecule count N0 in voxel v0."""
function set_voxel_ic!(state::SpatialCMEState{K}, v0::Int, N0::Int) where K
    sp = state.sp
    empty!(sp.states); empty!(sp.probs); empty!(sp.index); empty!(sp.global_id)
    add_state!(sp, CartesianIndex(ntuple(k -> k == v0 ? N0 : 0, Val(K))), 1.0)
end

# ── V-cycle step ───────────────────────────────────────────────────────────────

"""
    spatial_vcycle!(state, dt; kwargs...) → SpatialCMEState{K}

One multigrid V-cycle on the spatial CME:

1. **Pre-smooth**:  evolve fine distribution for τ_pre
2. **Coarsen**:     greedy flux-based pairing (make_coarsening_plan)
3. **Restrict**:    exact marginal onto coarse state space
4. **Coarse solve**: build coarse RDME, evolve for dt - τ_pre - τ_post
5. **Prolong**:     midpoint-split back to fine
6. **Post-smooth**: evolve fine distribution for τ_post
7. **Expand**:      add boundary states reachable from the new distribution
"""
function spatial_vcycle!(state::SpatialCMEState{K},
                          dt::Float64;
                          τ_pre::Float64     = 0.1 * dt,
                          τ_post::Float64    = 0.1 * dt,
                          threshold::Float64 = 0.0,
                          n_max::Int         = 40,
                          krylov_m::Int      = 30,
                          prune_tol::Float64 = 1e-14) where K

    sp   = state.sp
    mod  = state.model
    mesh = state.mesh
    RATES = Float64[]
    bc    = x -> all(c -> 0 ≤ c ≤ n_max, Tuple(x))

    # ── 1. Pre-smooth ──────────────────────────────────────────────────────────
    if τ_pre > 0 && length(sp) > 0
        sys = build_graph_rdme_system(mod, mesh)
        A,  = build_generator(sp, sys, RATES, 0.0)
        sp.probs .= expv(τ_pre, A, sp.probs; m=krylov_m)
        sp.probs .= max.(0.0, sp.probs)
        s = sum(sp.probs); s > 0 && (sp.probs ./= s)
        prune_threshold!(sp, prune_tol)
    end

    # ── 2. Coarsen ─────────────────────────────────────────────────────────────
    μ          = mean_counts(state)
    plan       = make_coarsening_plan(mesh, μ; threshold)
    mod_c      = coarsen_model(mod,  plan)
    mesh_c     = coarsen_mesh(mesh, plan)

    # ── 3. Restrict ────────────────────────────────────────────────────────────
    coarse_sp  = restrict(sp, plan)

    # ── 4. Coarse solve ────────────────────────────────────────────────────────
    dt_c = dt - τ_pre - τ_post
    if dt_c > 0 && length(coarse_sp) > 0
        _coarse_solve!(coarse_sp, mod_c, mesh_c, dt_c, krylov_m, RATES, prune_tol)
    end

    # ── 5. Prolong ─────────────────────────────────────────────────────────────
    fine_sp = prolong(coarse_sp, plan, K; split=:equal)
    # keep any fine states from before that weren't in the coarse solution
    for s in sp.states
        get(fine_sp.index, s, 0) == 0 && add_state!(fine_sp, s, 0.0)
    end

    # ── 6. Post-smooth ─────────────────────────────────────────────────────────
    if τ_post > 0 && length(fine_sp) > 0
        sys = build_graph_rdme_system(mod, mesh)
        A,  = build_generator(fine_sp, sys, RATES, 0.0)
        fine_sp.probs .= expv(τ_post, A, fine_sp.probs; m=krylov_m)
        fine_sp.probs .= max.(0.0, fine_sp.probs)
        s = sum(fine_sp.probs); s > 0 && (fine_sp.probs ./= s)
        prune_threshold!(fine_sp, prune_tol)
    end

    # ── 7. Expand ───────────────────────────────────────────────────────────────
    sys2 = build_graph_rdme_system(mod, mesh)
    expand!(fine_sp, sys2, bc; depth=1)

    SpatialCMEState{K}(fine_sp, mod, mesh, state.level)
end

# Function barrier: coarse solve with K_c determined at runtime
function _coarse_solve!(sp::StateSpace{CartesianIndex{K_c},Float64},
                        model::GraphRDMEModel,
                        mesh::VoxelMesh,
                        dt::Float64,
                        m::Int,
                        rates::Vector{Float64},
                        prune_tol::Float64) where K_c
    length(sp) == 0 && return
    sys = build_graph_rdme_system(model, mesh)
    A,  = build_generator(sp, sys, rates, 0.0)
    sp.probs .= expv(dt, A, sp.probs; m=m)
    sp.probs .= max.(0.0, sp.probs)
    s = sum(sp.probs); s > 0 && (sp.probs ./= s)
    prune_threshold!(sp, prune_tol)
end

# ── Mean counts (exact FSP, no mean-field) ────────────────────────────────────

"""
    mean_counts(state::SpatialCMEState{K}) → Vector{Float64}

Exact FSP mean: ⟨nₖ⟩ = Σ_{configurations n} nₖ · p(n).
"""
function mean_counts(state::SpatialCMEState{K}) where K
    μ = zeros(K)
    for (ci, p) in zip(state.sp.states, state.sp.probs)
        t = Tuple(ci)
        @inbounds for v in 1:K; μ[v] += t[v] * p; end
    end
    μ
end

"""
    mean_counts(sp::StateSpace{CartesianIndex{K}}) → Vector{Float64}

Overload accepting a `StateSpace` directly.
"""
function mean_counts(sp::StateSpace{CartesianIndex{K}, Float64}) where K
    μ = zeros(K)
    for (ci, p) in zip(sp.states, sp.probs)
        t = Tuple(ci)
        @inbounds for v in 1:K; μ[v] += t[v] * p; end
    end
    μ
end
