# Mixed fine/lumped joint FSP over an arbitrary VoxelGraph.
#
# Generalises the proven 1-D windowed chain stepper (examples/schlogl_chain_basin.jl)
# to any geometry. The genuine joint distribution is carried ONLY over the set of
# ACTIVE voxels (the moving front); every settled voxel is collapsed to a single
# reservoir basin label (1=lo, 2=hi) and feeds its active neighbours at its typical
# count D·n_basin. Collapsing a settled voxel = the multigrid restriction; promoting
# a reservoir back to counts = the prolongation. Because the active window stays a
# fixed width as the front travels, |S| = n_max^(#active) stays bounded regardless
# of how large the graph is — the plain full joint explodes.
#
# Geometry enters ONLY through graph.edges / adjacency_list, so the same stepper
# runs a 1-D chain, a 2-D grid droplet, or any graph. Colour / occupancy is the
# exact hi-probability P(n>n_un) read straight off the joint (never a mean field).

using ExponentialUtilities: expv

"""
    MixedFSP{K}

Mixed fine/lumped joint FSP for Schlögl kinetics on a VoxelGraph. Active voxels
carry the genuine joint over `sp`; reservoir voxels are pinned to a basin label.
"""
mutable struct MixedFSP{K}
    graph    :: VoxelGraph
    nb       :: Vector{Vector{Int}}        # adjacency list (cached)
    c1 :: Float64; c2 :: Float64; c3 :: Float64; c4 :: Float64; D :: Float64
    n_lo :: Int; n_un :: Int; n_hi :: Int  # lo basin / separatrix / hi basin
    Vact :: Set{Int}                       # active (genuine-joint) voxels
    res  :: Vector{Int}                    # basin label of every voxel (1=lo,2=hi)
    sp   :: StateSpace{CartesianIndex{K}, Float64}
    ε_prune :: Float64
    krylov_m :: Int
end

nrep(f::MixedFSP, b::Int) = b == 1 ? f.n_lo : f.n_hi

"""
    MixedFSP(model, graph; n_lo, n_un, n_hi, Vact, res, ε_prune, krylov_m)

Build a mixed stepper. `Vact` is the initial active set, `res` the initial basin
labels (defaults: all lo). The caller seeds `sp` afterwards via `set_seed!`.
"""
function MixedFSP(model::SchloglModel1D, graph::VoxelGraph;
                  n_lo::Int, n_un::Int, n_hi::Int,
                  Vact::Set{Int}, res::Union{Nothing,Vector{Int}} = nothing,
                  ε_prune::Float64 = 1e-6, krylov_m::Int = 30)
    K = graph.n_voxels
    r = res === nothing ? fill(1, K) : copy(res)
    MixedFSP{K}(graph, adjacency_list(graph),
                model.c1, model.c2, model.c3, model.c4, model.D,
                n_lo, n_un, n_hi, copy(Vact), r,
                StateSpace{CartesianIndex{K}, Float64}(), ε_prune, krylov_m)
end

"""Seed the joint with `prob` at configuration `u0` (a CartesianIndex{K})."""
function set_seed!(f::MixedFSP{K}, u0::CartesianIndex{K}, prob::Float64 = 1.0) where {K}
    add_state!(f.sp, u0, prob)
    f
end

# ── mixed generator: full Schlögl+diffusion on Vact, reservoir flux elsewhere ──
function build_mixed_generator(f::MixedFSP{K}) where {K}
    c1, c2, c3, c4, D = f.c1, f.c2, f.c3, f.c4, f.D
    Vact = f.Vact; res = f.res
    stoichs = CartesianIndex{K}[]; props = Function[]
    eu(k) = CartesianIndex(ntuple(i -> i == k ? 1 : 0, K))
    for k in Vact
        ek = eu(k)
        push!(stoichs,  ek); push!(props, (x,rv,t) -> c3)
        push!(stoichs, -ek); push!(props, (x,rv,t) -> c4*max(0,Tuple(x)[k]))
        push!(stoichs,  ek); push!(props, (x,rv,t) -> (n=Tuple(x)[k]; c1*max(0,n)*max(0,n-1)/2))
        push!(stoichs, -ek); push!(props, (x,rv,t) -> (n=Tuple(x)[k]; c2*max(0,n)*max(0,n-1)*max(0,n-2)/6))
    end
    for (i,j) in f.graph.edges
        ei = eu(i); ej = eu(j)
        ia = i in Vact; ja = j in Vact
        if ia && ja                       # active–active: ordinary hop
            push!(stoichs, -ei+ej); push!(props, (x,rv,t) -> D*max(0,Tuple(x)[i]))
            push!(stoichs,  ei-ej); push!(props, (x,rv,t) -> D*max(0,Tuple(x)[j]))
        elseif ia && !ja                  # j reservoir: out vanishes, in at D·n_res
            push!(stoichs, -ei); push!(props, (x,rv,t) -> D*max(0,Tuple(x)[i]))
            push!(stoichs,  ei); push!(props, (x,rv,t) -> D*nrep(f, res[j]))
        elseif ja && !ia                  # i reservoir
            push!(stoichs, -ej); push!(props, (x,rv,t) -> D*max(0,Tuple(x)[j]))
            push!(stoichs,  ej); push!(props, (x,rv,t) -> D*nrep(f, res[i]))
        end                               # reservoir–reservoir: nothing
    end
    DiscreteStochasticSystem{CartesianIndex{K}}(stoichs, props)
end

# ── per-voxel hi-probability P(n > n_un) read straight from the joint ─────────
function hi_prob(f::MixedFSP{K}) where {K}
    ϕ = zeros(K)
    for v in 1:K
        v in f.Vact || (ϕ[v] = f.res[v] == 2 ? 1.0 : 0.0)
    end
    for (s,p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(s)
        for v in f.Vact
            t[v] > f.n_un && (ϕ[v] += p)
        end
    end
    ϕ
end

# ── slide the window: promote lo reservoirs the front reached, demote settled hi
function reclassify!(f::MixedFSP{K}, ϕ::Vector{Float64};
                     promote_thresh::Float64 = 0.5,
                     demote_thresh::Float64 = 0.9) where {K}
    Vact, res, nb = f.Vact, f.res, f.nb
    # PROMOTE a lo reservoir to keep a count-resolved low margin ahead of the front.
    # The leading active voxel is always partly drained by the frozen lo reservoir
    # behind it, so its ϕ stalls and can never itself trigger promotion. Instead fire
    # off the SECOND active voxel from the edge: promote reservoir r if its active
    # neighbour u has ANOTHER active neighbour w≠r that is solidly rising (ϕ>thresh).
    for r in 1:K
        (r in Vact) && continue
        res[r] == 1 || continue
        promote = false
        for u in nb[r]
            (u in Vact) || continue
            for w in nb[u]
                if w != r && (w in Vact) && ϕ[w] > promote_thresh
                    promote = true; break
                end
            end
            promote && break
        end
        promote && push!(Vact, r)
    end
    # DEMOTE a settled-high active voxel whose ALL graph-neighbours are also high →
    # hi reservoir (the high bulk behind the front collapses to one frozen value).
    is_hi(u) = (u in Vact) ? ϕ[u] > 0.8 : res[u] == 2
    for v in collect(Vact)
        if ϕ[v] > demote_thresh && all(is_hi, nb[v])
            delete!(Vact, v); res[v] = 2
        end
    end
    f
end

# ── collapse demoted coords to their reservoir value, merge duplicates ────────
function project_to_window!(f::MixedFSP{K}) where {K}
    Vact, res = f.Vact, f.res
    acc = Dict{CartesianIndex{K}, Float64}()
    for (s,p) in zip(f.sp.states, f.sp.probs)
        c = CartesianIndex(ntuple(v -> v in Vact ? Tuple(s)[v] : nrep(f, res[v]), K))
        acc[c] = get(acc, c, 0.0) + p
    end
    out = StateSpace{CartesianIndex{K}, Float64}()
    for (c,p) in acc
        add_state!(out, c, p)
    end
    renormalize!(out)
    f.sp = out
    f
end

"""
    step!(f, dt; reclassify=true)

Advance the mixed joint by `dt`: expand support, build the mixed generator, evolve
with expv, prune, then (if `reclassify`) slide the active window and re-project the
settled voxels onto reservoir labels.
"""
function step!(f::MixedFSP{K}, dt::Float64; reclassify::Bool = true) where {K}
    sys = build_mixed_generator(f)
    expand!(f.sp, sys, x -> all(>=(0), Tuple(x)); depth = 1)
    A, = build_generator(f.sp, sys, nothing, 0.0; absorbing=false)
    f.sp.probs .= expv(dt, A, f.sp.probs; m = clamp(f.krylov_m, 4, length(f.sp)))
    f.sp.probs .= max.(0.0, f.sp.probs)
    prune_threshold!(f.sp, f.ε_prune); renormalize!(f.sp)
    if reclassify
        ϕ = hi_prob(f)
        reclassify!(f, ϕ)
        project_to_window!(f)
    end
    f
end

# ── per-voxel marginal P_v(n) read straight off the mixed joint ───────────────
function mixed_voxel_marginal(f::MixedFSP{K}, v::Int; n_max::Int) where {K}
    m = zeros(n_max + 1)
    if v in f.Vact
        for (s,p) in zip(f.sp.states, f.sp.probs)
            n = Tuple(s)[v]
            0 <= n <= n_max && (m[n + 1] += p)
        end
    else
        nr = nrep(f, f.res[v])
        nr <= n_max && (m[nr + 1] = 1.0)
    end
    m
end

mixed_voxel_marginals(f::MixedFSP{K}; n_max::Int) where {K} =
    [mixed_voxel_marginal(f, v; n_max=n_max) for v in 1:f.graph.n_voxels]
