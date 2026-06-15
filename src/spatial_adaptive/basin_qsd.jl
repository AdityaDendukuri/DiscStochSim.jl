"""
basin_qsd.jl

Basin-resolved quasi-stationary lumping for bistable 1-species birth–death
nodes (the Phase-2 reservoir of the graph-multigrid hybrid).

A bistable voxel relaxes fast *within* a basin and switches *slowly* between
basins (NCD / timescale separation). Strong (total-count) lumping collapses such
a voxel to a single mean that sits on the unstable separatrix — destroying
bistability. Instead we represent a lumped voxel by a discrete basin label
`b ∈ {lo, hi}`, with:

  • `qsd[b]`     : the quasi-stationary distribution inside basin b
                  (lead eigenvector of the generator restricted/killed at the
                  separatrix — same kernel as `_build_patch_qsd_2s`),
  • `escape[b]`  : the basin-escape rate (= −λ_lead), i.e. the slow switching
                  channel rate lo↔hi,
  • `mean[b]`    : ⟨n⟩ under qsd[b]  (the value the lump emits as diffusion flux),
  • `birth/death`: nonlinear propensities AVERAGED under qsd[b] → constants
                  (weak lumpability: nonlinear rates → effective constants under
                  the conditional stationary law).

Rates are read straight off a column-convention birth–death generator
(`L[j,i]` = rate i→j), so this is model-agnostic: any 1-species bistable system
whose single-voxel generator is available (e.g. `schlogl_generator`).
"""

struct BasinQSD
    n_un   :: Int                         # separatrix (unstable fixed point)
    ns     :: NTuple{2, UnitRange{Int}}   # (lo range, hi range) of counts
    qsd    :: NTuple{2, Vector{Float64}}  # QSD over ns[b]
    mean   :: NTuple{2, Float64}          # ⟨n⟩ under qsd[b]
    var    :: NTuple{2, Float64}          # Var(n) under qsd[b]
    escape :: NTuple{2, Float64}          # escape rate −λ_lead  (lo→hi, hi→lo)
    birth  :: NTuple{2, Float64}          # QSD-averaged birth propensity
    death  :: NTuple{2, Float64}          # QSD-averaged death propensity
end

"""
    build_basin_qsd(L, n_un; n_max=size(L,1)-1) -> BasinQSD

`L` is a column-convention single-voxel birth–death generator
(`L[n+2,n+1]` = birth rate from n, `L[n,n+1]` = death rate from n). `n_un` is the
separatrix count. Basins are a DISJOINT partition consistent with `basin_of`:
lo = counts `[0, n_un]`, hi = `[n_un+1, n_max]`. Each restricted block is killed at
its boundary; the lead eigenvector is that basin's QSD.
"""
function build_basin_qsd(L::AbstractMatrix, n_un::Int; n_max::Int = size(L, 1) - 1)
    ranges = (0:n_un, (n_un + 1):n_max)
    qsd  = Vector{Vector{Float64}}(undef, 2)
    mean = zeros(2); var = zeros(2); esc = zeros(2); bir = zeros(2); dea = zeros(2)
    for (b, ns) in enumerate(ranges)
        idx = collect(ns) .+ 1
        B   = L[idx, idx]                         # restricted (killed at separatrix)
        vals, vecs = eigen(Matrix(B))
        lead = argmax(real.(vals))                # eigenvalue nearest 0
        π    = abs.(real.(vecs[:, lead])); π ./= sum(π)
        qsd[b]  = π
        m       = sum(collect(ns) .* π)
        mean[b] = m
        var[b]  = sum((collect(ns) .- m).^2 .* π)
        esc[b]  = -real(vals[lead])               # > 0 : basin escape
        # birth/death read off the generator off-diagonals, averaged under QSD
        bsum = 0.0; dsum = 0.0
        for (k, n) in enumerate(ns)
            b_n = n < n_max ? L[n + 2, n + 1] : 0.0   # rate n → n+1
            d_n = n > 0     ? L[n,     n + 1] : 0.0   # rate n → n−1
            bsum += π[k] * b_n
            dsum += π[k] * d_n
        end
        bir[b] = bsum; dea[b] = dsum
    end
    BasinQSD(n_un, (ranges[1], ranges[2]),
             (qsd[1], qsd[2]),
             (mean[1], mean[2]), (var[1], var[2]),
             (esc[1], esc[2]), (bir[1], bir[2]), (dea[1], dea[2]))
end

"""
    basin_of(bq::BasinQSD, n) -> 1 (lo) or 2 (hi)

Classify a count `n` into its basin by the separatrix.
"""
basin_of(bq::BasinQSD, n::Integer) = n <= bq.n_un ? 1 : 2

# ─── restrict / prolong one voxel axis of a joint StateSpace ────────────────────
#
# A lumped axis stores the basin LABEL (0 = lo, 1 = hi) in place of the count.
# `restrict_voxel` marginalizes voxel v's count into its basin label (lossy only
# WITHIN a basin — the basin masses P(lo)/P(hi) are preserved exactly). The full
# joint over all OTHER voxels, and the correlation between v's basin and the rest,
# are retained. `prolong_voxel` reconstructs the count axis by spreading each
# label's mass over that basin's QSD — restoring within-basin variance the mean
# representative would drop, which is what keeps the distribution bimodal.

@inline function _setaxis(s::CartesianIndex{K}, v::Int, val::Int) where {K}
    CartesianIndex(ntuple(i -> i == v ? val : s[i], K))
end

"""
    within_basin_coupling(sp, v, w, bq) -> Float64

The within-basin covariance `Cov(n_v, n_w)` that lumping voxel `v` would DESTROY:
the total `Cov(n_v,n_w)` minus the basin-level part (the across-basin correlation
that the basin label preserves). This is exactly the error a lump incurs at the
`v↔w` edge, so it is the quantity the restrict criterion must bound.
"""
function within_basin_coupling(sp::StateSpace{CartesianIndex{K}}, v::Int, w::Int,
                               bq::BasinQSD) where {K}
    μv = μw = 0.0; Cvw = 0.0
    # basin-level moments: replace each count by its within-basin mean given its label
    bμ = (0.0, 0.0); bP = (0.0, 0.0)   # placeholder, recomputed below
    sv = zeros(2); sw = zeros(2); Pb = zeros(2)          # per-basin sums for v, w
    for (s, p) in zip(sp.states, sp.probs)
        nv = s[v]; nw = s[w]
        μv += nv*p; μw += nw*p; Cvw += nv*nw*p
        b = basin_of(bq, nv); sv[b] += nv*p; sw[b] += nw*p; Pb[b] += p
    end
    Cvw -= μv*μw                                          # total Cov
    # basin-level Cov: use v's basin-mean as its coarse value
    mb = (Pb[1] > 0 ? sv[1]/Pb[1] : 0.0, Pb[2] > 0 ? sv[2]/Pb[2] : 0.0)
    Cbasin = 0.0
    Cbasin += Pb[1]*mb[1]*(Pb[1] > 0 ? sw[1]/Pb[1] : 0.0)
    Cbasin += Pb[2]*mb[2]*(Pb[2] > 0 ? sw[2]/Pb[2] : 0.0)
    Cbasin -= (Pb[1]*mb[1] + Pb[2]*mb[2]) * μw
    Cvw - Cbasin                                          # within-basin part (lost)
end

"""
    lumpable_voxels(sp, graph, bq; ε_cov=1e-2, skip=()) -> Vector{Int}

Voxels safe to lump: those whose within-basin coupling to every graph NEIGHBOUR is
below `ε_cov` (so the covariance lost at each incident edge is bounded). This is the
NCD/weak-lumpability condition measured directly, not a QSD-shape assumption. Pass
`skip` to exclude voxels (e.g. the active front) unconditionally.
"""
function lumpable_voxels(sp::StateSpace{CartesianIndex{K}}, graph, bq::BasinQSD;
                         ε_cov::Float64 = 1e-2, skip = ()) where {K}
    nbrs = [Int[] for _ in 1:K]
    for (a, b) in graph.edges; push!(nbrs[a], b); push!(nbrs[b], a); end
    out = Int[]
    for v in 1:K
        v in skip && continue
        all(abs(within_basin_coupling(sp, v, w, bq)) <= ε_cov for w in nbrs[v]) &&
            push!(out, v)
    end
    out
end

# ─── LumpInfo: the smart-prolongation operator (carries pre-coarsening info) ────
#
# When we coarsen voxel v we already SEE its exact within-basin shape in the joint
# — so we store it and replay it on prolongation. The prolongation operator IS the
# pre-coarsening fine structure: no assumed/QSD shape, means & variance restored
# exactly. The basin label carries the bistability (which attractor); the stored
# conditional carries the within-basin detail; `escape` (the one QSD-derived
# scalar) carries the slow basin-switch for the coarse evolve. For a basin the
# voxel had not visited at coarsening (`P==0`), we fall back to the isolated QSD.

struct LumpInfo
    v      :: Int
    ns     :: NTuple{2, UnitRange{Int}}
    cond   :: NTuple{2, Vector{Float64}}   # per-basin conditional fine shape (normalized)
    P      :: NTuple{2, Float64}           # basin masses at coarsening time
    mean   :: NTuple{2, Float64}           # ⟨n⟩ within each basin (for diffusion flux)
    escape :: NTuple{2, Float64}           # basin-switch rates lo→hi, hi→lo
end

"""
    restrict_voxel(sp, v, bq) -> (sp_reduced, ::LumpInfo)

Collapse voxel `v`'s count axis to its basin label (0=lo, 1=hi). Mass is conserved
and basin masses / cross-voxel correlations are exact. The returned `LumpInfo` is
the smart-prolongation operator: it records the MEASURED within-basin conditional
shape so `prolong_voxel` can restore it exactly.
"""
function restrict_voxel(sp::StateSpace{CartesianIndex{K}, T}, v::Int,
                        bq::BasinQSD) where {K, T}
    out  = StateSpace{CartesianIndex{K}, T}()
    cond = (zeros(length(bq.ns[1])), zeros(length(bq.ns[2])))
    P    = zeros(2)
    for (s, p) in zip(sp.states, sp.probs)
        n = s[v]; b = basin_of(bq, n)
        ℓ  = b - 1
        sr = _setaxis(s, v, ℓ)
        idx = get(out.index, sr, 0)
        idx == 0 ? add_state!(out, sr, p) : (out.probs[idx] += p)
        cond[b][n - first(bq.ns[b]) + 1] += p
        P[b] += p
    end
    c  = Vector{Vector{Float64}}(undef, 2)
    mn = zeros(2)
    for b in 1:2
        if P[b] > 0
            c[b]  = cond[b] ./ P[b]                 # measured conditional shape
            mn[b] = sum(collect(bq.ns[b]) .* c[b])
        else
            c[b]  = copy(bq.qsd[b])                 # unseen basin → QSD fallback
            mn[b] = bq.mean[b]
        end
    end
    li = LumpInfo(v, bq.ns, (c[1], c[2]), (P[1], P[2]),
                  (mn[1], mn[2]), bq.escape)
    (out, li)
end

"""
    prolong_voxel(sp_red, li::LumpInfo; wtol=1e-14) -> StateSpace

Inverse of `restrict_voxel`: expand voxel `li.v`'s basin label back to a count axis
by replaying the stored per-basin conditional shape (the pre-coarsening info).
"""
function prolong_voxel(sp_red::StateSpace{CartesianIndex{K}, T}, li::LumpInfo;
                       wtol::Float64 = 1e-14) where {K, T}
    out = StateSpace{CartesianIndex{K}, T}()
    for (s, p) in zip(sp_red.states, sp_red.probs)
        b  = s[li.v] + 1                     # label 0/1 → basin 1/2
        ns = li.ns[b]; shape = li.cond[b]
        for (k, n) in enumerate(ns)
            w = p * shape[k]
            w > wtol || continue
            sf  = _setaxis(s, li.v, n)
            idx = get(out.index, sf, 0)
            idx == 0 ? add_state!(out, sf, w) : (out.probs[idx] += w)
        end
    end
    out
end
