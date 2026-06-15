"""
polarisome_fsp.jl

Spatial FSP for the Lawson et al. (2013) polarisome model — 2-species RDME on a
1-D periodic ring membrane.

Species:
  B = Bni1 (formin, membrane-bound)
  S = Spa2 (scaffold, membrane-bound)

Reactions per voxel k (S-model approximation):
  ∅ → B   rate = B_on*(B_pool - nB) + B_fb*nS   [Cdc42a + Spa2 feedback, finite pool]
  B → ∅   rate = B_off * nB
  ∅ → S   rate = S_on * nB                        [actin-cable transport, prop. to Bni1]
  S → ∅   rate = A_off * Km/(Km + nS) * nS        [Michaelis-Menten inhibition by Spa2]

Diffusion (both species, same D):
  hop rate D/dx² per neighbor (mean-field in-flux from neighbor means)

State: (nB, nS), flat index s = nB*(NS+1) + nS + 1 (1-indexed, Julia convention)

Bistable parameters give:
  lo ≈ (nB≈2, nS≈2)   dispersed background
  hi ≈ (nB≈10, nS≈30)  polarisome cluster
"""

using SparseArrays, ExponentialUtilities

# ─── Model ────────────────────────────────────────────────────────────────────

struct PolarisomeModel
    B_on   :: Float64   # s⁻¹  base Bni1 on-rate (Cdc42a-driven, per pool unit)
    B_fb   :: Float64   # s⁻¹/molecule  Spa2 → Bni1 feedback
    B_off  :: Float64   # s⁻¹/molecule  Bni1 off-rate
    B_pool :: Int       # max Bni1 from pool per voxel (finite cytoplasmic pool)
    S_on   :: Float64   # s⁻¹/molecule  Bni1 → Spa2 actin-transport rate
    A_off  :: Float64   # s⁻¹  basal Spa2 removal rate
    Km     :: Float64   # molecules  Michaelis-Menten constant (Spa2 self-inhibits removal)
    D      :: Float64   # s⁻¹  membrane diffusion D/dx² (same for B and S)
    n_B_max :: Int
    n_S_max :: Int
end

# Convenient defaults that give lo≈(2,2), hi≈(10,30)
function PolarisomeModel(; B_on=0.05, B_fb=0.15, B_off=0.5, B_pool=20,
                           S_on=0.5,  A_off=0.8, Km=8.0, D=1.0,
                           n_B_max=25, n_S_max=45)
    PolarisomeModel(B_on, B_fb, B_off, B_pool, S_on, A_off, Km, D, n_B_max, n_S_max)
end

# ─── State indexing ───────────────────────────────────────────────────────────

@inline state_dim(m::PolarisomeModel)        = (m.n_B_max+1)*(m.n_S_max+1)
@inline state_idx(m::PolarisomeModel, nB, nS) = nB*(m.n_S_max+1) + nS + 1

function state_means(m::PolarisomeModel, p::AbstractVector)
    NS1 = m.n_S_max + 1
    μB = 0.0; μS = 0.0
    @inbounds for nB in 0:m.n_B_max, nS in 0:m.n_S_max
        pk = p[nB*NS1 + nS + 1]
        μB += nB * pk
        μS += nS * pk
    end
    μB, μS
end

# ─── Single-voxel CME generator ───────────────────────────────────────────────

"""
    polarisome_generator(m, in_B, in_S, n_nb)

Build the (dim × dim) CME generator for one voxel.

  in_B  = D * Σ_{j∈neighbors} ⟨nB⟩_j   (mean-field Bni1 influx)
  in_S  = D * Σ_{j∈neighbors} ⟨nS⟩_j   (mean-field Spa2 influx)
  n_nb  = number of neighbors (outflux scales with n_nb)
"""
function polarisome_generator(m::PolarisomeModel,
                               in_B::Float64, in_S::Float64, n_nb::Int)
    NB  = m.n_B_max;  NS  = m.n_S_max
    NS1 = NS + 1
    dim = (NB+1)*NS1
    Q   = zeros(Float64, dim, dim)

    @inbounds for nB in 0:NB, nS in 0:NS
        s = nB*NS1 + nS + 1    # source state (1-indexed)

        # 1. Bni1 birth: (nB,nS) → (nB+1,nS)
        pool_avail = max(m.B_pool - nB, 0)
        λB = m.B_on * pool_avail + m.B_fb * nS + in_B
        if nB < NB && λB > 0.0
            t = s + NS1
            Q[t, s] += λB;  Q[s, s] -= λB
        end

        # 2. Bni1 death: (nB,nS) → (nB-1,nS)  [removal + n_nb diffusion-out channels]
        μB = (m.B_off + n_nb * m.D) * nB
        if nB > 0 && μB > 0.0
            t = s - NS1
            Q[t, s] += μB;  Q[s, s] -= μB
        end

        # 3. Spa2 birth: (nB,nS) → (nB,nS+1)
        λS = m.S_on * nB + in_S
        if nS < NS && λS > 0.0
            t = s + 1
            Q[t, s] += λS;  Q[s, s] -= λS
        end

        # 4. Spa2 death: (nB,nS) → (nB,nS-1)  [Michaelis-Menten + diffusion-out]
        menten = m.A_off * m.Km / (m.Km + nS)
        μS = (menten + n_nb * m.D) * nS
        if nS > 0 && μS > 0.0
            t = s - 1
            Q[t, s] += μS;  Q[s, s] -= μS
        end
    end
    Q
end

# ─── Fixed-point finder ────────────────────────────────────────────────────────

"""
Find the two stable fixed points of the deterministic rate equations (lo and hi basin).
Returns (nB_lo, nS_lo, nB_hi, nS_hi) as Float64.
"""
function polarisome_fixed_points(m::PolarisomeModel)
    # Nullcline intersection:
    # B-nullcline:  nB = (B_on*(B_pool-nB) + B_fb*nS) / B_off
    #             → nB*(B_off+B_on) = B_on*B_pool + B_fb*nS
    #             → nB = (B_on*B_pool + B_fb*nS) / (B_off+B_on)   [for nB < B_pool]
    # S-nullcline:  S_on*nB = A_off*Km/(Km+nS)*nS
    #             → nB = A_off*Km*nS / (S_on*(Km+nS))

    fps = Tuple{Float64,Float64}[]
    α = m.B_off + m.B_on
    for nS in LinRange(1e-6, Float64(m.n_S_max), 10_000)
        nB_Bnull = (m.B_on * m.B_pool + m.B_fb * nS) / α
        nB_Snull = m.A_off * m.Km * nS / (m.S_on * (m.Km + nS))
        push!(fps, (nB_Bnull - nB_Snull, nS))
    end
    # Find sign changes (zeros of nB_Bnull - nB_Snull)
    zeros_nS = Float64[]
    for i in 1:length(fps)-1
        if fps[i][1] * fps[i+1][1] < 0
            # Linear interpolation
            s1, s2 = fps[i][2], fps[i+1][2]
            v1, v2 = fps[i][1], fps[i+1][1]
            push!(zeros_nS, s1 - v1*(s2-s1)/(v2-v1))
        end
    end

    isempty(zeros_nS) && error("No fixed points found — check parameters")

    results = map(zeros_nS) do nS
        nB = (m.B_on * m.B_pool + m.B_fb * nS) / α
        (nB, nS)
    end

    # Sort by nS → first = lo, last = hi
    sort!(results; by = x -> x[2])
    nB_lo, nS_lo = results[1]
    nB_hi, nS_hi = results[end]
    nB_lo, nS_lo, nB_hi, nS_hi
end

# ─── Equilibrium distributions ────────────────────────────────────────────────

"""Compute equilibrium distribution by running a single-voxel CME from a delta IC at (nB0,nS0)."""
function polarisome_equil_dist(m::PolarisomeModel, nB0::Int, nS0::Int;
                                dt=0.5, T=200.0, krylov_m=30)
    dim = state_dim(m)
    p   = zeros(dim)
    p[state_idx(m, nB0, nS0)] = 1.0
    Q = polarisome_generator(m, 0.0, 0.0, 0)   # isolated voxel
    for _ in 1:round(Int, T/dt)
        p = expv(dt, Q, p; m=min(krylov_m, dim))
    end
    p ./= sum(p)
    p
end

# ─── Spatial FSP ──────────────────────────────────────────────────────────────

mutable struct PolarisomeFSP
    model     :: PolarisomeModel
    K         :: Int                          # ring voxels (periodic)
    dists     :: Dict{Int, Vector{Float64}}  # active voxels: full (nB,nS) CME dist
    hi_voxels :: Set{Int}                    # voxels in hi equilibrium
    lo_voxels :: Set{Int}                    # voxels in lo equilibrium
    equil_lo  :: Vector{Float64}             # precomputed lo equilibrium dist
    equil_hi  :: Vector{Float64}             # precomputed hi equilibrium dist
    nB_lo     :: Float64;  nS_lo :: Float64
    nB_hi     :: Float64;  nS_hi :: Float64
    nB_thresh :: Int;      nS_thresh :: Int  # basin midpoints for phi_hi
    t         :: Float64
    ε_expand  :: Float64   # excess-flux threshold to activate a lo voxel
    ε_phi_hi  :: Float64   # phi_hi threshold → deactivate to hi  (e.g. 0.65)
    ε_phi_lo  :: Float64   # phi_hi threshold → deactivate to lo  (e.g. 0.05)
    krylov_m  :: Int
end

function PolarisomeFSP(model::PolarisomeModel, K::Int;
                        ε_expand::Float64=0.1,
                        ε_phi_hi::Float64=0.65, ε_phi_lo::Float64=0.05,
                        krylov_m::Int=30)
    nB_lo, nS_lo, nB_hi, nS_hi = polarisome_fixed_points(model)
    @printf("  Fixed points: lo=(%.1f, %.1f)  hi=(%.1f, %.1f)\n",
            nB_lo, nS_lo, nB_hi, nS_hi)

    equil_lo = polarisome_equil_dist(model, round(Int,nB_lo), round(Int,nS_lo))
    equil_hi = polarisome_equil_dist(model, round(Int,nB_hi), round(Int,nS_hi))

    nB_thresh = round(Int, (nB_lo + nB_hi) / 2)
    nS_thresh = round(Int, (nS_lo + nS_hi) / 2)

    lo_voxels = Set(1:K)
    PolarisomeFSP(model, K,
                  Dict{Int, Vector{Float64}}(),
                  Set{Int}(),
                  lo_voxels,
                  equil_lo, equil_hi,
                  nB_lo, nS_lo, nB_hi, nS_hi,
                  nB_thresh, nS_thresh,
                  0.0, ε_expand, ε_phi_hi, ε_phi_lo, krylov_m)
end

function Base.setindex!(s::PolarisomeFSP, p::Vector{Float64}, k::Int)
    s.dists[k] = p
    delete!(s.lo_voxels, k)
    delete!(s.hi_voxels, k)
end

# Voxel mean (B,S) — from active dist, equil_hi, or equil_lo
function voxel_mean_BS(s::PolarisomeFSP, k::Int)
    if haskey(s.dists, k)
        return state_means(s.model, s.dists[k])
    elseif k in s.hi_voxels
        return s.nB_hi, s.nS_hi
    else
        return s.nB_lo, s.nS_lo
    end
end

# Ring neighbors (periodic)
@inline ring_neighbors(k::Int, K::Int) = (mod1(k-1, K), mod1(k+1, K))

# ─── step! ────────────────────────────────────────────────────────────────────

function step!(s::PolarisomeFSP, dt::Float64)
    m = s.model

    # ── Cache neighbor means ─────────────────────────────────────────────────
    mean_B = zeros(s.K)
    mean_S = zeros(s.K)
    for k in 1:s.K
        mean_B[k], mean_S[k] = voxel_mean_BS(s, k)
    end

    # ── Step 1: evolve active voxels ─────────────────────────────────────────
    to_deactivate_hi = Int[]
    to_deactivate_lo = Int[]

    for (k, p) in s.dists
        nb1, nb2 = ring_neighbors(k, s.K)
        in_B = m.D * (mean_B[nb1] + mean_B[nb2])
        in_S = m.D * (mean_S[nb1] + mean_S[nb2])

        Q = polarisome_generator(m, in_B, in_S, 2)
        p_new = expv(dt, Q, p; m=min(s.krylov_m, state_dim(m)))
        p_new ./= sum(p_new)
        s.dists[k] = p_new

        # Check convergence via phi_hi = P(nB≥nB_thresh AND nS≥nS_thresh)
        NS1 = m.n_S_max + 1
        phi_k = 0.0
        @inbounds for nb in s.nB_thresh:m.n_B_max, ns in s.nS_thresh:m.n_S_max
            phi_k += p_new[nb*NS1 + ns + 1]
        end
        if phi_k > s.ε_phi_hi
            push!(to_deactivate_hi, k)
        elseif phi_k < s.ε_phi_lo
            push!(to_deactivate_lo, k)
        end
    end

    for k in to_deactivate_hi
        delete!(s.dists, k);  push!(s.hi_voxels, k)
    end
    for k in to_deactivate_lo
        delete!(s.dists, k);  push!(s.lo_voxels, k)
    end

    # ── Step 2: activate lo voxels adjacent to hi voxels ────────────────────
    # Use excess flux above the lo-equilibrium baseline so that two lo neighbors
    # (baseline = D*2*nB_lo) never trigger activation; only a hi neighbor does.
    baseline_flux = m.D * 2 * s.nB_lo * dt
    to_activate = Int[]
    for k in s.lo_voxels
        nb1, nb2 = ring_neighbors(k, s.K)
        flux_in = m.D * (mean_B[nb1] + mean_B[nb2]) * dt
        flux_in - baseline_flux > s.ε_expand && push!(to_activate, k)
    end

    for k in to_activate
        delete!(s.lo_voxels, k)
        # Start from equil-lo distribution
        s.dists[k] = copy(s.equil_lo)
    end

    s.t += dt
    nothing
end

# ─── Snapshot helpers ─────────────────────────────────────────────────────────

function snap_means(s::PolarisomeFSP)
    B = zeros(s.K)
    S = zeros(s.K)
    for k in 1:s.K
        B[k], S[k] = voxel_mean_BS(s, k)
    end
    B, S
end

function snap_phi_hi(s::PolarisomeFSP)
    # P(nB > nB_lo + δ AND nS > nS_lo + δ) — fraction in hi basin
    m  = s.model
    NS1 = m.n_S_max + 1
    nB_thresh = round(Int, (s.nB_lo + s.nB_hi) / 2)
    nS_thresh = round(Int, (s.nS_lo + s.nS_hi) / 2)

    phi = zeros(s.K)
    for k in 1:s.K
        if k in s.hi_voxels
            phi[k] = 1.0
        elseif k in s.lo_voxels
            phi[k] = 0.0
        else
            p = s.dists[k]
            for nB in nB_thresh:m.n_B_max, nS in nS_thresh:m.n_S_max
                phi[k] += p[nB*NS1 + nS + 1]
            end
        end
    end
    phi
end

n_active(s::PolarisomeFSP) = length(s.dists)
