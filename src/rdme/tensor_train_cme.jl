"""
Tensor Train CME for the Spatial RDME
======================================

Represents the joint probability distribution P(n¹,...,nᴹ,t) over M voxels
as a Matrix Product State (MPS) / Tensor Train (TT):

  P(n¹,...,nᴹ) = G₁[1,n¹,α₁] × G₂[α₁,n²,α₂] × ... × Gₘ[α_{M-1},nᴹ,1]

Each Gₖ has shape [r_{k-1}, n_max+1, r_k].  Bond dimension r controls the
accuracy/cost tradeoff:

  r = 1 : product measure — independent per-voxel marginals (= SpatialFSP)
  r > 1 : captures inter-voxel correlations up to range ~ log(r)
  r → ∞ : exact joint distribution

Cost per time step: O(M × r² × n_max²) — LINEAR in M for fixed r.
Compare to direct FSP: O((n_max+1)^M) — exponential in M.

Time evolution via Strang splitting:
  1. exp(dt/2 × L_local):  per-voxel CME (reactions), no bond growth
  2. exp(dt   × L_diff):   nearest-neighbor diffusion + SVD truncation
  3. exp(dt/2 × L_local):  per-voxel CME again

For 2D grids (Kx × Ky): use alternating-direction operator splitting —
row-diffusion sweeps then column-diffusion sweeps.  Column neighbors are
non-adjacent in row-major ordering; bring them adjacent via pairwise SWAPs
(TEBD approach) at cost O(Ky) per column pair.

References:
  - Vidal (2003): TEBD algorithm for many-body quantum systems
  - White (1992): DMRG for 1D quantum chains
  - Haegeman et al. (2011): TDVP for MPS
"""

using LinearAlgebra, ExponentialUtilities

# ── TensorTrain type ──────────────────────────────────────────────────────────

"""
    CMETensorTrain

Tensor train representation of a joint RDME probability distribution.

Fields
------
- `tensors`  : Vector of M tensors, each of shape [r_left, n_max+1, r_right]
- `n_max`    : maximum molecule count per voxel (FSP truncation)
- `max_r`    : maximum bond dimension (controls accuracy)

Indexing: site k has shape tensors[k][α_left, n_k, α_right],
          where n_k ∈ 0:n_max and α_left/right ∈ 1:r.
"""
struct CMETensorTrain
    tensors :: Vector{Array{Float64, 3}}   # [r, n_max+1, r] per site
    n_max   :: Int
    max_r   :: Int
end

Base.length(tt::CMETensorTrain) = length(tt.tensors)

"""
    CMETensorTrain(M, n_max, max_r; ic)

Initialise a TT with `M` sites.  `ic` can be:
- `CartesianIndex{M}`: delta distribution at a specific state
- `Vector{Vector{Float64}}`: product of per-site distributions
"""
function CMETensorTrain(M::Int, n_max::Int, max_r::Int;
                         ic::Union{CartesianIndex, Vector{Vector{Float64}}, Nothing} = nothing)
    tensors = [zeros(Float64, 1, n_max+1, 1) for _ in 1:M]

    if ic === nothing
        # All molecules at n=0
        for k in 1:M; tensors[k][1, 1, 1] = 1.0; end

    elseif ic isa CartesianIndex
        # Delta distribution at ic
        ns = Tuple(ic)
        M == length(ns) || error("IC dimension $(length(ns)) ≠ M=$M")
        for k in 1:M
            n = ns[k]
            0 ≤ n ≤ n_max || error("IC molecule count $n at voxel $k exceeds n_max=$n_max")
            tensors[k][1, n+1, 1] = 1.0
        end

    elseif ic isa Vector{Vector{Float64}}
        # Product state from per-site marginals
        length(ic) == M || error("IC length $(length(ic)) ≠ M=$M")
        for k in 1:M
            p = ic[k]
            n = min(length(p), n_max+1)
            tensors[k][1, 1:n, 1] = p[1:n]
            tensors[k][1, :, 1] ./= sum(tensors[k][1, :, 1])
        end
    end

    CMETensorTrain(tensors, n_max, max_r)
end

# ── Norm and marginals ────────────────────────────────────────────────────────

"""
    tt_norm(tt) → Float64

Total probability mass Σ_{n¹,...,nᴹ} P(n¹,...,nᴹ).
Computed via left-to-right contraction of (Σ_n G_k) tensors.
"""
function tt_norm(tt::CMETensorTrain)
    # Accumulate left boundary vector: shape [r_right]
    v = ones(Float64, 1)   # r_left = 1 at left boundary
    for k in 1:length(tt)
        G = tt.tensors[k]   # [r_left, n, r_right]
        r_l, _, r_r = size(G)
        # Contract over n: Σ_n G[α,n,β]
        S = reshape(sum(G; dims=2), r_l, r_r)   # [r_left, r_right]
        v = v' * S                                # [1, r_right] → [r_right]
        v = vec(v)
    end
    v[1]
end

"""
    tt_marginals(tt) → Vector{Vector{Float64}}

Per-voxel marginal distributions: p_k(n) = Σ_{n^{-k}} P(n¹,...,nᴹ).
Returns a length-M vector, each element a length-(n_max+1) probability vector.
"""
function tt_marginals(tt::CMETensorTrain)
    M = length(tt)
    n = tt.n_max + 1

    # Left-boundary vectors: L[k][α] = Σ_{n¹,...,n^{k-1}} P up to site k-1
    L = Vector{Vector{Float64}}(undef, M+1)
    L[1] = ones(Float64, 1)
    for k in 1:M
        G = tt.tensors[k]   # [r_l, n, r_r]
        r_l, _, r_r = size(G)
        S = reshape(sum(G; dims=2), r_l, r_r)
        L[k+1] = S' * L[k]   # [r_r]
    end

    # Right-boundary vectors: R[k][α] = Σ_{n^{k+1},...,nᴹ} P from site k+1
    R = Vector{Vector{Float64}}(undef, M+1)
    R[M+1] = ones(Float64, 1)
    for k in M:-1:1
        G = tt.tensors[k]
        r_l, _, r_r = size(G)
        S = reshape(sum(G; dims=2), r_l, r_r)
        R[k] = S * R[k+1]   # [r_l]
    end

    # Marginal at each site: p_k(n) = L[k]ᵀ × G_k[:,n,:] × R[k+1]
    margs = Vector{Vector{Float64}}(undef, M)
    for k in 1:M
        G = tt.tensors[k]   # [r_l, n, r_r]
        r_l, nv, r_r = size(G)
        pk = zeros(nv)
        for ni in 1:nv
            pk[ni] = dot(L[k], G[:, ni, :] * R[k+1])
        end
        margs[k] = pk
    end
    margs
end

"""
    tt_means(tt) → Vector{Float64}

Exact FSP mean counts ⟨n_k⟩ = Σ_n n × p_k(n) for each voxel.
"""
function tt_means(tt::CMETensorTrain)
    margs = tt_marginals(tt)
    [sum((0:tt.n_max) .* m) for m in margs]
end

# ── SVD truncation ────────────────────────────────────────────────────────────

"""
    tt_compress!(tt; tol=1e-14)

Sweep right-to-left, performing QR/SVD at each bond and truncating
singular values below `tol` (while also respecting `tt.max_r`).
Left-normalises the result.
"""
function tt_compress!(tt::CMETensorTrain; tol::Float64 = 1e-14)
    M = length(tt)
    # Right-to-left sweep: SVD at each bond, truncate
    for k in M:-1:2
        G = tt.tensors[k]   # [r_l, n, r_r]
        r_l, nv, r_r = size(G)
        # Reshape to matrix [r_l, n×r_r]
        mat = reshape(G, r_l, nv * r_r)
        U, S, V = svd(mat)
        # Truncate
        r_keep = min(tt.max_r, sum(S .> tol * S[1]), length(S))
        r_keep = max(r_keep, 1)
        # Update G_k ← V'[:r_keep, :] reshaped
        tt.tensors[k] = reshape(V[:, 1:r_keep]', r_keep, nv, r_r)
        # Absorb U × diag(S[:r_keep]) into G_{k-1}
        G_prev = tt.tensors[k-1]   # [r_ll, n, r_l]
        r_ll, nv_prev, _ = size(G_prev)
        # Reshape G_prev to [r_ll × n, r_l], multiply with US
        US = U[:, 1:r_keep] * Diagonal(S[1:r_keep])   # [r_l, r_keep]
        G_prev_mat = reshape(G_prev, r_ll * nv_prev, r_l)
        new_mat = G_prev_mat * US   # [r_ll × n, r_keep]
        tt.tensors[k-1] = reshape(new_mat, r_ll, nv_prev, r_keep)
    end
    tt
end

# ── Local evolution (reactions — no bond dimension growth) ────────────────────

"""
    tt_evolve_local!(tt, k, L_local, dt)

Apply the per-voxel CME matrix exponential exp(L_local × dt) to site k.
`L_local` is the (n_max+1) × (n_max+1) generator matrix for one voxel.

Gₖ'[α,n',β] = Σ_n exp(L×dt)[n',n] × Gₖ[α,n,β]

Cost: O(r² × n_max²) — no bond dimension change.
"""
function tt_evolve_local!(tt::CMETensorTrain, k::Int,
                           L_local::AbstractMatrix{Float64}, dt::Float64;
                           krylov_m::Int = 20)
    G = tt.tensors[k]   # [r_l, n, r_r]
    r_l, nv, r_r = size(G)
    G_new = similar(G)
    # Apply exp(L dt) to each (α, β) slice along the n axis
    for α in 1:r_l, β in 1:r_r
        v = @view G[α, :, β]
        G_new[α, :, β] = expv(dt, L_local, v; m = min(krylov_m, nv))
    end
    tt.tensors[k] = G_new
    tt
end

"""
    tt_evolve_all_local!(tt, L_locals, dt)

Apply local CME evolution to ALL sites simultaneously.
`L_locals` is either a single matrix (same at every site) or a Vector of M matrices.
"""
function tt_evolve_all_local!(tt::CMETensorTrain,
                               L_locals::Union{Matrix{Float64}, Vector{<:AbstractMatrix{Float64}}},
                               dt::Float64; krylov_m::Int = 20)
    M = length(tt)
    for k in 1:M
        L = L_locals isa Matrix ? L_locals : L_locals[k]
        tt_evolve_local!(tt, k, L, dt; krylov_m)
    end
    tt
end

# ── Pair evolution (diffusion — grows then truncates bond dimension) ──────────

"""
    _diffusion_generator(D, n_max) → Matrix{Float64}

Build the (n_max+1)² × (n_max+1)² generator for one pair of voxels
coupled by diffusion rate D (= D/dx²):

  L_diff[n₁',n₂'; n₁,n₂] = D×n₁×δ(n₁'-n₁+1)×δ(n₂'-n₂-1)   [hop k→k']
                           + D×n₂×δ(n₁'-n₁-1)×δ(n₂'-n₂+1)   [hop k'→k]
                           - D×(n₁+n₂)×δ(n₁',n₁)×δ(n₂',n₂)  [diagonal]

Indexing: combined index (n₁,n₂) → n₁*(n_max+1) + n₂ (row-major).
"""
function _diffusion_generator(D::Float64, n_max::Int)
    nv = n_max + 1
    N  = nv * nv
    L  = zeros(Float64, N, N)
    for n1 in 0:n_max, n2 in 0:n_max
        col = n1 * nv + n2 + 1   # 1-indexed column
        # Diagonal: total outflow rate
        L[col, col] -= D * (n1 + n2)
        # k → k' hop: n1 decreases by 1, n2 increases by 1
        if n1 > 0 && n2 < n_max
            row = (n1-1) * nv + (n2+1) + 1
            L[row, col] += D * n1
        end
        # k' → k hop: n2 decreases by 1, n1 increases by 1
        if n2 > 0 && n1 < n_max
            row = (n1+1) * nv + (n2-1) + 1
            L[row, col] += D * n2
        end
    end
    L
end

"""
    tt_evolve_pair!(tt, k, L_diff_pair, dt, max_r; tol)

Apply the two-site diffusion evolution to the bond between sites k and k+1.

Steps:
  1. Contract:  M[α,n₁,n₂,β] = Σ_γ Gₖ[α,n₁,γ] × G_{k+1}[γ,n₂,β]
  2. Evolve:    apply exp(L_diff × dt) to the (n₁,n₂) pair index
  3. Truncate:  SVD M_new → Gₖ' × G'_{k+1}, keep r ≤ max_r

Cost: O(r_left × n_max² × r_right) for contraction,
      O(n_max⁴ × r²) for the expv,
      O(r × n_max² × r × min(r×n_max, r×n_max)) for SVD.

The SVD truncation is where accuracy is traded for tractability.
"""
function tt_evolve_pair!(tt::CMETensorTrain, k::Int,
                          L_diff_pair::AbstractMatrix{Float64},
                          dt::Float64;
                          tol::Float64 = 1e-14,
                          krylov_m::Int = 30)
    Gk  = tt.tensors[k]      # [r_l, nv, r_m]
    Gk1 = tt.tensors[k+1]    # [r_m, nv, r_r]
    r_l, nv, r_m = size(Gk)
    _,   _,  r_r = size(Gk1)

    # 1. Contract the two-site tensor M[α,n₁,n₂,β]
    # Efficient: reshape and matrix multiply
    Gk_mat  = reshape(Gk,  r_l * nv, r_m)            # [r_l×n, r_m]
    Gk1_mat = reshape(Gk1, r_m, nv * r_r)             # [r_m, n×r_r]
    M_mat   = Gk_mat * Gk1_mat                         # [r_l×n, n×r_r]
    M       = reshape(M_mat, r_l, nv, nv, r_r)         # [r_l, n₁, n₂, r_r]

    # 2. Apply exp(L_diff dt) to the (n₁, n₂) combined index
    # Reshape M to [r_l × r_r, n₁ × n₂] (pair index last)
    # and apply exp to each (α,β) "slice"
    nv2 = nv * nv
    M_reshaped = reshape(permutedims(M, (1,4,2,3)), r_l * r_r, nv2)  # [r_l×r_r, n²]
    M_new = similar(M_reshaped)
    for ab in 1:r_l * r_r
        M_new[ab, :] = expv(dt, L_diff_pair, M_reshaped[ab, :]; m = min(krylov_m, nv2))
    end
    # Reshape back to [r_l, r_r, n₁, n₂] then [r_l, n₁, n₂, r_r]
    M_new_4d = permutedims(reshape(M_new, r_l, r_r, nv, nv), (1, 3, 4, 2))

    # 3. SVD truncate the (k, k+1) bond
    # Reshape to [r_l × n₁, n₂ × r_r] for SVD
    M_svd = reshape(M_new_4d, r_l * nv, nv * r_r)
    U, S, V = svd(M_svd)

    # Truncate to max_r (also drop negligible singular values)
    r_max_effective = min(tt.max_r, length(S))
    r_keep = max(1, sum(S .> tol * max(S[1], 1e-100)))
    r_keep = min(r_keep, r_max_effective)

    # Reconstruct the two site tensors
    tt.tensors[k]   = reshape(U[:, 1:r_keep], r_l, nv, r_keep)
    Vt_scaled       = Diagonal(S[1:r_keep]) * V[:, 1:r_keep]'
    tt.tensors[k+1] = reshape(Vt_scaled, r_keep, nv, r_r)

    tt
end

# ── One time step: Strang splitting for 1D RDME ───────────────────────────────

"""
    tt_step_1d!(tt, L_local, D, dt; krylov_m, tol)

One time step for a 1D RDME chain using Strang splitting:

  1. exp(dt/2 × Σ_k L_k^local)  — per-voxel reactions (left-to-right)
  2. exp(dt × L_diffusion)       — nearest-neighbor diffusion + SVD truncation
  3. exp(dt/2 × Σ_k L_k^local)  — per-voxel reactions (right-to-left)

`L_local` is the (n_max+1) × (n_max+1) CME generator for reactions at one voxel.
`D` is the diffusion rate D/dx² (same for all bonds).
"""
function tt_step_1d!(tt::CMETensorTrain,
                      L_local::Matrix{Float64},
                      D::Float64,
                      dt::Float64;
                      krylov_m::Int  = 20,
                      tol::Float64   = 1e-14)
    M       = length(tt)
    L_diff  = _diffusion_generator(D, tt.n_max)

    # ── 1. Half-step local evolution (left → right) ───────────────────────────
    for k in 1:M
        tt_evolve_local!(tt, k, L_local, dt/2; krylov_m)
    end

    # ── 2. Full-step diffusion (left → right sweep, nearest-neighbor pairs) ───
    for k in 1:M-1
        tt_evolve_pair!(tt, k, L_diff, dt; tol, krylov_m)
    end

    # ── 3. Half-step local evolution (right → left) ───────────────────────────
    for k in M:-1:1
        tt_evolve_local!(tt, k, L_local, dt/2; krylov_m)
    end

    # Compress to enforce max_r globally, then renormalize.
    # SVD truncation loses probability mass (analogous to FSP probability loss).
    # Renormalizing keeps the distribution properly normalised; the truncation
    # error is controlled by max_r and tol.
    tt_compress!(tt; tol)
    n = tt_norm(tt)
    n > 1e-14 && (tt.tensors[1] ./= n)
    tt
end

# ── 2D RDME via operator splitting ────────────────────────────────────────────

"""
    tt_step_2d!(tt, L_local, D, dt, Kx, Ky; krylov_m, tol)

One time step for a 2D RDME grid (Kx × Ky voxels) using ADI-style splitting:

  1. exp(dt/4 × L_reactions)       — reactions
  2. exp(dt/2 × L_row_diffusion)   — row-direction diffusion
  3. exp(dt/4 × L_reactions)       — reactions
  4. exp(dt/2 × L_col_diffusion)   — column-direction diffusion
  5. exp(dt/4 × L_reactions)       — reactions
  6. exp(dt/2 × L_row_diffusion)   — row-direction diffusion (symmetry)
  7. exp(dt/4 × L_reactions)       — reactions

Voxels are ordered row-major: site k = (i-1)*Ky + j for row i, column j.

Row diffusion (j ↔ j+1 within row i): sites k = (i-1)*Ky+j and (i-1)*Ky+j+1
  → NEAREST NEIGHBOR in the 1D ordering → direct pair evolution.

Column diffusion (i ↔ i+1 within column j): sites k = (i-1)*Ky+j and i*Ky+j
  → Ky-1 sites apart in the 1D ordering.
  Handled via SWAP gates: bring k and k+Ky adjacent, evolve pair, SWAP back.
  Cost: O(Ky) SWAP operations per column pair.
"""
function tt_step_2d!(tt::CMETensorTrain,
                      L_local::Matrix{Float64},
                      D::Float64,
                      dt::Float64,
                      Kx::Int, Ky::Int;
                      krylov_m::Int  = 20,
                      tol::Float64   = 1e-14)
    M = Kx * Ky
    length(tt) == M || error("TT length $(length(tt)) ≠ Kx*Ky=$M")
    L_diff = _diffusion_generator(D, tt.n_max)

    # ── Helper: swap sites k and k+1 ──────────────────────────────────────────
    function swap_sites!(k::Int)
        Gk  = tt.tensors[k]      # [r_l, nv, r_m]
        Gk1 = tt.tensors[k+1]    # [r_m, nv, r_r]
        r_l, nv, r_m = size(Gk)
        _,   _,  r_r = size(Gk1)

        # Contract, permute n₁↔n₂, SVD
        M_mat = reshape(Gk, r_l*nv, r_m) * reshape(Gk1, r_m, nv*r_r)
        M4d   = reshape(M_mat, r_l, nv, nv, r_r)
        M4d_T = permutedims(M4d, (1, 3, 2, 4))   # swap n₁↔n₂
        M_svd = reshape(M4d_T, r_l*nv, nv*r_r)
        U, S, V = svd(M_svd)
        r_keep = min(tt.max_r, sum(S .> tol * max(S[1], 1e-100)))
        r_keep = max(1, r_keep)
        tt.tensors[k]   = reshape(U[:, 1:r_keep], r_l, nv, r_keep)
        tt.tensors[k+1] = reshape((Diagonal(S[1:r_keep]) * V[:,1:r_keep]'), r_keep, nv, r_r)
    end

    # ── Helper: row diffusion pass ────────────────────────────────────────────
    function row_diffusion!(dt_sub)
        for i in 1:Kx
            for j in 1:Ky-1
                k = (i-1)*Ky + j
                tt_evolve_pair!(tt, k, L_diff, dt_sub; tol, krylov_m)
            end
        end
    end

    # ── Helper: column diffusion via SWAP gates ───────────────────────────────
    # Bring column neighbor (i,j) and (i+1,j) adjacent by SWAPping through
    # the Ky-1 sites between them, evolve pair, then SWAP back.
    function col_diffusion!(dt_sub)
        for j in 1:Ky         # each column
            for i in 1:Kx-1   # each row interface in this column
                # Site of (i,j) in row-major = (i-1)*Ky+j (1-indexed → base (i-1)*Ky+j)
                ka = (i-1)*Ky + j   # (i,j)
                kb = i*Ky + j       # (i+1,j)
                # Bubble ka rightward to be adjacent to kb
                for s in ka:kb-2
                    swap_sites!(s)
                end
                # Now ka is at position kb-1, kb is at kb: evolve pair
                tt_evolve_pair!(tt, kb-1, L_diff, dt_sub; tol, krylov_m)
                # Bubble back
                for s in kb-2:-1:ka
                    swap_sites!(s)
                end
            end
        end
    end

    # ── ADI Strang splitting ──────────────────────────────────────────────────
    # Reaction quarter-steps interspersed with full diffusion half-steps
    for k in 1:M; tt_evolve_local!(tt, k, L_local, dt/4; krylov_m); end
    row_diffusion!(dt/2)
    for k in 1:M; tt_evolve_local!(tt, k, L_local, dt/4; krylov_m); end
    col_diffusion!(dt/2)
    for k in 1:M; tt_evolve_local!(tt, k, L_local, dt/4; krylov_m); end
    row_diffusion!(dt/2)
    for k in 1:M; tt_evolve_local!(tt, k, L_local, dt/4; krylov_m); end

    tt_compress!(tt; tol)
    tt
end

# ── Local CME generator builders ─────────────────────────────────────────────

"""
    birth_death_generator(k_b, k_d, n_max) → Matrix{Float64}

CME generator for ∅ →^{k_b} A →^{k_d} ∅ at a single voxel.
Shape: (n_max+1) × (n_max+1).
"""
function birth_death_generator(k_b::Float64, k_d::Float64, n_max::Int)
    nv = n_max + 1
    L  = zeros(Float64, nv, nv)
    for n in 0:n_max
        L[n+1, n+1] -= k_b                             # birth outflow
        n > 0 && (L[n+1, n+1] -= k_d * n)              # death outflow
        n < n_max && (L[n+2, n+1]  += k_b)             # birth inflow
        n > 0     && (L[n,   n+1]  += k_d * n)         # death inflow
    end
    L
end

"""
    schlogl_generator(model, n_max) → Matrix{Float64}

CME generator for the single-voxel Schlögl model.
"""
function schlogl_generator(model::SchloglModel1D, n_max::Int)
    c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
    nv = n_max + 1
    L  = zeros(Float64, nv, nv)
    for n in 0:n_max
        # Outflows
        α_birth  = c3
        α_death  = c4 * n
        α_auto   = n > 1 ? c1 * n * (n-1) / 2 : 0.0
        α_rev    = n > 2 ? c2 * n * (n-1) * (n-2) / 6 : 0.0
        L[n+1, n+1] -= (α_birth + α_death + α_auto + α_rev)
        # Inflows
        n < n_max && (L[n+2, n+1] += α_birth + α_auto)
        n > 0     && (L[n,   n+1] += α_death + α_rev)
    end
    L
end

# ── Diagnostics ───────────────────────────────────────────────────────────────

"""
    tt_bond_dims(tt) → Vector{Int}

Return the bond dimensions (r₁, r₂, ..., r_{M-1}) between adjacent sites.
"""
function tt_bond_dims(tt::CMETensorTrain)
    [size(tt.tensors[k], 3) for k in 1:length(tt)-1]
end

"""
    tt_max_bond(tt) → Int

Maximum bond dimension currently in the TT.
"""
tt_max_bond(tt::CMETensorTrain) = maximum(tt_bond_dims(tt); init=1)

function Base.show(io::IO, tt::CMETensorTrain)
    print(io, "CMETensorTrain(M=$(length(tt)), n_max=$(tt.n_max), " *
              "max_r=$(tt.max_r), r_current=$(tt_max_bond(tt)), " *
              "norm=$(round(tt_norm(tt), digits=4)))")
end
