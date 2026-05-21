"""
Row-TT Solver for 2D RDME
==========================

Represents the joint RDME distribution as a product of ROW tensor trains:

  P(n^{1,:}, ..., n^{Kx,:}) ≈ Π_k TT_k(n^{k,1}, ..., n^{k,Ky})

Each row k maintains an EXACT 1D joint distribution P(n^{k,1},...,n^{k,Ky})
as a CMETensorTrain of length Ky.  The row-direction (horizontal) diffusion
is handled exactly within each row TT.  The column-direction (vertical)
diffusion between rows is handled via mean-field: row k sees the EXACT
marginal means ⟨n^{k±1,j}⟩ from neighboring row TTs as effective birth/death.

Why this is better than per-voxel (SpatialFSP):
  - SpatialFSP: P ≈ Π_{k,j} p_{k,j}(n)  — zero inter-voxel correlations
  - Row-TT:     P ≈ Π_k TT_k              — exact within-row correlations

For birth-death RDME the inter-row correlations are small (correlation length
ξ ≈ √(D/k_d) ≈ 2 voxels), so mean-field between rows is highly accurate.

The row-TT gives access to:
  - Per-voxel means:       ⟨n^{k,j}⟩  (same as SpatialFSP, exact for linear)
  - Within-row marginals:  P(n^{k,j})  (full distribution, not just mean)
  - Within-row joint:      P(n^{k,j}, n^{k,j'})  (joint over any row pair)
  - Within-row mutual information, correlations, etc.

Cost per time step:  O(Kx × Ky × r² × n_max²)  — linear in both Kx and Ky
Storage:             O(Kx × Ky × r² × n_max)
"""

# ── Site-dependent stepping ────────────────────────────────────────────────────

"""
    tt_step_row!(tt, L_locals, D_row, dt; krylov_m, tol)

One time step for a 1D row TT with SITE-DEPENDENT local generators.

`L_locals[j]` is the (n_max+1)×(n_max+1) CME generator for site j,
including both reactions AND the mean-field inter-row column coupling.

Strang splitting: local/2 → row-diffusion → local/2 + compress + renorm.
"""
function tt_step_row!(tt::CMETensorTrain,
                       L_locals::Vector{<:AbstractMatrix{Float64}},
                       D_row::Float64,
                       dt::Float64;
                       krylov_m::Int  = 20,
                       tol::Float64   = 1e-14)
    M       = length(tt)
    L_diff  = _diffusion_generator(D_row, tt.n_max)

    for k in 1:M; tt_evolve_local!(tt, k, L_locals[k], dt/2; krylov_m); end
    for k in 1:M-1; tt_evolve_pair!(tt, k, L_diff, dt; tol, krylov_m); end
    for k in M:-1:1; tt_evolve_local!(tt, k, L_locals[k], dt/2; krylov_m); end

    tt_compress!(tt; tol)
    n = tt_norm(tt); n > 1e-14 && (tt.tensors[1] ./= n)
    tt
end

# ── RowTTSolver ────────────────────────────────────────────────────────────────

"""
    RowTTSolver

2D RDME solver: Kx rows, each a 1D TT of Ky sites.

Fields
------
- `row_tts`  : Kx CMETensorTrains, each length Ky
- `means`    : Kx × Ky matrix of current ⟨n^{k,j}⟩ (exact from each row TT)
- `Kx`, `Ky` : grid dimensions
- `k_b`, `k_d` : base birth/death rates
- `D_row`, `D_col` : horizontal / vertical diffusion rates D/dx²
- `n_max`    : FSP truncation per voxel
- `max_r`    : maximum TT bond dimension
- `t`        : current time
"""
mutable struct RowTTSolver
    row_tts  :: Vector{CMETensorTrain}
    means    :: Matrix{Float64}          # [Kx, Ky] exact row-TT means
    Kx       :: Int
    Ky       :: Int
    k_b      :: Float64
    k_d      :: Float64
    D_row    :: Float64
    D_col    :: Float64
    n_max    :: Int
    max_r    :: Int
    t        :: Float64
end

"""
    RowTTSolver(Kx, Ky, k_b, k_d, D; n_max, max_r)

Initialise with all voxels empty (n=0 everywhere).
"""
function RowTTSolver(Kx::Int, Ky::Int,
                      k_b::Float64, k_d::Float64,
                      D::Float64;
                      D_row::Float64 = D,
                      D_col::Float64 = D,
                      n_max::Int     = 12,
                      max_r::Int     = 4)
    row_tts = [CMETensorTrain(Ky, n_max, max_r) for _ in 1:Kx]
    means   = zeros(Float64, Kx, Ky)
    RowTTSolver(row_tts, means, Kx, Ky, k_b, k_d, D_row, D_col, n_max, max_r, 0.0)
end

"""
    set_ic!(solver, ki, ji, n0)

Place `n0` molecules at grid position (row ki, column ji).
Resets the TT for row ki with a delta IC at column ji = n0 molecules.
"""
function set_ic!(solver::RowTTSolver, ki::Int, ji::Int, n0::Int)
    Ky = solver.Ky; nv = solver.n_max + 1
    tt = solver.row_tts[ki]
    for j in 1:Ky
        r_l = size(tt.tensors[j], 1)
        r_r = size(tt.tensors[j], 3)
        tt.tensors[j] = zeros(r_l, nv, r_r)
        n = (j == ji) ? min(n0, solver.n_max) : 0
        tt.tensors[j][1, n+1, 1] = 1.0
    end
    solver.means[ki, :] = tt_means(tt)
end

# ── Per-row generator with inter-row coupling ─────────────────────────────────

"""
    _row_generators(solver, k) → Vector{Matrix{Float64}}

Build site-dependent CME generators for row k, incorporating the mean-field
column coupling to rows k-1 and k+1.

For voxel (k,j), the column diffusion contributes:
  - Gain from row k-1:  rate D_col × ⟨n^{k-1,j}⟩  (effective zeroth-order birth)
  - Loss to row k-1:    rate D_col × n^{k,j}          (effective first-order death)
  - Gain from row k+1:  rate D_col × ⟨n^{k+1,j}⟩
  - Loss to row k+1:    rate D_col × n^{k,j}

This modifies the effective birth/death rates:
  k_b_eff[j] = k_b + D_col × (μ_above[j] + μ_below[j])
  k_d_eff[j] = k_d + D_col × (n_col_neighbors)
"""
function _row_generators(solver::RowTTSolver, k::Int)
    Ky = solver.Ky; D_col = solver.D_col; n_max = solver.n_max
    k_b = solver.k_b; k_d = solver.k_d

    has_above = k > 1; has_below = k < solver.Kx
    n_col_nbrs = Float64(has_above + has_below)

    L_vec = Vector{Matrix{Float64}}(undef, Ky)
    for j in 1:Ky
        μ_above = has_above ? solver.means[k-1, j] : 0.0
        μ_below = has_below ? solver.means[k+1, j] : 0.0
        k_b_eff = k_b + D_col * (μ_above + μ_below)
        k_d_eff = k_d + D_col * n_col_nbrs
        L_vec[j] = birth_death_generator(k_b_eff, k_d_eff, n_max)
    end
    L_vec
end

# ── Time step ──────────────────────────────────────────────────────────────────

"""
    step!(solver, dt; krylov_m, tol)

Advance the RowTTSolver by one time step `dt`.

Order:
  1. Snapshot current means from all row TTs
  2. Build site-dependent generators for each row (includes inter-row coupling)
  3. Update each row TT with `tt_step_row!`
  4. Refresh the means matrix
"""
function step!(solver::RowTTSolver, dt::Float64;
               krylov_m::Int  = 20,
               tol::Float64   = 1e-14)
    # Step 1: refresh means (needed for inter-row coupling)
    for k in 1:solver.Kx
        solver.means[k, :] = tt_means(solver.row_tts[k])
    end

    # Steps 2 & 3: update each row TT
    for k in 1:solver.Kx
        L_row = _row_generators(solver, k)
        solver.row_tts[k] = tt_step_row!(solver.row_tts[k], L_row,
                                          solver.D_row, dt; krylov_m, tol)
    end

    # Step 4: refresh means for output
    for k in 1:solver.Kx
        solver.means[k, :] = tt_means(solver.row_tts[k])
    end

    solver.t += dt
    solver
end

# ── Diagnostics and observables ───────────────────────────────────────────────

"""
    mean_grid(solver) → Matrix{Float64}

Return the Kx × Ky matrix of exact FSP means ⟨n^{k,j}⟩.
These are computed from the full joint row distributions.
"""
mean_grid(solver::RowTTSolver) = copy(solver.means)

"""
    max_bond_grid(solver) → Matrix{Float64}

Return a Kx × (Ky-1) matrix of maximum bond dimensions per row.
High bond dimension = high within-row entanglement = strong correlations.
"""
function max_bond_grid(solver::RowTTSolver)
    bonds = zeros(Float64, solver.Kx, solver.Ky - 1)
    for k in 1:solver.Kx
        bonds[k, :] = Float64.(tt_bond_dims(solver.row_tts[k]))
    end
    bonds
end

"""
    row_joint_marginal(solver, k, j1, j2) → Matrix{Float64}

Compute the 2-voxel joint marginal P(n^{k,j1}, n^{k,j2}) from row k's TT.
Returns a (n_max+1) × (n_max+1) matrix, normalised to sum = 1.

This captures WITHIN-ROW CORRELATIONS — information per-voxel methods cannot provide.
"""
function row_joint_marginal(solver::RowTTSolver, k::Int, j1::Int, j2::Int)
    j1, j2 = min(j1, j2), max(j1, j2)
    tt   = solver.row_tts[k]
    nv   = solver.n_max + 1
    M    = length(tt)

    # Left environment: env_L[α] = Σ_{n¹,...,n^{j1-1}} P(...) (row vector of length r_left_j1)
    env_L = ones(Float64, 1, 1)   # [1, r_left=1] at the left boundary
    for s in 1:j1-1
        G  = tt.tensors[s]; r_l, _, r_r = size(G)
        S  = reshape(sum(G; dims=2), r_l, r_r)   # Σ_n G[α,n,β]
        env_L = env_L * S                          # accumulate: [1, r_r]
    end   # env_L: [1, r_left_j1]

    # Right environment: env_R[β] = Σ_{n^{j2+1},...,n^M} P(...) (column vector)
    env_R = ones(Float64, 1, 1)   # [r_right=1, 1] at the right boundary
    for s in M:-1:j2+1
        G  = tt.tensors[s]; r_l, _, r_r = size(G)
        S  = reshape(sum(G; dims=2), r_l, r_r)
        env_R = S * env_R                          # accumulate: [r_l, 1]
    end   # env_R: [r_right_j2, 1]

    # Middle environment: contract tensors j1+1 .. j2-1 (trace over n)
    G1 = tt.tensors[j1]; r_l1, _, r_r1 = size(G1)
    mid = Matrix{Float64}(I, r_r1, r_r1)
    for s in j1+1:j2-1
        G  = tt.tensors[s]; r_l, _, r_r = size(G)
        S  = reshape(sum(G; dims=2), r_l, r_r)
        mid = mid * S
    end   # mid: [r_r1, r_l2]

    G2 = tt.tensors[j2]; r_l2, _, r_r2 = size(G2)

    # P(n1,n2) = env_L × G1[:,n1,:] × mid × G2[:,n2,:] × env_R
    P = zeros(Float64, nv, nv)
    for n1 in 1:nv
        lG1M = env_L * G1[:, n1, :] * mid    # [1, r_l2]
        for n2 in 1:nv
            P[n1, n2] = (lG1M * G2[:, n2, :] * env_R)[1, 1]
        end
    end
    P .= max.(P, 0.0)
    s = sum(P); s > 0.0 && (P ./= s)
    P
end

"""
    tt_norm_all(solver) → Vector{Float64}

Return the norm of each row TT (should all be ≈1).
"""
tt_norm_all(solver::RowTTSolver) = [tt_norm(tt) for tt in solver.row_tts]

function Base.show(io::IO, s::RowTTSolver)
    r_max = maximum(tt_max_bond(tt) for tt in s.row_tts)
    print(io, "RowTTSolver($(s.Kx)×$(s.Ky), r_max=$r_max, t=$(round(s.t,digits=2)))")
end
