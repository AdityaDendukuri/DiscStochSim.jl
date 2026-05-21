using LinearAlgebra: Tridiagonal

"""
Schlögl Model as AbstractSpatialRDMEModel
==========================================

Adds 2D spatial FSP support to SchloglModel1D (defined in rdme_model.jl).
No hardcoded bistability logic — the algorithm handles both basins automatically:
  - Activation:  flux-based (D×⟨n_k⟩ > ε_expand) — same for any model
  - Evolution:   per-voxel Schlögl CME via _single_voxel_generator (overridden here)
  - Equil:       two-basin criterion |μ - n_lo| < ε OR |μ - n_hi| < ε
  - Equil update: expv on single-voxel Schlögl CME — same call path as BranchingDeathRDME

Key: `_single_voxel_generator` dispatch on SchloglModel1D handles all four
Schlögl reactions (c1–c4) plus diffusion boundary — no linear birth-death assumption.
"""

# ── Interface methods for SpatialFSP ──────────────────────────────────────────

jump_rate(m::SchloglModel1D) = m.D   # D/dx² assuming dx=1

# Downward rate for _single_voxel_generator_range dispatch
# (degradation c4·n + reverse autocatalysis c2·n(n-1)(n-2)/6 + diffusion out)
function _voxel_down_rate(m::SchloglModel1D, n::Int, out_coeff::Float64)
    m.c4*n + (n >= 3 ? m.c2*n*(n-1)*(n-2)/6 : 0.0) + out_coeff*n
end

# "Birth-like" = all upward propensities at count n
local_birth_rate(m::SchloglModel1D, n::Int) =
    m.c3 + (n >= 2 ? m.c1 * n * (n-1) / 2 : 0.0)

# For multi-voxel coarse block (n_sites voxels, total count nm)
local_total_birth_rate(m::SchloglModel1D, nm::Int, n_sites::Int=1) =
    n_sites * m.c3 + (nm >= 2 ? m.c1 * nm * (nm-1) / 2 : 0.0)

has_finite_local_ss(::SchloglModel1D) = false   # bistable: two stable states
ss_mean(::SchloglModel1D) = error("Schlögl is bistable — no unique steady-state mean")

# ── Single-voxel CME generator: full nonlinear Schlögl propensities ──────────

"""
    _single_voxel_generator(model::SchloglModel1D, in_rate, out_coeff, n_max)

Override for Schlögl: includes all four reactions (c1–c4) plus mean-field
diffusion boundary.  Replaces the linear birth-death formula used for
BirthDeathRDME / BranchingDeathRDME.

Called by the equil update (expv) and the equil detection (stationary check).
"""
function _single_voxel_generator(model::SchloglModel1D,
                                   in_rate  :: Float64,
                                   out_coeff :: Float64,
                                   n_max    :: Int;
                                   n_sites  :: Int = 1)
    c1,c2,c3,c4 = model.c1, model.c2, model.c3, model.c4
    nv = n_max + 1
    Q  = zeros(Float64, nv, nv)
    for col in 1:nv
        n = col - 1    # molecule count (0-indexed)

        # ── Upward: production + autocatalysis + diffusion in ─────────────────
        # 2X→3X autocatalysis propensity = c1·n(n-1)/2
        up = in_rate + n_sites*c3 + (n >= 2 ? c1*n*(n-1)/2 : 0.0)
        if col < nv && up > 0
            Q[col+1, col] += up;   Q[col, col] -= up
        end

        # ── Downward: degradation + reverse autocatalysis + diffusion out ──────
        # 3X→2X propensity = c2·n(n-1)(n-2)/6
        down = c4*n + (n >= 3 ? c2*n*(n-1)*(n-2)/6 : 0.0) + out_coeff*n
        if col > 1 && down > 0
            Q[col-1, col] += down;  Q[col, col] -= down
        end
    end
    Q
end

"""
    _single_voxel_generator_range(model::SchloglModel1D, ...)

Override returning a sparse `Tridiagonal` instead of a dense matrix.
The Schlögl CME is a birth-death chain — only O(n) non-zeros vs O(n²) dense.
With `expv` from ExponentialUtilities, O(n) mv products replace O(n²),
giving ~100× speedup for n≈100 states at large diffusion rates.
"""
function _single_voxel_generator_range(model::SchloglModel1D,
                                        in_rate::Float64, out_coeff::Float64,
                                        n_lo::Int, n_hi::Int;
                                        n_sites::Int = 1)
    nv = n_hi - n_lo + 1
    dl = zeros(Float64, nv - 1)   # sub-diagonal:   downward rates (col → col-1)
    d  = zeros(Float64, nv)        # main diagonal:  -sum of exit rates
    du = zeros(Float64, nv - 1)   # super-diagonal: upward rates   (col → col+1)
    for col in 1:nv
        n    = col - 1 + n_lo
        up   = in_rate + local_total_birth_rate(model, n, n_sites)
        down = _voxel_down_rate(model, n, out_coeff)
        if col < nv && up > 0
            dl[col]  = up;  d[col] -= up
        end
        if col > 1 && down > 0
            du[col-1] = down;  d[col] -= down
        end
    end
    LinearAlgebra.Tridiagonal(dl, d, du)
end
