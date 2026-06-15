# ─── Smoothed-aggregation algebraic multigrid (SA-AMG) ───────────────────────
#
# The single, front-and-center multigrid solver for the package.  It builds the
# multigrid hierarchy *purely from the matrix entries* of a sparse generator —
# no geometry, no per-system prolongation.  The same code therefore drives every
# example (annulus annihilation, Schlögl ring, 2-D birth–death) once each hands
# over its operator.
#
# Used as an implicit linear solver for backward Euler on a CME generator A:
# solve  (I − dt·A) xⁿ⁺¹ = xⁿ  with an SA-AMG V-cycle.
#
# Two design choices are forced by the structure of CME generators (which are
# nonsymmetric and *column*-diagonally-dominant — conservative, columns of A
# sum to ~0):
#
#   1. **Gauss–Seidel smoother**, not damped Jacobi.  Column-dominance means
#      ρ(I − ω D⁻¹S) can exceed 1, so Jacobi diverges; GS converges for any
#      M-matrix regardless of symmetry.
#   2. **Petrov–Galerkin** transfer: the prolongator P is smoothed with S and
#      the restriction R is built from a separate prolongator smoothed with Sᵀ,
#      so restriction tracks the *left* near-null mode.  Sc = R·S·P.
#
# The near-null-space is the plain constant (normalised piecewise-constant over
# each aggregate).  Validated to reproduce the bimodal quasi-stationary
# distribution at a bistable Schlögl interface without a QSD seed.

"""
    AMGLevel

One level of an SA-AMG hierarchy.  `P`/`R`/`DL` are empty on the coarsest level,
where `fact` holds a direct LU factorisation instead.
"""
struct AMGLevel
    S    :: SparseMatrixCSC{Float64,Int}
    P    :: SparseMatrixCSC{Float64,Int}   # prolongation   (empty at coarsest)
    R    :: SparseMatrixCSC{Float64,Int}   # restriction    (empty at coarsest)
    DL   :: SparseMatrixCSC{Float64,Int}   # lower-tri D+L of S, for Gauss-Seidel
    fact :: Union{Nothing,SparseArrays.UMFPACK.UmfpackLU{Float64,Int}}
end

# Strength of connection: i~j strong if |S_ij| ≥ θ·√(|S_ii·S_jj|).  Symmetrised.
function sa_strength(S::SparseMatrixCSC; θ::Float64=0.25)
    n = size(S, 1)
    d = abs.(diag(S))
    rows = rowvals(S); vals = nonzeros(S)
    I = Int[]; J = Int[]
    for j in 1:n, p in nzrange(S, j)
        i = rows[p]
        i == j && continue
        if abs(vals[p]) ≥ θ * sqrt(d[i] * d[j])
            push!(I, i); push!(J, j)
        end
    end
    G = sparse(I, J, trues(length(I)), n, n)
    G = G .| permutedims(G)
    dropzeros!(G)
    G
end

# Greedy 3-pass aggregation (standard SA, Vaněk et al.).
function sa_aggregate(G::SparseMatrixCSC{Bool,Int})
    n = size(G, 1)
    agg = zeros(Int, n)
    rows = rowvals(G)
    nbrs(i) = @view rows[nzrange(G, i)]
    nagg = 0
    for i in 1:n                                   # pass 1: seed from free neighbourhoods
        agg[i] != 0 && continue
        if all(agg[j] == 0 for j in nbrs(i))
            nagg += 1; agg[i] = nagg
            for j in nbrs(i); agg[j] = nagg; end
        end
    end
    for i in 1:n                                   # pass 2: attach to adjacent aggregate
        agg[i] != 0 && continue
        for j in nbrs(i)
            if agg[j] != 0; agg[i] = agg[j]; break; end
        end
    end
    for i in 1:n                                   # pass 3: leftovers self-aggregate
        if agg[i] == 0; nagg += 1; agg[i] = nagg; end
    end
    agg, nagg
end

# Tentative prolongator: normalised piecewise-constant (orthonormal columns).
function sa_tentative(agg::Vector{Int}, nagg::Int)
    n = length(agg)
    cnt = zeros(Int, nagg)
    for a in agg; cnt[a] += 1; end
    sparse(collect(1:n), agg, [1.0 / sqrt(cnt[agg[i]]) for i in 1:n], n, nagg)
end

# Jacobi-smoothed prolongator P = (I − ω D⁻¹ S) P₀.  Pass S for the right
# near-null mode, Sᵀ for the left (restriction) mode.
function sa_smooth_P(S::SparseMatrixCSC, P0::SparseMatrixCSC; ω::Float64=2/3)
    Dinv = Diagonal(1.0 ./ diag(S))
    P0 - ω * (Dinv * (S * P0))
end

"""
    build_sa_amg(S; strength_from=nothing, θ_str=0.25, ω=2/3,
                 max_coarse=60, max_levels=12, verbose=false) -> Vector{AMGLevel}

Build an SA-AMG hierarchy (finest first) for the sparse matrix `S`.  `S` is the
operator the V-cycle will solve against — e.g. the backward-Euler matrix
`I − dt·A` of a CME generator `A`.

`strength_from` is the matrix whose entries define the strength-of-connection
graph (and hence the aggregation).  Default `nothing` uses `S` itself.  For a
backward-Euler operator `I − dt·A` pass the generator `A`: the identity shift
makes `S`'s diagonal dominate, so strength measured on `S` sees nearly every
off-diagonal as weak and the finest level barely coarsens.  Measuring strength
on `A` recovers the true connection structure.  The strength operator is
coarsened by the same Galerkin product as `S`.
"""
function build_sa_amg(S0::SparseMatrixCSC{Float64,Int};
                      strength_from::Union{Nothing,SparseMatrixCSC{Float64,Int}}=nothing,
                      θ_str::Float64=0.25, ω::Float64=2/3,
                      max_coarse::Int=60, max_levels::Int=12, verbose::Bool=false)
    levels = AMGLevel[]
    S = S0
    Astr = strength_from === nothing ? S0 : strength_from
    for lvl in 1:max_levels
        n = size(S, 1)
        if n ≤ max_coarse || lvl == max_levels
            push!(levels, AMGLevel(S, spzeros(0,0), spzeros(0,0), spzeros(0,0), lu(S)))
            verbose && println("  SA-AMG level $lvl (coarsest): n=$n  direct solve")
            break
        end
        G          = sa_strength(Astr; θ=θ_str)
        agg, nagg  = sa_aggregate(G)
        P0         = sa_tentative(agg, nagg)
        P  = sa_smooth_P(S,           P0; ω=ω); dropzeros!(P)
        Pr = sa_smooth_P(sparse(S'),  P0; ω=ω); dropzeros!(Pr)
        R  = sparse(Pr')
        Sc = sparse(R * S * P)
        push!(levels, AMGLevel(S, P, R, tril(S), nothing))
        verbose && println("  SA-AMG level $lvl: n=$n → $nagg  (agg ratio ",
                           round(n/nagg; digits=2), ", P nnz/row ",
                           round(nnz(P)/n; digits=1), ")")
        S    = Sc
        Astr = sparse(R * Astr * P)
    end
    levels
end

# Forward Gauss–Seidel smoother:  x ← x + (D+L)⁻¹ (b − S x).
function amg_gauss_seidel!(x, b, lvl::AMGLevel, ν::Int)
    for _ in 1:ν
        x .+= lvl.DL \ (b .- lvl.S * x)
    end
end

"""
    amg_vcycle!(x, b, levels, l=1, ν1=2, ν2=2)

One SA-AMG V-cycle on level `l`, updating `x` in place toward solving `S x = b`.
"""
function amg_vcycle!(x, b, levels::Vector{AMGLevel}, l::Int=1, ν1::Int=2, ν2::Int=2)
    lvl = levels[l]
    if lvl.fact !== nothing
        x .= lvl.fact \ b
        return x
    end
    amg_gauss_seidel!(x, b, lvl, ν1)
    rc = lvl.R * (b .- lvl.S * x)
    ec = zeros(length(rc))
    amg_vcycle!(ec, rc, levels, l+1, ν1, ν2)
    x .+= lvl.P * ec
    amg_gauss_seidel!(x, b, lvl, ν2)
    x
end

"""
    amg_solve(levels, b; n_cycles=20, ν1=2, ν2=2, rtol=1e-10, x0=nothing,
              fallback=true, verbose=false) -> x

Solve `S x = b` (where `S = levels[1].S`) by repeated V-cycles until the relative
residual drops below `rtol` or `n_cycles` is reached.

If the V-cycle diverges (non-finite or increasing residual) or fails to reach
`rtol` within `n_cycles`, and `fallback` is true, solve directly with `S \\ b`
instead of returning a garbage iterate.  This guard is essential when chaining
backward-Euler steps: an unguarded diverging solve poisons every later step.
"""
function amg_solve(levels::Vector{AMGLevel}, b::AbstractVector;
                   n_cycles::Int=20, ν1::Int=2, ν2::Int=2, rtol::Float64=1e-10,
                   x0=nothing, fallback::Bool=true, verbose::Bool=false)
    S = levels[1].S
    x = x0 === nothing ? zeros(length(b)) : copy(x0)
    r0 = norm(b)
    r0 == 0 && return x
    converged = false
    prev = Inf
    for c in 1:n_cycles
        amg_vcycle!(x, b, levels, 1, ν1, ν2)
        rc = norm(b .- S*x)
        verbose && println("    SA-AMG cycle $c  rel.resid = ", rc/r0)
        if !isfinite(rc) || rc > 2 * prev          # divergence
            verbose && println("    SA-AMG diverged — falling back to direct solve")
            break
        end
        prev = rc
        if rc ≤ rtol*r0; converged = true; break; end
    end
    if !converged && fallback
        x .= S \ b
    end
    x
end

"""
    amg_be_step(A, p, dt; n_cycles=20, rtol=1e-10, kwargs...) -> p_next

One backward-Euler step of  dp/dt = A·p  via SA-AMG: builds the hierarchy for
`I − dt·A` and solves `(I − dt·A) p_next = p`.  Convenience wrapper; when taking
many steps with a fixed `A`/`dt`, call `build_sa_amg` once and reuse `amg_solve`.
"""
function amg_be_step(A::SparseMatrixCSC{Float64,Int}, p::AbstractVector, dt::Real;
                     n_cycles::Int=20, rtol::Float64=1e-10, kwargs...)
    n = size(A, 1)
    S = sparse(I, n, n) - dt * A
    levels = build_sa_amg(S; strength_from=A, kwargs...)
    amg_solve(levels, p; n_cycles, rtol, x0=p)
end
