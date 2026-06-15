"""
graph_vcycle.jl

Basin-aware two-grid V-cycle for the implicit (backward-Euler) CME step on a
`GraphJointFSP` joint state space.

This is a genuine algebraic multigrid cycle whose grid-transfer operators are the
bistability-preserving basin lumping:

  • Restriction  R : fine → coarse   (aggregate each lumped voxel's counts into its
                                       basin label; a 0/1 summation matrix)
  • Prolongation P : coarse → fine   (smart prolong — redistribute coarse mass over
                                       the CURRENT fine conditional shape; `R·P = I`)
  • Coarse operator  A_c = R·A·P     (Galerkin — no hand-derived lumped generator)

The coarse DOF distinguishes basins, so the slow bistable-switching mode lives on
the coarse grid (cheap to solve); the fine smoother recovers the local within-basin
diffusive correlation the lumping drops. Solving `(I − dt·A) x = p_old` by this cycle
returns one backward-Euler step.
"""

@inline function _coarse_image(s::CartesianIndex{K}, Vset, bq::BasinQSD) where {K}
    CartesianIndex(ntuple(i -> i in Vset ? basin_of(bq, s[i]) - 1 : s[i], K))
end

"""
    build_basin_transfer(sp_f, p_ref, Vset, bq) -> (R, P, n_c, fc)

Sparse grid-transfer operators for lumping the voxels in `Vset`. `p_ref` is the
reference fine distribution defining the smart-prolongation conditional. Guarantees
`R·P = I_c` exactly. `fc[i]` is the coarse index of fine state `i`.
"""
function build_basin_transfer(sp_f::StateSpace{CartesianIndex{K}, T},
                              p_ref::Vector{Float64}, Vset, bq::BasinQSD) where {K, T}
    n_f  = length(sp_f)
    cidx = Dict{CartesianIndex{K}, Int}()
    cprob = Float64[]
    fc   = Vector{Int}(undef, n_f)
    for (i, s) in enumerate(sp_f.states)
        c = _coarse_image(s, Vset, bq)
        j = get(cidx, c, 0)
        if j == 0
            push!(cprob, 0.0); j = length(cprob); cidx[c] = j
        end
        fc[i] = j
        cprob[j] += p_ref[i]
    end
    n_c = length(cprob)
    R = sparse(fc, 1:n_f, ones(Float64, n_f), n_c, n_f)
    pv = Vector{Float64}(undef, n_f)
    @inbounds for i in 1:n_f
        cp = cprob[fc[i]]
        pv[i] = cp > 0 ? p_ref[i] / cp : 0.0
    end
    P = sparse(1:n_f, fc, pv, n_f, n_c)
    (R, P, n_c, fc)
end

"""
    basin_vcycle_be(A, p_old, sp_f, Vset, bq; dt, ncycle=20, ν=2,
                    smoother=:sgs, ω=0.6, tol=1e-10) -> (x, residuals, n_c)

Solve `(I − dt·A) x = p_old` (one backward-Euler CME step) with a basin-aware
two-grid V-cycle. `A` is the column-convention generator over `sp_f`. Direct coarse
solve; the smoother is selectable (see `student_project.md`):

  • `:jacobi` — weighted Jacobi, relaxation `ω` (pre and post)
  • `:gs`     — Gauss-Seidel, forward sweep both pre and post
  • `:sgs`    — symmetric Gauss-Seidel: forward pre-smooth, backward post-smooth

Returns the step `x`, the per-cycle residual norms (convergence history), and `n_c`.
"""
function basin_vcycle_be(A::AbstractMatrix, p_old::Vector{Float64},
                         sp_f::StateSpace{CartesianIndex{K}}, Vset, bq::BasinQSD;
                         dt::Float64, ncycle::Int = 20, ν::Int = 2,
                         smoother::Symbol = :sgs, ω::Float64 = 0.6,
                         tol::Float64 = 1e-10) where {K}
    n  = length(p_old)
    M  = sparse(1.0I, n, n) - dt * A
    b  = copy(p_old)
    R, P, n_c, _ = build_basin_transfer(sp_f, p_old, Vset, bq)
    A_c = R * A * P
    M_c = sparse(1.0I, n_c, n_c) - dt * A_c
    Mc  = lu(M_c)
    D   = Vector(diag(M))
    DL  = LowerTriangular(tril(M));  U  = triu(M, 1)    # forward GS split
    DU  = UpperTriangular(triu(M));  Lo = tril(M, -1)   # backward GS split
    jac!(y)  = (for _ in 1:ν; y .+= ω .* (b .- M * y) ./ D; end)
    fwd!(y)  = (for _ in 1:ν; y .= DL \ (b - U * y);  end)
    bwd!(y)  = (for _ in 1:ν; y .= DU \ (b - Lo * y); end)
    pre!  = smoother === :jacobi ? jac! : fwd!
    post! = smoother === :jacobi ? jac! : (smoother === :sgs ? bwd! : fwd!)
    x   = copy(b)
    res = Float64[]
    for _ in 1:ncycle
        pre!(x)
        r = b - M * x
        push!(res, norm(r))
        norm(r) < tol && break
        x .+= P * (Mc \ (R * r))
        post!(x)
    end
    (x, res, n_c)
end
