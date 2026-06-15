# ─── Unified SA-AMG FSP stepper ──────────────────────────────────────────────
#
# The single time-stepping driver for every spatial example in the package.  It
# is deliberately geometry- and chemistry-agnostic: it operates only on a
# `StateSpace` of realised configurations and a `DiscreteStochasticSystem`
# (stoichiometry + propensities).  Each system — 2-D birth–death, Schlögl ring,
# Schlögl 2-D, annihilation — supplies its own `DiscreteStochasticSystem` via the
# existing builders (`build_graph_rdme_system`, `build_schlogl_graph_system`, …);
# this stepper is identical for all of them.
#
# One step is a θ-method implicit move on the CME generator A:
#
#     (I − θ·dt·A) pⁿ⁺¹ = (I + (1−θ)·dt·A) pⁿ
#
# solved with an SA-AMG V-cycle (build_sa_amg / amg_solve).  θ = 1 is backward
# Euler (monotone, no negative probabilities — the default); θ = ½ is
# Crank–Nicolson (2nd-order, more accurate but may produce small negatives that
# are clamped).  The generator is reassembled each step from the *realised*
# support, which is cheap (O(|S|·n_rxn)); the support is grown by `expand!` and
# trimmed by threshold pruning.
#
# There is no geometric coarsening anywhere: SA-AMG *is* the multigrid, built
# purely from the matrix entries of A.  This replaces all the bespoke per-system
# V-cycles.

"""
    UnifiedFSP{K}

State + configuration for the unified SA-AMG FSP stepper over `K`-component
`CartesianIndex` configurations.

Fields
------
- `sp`           : StateSpace of realised configurations + probabilities
- `sys`          : DiscreteStochasticSystem (stoichiometry + propensities)
- `bc`           : predicate `x -> Bool`, the FSP truncation (e.g. counts ≤ n_max)
- `dt`           : time step
- `θ`            : θ-method parameter (1.0 = backward Euler, 0.5 = Crank–Nicolson)
- `expand_depth` : rstep expansion depth applied before each solve
- `prune_tol`    : drop configurations with probability below this after each step
- `n_cycles`     : max SA-AMG V-cycles per solve
- `rtol`         : SA-AMG relative-residual tolerance
- `t`            : current time (mutable)
"""
mutable struct UnifiedFSP{K}
    sp           :: StateSpace{CartesianIndex{K}, Float64}
    sys          :: DiscreteStochasticSystem{CartesianIndex{K}}
    bc           :: Function
    dt           :: Float64
    θ            :: Float64
    expand_depth :: Int
    prune_tol    :: Float64
    n_cycles     :: Int
    rtol         :: Float64
    t            :: Float64
end

"""
    UnifiedFSP(sys, bc; dt, θ=1.0, expand_depth=1, prune_tol=1e-10,
               n_cycles=30, rtol=1e-10, ic=nothing) -> UnifiedFSP{K}

Construct a stepper for `sys`.  `K` is inferred from the system's element type.
`ic` is an optional `Vector{Pair{CartesianIndex{K},Float64}}` of initial
(configuration ⇒ probability); if omitted the all-zero configuration gets
probability 1.
"""
function UnifiedFSP(sys::DiscreteStochasticSystem{CartesianIndex{K}}, bc::Function;
                    dt::Float64, θ::Float64=1.0, expand_depth::Int=1,
                    prune_tol::Float64=1e-10, n_cycles::Int=30, rtol::Float64=1e-10,
                    ic=nothing) where K
    sp = StateSpace{CartesianIndex{K}, Float64}()
    if ic === nothing
        add_state!(sp, CartesianIndex(ntuple(_ -> 0, Val(K))), 1.0)
    else
        for (x, p) in ic; add_state!(sp, x, Float64(p)); end
    end
    UnifiedFSP{K}(sp, sys, bc, dt, θ, expand_depth, prune_tol, n_cycles, rtol, 0.0)
end

n_active(f::UnifiedFSP) = length(f.sp)

"""
    step!(f::UnifiedFSP) -> f

Advance one `θ`-method step of `dt`: expand the support, assemble the generator,
solve the implicit system with SA-AMG, clamp/renormalise, and prune.
"""
function step!(f::UnifiedFSP{K}) where K
    expand!(f.sp, f.sys, f.bc; depth=f.expand_depth)
    A, = build_generator(f.sp, f.sys, Float64[], f.t)
    n  = size(A, 1)
    p  = f.sp.probs

    if f.θ == 1.0
        p_new = amg_be_step(A, p, f.dt; n_cycles=f.n_cycles, rtol=f.rtol)
    else
        rhs    = p .+ ((1.0 - f.θ) * f.dt) .* (A * p)
        S      = sparse(I, n, n) - (f.θ * f.dt) * A
        levels = build_sa_amg(S; strength_from=A)
        p_new  = amg_solve(levels, rhs; n_cycles=f.n_cycles, rtol=f.rtol, x0=p)
    end

    p_new .= max.(0.0, p_new)          # FSP probabilities are non-negative
    f.sp.probs .= p_new
    renormalize!(f.sp)
    prune_threshold!(f.sp, f.prune_tol)
    f.t += f.dt
    f
end

"""
    mean_counts(f::UnifiedFSP) -> Vector{Float64}

Exact per-voxel FSP mean ⟨n_k⟩ = Σ n_k · p(n), from the realised distribution
(no mean-field).
"""
function mean_counts(f::UnifiedFSP{K}) where K
    μ = zeros(K)
    for (ci, p) in zip(f.sp.states, f.sp.probs)
        t = Tuple(ci)
        @inbounds for v in 1:K; μ[v] += t[v] * p; end
    end
    μ
end
