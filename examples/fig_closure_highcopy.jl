#=
fig_closure_highcopy.jl

The abstract's coarse-generator theorem-corollary (gap 2):
  "prolongation given by closed-form binomial/multinomial stationary distributions
   arising from fast intra-voxel diffusion. The resulting coarse-grid generator is
   exact for linear reaction networks, with closure error that vanishes in the
   high-copy-number limit."

Spatial aggregation: merge two adjacent voxels (fast diffusion between them) into a
super-voxel carrying the TOTAL count Y = n1 + n2. The fast intra-voxel diffusion
equilibrium splits Y between the two halves as Binomial(Y, 1/2) — this is the
closed-form prolongation P. Restriction R sums the (n1,n2) anti-diagonals onto Y.
The coarse generator is the Galerkin product  A_c = R A P  (R P = I).

We compare the coarse-model stationary distribution π(A_c) against the EXACT coarse
stationary R·π(A) (project the true fine 2-voxel stationary onto Y), as a function
of copy number (set by the volume Ω). Two networks:

  LINEAR    : ∅ →(β) A ,  A →(δ) ∅            — Y is autonomously Markov ⇒ A_c is
              EXACT at every Ω and every D (R A = A_Y R, so A_c = A_Y).
  NONLINEAR : ∅ →(β) A ,  2A →(c) ∅           — the death propensity c·n(n-1)/2 is
              nonlinear, so the binomial closure carries an O(1/⟨n⟩) error that
              vanishes as the copy number grows.

Copy-number scaling: β = β0·Ω, c = c0/Ω ⇒ ⟨n⟩ ∝ Ω (linear: δ fixed ⇒ ⟨n⟩ = β/δ ∝ Ω).

Output: paper/figures/fig_closure_highcopy.{png,pdf}
=#
using Printf, LinearAlgebra, SparseArrays, CairoMakie

const D   = 10.0     # fast intra-voxel diffusion (sets the binomial equilibrium)
const β0  = 1.0
const δ   = 0.2      # linear death
const c0  = 0.5      # nonlinear (annihilation) constant

# per-voxel propensities
birth_rate(Ω)            = β0 * Ω
death_lin(n)             = δ * n
death_nl(n, Ω)           = (c0 / Ω) * n * (n - 1) / 2

# choose a per-voxel truncation comfortably above the mean
function nmax_for(Ω, nonlinear)
    μ = nonlinear ? Ω * sqrt(2 * β0 / c0) : (β0 * Ω) / δ
    max(12, ceil(Int, 4.5 * μ) + 8)
end

# ── two-voxel fine generator over (n1,n2), col-major, voxel 1 fastest axis ───────
function build_2voxel(Ω, nonlinear)
    nmax = nmax_for(Ω, nonlinear)
    L = nmax + 1
    N = L * L
    idx(n1, n2) = n1 + n2 * L + 1
    death(n) = nonlinear ? death_nl(n, Ω) : death_lin(n)
    β = birth_rate(Ω)
    I = Int[]; J = Int[]; V = Float64[]
    push3!(i, j, r) = (push!(I, i); push!(J, j); push!(V, r))
    for n1 in 0:nmax, n2 in 0:nmax
        col = idx(n1, n2)
        # reactions in voxel 1
        n1 < nmax && (push3!(idx(n1+1, n2), col, β);          push3!(col, col, -β))
        n1 > 0    && (r = death(n1); push3!(idx(n1-1, n2), col, r); push3!(col, col, -r))
        # reactions in voxel 2
        n2 < nmax && (push3!(idx(n1, n2+1), col, β);          push3!(col, col, -β))
        n2 > 0    && (r = death(n2); push3!(idx(n1, n2-1), col, r); push3!(col, col, -r))
        # diffusion 1→2 and 2→1 (conserves Y)
        if n1 > 0 && n2 < nmax; r = D*n1; push3!(idx(n1-1, n2+1), col, r); push3!(col, col, -r); end
        if n2 > 0 && n1 < nmax; r = D*n2; push3!(idx(n1+1, n2-1), col, r); push3!(col, col, -r); end
    end
    sparse(I, J, V, N, N), nmax
end

# ── restriction R (sum anti-diagonals → Y) and binomial prolongation P ───────────
function build_RP(nmax)
    L = nmax + 1
    N = L * L
    Ymax = 2 * nmax
    nc = Ymax + 1
    idx(n1, n2) = n1 + n2 * L + 1
    # R: nc × N  (R[Y+1, fine] = 1 if n1+n2 == Y)
    Ri = Int[]; Rj = Int[]; Rv = Float64[]
    for n1 in 0:nmax, n2 in 0:nmax
        push!(Ri, (n1+n2)+1); push!(Rj, idx(n1,n2)); push!(Rv, 1.0)
    end
    R = sparse(Ri, Rj, Rv, nc, N)
    # P: N × nc  (binomial split of Y over the two halves, truncated+renormalized)
    logfact = [sum(log, 1:k; init=0.0) for k in 0:Ymax]
    binom(Y, k) = exp(logfact[Y+1] - logfact[k+1] - logfact[Y-k+1] - Y*log(2))
    Pi = Int[]; Pj = Int[]; Pv = Float64[]
    for Y in 0:Ymax
        # valid splits with n1,n2 ≤ nmax
        n1lo = max(0, Y - nmax); n1hi = min(Y, nmax)
        w = [binom(Y, n1) for n1 in n1lo:n1hi]; s = sum(w)
        for (k, n1) in enumerate(n1lo:n1hi)
            push!(Pi, idx(n1, Y-n1)); push!(Pj, Y+1); push!(Pv, w[k]/s)
        end
    end
    P = sparse(Pi, Pj, Pv, N, nc)
    R, P
end

# ── stationary distribution of a column-sum-zero generator G (G π = 0, Σπ = 1) ────
function stationary(G)
    n = size(G, 1)
    A = vcat(G[1:n-1, :], sparse(ones(1, n)))
    b = zeros(n); b[n] = 1.0
    π = A \ b
    π ./= sum(π)
    π
end

mean_Y(πc) = sum((0:length(πc)-1) .* πc)

# ── sweep copy number ────────────────────────────────────────────────────────────
Ωs = [1, 2, 4, 8, 16]
results = Dict{Bool, NamedTuple}()
for nonlinear in (false, true)
    μs = Float64[]; l1 = Float64[]; emean = Float64[]
    for Ω in Ωs
        A, nmax = build_2voxel(Ω, nonlinear)
        R, P = build_RP(nmax)
        π_fine = stationary(A)
        πc_true = R * π_fine; πc_true ./= sum(πc_true)   # exact coarse stationary
        A_c = R * A * P
        πc_mod = stationary(A_c)                          # binomial-closure coarse model
        push!(μs, mean_Y(πc_true))
        push!(l1, sum(abs.(πc_mod .- πc_true)))
        push!(emean, abs(mean_Y(πc_mod) - mean_Y(πc_true)) / mean_Y(πc_true))
        @printf("%-9s Ω=%-3d ⟨Y⟩=%6.2f  L1=%.3e  rel⟨Y⟩err=%.3e\n",
                nonlinear ? "nonlin" : "linear", Ω, mean_Y(πc_true), l1[end], emean[end])
    end
    results[nonlinear] = (μ=μs, l1=l1, emean=emean)
end

# distributions at low vs high copy for the nonlinear case (panel b)
function coarse_pair(Ω)
    A, nmax = build_2voxel(Ω, true); R, P = build_RP(nmax)
    π_fine = stationary(A); πc_true = R*π_fine; πc_true ./= sum(πc_true)
    πc_mod = stationary(R*A*P)
    πc_true, πc_mod
end
lowΩ, highΩ = 1, 16
t_lo, m_lo = coarse_pair(lowΩ)
t_hi, m_hi = coarse_pair(highΩ)

# ── figure ───────────────────────────────────────────────────────────────────────
set_theme!(theme_minimal())
fig = Figure(size=(1180, 380), fontsize=15)

# (a) closure error vs copy number
ax1 = Axis(fig[1,1], title="(a) closure error vs copy number",
           xlabel="copy number  ⟨Y⟩", ylabel="rel. error in stationary ⟨Y⟩",
           xscale=log10, yscale=log10)
rl = results[false]; rn = results[true]
lines!(ax1, rl.μ, max.(rl.emean, 1e-16), color=:steelblue, linewidth=2.5,
       label="linear network (exact)")
scatter!(ax1, rl.μ, max.(rl.emean, 1e-16), color=:steelblue, markersize=9)
lines!(ax1, rn.μ, rn.emean, color=:crimson, linewidth=2.5,
       label="nonlinear (binomial closure)")
scatter!(ax1, rn.μ, rn.emean, color=:crimson, markersize=9)
# reference 1/⟨n⟩ slope guide
g0 = rn.emean[1] * rn.μ[1]
lines!(ax1, rn.μ, g0 ./ rn.μ, color=:gray, linestyle=:dash, label="∝ 1/⟨Y⟩")
axislegend(ax1, position=:lb, framevisible=false, labelsize=10)

# (b) coarse stationary: exact vs binomial-closure model, low vs high copy
ax2 = Axis(fig[1,2], title="(b) coarse stationary: exact vs closure",
           xlabel="Y / ⟨Y⟩  (rescaled)", ylabel="P(Y) · ⟨Y⟩")
rescale!(p) = (μ = mean_Y(p); ((0:length(p)-1) ./ μ, p .* μ))
xl, yl = rescale!(t_lo); xlm, ylm = rescale!(m_lo)
xh, yh = rescale!(t_hi); xhm, yhm = rescale!(m_hi)
lines!(ax2, xl, yl, color=(:crimson,0.5), linewidth=2.0, label="exact, ⟨Y⟩≈$(round(Int,mean_Y(t_lo)))")
lines!(ax2, xlm, ylm, color=:crimson, linewidth=2.0, linestyle=:dash, label="closure, low copy")
lines!(ax2, xh, yh, color=(:navy,0.5), linewidth=2.0, label="exact, ⟨Y⟩≈$(round(Int,mean_Y(t_hi)))")
lines!(ax2, xhm, yhm, color=:navy, linewidth=2.0, linestyle=:dash, label="closure, high copy")
axislegend(ax2, position=:rt, framevisible=false, labelsize=9)

# (c) full-distribution L1 error vs copy number
ax3 = Axis(fig[1,3], title="(c) distribution error vs copy number",
           xlabel="copy number  ⟨Y⟩", ylabel="L1‖π_coarse − R·π_fine‖",
           xscale=log10, yscale=log10)
lines!(ax3, rl.μ, max.(rl.l1, 1e-16), color=:steelblue, linewidth=2.5, label="linear (exact)")
scatter!(ax3, rl.μ, max.(rl.l1, 1e-16), color=:steelblue, markersize=9)
lines!(ax3, rn.μ, rn.l1, color=:crimson, linewidth=2.5, label="nonlinear")
scatter!(ax3, rn.μ, rn.l1, color=:crimson, markersize=9)
axislegend(ax3, position=:lb, framevisible=false, labelsize=10)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
save(joinpath(outdir, "fig_closure_highcopy.png"), fig)
save(joinpath(outdir, "fig_closure_highcopy.pdf"), fig)
println("saved paper/figures/fig_closure_highcopy.{png,pdf}")
