# -----------------------------------------------------------------------------
# ELSE escape law via the RESOLVENT (Laplace transform) — scalable sampler.
#
# The eigendecomposition of A_JJ doesn't scale (dense). Mark's route: the escape
# law lives in the resolvent (sI - A_JJ)⁻¹, the Laplace transform of the semigroup
#     (sI - A_JJ)⁻¹ p0 = ∫₀^∞ e^{-st} e^{tA_JJ} p0 dt.
# The escape-density transform is  f̂(s) = Δᵀ (sI - A_JJ)⁻¹ p0.  We evaluate it at
# a FIXED set of Bromwich nodes s_k = a + i kπ/T (one sparse solve each, shared
# across all times and all samples), invert by the Fourier-series method to get
# f(t) on a grid, integrate to the CDF, and sample by inverse-CDF.
#
# Why the density (not the survival): started deep in the well, Δᵀp0 = 0 and the
# first many derivatives vanish, so f is smooth with f(0)=0 → f̂(s) decays fast →
# a few dozen nodes suffice. (Inverting the survival Γ_J has a 1/s tail → slow.)
#
# Validation on the toy toggle: f(t) vs the exact eigendecomposition, and the
# sampled switch-time CDF vs plain Gillespie.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, SparseArrays, Statistics, Random, Plots
using DiscStochSim: build_generator, DiscreteStochasticSystem, StateSpace, add_state!
gr()

rn = @reaction_network begin
    hillr(Y, b, K, n), ∅ --> X
    hillr(X, b, K, n), ∅ --> Y
    d, X --> ∅
    d, Y --> ∅
end
model = DiscreteStochasticSystem(rn)
hillr(z, v, K, n) = v * K^n / (K^n + z^n)

basin(G) = begin
    sp = StateSpace{CartesianIndex, Float64}()
    for x in 0:G, y in 0:G
        x >= y && add_state!(sp, CartesianIndex(x, y), 0.0)
    end
    sp
end
function xhigh_mode(b, K, n, d)
    x, y = float(b), 0.0
    for _ in 1:2000; x = hillr(y,b,K,n)/d; y = hillr(x,b,K,n)/d; end
    (round(Int, x), round(Int, y))
end

function subnetwork(G, b, K, n, d)
    sp = basin(G)
    A, _, _, Δ = build_generator(sp, model, [b,K,n,d], 0.0; absorbing=true)  # SPARSE
    s0 = sp.index[CartesianIndex(xhigh_mode(b,K,n,d)...)]
    p0 = zeros(length(sp)); p0[s0] = 1.0
    (sp=sp, A=A, Δ=Δ, p0=p0)
end

# ---- Laplace / resolvent sampler -------------------------------------------
# f̂(s) = Δᵀ (sI - A)⁻¹ p0  at fixed nodes s_k = a + i kπ/T,  k = 0..N.
function resolvent_nodes(sub; T, a, N)
    n = length(sub.p0)
    Ip0 = ComplexF64.(sub.p0)
    fhat = Vector{ComplexF64}(undef, N+1)
    for k in 0:N
        s = a + im*(k*π/T)
        v = (s*I - sub.A) \ Ip0            # one sparse (complex) solve — shared
        fhat[k+1] = dot(sub.Δ, v)          # Δᵀ v  (Δ real)
    end
    fhat                                    # N+1 numbers = the whole escape law
end

# Fourier-series inversion:  f(t) ≈ (e^{at}/T)[ ½Re f̂(a) + Σ Re(f̂_k e^{i kπt/T}) ].
function density(fhat, t; T, a)
    N = length(fhat) - 1
    s = 0.5*real(fhat[1])
    for k in 1:N
        s += real(fhat[k+1] * cis(k*π*t/T))
    end
    max(0.0, (exp(a*t)/T) * s)              # clamp tiny negative ripple
end

# tabulate density → CDF → inverse-CDF sampler (nodes fhat precomputed once)
function laplace_sampler(fhat; T, a, tmax, ngrid=4000)
    tg = range(0, tmax; length=ngrid)                       # tabulate over mass region
    fg = [density(fhat, t; T=T, a=a) for t in tg]
    dt = step(tg)
    Fg = cumsum(0.5 .* (fg .+ [0.0; fg[1:end-1]]) .* dt)     # cumulative ∫f
    Fg ./= Fg[end]                                          # normalize to 1
    (tg=collect(tg), fg=fg, Fg=Fg)
end

sample_time(S, rng) = S.tg[searchsortedfirst(S.Fg, rand(rng))]

# ---- exact eigendecomposition reference (toy only) -------------------------
function eig_density(sub)
    F = eigen(Matrix(sub.A))
    a = vec(sum(F.vectors, dims=1)) .* (F.vectors \ sub.p0)  # Γ_J = Σ a_k e^{λ_k t}
    t -> real(-sum(a .* F.values .* exp.(F.values .* t)))    # f = -Γ_J'
end

# ---- plain Gillespie reference ---------------------------------------------
function gillespie_switch(b, K, n, d, rng)
    x, y = xhigh_mode(b,K,n,d); t = 0.0
    while true
        r1=hillr(y,b,K,n); r2=hillr(x,b,K,n); r3=d*x; r4=d*y; tot=r1+r2+r3+r4
        t += randexp(rng)/tot
        u = rand(rng)*tot
        if     u<r1;      x+=1
        elseif u<r1+r2;   y+=1
        elseif u<r1+r2+r3;x-=1
        else              y-=1 end
        y > x && return t
    end
end

# ---- run --------------------------------------------------------------------
b, K, d, G, n = 20.0, 10.0, 1.0, 45, 3.0
sub = subnetwork(G, b, K, n, d)
mean_exact = sum(-(sub.A \ sub.p0))                          # s=0 resolvent = Z p0
println("|J| = ", length(sub.p0), "   mean escape (s=0 solve) = ", round(mean_exact; sigdigits=6))

# fixed Bromwich nodes: large T + small a keeps e^{at} amplification mild over [0,tmax]
Tlap, alap, Nlap, tmax = 500.0, 0.015, 500, 220.0
fhat = resolvent_nodes(sub; T=Tlap, a=alap, N=Nlap)          # Nlap+1 solves, done ONCE
S = laplace_sampler(fhat; T=Tlap, a=alap, tmax=tmax)
println("Laplace sampler: ", Nlap+1, " resolvent solves (shared across all t & samples)")

fexact = eig_density(sub)
M = 5000
rng = MersenneTwister(1); g_t = [gillespie_switch(b,K,n,d,rng) for _ in 1:M]
rng2 = MersenneTwister(2); l_t = [sample_time(S, rng2) for _ in 1:M]
println("mean switch time  Gillespie = ", round(mean(g_t);sigdigits=5),
        "   Laplace-sampled = ", round(mean(l_t);sigdigits=5),
        "   (exact = ", round(mean_exact;sigdigits=5), ")")

# ---- figures ---------------------------------------------------------------
mkpath(joinpath(@__DIR__, "..", "figs"))
tt = range(0, 160; length=300)
f1 = plot(tt, [fexact(t) for t in tt]; lw=4, color=1, label="exact (eigendecomp)",
          xlabel="escape time t", ylabel="density f(t)",
          title="ELSE escape density via resolvent ($(Nlap+1) solves)", size=(760,450))
plot!(f1, tt, [density(fhat, t; T=Tlap, a=alap) for t in tt];
      lw=2, ls=:dash, color=3, label="Laplace (resolvent)")
savefig(f1, joinpath(@__DIR__, "..", "figs", "toggle_laplace_density.png"))

ts = range(0, quantile(g_t,0.995); length=200)
cdf(v) = [count(<=(t),v)/length(v) for t in ts]
f2 = plot(ts, cdf(g_t); lw=3, color=1, label="Gillespie", xlabel="switch time", ylabel="CDF",
          title="Laplace-sampled switch times = Gillespie", legend=:bottomright, size=(760,450))
plot!(f2, ts, cdf(l_t); lw=2, ls=:dash, color=3, label="Laplace sampler")
savefig(f2, joinpath(@__DIR__, "..", "figs", "toggle_laplace_cdf.png"))

println("\nfigures: figs/toggle_laplace_density.png, figs/toggle_laplace_cdf.png")
