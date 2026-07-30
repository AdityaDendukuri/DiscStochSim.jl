# -----------------------------------------------------------------------------
# Boundary closures for a truncated CME: leak vs reflect vs ELSE reinjection.
#
# Every FSP truncates the state space; the open question is what to do with the
# probability flux that crosses the truncation boundary. Three closures on the
# interior J = {0..M}:
#
#   leak    : evolve the absorbing block A_JJ (classic FSP). Mass that exits J
#             is lost -- 1^T e^{tA_JJ} p < 1.
#   reflect : evolve the compressed generator A~_JJ = A_JJ + diag(1^T A_{J'J}).
#             Boundary outflow is folded back onto the diagonal: mass conserved,
#             but returned to the state it LEFT from (the wall).
#   ELSE(s) : A_JJ + A_{JJ'} (-A_{J'J'})^{-1} A_{J'J} over an exterior shell
#             J' = {M+1..M+s}. The fundamental matrix Z = -A_{J'J'}^{-1} routes
#             escaping mass through the true exterior and reinjects it where the
#             exterior actually delivers it. For s = full exterior this is the
#             exact Schur-complement marginalization of the stationary law.
#
# System: immigration-death-BURST, a NON-detailed-balance CME
#     ∅ --b--> X,   X --mu--> ∅,   Jc*X --c--> ∅   (X -> X - Jc)
# Jump sizes 1 (death) and Jc (burst) break detailed balance, so exterior
# excursions return to a SPREAD of interior states (a burst lands Jc below where
# it fires) rather than only at the wall. This is exactly where "reflect at the
# wall" mis-places mass and ELSE reinjection is needed. (In a reversible 1D
# birth-death chain, reflect already reproduces the exact stationary law, so the
# non-DB burst is essential to the comparison.)
#
# Ground truth: stationary law of the full chain on {0..Mbig}, conditioned on J.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, Plots
using DiscStochSim: build_generator, DiscreteStochasticSystem, StateSpace, add_state!
gr()

# --- model -----------------------------------------------------------------
Jc = 4                        # burst jump size
b, mu, c = 80.0, 1.0, 0.08    # immigration / death / burst rates
rn = @reaction_network begin
    b,  0 --> X
    mu, X --> 0
    c,  $Jc*X --> 0           # Jc-th order burst: X -> X - Jc
end
model = DiscreteStochasticSystem(rn)
rates = [b, mu, c]

chain_space(lo, hi) = begin
    sp = StateSpace{CartesianIndex, Float64}()
    for n in lo:hi; add_state!(sp, CartesianIndex(n), 0.0); end
    sp
end

# leading (Perron) eigenvector of a generator, normalized to a probability.
function stationary(A)
    F = eigen(Matrix(A))
    v = real.(F.vectors[:, argmax(real.(F.values))])
    v ./= sum(v)
    all(v .>= -1e-10) || (v = abs.(v); v ./= sum(v))
    v
end
tv(p, q) = 0.5 * sum(abs.(p .- q))

M, Mbig = 10, 90

# ground truth
Afull, = build_generator(chain_space(0, Mbig), model, rates, 0.0; absorbing=false)
pi_full = stationary(Afull)
pi_true = pi_full[1:M+1] ./ sum(pi_full[1:M+1])

# leak + reflect on J
p_leak    = stationary(build_generator(chain_space(0, M), model, rates, 0.0; absorbing=true)[1])
p_reflect = stationary(build_generator(chain_space(0, M), model, rates, 0.0; absorbing=false)[1])
tv_leak, tv_reflect = tv(p_leak, pi_true), tv(p_reflect, pi_true)

# ELSE reinjection over shells of increasing thickness
shells = 1:(Mbig - M)
p_else = Dict{Int,Vector{Float64}}()
tv_else = Float64[]
for s in shells
    A = Matrix(build_generator(chain_space(0, M + s), model, rates, 0.0; absorbing=true)[1])
    J, Pp = 1:(M+1), (M+2):(M+1+s)
    Gelse = A[J,J] - A[J,Pp] * (A[Pp,Pp] \ A[Pp,J])      # Schur complement
    ps = stationary(Gelse)
    p_else[s] = ps
    push!(tv_else, tv(ps, pi_true))
end

# --- report table ----------------------------------------------------------
println("immigration-death-burst (non-DB): Jc=$Jc, b=$b mu=$mu c=$c")
println("interior M=$M, exterior Mbig=$Mbig | mean(X)=",
        round(sum((0:Mbig).*pi_full); digits=2),
        "  P(X>M)=", round(1 - sum(pi_full[1:M+1]); sigdigits=3))
println("-"^52)
println("closure      TV to exact conditional")
println("leak         ", round(tv_leak;    sigdigits=4))
println("reflect      ", round(tv_reflect; sigdigits=4))
for s in shells
    (s <= 6 || s == shells[end]) && println("ELSE s=$s      ", round(tv_else[s]; sigdigits=4))
end

# --- figure 1: stationary laws ---------------------------------------------
xs = 0:M
f1 = plot(xs, pi_true; seriestype=:path, marker=:circle, ms=5, lw=2, color=:black,
          label="exact (conditional)", xlabel="X", ylabel="stationary probability",
          title="Closures at a deliberately TIGHT boundary  (M=$M, P(X>M)=$(round(1-sum(pi_full[1:M+1]);sigdigits=2)))",
          legend=:topleft, size=(760,460))
plot!(f1, xs, p_leak;      marker=:utriangle, ms=4, lw=1.5, ls=:dash, color=1, label="leak  (TV=$(round(tv_leak;sigdigits=2)))")
plot!(f1, xs, p_reflect;   marker=:diamond,   ms=4, lw=1.5, ls=:dash, color=2, label="reflect  (TV=$(round(tv_reflect;sigdigits=2)))")
plot!(f1, xs, p_else[2];   marker=:star5,     ms=5, lw=1.5, color=3,          label="ELSE shell s=2  (TV=$(round(tv_else[2];sigdigits=2)))")
mkpath(joinpath(@__DIR__, "..", "figs"))
savefig(f1, joinpath(@__DIR__, "..", "figs", "closure_distributions.png"))

# --- figure 2: TV vs shell thickness ---------------------------------------
f2 = plot(collect(shells), max.(tv_else, 1e-16); marker=:circle, ms=4, lw=2, color=3,
          yscale=:log10, xlabel="ELSE exterior shell thickness  s",
          ylabel="TV to exact conditional", legend=:topright, size=(760,460),
          title="ELSE reinjection converges to the exact marginal", label="ELSE(s)")
hline!(f2, [tv_leak];    ls=:dash, color=1, lw=2, label="leak")
hline!(f2, [tv_reflect]; ls=:dash, color=2, lw=2, label="reflect")
savefig(f2, joinpath(@__DIR__, "..", "figs", "closure_convergence.png"))

# --- figure 3: FAIR comparison at a FIXED quantity of interest -------------
# We want the stationary law on {0..M0}. Given s EXTRA states of budget, is it
# better to (a) spend them on an ELSE reinjection shell at the M0 boundary, or
# (b) simply ENLARGE the interior to M0+s (leak or reflect) and marginalize back
# to {0..M0}? All three use M0+1+s states and produce a law on {0..M0}, so this
# is a genuine apples-to-apples test of whether exact reinjection is worth it.
M0 = 10
pi0 = pi_full[1:M0+1] ./ sum(pi_full[1:M0+1])
restrict_tv(pfull) = tv(pfull[1:M0+1] ./ sum(pfull[1:M0+1]), pi0)

Ssweep = 0:2:24
else_tv = Float64[]; leakR_tv = Float64[]; reflR_tv = Float64[]
for s in Ssweep
    # (a) ELSE: interior M0 + shell s, exact reinjection at the M0 boundary
    if s == 0
        push!(else_tv, tv(stationary(build_generator(chain_space(0,M0), model, rates, 0.0; absorbing=true)[1]), pi0))
    else
        A = Matrix(build_generator(chain_space(0,M0+s), model, rates, 0.0; absorbing=true)[1])
        J, Pp = 1:(M0+1), (M0+2):(M0+1+s)
        G = A[J,J] - A[J,Pp]*(A[Pp,Pp]\A[Pp,J])
        push!(else_tv, tv(stationary(G), pi0))
    end
    # (b) enlarge interior to M0+s, then marginalize back to {0..M0}
    push!(leakR_tv, restrict_tv(stationary(build_generator(chain_space(0,M0+s), model, rates, 0.0; absorbing=true)[1])))
    push!(reflR_tv, restrict_tv(stationary(build_generator(chain_space(0,M0+s), model, rates, 0.0; absorbing=false)[1])))
end
budget = collect(Ssweep) .+ (M0 + 1)
f3 = plot(budget, max.(leakR_tv,1e-16); marker=:utriangle, ms=5, lw=2, color=1,
          yscale=:log10, xlabel="total states used (interior + shell)",
          ylabel="TV to exact law on {0..$M0}", legend=:topright, size=(780,470),
          title="Extra states: exact reinjection shell vs bigger interior",
          label="enlarge interior + leak, marginalize")
plot!(f3, budget, max.(reflR_tv,1e-16); marker=:diamond, ms=5, lw=2, color=2, label="enlarge interior + reflect, marginalize")
plot!(f3, budget, max.(else_tv,1e-16);  marker=:star5,   ms=6, lw=2, color=3, label="ELSE reinjection shell at M0")
savefig(f3, joinpath(@__DIR__, "..", "figs", "closure_fair.png"))

println("\nFAIR test -- law on {0..$M0}, s extra states spent as shell vs interior:")
println(" states | ELSE-shell | leak+enlarge | reflect+enlarge")
for (i, s) in enumerate(Ssweep)
    println("  $(budget[i])    | ", rpad(round(else_tv[i];sigdigits=3),10), " | ",
            rpad(round(leakR_tv[i];sigdigits=3),12), " | ", round(reflR_tv[i];sigdigits=3))
end
println("\nfigures: figs/closure_distributions.png, figs/closure_convergence.png, figs/closure_fair.png")
