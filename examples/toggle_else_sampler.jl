# -----------------------------------------------------------------------------
# Simplest ELSE-accelerated sampler:  a Gillespie whose well sojourn is replaced
# by one exact ELSE draw.   [toggle switch]
#
#   PLAIN GILLESPIE : simulate every reaction until the trajectory first crosses
#                     the separatrix x=y (the switch).  Cost = every fast in-well
#                     flip, thousands of events per switch.
#   ELSE SAMPLER    : the well J={x≥y} is a fast subnetwork with the crossing
#                     absorbed, R = A_JJ.  Its first-passage (exit) time is a
#                     phase-type law with survival S(t)=1ᵀe^{tR}p0; the exit state
#                     is (cross-rate)⊙(Ze_init).  We eigendecompose R ONCE (shared
#                     across the whole batch) and then draw each switch in one
#                     shot -- exact time + exit state -- with no in-well stepping.
#
# Check: ELSE-sampled switch times match plain Gillespie (same distribution), at
# a fraction of the per-sample work.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, Statistics, Random, Plots
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
    for _ in 1:2000
        x = hillr(y, b, K, n) / d;  y = hillr(x, b, K, n) / d
    end
    (round(Int, x), round(Int, y))
end

# ---- plain Gillespie: every reaction, to first separatrix crossing ----------
function gillespie_switch(b, K, n, d, rng)
    x, y = xhigh_mode(b, K, n, d);  t, ev = 0.0, 0
    while true
        r1 = hillr(y,b,K,n); r2 = hillr(x,b,K,n); r3 = d*x; r4 = d*y
        tot = r1+r2+r3+r4
        t += randexp(rng)/tot; ev += 1
        xprev = x
        u = rand(rng)*tot
        if     u < r1;        x += 1
        elseif u < r1+r2;     y += 1
        elseif u < r1+r2+r3;  x -= 1
        else                  y -= 1
        end
        y > x && return (t=t, events=ev, x=xprev)   # crossed from (xprev,xprev)
    end
end

# ---- ELSE well operator: built ONCE, shared across the whole batch ----------
function well_else(G, b, K, n, d)
    sp = basin(G)
    A, _, _, exitrate = build_generator(sp, model, [b,K,n,d], 0.0; absorbing=true)
    R  = Matrix(A)
    s0 = sp.index[CartesianIndex(xhigh_mode(b,K,n,d)...)]
    # exit-TIME law:  S(t) = 1ᵀ e^{tR} p0 = Σ_k a_k e^{λ_k t}   via R = V Λ V⁻¹
    F    = eigen(R)
    Vinv = inv(F.vectors)
    a    = vec(sum(F.vectors, dims=1)) .* Vinv[:, s0]     # (1ᵀv_k)·(V⁻¹p0)_k
    # exit-STATE law:  (cross-rate) ⊙ (Z p0),  Z p0 = -R⁻¹ p0 = occupation times
    e0 = zeros(size(R,1)); e0[s0] = 1.0
    occ = -(R \ e0)
    (sp=sp, λ=F.values, a=a, exitrate=exitrate, occ=occ, mfpt=sum(occ))
end

Ssurv(w, t) = real(sum(w.a .* exp.(w.λ .* t)))          # survival at time t

# one exact ELSE switch: inverse-CDF for the time, weighted draw for the state.
function else_switch(w, rng)
    U = rand(rng)                                        # sample S(τ)=U, S decreasing
    hi = w.mfpt
    while Ssurv(w, hi) > U; hi *= 2; end                 # bracket [0, hi]
    lo = 0.0
    for _ in 1:60
        mid = 0.5*(lo+hi)
        Ssurv(w, mid) > U ? (lo = mid) : (hi = mid)      # bisection
    end
    τ = 0.5*(lo+hi)
    wts = w.exitrate .* w.occ                            # exit-state weights
    r = rand(rng)*sum(wts); c = 0.0; xexit = w.sp.states[end][1]
    for j in eachindex(wts)
        c += wts[j]
        if c >= r; xexit = w.sp.states[j][1]; break; end
    end
    (t=τ, x=xexit)
end

# ---- run both samplers ------------------------------------------------------
b, K, d, G, n = 20.0, 10.0, 1.0, 45, 3.0
M = 5000
rng = MersenneTwister(1)

w = well_else(G, b, K, n, d)                             # one shared solve
println("ELSE well: |J| = ", length(w.sp), ",  MFPT = ", round(w.mfpt; sigdigits=6))

g_t = Float64[]; g_x = Int[]; g_ev = Int[]
for _ in 1:M; s = gillespie_switch(b,K,n,d,rng); push!(g_t,s.t); push!(g_x,s.x); push!(g_ev,s.events); end

e_t = Float64[]; e_x = Int[]
for _ in 1:M; s = else_switch(w, rng); push!(e_t,s.t); push!(e_x,s.x); end

println("\nswitch time   Gillespie = ", round(mean(g_t);sigdigits=5), " ±", round(std(g_t)/sqrt(M);sigdigits=2),
        "   ELSE = ", round(mean(e_t);sigdigits=5), " ±", round(std(e_t)/sqrt(M);sigdigits=2))
println("work/switch   Gillespie = ", round(Int, mean(g_ev)), " reaction events",
        "     ELSE = 1 shared eigendecomp + ", M, " cheap draws")

# ---- figure: empirical CDFs overlay + exit-state histograms ----------------
mkpath(joinpath(@__DIR__, "..", "figs"))
ts = range(0, quantile(g_t, 0.995); length=200)
cdf(v, ts) = [count(<=(t), v)/length(v) for t in ts]
f1 = plot(ts, cdf(g_t, ts); lw=3, color=1, label="Gillespie (every reaction)",
          xlabel="first-switch time", ylabel="CDF",
          title="ELSE-sampled switch times = Gillespie (n=$n)", legend=:bottomright, size=(760,460))
plot!(f1, ts, cdf(e_t, ts); lw=2, ls=:dash, color=3, label="ELSE (one draw / switch)")
savefig(f1, joinpath(@__DIR__, "..", "figs", "toggle_sampler_cdf.png"))

xr = minimum(min.(minimum(g_x),minimum(e_x))):maximum(max.(maximum(g_x),maximum(e_x)))
hist(v, xr) = [count(==(x), v)/length(v) for x in xr]
f2 = bar(collect(xr), hist(g_x, xr); bar_width=0.6, color=1, alpha=0.5, label="Gillespie",
         xlabel="crossing coordinate x", ylabel="exit probability",
         title="ELSE exit state = Gillespie", size=(760,460))
scatter!(f2, collect(xr), hist(e_x, xr); marker=:circle, ms=5, color=3, label="ELSE")
savefig(f2, joinpath(@__DIR__, "..", "figs", "toggle_sampler_exit.png"))

println("\nfigures: figs/toggle_sampler_cdf.png, figs/toggle_sampler_exit.png")
