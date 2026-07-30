# -----------------------------------------------------------------------------
# ELSE on a multistable CME:  genetic toggle switch.   [scaling side of §4.3]
#
#   X production repressed by Y,  Y production repressed by X  (mutual Hill
#   repression), plus first-order degradation.  Two metastable basins:
#   X-high/Y-low and X-low/Y-high, separated by the diagonal x=y.  Within a
#   basin the 2-D fluctuations equilibrate fast; a rare noise excursion crosses
#   the separatrix and commits to the other basin.
#
# This lifts the 1-D fast/slow chain (fast_slow_else.jl) to a genuine 2-D basin
# and a genuine escape SURFACE, so ELSE's identity is exercised on both of its
# outputs:
#
#   ELSE : the fundamental matrix Z = -R^{-1} of the basin R = A_JJ (crossing
#          absorbed) gives in ONE linear solve
#            - the EXACT first-crossing MFPT      1^T Z e_init,
#            - the EXACT exit distribution        (cross-rate) ⊙ (Z e_init),
#          i.e. WHEN and WHERE you leave the basin.  Cost is one |J|-state solve,
#          independent of how rare the crossing is.
#   SSA  : exact first-passage reference.  It must simulate every fast in-basin
#          fluctuation, so event count to one crossing grows as the barrier
#          deepens; ELSE stays at one solve.
#
# NOTE (transient vs recurrent).  Unlike the fast/slow example -- where C is a
# true absorbing sink and the FSP survival integral ∫P(C=0,t)dt is the exact
# MFPT -- the toggle is RECURRENT: on the closed CME P(x≥y) relaxes to its
# stationary plateau, not to 0, because trajectories re-cross.  The separatrix
# crossing is a first-passage on a recurrent chain, so the exact reference is
# SSA first-passage (equivalently, the FSP of the basin with the crossing made
# absorbing -- which is exactly the operator R that ELSE inverts, so it cannot
# be an independent check).  This is the transient/recurrent split of the
# resolvent theory (else_theory.md).
#
# Hard axis: the Hill cooperativity n deepens the barrier (rarer crossing) while
# the grid and basin stay FIXED, so ELSE's work is flat exactly where SSA's grows.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, Plots, Statistics, Random
using DiscStochSim: build_generator, DiscreteStochasticSystem, StateSpace, add_state!
gr()

# X <- repressed by Y,  Y <- repressed by X;  hillr(z,v,K,n) = v K^n/(K^n+z^n).
rn = @reaction_network begin
    hillr(Y, b, K, n), ∅ --> X
    hillr(X, b, K, n), ∅ --> Y
    d, X --> ∅
    d, Y --> ∅
end
model = DiscreteStochasticSystem(rn)           # rates order: [b, K, n, d]

hillr(z, v, K, n) = v * K^n / (K^n + z^n)

# X-high basin J = {x ≥ y} inside the grid 0:G.  Crossing to x<y = escape.
basin(G) = begin
    sp = StateSpace{CartesianIndex, Float64}()
    for x in 0:G, y in 0:G
        x >= y && add_state!(sp, CartesianIndex(x, y), 0.0)
    end
    sp
end

# Deterministic X-high fixed point (mode), used as the initial state.
function xhigh_mode(b, K, n, d)
    x, y = float(b), 0.0
    for _ in 1:2000
        x = hillr(y, b, K, n) / d
        y = hillr(x, b, K, n) / d
    end
    (round(Int, x), round(Int, y))
end

# ELSE: exact first-crossing MFPT + exit distribution via the fundamental matrix.
# Returns (mfpt, exit_x, exit_p): the exit distribution indexed by the separatrix
# state (x=y) the walk crosses FROM.
function else_escape(G, b, K, n, d)
    sp = basin(G)
    rates = [b, K, n, d]
    R = build_generator(sp, model, rates, 0.0; absorbing=true)[1]
    x0 = xhigh_mode(b, K, n, d)
    e_init = zeros(length(sp)); e_init[sp.index[CartesianIndex(x0...)]] = 1.0
    occ = -(Matrix(R) \ e_init)                # Z e_init = occupation times
    mfpt = sum(occ)                            # 1^T Z e_init
    # crossing happens only from the diagonal x=y (Y-birth y→y+1 or X-death x→x-1)
    xs, ws = Int[], Float64[]
    for (j, s) in enumerate(sp.states)
        x, y = s[1], s[2]
        x == y || continue
        rate = hillr(x, b, K, n) + d * x       # prod_Y (crosses up) + deg_X (crosses left)
        rate > 0 || continue
        push!(xs, x); push!(ws, rate * occ[j])
    end
    p = sortperm(xs)
    (mfpt=mfpt, exit_x=xs[p], exit_p=ws[p] ./ sum(ws))
end

# SSA: Gillespie to first separatrix crossing.
# Returns (time, n_events, m) where m = separatrix state (x=y) crossed FROM.
function ssa_cross(G, b, K, n, d, rng)
    x0 = xhigh_mode(b, K, n, d)
    x, y, t, ev = x0[1], x0[2], 0.0, 0
    while true
        r1 = hillr(y, b, K, n)          # ∅->X
        r2 = hillr(x, b, K, n)          # ∅->Y
        r3 = d * x                      # X->∅
        r4 = d * y                      # Y->∅
        tot = r1 + r2 + r3 + r4
        t += randexp(rng) / tot; ev += 1
        xprev = x                       # pre-step diagonal value (x==y at a crossing step)
        u = rand(rng) * tot
        if     u < r1;            x += 1
        elseif u < r1+r2;         y += 1
        elseif u < r1+r2+r3;      x -= 1
        else                      y -= 1
        end
        y > x && return (t, ev, xprev)   # crossed the separatrix from state (xprev,xprev)
    end
end
function ssa_stats(G, b, K, n, d; nrun=1500, seed=1, exit_hist=false)
    rng = MersenneTwister(seed)
    ts = Vector{Float64}(undef, nrun); evs = Vector{Int}(undef, nrun)
    hist = Dict{Int,Int}()
    for i in 1:nrun
        t, ev, m = ssa_cross(G, b, K, n, d, rng)
        ts[i] = t; evs[i] = ev
        exit_hist && (hist[m] = get(hist, m, 0) + 1)
    end
    xs = sort(collect(keys(hist)))
    (mfpt=mean(ts), mfpt_se=std(ts)/sqrt(nrun), events=mean(evs),
     exit_x=xs, exit_p=[hist[x]/nrun for x in xs])
end

# ---- parameters ------------------------------------------------------------
b, K, d, G = 20.0, 10.0, 1.0, 45          # grid & basin fixed across the sweep
n0 = 3.0

# ---- exactness certificate (representative n):  ELSE vs SSA first-passage ---
println("EXACTNESS  (toggle b=$b K=$K d=$d n=$n0, grid 0:$G, |J|=$(length(basin(G))))")
E  = else_escape(G, b, K, n0, d)
S0 = ssa_stats(G, b, K, n0, d; nrun=20000, exit_hist=true)
println("  init X-high mode = ", xhigh_mode(b, K, n0, d))
println("  ELSE MFPT = ", round(E.mfpt; sigdigits=6),
        "   SSA MFPT = ", round(S0.mfpt; sigdigits=6), " ±", round(S0.mfpt_se; sigdigits=2),
        "   (z = ", round(abs(E.mfpt-S0.mfpt)/S0.mfpt_se; sigdigits=2), " se)")
tv = 0.5*sum(abs, [ get(Dict(zip(E.exit_x,E.exit_p)), x, 0.0) -
                    get(Dict(zip(S0.exit_x,S0.exit_p)), x, 0.0)
                    for x in union(E.exit_x, S0.exit_x) ])
println("  exit-dist total-variation (ELSE vs SSA, ", 20000, " samples) = ", round(tv; sigdigits=3))

# ---- barrier sweep:  ELSE work flat, SSA events grow -----------------------
ns = [2.0, 2.5, 3.0, 3.5, 4.0]
me, ssa_m, ssa_s, ssa_e = Float64[], Float64[], Float64[], Float64[]
println("\nBARRIER SWEEP  (grid/basin fixed at |J|=$(length(basin(G))))")
println("   n  | ELSE MFPT | SSA MFPT (±se)      | SSA events | ELSE work")
for n in ns
    push!(me, else_escape(G, b, K, n, d).mfpt)
    s = ssa_stats(G, b, K, n, d; nrun=1500)
    push!(ssa_m, s.mfpt); push!(ssa_s, s.mfpt_se); push!(ssa_e, s.events)
    println("  ", rpad(n,3), " | ", rpad(round(me[end];sigdigits=5),9), " | ",
            rpad(string(round(s.mfpt;sigdigits=5))*" ±"*string(round(s.mfpt_se;sigdigits=2)),18),
            " | ", rpad(round(Int,s.events),10), " | ", length(basin(G)))
end

# ---- figure 1: exit distribution over the separatrix (ELSE vs SSA) ---------
mkpath(joinpath(@__DIR__, "..", "figs"))
f1 = bar(E.exit_x, E.exit_p; bar_width=0.6, color=3, alpha=0.6, label="ELSE (Z e_init, exact)",
         xlabel="crossing coordinate  x (= y on separatrix)", ylabel="exit probability",
         title="ELSE gives WHERE you leave the basin (n=$n0)", size=(760,460))
scatter!(f1, S0.exit_x, S0.exit_p; marker=:circle, ms=5, color=1, label="SSA (20k crossings)")
savefig(f1, joinpath(@__DIR__, "..", "figs", "toggle_exit.png"))

# ---- figure 2: cost vs barrier -- ELSE flat, SSA grows ---------------------
f2 = plot(ns, ssa_e; marker=:circle, ms=5, lw=2, color=1, yscale=:log10,
          xlabel="Hill cooperativity  n  (deeper barrier →)",
          ylabel="work to one exact crossing",
          title="ELSE cost is flat in the barrier SSA chokes on",
          label="SSA (reaction events)", legend=:left, size=(760,460))
hline!(f2, [length(basin(G))]; ls=:dash, lw=2, color=3,
       label="ELSE (one $(length(basin(G)))-state solve)")
savefig(f2, joinpath(@__DIR__, "..", "figs", "toggle_cost.png"))

println("\nfigures: figs/toggle_exit.png, figs/toggle_cost.png")
