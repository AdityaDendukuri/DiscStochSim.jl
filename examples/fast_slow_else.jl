# -----------------------------------------------------------------------------
# ELSE on a stiff CME:  A <=> B (fast),  B -> C (slow).   [prototype for §4.3]
#
# The (A,B) manifold (A+B = N, C=0) equilibrates over many fast flips, then
# rarely drains one molecule to C. First-escape kinetics illustrate ELSE's core
# identity on a chemical master equation:
#
#   ELSE : the fundamental matrix Z = -R^{-1} of the fast subnetwork R = A_JJ
#          (escape absorbed) gives the EXACT first-escape MFPT 1^T Z e_init and
#          the exact exit-state distribution -- ONE linear solve on the N+1
#          manifold states, cost independent of the stiffness kf/ks.
#   FSP  : flux-FSP on the full CME certifies exactness (survival integral).
#   SSA  : exact but must simulate every fast A<=>B flip -> cost grows ~ kf/ks.
#   QSSA : fast-equilibrium reduction -> cheap but biased at finite kf/ks.
#
# Result: ELSE matches FSP to ~1e-8, is unbiased where QSSA is not, and its cost
# is flat in the stiffness that makes SSA blow up.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, Plots, Statistics, Random
using DiscStochSim: build_generator, DiscreteStochasticSystem, StateSpace, add_state!
gr()

rn = @reaction_network begin
    kf, A --> B
    kr, B --> A
    ks, B --> C
end
model = DiscreteStochasticSystem(rn)

submanifold(N) = begin
    sp = StateSpace{CartesianIndex, Float64}()
    for a in 0:N; add_state!(sp, CartesianIndex(a, N-a, 0), 0.0); end
    sp
end

# ELSE: exact first-escape MFPT via the fundamental matrix of the fast manifold.
function else_mfpt(N, kf, kr, ks)
    R = Matrix(build_generator(submanifold(N), model, [kf,kr,ks], 0.0; absorbing=true)[1])
    e_init = zeros(N+1); e_init[N+1] = 1.0        # start at (a=N, b=0, c=0)
    sum(-(R \ e_init))                            # 1^T Z e_init,  Z = -R^{-1}
end

# FSP ground truth: integrate survival P(C=0, t) on the full CME.
function fsp_mfpt(N, kf, kr, ks; tf=60.0)
    prob = FSPProblem(rn, CartesianIndex(N,0,0), (0.0, tf), [kf,kr,ks]; bounds=(0, N))
    alg  = AdaptiveFSP(ε_dt=0.2, prob_quantile=0.0, expand_method=:stoich,
                       expansion_depth=2N, save_interval=1)
    sol, = solve(prob, alg)
    surv = [ (m = marginal(sol, i, 3); j = findfirst(==(0), m.values);
              j === nothing ? 0.0 : m.probs[j]) for i in 1:length(sol) ]
    sum(0.5 .* (surv[1:end-1] .+ surv[2:end]) .* diff(sol.t))
end

# SSA: Gillespie to first escape; returns (time_to_escape, n_events).
function ssa_first_escape(N, kf, kr, ks, rng)
    a, b, t, ev = N, 0, 0.0, 0
    while true
        r1, r2, r3 = kf*a, kr*b, ks*b          # A->B, B->A, B->C
        tot = r1 + r2 + r3
        t += randexp(rng) / tot
        ev += 1
        u = rand(rng) * tot
        if u < r1;        a -= 1; b += 1        # A->B
        elseif u < r1+r2; a += 1; b -= 1        # B->A
        else              return (t, ev)         # B->C : first escape
        end
    end
end
function ssa_stats(N, kf, kr, ks; nrun=3000, seed=1)
    rng = MersenneTwister(seed)
    ts = Vector{Float64}(undef, nrun); evs = Vector{Int}(undef, nrun)
    for i in 1:nrun; ts[i], evs[i] = ssa_first_escape(N, kf, kr, ks, rng); end
    (mfpt=mean(ts), mfpt_se=std(ts)/sqrt(nrun), events=mean(evs))
end

# --- exactness certificate --------------------------------------------------
N, ks = 20, 0.05
kf0 = kr0 = 5.0
println("EXACTNESS  (N=$N, kf=kr=$kf0, ks=$ks,  stiffness kf/ks=$(Int(kf0/ks)))")
me, mf = else_mfpt(N,kf0,kr0,ks), fsp_mfpt(N,kf0,kr0,ks)
println("  ELSE MFPT = ", round(me;sigdigits=7), "   FSP MFPT = ", round(mf;sigdigits=7),
        "   rel.diff = ", round(abs(me-mf)/mf;sigdigits=3))
println("  QSSA MFPT = ", round(1/(ks*N*kf0/(kf0+kr0)); sigdigits=6), "  (biased)")

# --- stiffness sweep --------------------------------------------------------
ratios = [2, 5, 10, 20, 50, 100, 200]          # stiffness kf/ks
mfpt_qssa = 1/(ks*N/2)                          # constant (equil fraction fixed)
me_r, ssa_mfpt, ssa_se, ssa_ev = Float64[], Float64[], Float64[], Float64[]
println("\nSTIFFNESS SWEEP  (QSSA MFPT = $(round(mfpt_qssa;sigdigits=4)), constant)")
println(" kf/ks | ELSE MFPT | SSA MFPT (±se)      | SSA events | ELSE work")
for r in ratios
    kf = kr = r*ks
    push!(me_r, else_mfpt(N,kf,kr,ks))
    s = ssa_stats(N,kf,kr,ks; nrun=3000)
    push!(ssa_mfpt, s.mfpt); push!(ssa_se, s.mfpt_se); push!(ssa_ev, s.events)
    println(" ", rpad(r,5), " | ", rpad(round(me_r[end];sigdigits=5),9), " | ",
            rpad(string(round(s.mfpt;sigdigits=5))*" ±"*string(round(s.mfpt_se;sigdigits=2)),18),
            " | ", rpad(round(Int,s.events),10), " | ", N+1)
end

# --- figure 1: MFPT -- ELSE exact, QSSA biased, SSA noisy ------------------
mkpath(joinpath(@__DIR__, "..", "figs"))
f1 = plot(ratios, me_r; marker=:star5, ms=7, lw=2, color=3, xscale=:log10,
          xlabel="stiffness  kf/ks", ylabel="first-escape MFPT", legend=:right,
          title="ELSE is exact; QSSA is biased", label="ELSE (exact)", size=(760,460))
scatter!(f1, ratios, ssa_mfpt; yerror=ssa_se, marker=:circle, ms=4, color=1, label="SSA (±se)")
hline!(f1, [mfpt_qssa]; ls=:dash, lw=2, color=2, label="QSSA (biased)")
savefig(f1, joinpath(@__DIR__, "..", "figs", "fastslow_mfpt.png"))

# --- figure 2: cost vs stiffness -- ELSE flat, SSA linear ------------------
f2 = plot(ratios, ssa_ev; marker=:circle, ms=5, lw=2, color=1, xscale=:log10, yscale=:log10,
          xlabel="stiffness  kf/ks", ylabel="work to one exact escape",
          title="ELSE cost is flat in the stiffness SSA chokes on",
          label="SSA (reaction events)", legend=:topleft, size=(760,460))
hline!(f2, [N+1]; ls=:dash, lw=2, color=3, label="ELSE (one $(N+1)-state solve)")
savefig(f2, joinpath(@__DIR__, "..", "figs", "fastslow_cost.png"))

println("\nfigures: figs/fastslow_mfpt.png, figs/fastslow_cost.png")
