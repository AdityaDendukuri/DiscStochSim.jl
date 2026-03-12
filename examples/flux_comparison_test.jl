# Quick numerical comparison of Φ_total, Φ_max, and Φ_out for Oregonator and Robertson
using DiscStochSim, Catalyst, LinearAlgebra

function analyze(label, sp, model, rates, t=0.0)
    A, _, _, erb = build_generator(sp, model, rates, t)
    d = diag(A)   # compressed diagonal (negative internal rates)
    Phi_int   = sum(-d .* sp.probs)          # what :total used (internal flux)
    Phi_max   = maximum(-d .* sp.probs)      # what :maximum used
    Phi_out   = dot(erb, sp.probs)           # true boundary outflux (new rule)
    Phi_total = Phi_int + Phi_out            # true total propensity × p
    println("=== $label ===")
    println("  |J|              = $(length(sp))")
    println("  Φ_int (old :tot) = $(round(Phi_int;   sigdigits=4))")
    println("  Φ_max (old :max) = $(round(Phi_max;   sigdigits=4))")
    println("  Φ_out (new rule) = $(round(Phi_out;   sigdigits=4))")
    println("  Φ_total (Φ_int+Φ_out) = $(round(Phi_total; sigdigits=4))")
    println("  Φ_out/Φ_int      = $(round(Phi_out/max(Phi_int,1e-300); sigdigits=3))")
    println("  dt :int  (ε=1)   = $(round(1/Phi_int;  sigdigits=3))")
    println("  dt :max  (ε=1)   = $(round(1/Phi_max;  sigdigits=3))")
    println("  dt :Φout (ε=1)   = $(Phi_out > 0 ? round(1/Phi_out; sigdigits=3) : Inf)")
    println()
end

# ─── Oregonator ──────────────────────────────────────────────────────────────
rn_o = @reaction_network begin
    @species X(t) Y(t) Z(t)
    @parameters k1 k2 k3 k4 k5
    k1, Y --> X
    k2, X + Y --> 0
    k3, X --> 2X + Z
    k4, 2X --> 0
    k5, Z --> Y
end
y1s, y2s, y3s = 500.0, 1000.0, 2000.0
mu1s, mu2s    = 2000.0, 50000.0
rates_o = [mu1s/y2s, mu2s/(y1s*y2s), (mu1s+mu2s)/y1s, 2*mu1s/y1s^2, (mu1s+mu2s)/y3s]

prob_o = FSPProblem(rn_o, CartesianIndex(500, 1000, 2000), (0.0, 1.0), rates_o; bounds=(0,50_000))
model_o = prob_o.model; bc_o = prob_o.bc

sp_o = StateSpace{CartesianIndex{3}, Float64}()
add_state!(sp_o, CartesianIndex(500, 1000, 2000), 1.0)
using Random; Random.seed!(42)
for _ in 1:100
    expand_ssa!(sp_o, model_o, rates_o, 0.0, bc_o; depth=1)
end
renormalize!(sp_o)
analyze("Oregonator (100 SSA expansions from IC)", sp_o, model_o, rates_o)

# ─── Robertson slow phase ─────────────────────────────────────────────────────
rn_r = @reaction_network begin
    k1, A --> B
    k2, 2B --> C + B
    k3, C + B --> A + C
end
rates_r = [0.04, 3e6, 1.0]
prob_r = FSPProblem(rn_r, CartesianIndex(10_000, 0, 0), (0.0, 1e4), rates_r; bounds=(0,100_001))
model_r = prob_r.model; bc_r = prob_r.bc

# Slow phase: spread p over 200 (a, 0, c) states with a+c=10000, plus B=1 neighbors in J
sp_r = StateSpace{CartesianIndex{3}, Float64}()
N, n_s = 10_000, 200
for i in 1:n_s
    a = round(Int, N * (i - 0.5) / n_s)
    add_state!(sp_r, CartesianIndex(a, 0, N-a), 1.0/n_s)
end
for i in 1:n_s   # add B=1 neighbors with 0 prob (well-expanded)
    a = round(Int, N * (i - 0.5) / n_s)
    a > 0 && add_state!(sp_r, CartesianIndex(a-1, 1, N-a), 0.0)
end
renormalize!(sp_r)
analyze("Robertson slow phase (spread A/C, B=0, B=1 neighbors in J)", sp_r, model_r, rates_r)

# Robertson slow phase WITHOUT B=1 neighbors (tight J)
sp_r2 = StateSpace{CartesianIndex{3}, Float64}()
for i in 1:n_s
    a = round(Int, N * (i - 0.5) / n_s)
    add_state!(sp_r2, CartesianIndex(a, 0, N-a), 1.0/n_s)
end
renormalize!(sp_r2)
analyze("Robertson slow phase (spread A/C, B=0 ONLY, tight J)", sp_r2, model_r, rates_r)

# Robertson burst phase
sp_r3 = StateSpace{CartesianIndex{3}, Float64}()
for da in -3:3
    a = max(0, 500+da)
    add_state!(sp_r3, CartesianIndex(a, 2, N-a-2), 1.0/7)
end
renormalize!(sp_r3)
analyze("Robertson burst phase (B=2 states, tight J)", sp_r3, model_r, rates_r)
