# Test flux-based expansion: examine propensity distribution for Oregonator states
using DiscStochSim, Catalyst, LinearAlgebra, Statistics, Random

rn = @reaction_network begin
    @species X(t) Y(t) Z(t); @parameters k1 k2 k3 k4 k5
    k1, Y --> X; k2, X + Y --> 0; k3, X --> 2X + Z; k4, 2X --> 0; k5, Z --> Y
end
y1s,y2s,y3s=500.,1000.,2000.; mu1s,mu2s=2000.,50000.
rates=[mu1s/y2s,mu2s/(y1s*y2s),(mu1s+mu2s)/y1s,2*mu1s/y1s^2,(mu1s+mu2s)/y3s]
prob_o = FSPProblem(rn, CartesianIndex(500,1000,2000), (0.0,1.0), rates; bounds=(0,50_000))
model = prob_o.model; bc = prob_o.bc
rxn_names = ["k1:Y→X", "k2:X+Y→0", "k3:X→2X+Z", "k4:2X→0", "k5:Z→Y"]
n_rxn = length(model.stoichvecs)
x0 = CartesianIndex(500, 1000, 2000)

# ─── 1. Propensity distribution at IC ────────────────────────────────────────
alphas = Float64[Float64(model.propensities[k](x0, rates, 0.0)) for k in 1:n_rxn]
w_tot = sum(alphas)
println("=== Propensity distribution at IC (500, 1000, 2000) ===")
for k in 1:n_rxn
    pct = round(100*alphas[k]/w_tot; digits=1)
    println("  $(rxn_names[k]): α=$(round(alphas[k];sigdigits=4)) ($pct%)$(pct > 10 ? " ← flux-expand" : "")")
end
println("  w_total = $(round(w_tot;sigdigits=4))\n")

# ─── 2. SSA selection frequency (10000 samples) ──────────────────────────────
println("=== SSA selection frequency at IC (10000 samples) ===")
Random.seed!(1)
counts = zeros(Int, n_rxn)
for _ in 1:10_000
    u = rand() * w_tot; csum = 0.0
    for k in 1:n_rxn
        csum += alphas[k]
        if u <= csum; counts[k] += 1; break; end
    end
end
for k in 1:n_rxn
    println("  $(rxn_names[k]): expected $(round(100*alphas[k]/w_tot;digits=1))%,  got $(round(100*counts[k]/10_000;digits=1))%")
end
println()

# ─── 3. Build mid-simulation SSA-expanded state space ────────────────────────
sp = StateSpace{CartesianIndex{3}, Float64}()
add_state!(sp, x0, 1.0)
Random.seed!(42)
for _ in 1:200
    expand_ssa!(sp, model, rates, 0.0, bc; depth=1)
end
renormalize!(sp)
println("=== After 200 SSA expansions: |J| = $(length(sp)) ===")
active_states = [(sp.states[j], sp.probs[j]) for j in 1:length(sp) if sp.probs[j] > 1e-6]
println("  States with p > 1e-6: $(length(active_states))")

# For each active state: how many reactions would flux-expand add? What fraction already in J?
n_flux_total = Int[]; n_in_J_total = Int[]
for (x, _) in active_states
    αs = Float64[Float64(model.propensities[k](x, rates, 0.0)) for k in 1:n_rxn]
    wx = sum(αs)
    wx > 0 || continue
    n_f = sum(αs[k]/wx > 0.1 ? 1 : 0 for k in 1:n_rxn)
    n_j = sum((αs[k]/wx > 0.1 && haskey(sp.index, x + model.stoichvecs[k])) ? 1 : 0 for k in 1:n_rxn)
    push!(n_flux_total, n_f); push!(n_in_J_total, n_j)
end
println("  Flux neighbors (threshold=0.1): mean=$(round(mean(n_flux_total);digits=2)), range=$(minimum(n_flux_total))-$(maximum(n_flux_total))")
frac_covered = mean(n_in_J_total[i] / max(n_flux_total[i], 1) for i in 1:length(n_flux_total))
println("  Already in J: mean=$(round(mean(n_in_J_total);digits=2)) of $(round(mean(n_flux_total);digits=2)) flux neighbors")
println("  Coverage: $(round(100*frac_covered;digits=1))% of flux neighbors already in J via SSA")
println()

# ─── 4. Φ comparison ─────────────────────────────────────────────────────────
A, _, _, erb = build_generator(sp, model, rates, 0.0)
Phi_int = dot(-LinearAlgebra.diag(A), sp.probs)
Phi_out = dot(erb, sp.probs)
Phi_tot = Phi_int + Phi_out
println("=== Flux quantities (SSA-expanded J) ===")
println("  Φ_int = $(round(Phi_int; sigdigits=4))  (in J)")
println("  Φ_out = $(round(Phi_out; sigdigits=4))  (exits J)")
println("  Φ_out/Φ_tot = $(round(Phi_out/Phi_tot; sigdigits=3))")
println("  dt(ε=1, Φ_int) = $(round(1/Phi_int; sigdigits=3))")
println("  Approx with flux expand (3/5 neighbors in J): Φ_int ×3")
println("  → dt ≈ $(round(1/(3*Phi_int); sigdigits=3))  vs current $(round(1/Phi_int; sigdigits=3))")
