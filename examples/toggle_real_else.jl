# -----------------------------------------------------------------------------
# ELSE on the REAL toggle (Tian & Burrage params) — s=0 resolvent: exact mean
# hop time U-high → V-high, in ONE sparse solve.  Your FSP evolves to tf=30 and
# never sees this hop if it is rare; ELSE gives its mean directly.
#
# Basin J = U-high well = {U ≥ V} on a grid large enough to CONTAIN the well
# (mode ≈ (378,27), above your NGRID=300), so the only escape is the U=V crossing.
# -----------------------------------------------------------------------------
using DiscStochSim, Catalyst, LinearAlgebra, SparseArrays
using DiscStochSim: build_generator, DiscreteStochasticSystem, StateSpace, add_state!

rn = @reaction_network begin
    @species U(t) V(t)
    @parameters η α1 β1 K1 d1 dU α2 β2 K2 d2
    η*(α1 + β1*K1^3/(K1^3 + V^3)),   0 --> U
    d1 + dU,                           U --> 0
    η*(α2 + β2*K2^3/(K2^3 + U^3)),   0 --> V
    d2,                                V --> 0
end
rates = [1.0, 20.0, 400.0, 100.0, 1.0, 0.1/1.1, 20.0, 400.0, 100.0, 1.0]
model = DiscreteStochasticSystem(rn)

# deterministic U-high fixed point (to size the grid / sanity)
bU(V) = 1.0*(20 + 400*100.0^3/(100.0^3 + V^3)) / (1.0 + 0.1/1.1)
bV(U) = 1.0*(20 + 400*100.0^3/(100.0^3 + U^3)) / 1.0
function uhigh_mode()
    U, V = 85.0, 5.0
    for _ in 1:5000; U = 0.5U + 0.5*bU(V); V = 0.5V + 0.5*bV(U); end
    (round(Int,U), round(Int,V))
end

GU = 680          # top boundary far above the mode (~378) so top-edge escape is negligible
function uhigh_basin()
    sp = StateSpace{CartesianIndex, Float64}()
    for u in 0:GU, v in 0:u                 # J = {0 ≤ V ≤ U ≤ GU};  only escape should be U=V
        add_state!(sp, CartesianIndex(u, v), 0.0)
    end
    sp
end

mode = uhigh_mode()
println("U-high mode ≈ ", mode, "   (grid U:0..$GU)")

sp = uhigh_basin()
println("|J| = ", length(sp), " states — building generator...")
@time A, _, _, Δ = build_generator(sp, model, rates, 0.0; absorbing=true)

p0 = zeros(length(sp)); p0[sp.index[CartesianIndex(85, 5)]] = 1.0   # your u0 = (85,5)

println("solving  Z p0 = -A_JJ \\ p0  (one sparse solve)...")
@time occ = -(A \ p0)                     # s = 0 resolvent = Z p0
mfpt = sum(occ)                           # mean hop time U-high → V-high

# escape flux ∝ Δ ∘ (Z p0);  check it concentrates on the separatrix U=V
wts = Δ .* occ
onsep = [sp.states[j][1] == sp.states[j][2] for j in eachindex(wts)]
sep_frac = sum(wts[onsep]) / sum(wts)                       # fraction escaping across U=V
exit_mean = sum(sp.states[j][1]*wts[j] for j in eachindex(wts) if onsep[j]) / sum(wts[onsep])

println("\nmean hop time  U-high → V-high  =  ", round(mfpt; sigdigits=6))
println("  (your FSP runs to tf = 30 — ratio hop/tf = ", round(mfpt/30; sigdigits=4), ")")
println("escape flux across separatrix U=V =  ", round(100*sep_frac; sigdigits=4), "%  (rest = spurious grid-edge)")
println("mean crossing coordinate U(=V)    =  ", round(exit_mean; sigdigits=4))
println("occupation check: max occ at state ", sp.states[argmax(occ)])
