"""
validate_unified_rdme.jl

Validates the unified RDME builder `build_rdme_joint` (any Catalyst network ×
VoxelGraph). Single-species: a hand-built Schlögl reference over CartesianIndex{3}.
Multi-species: the well-mixed (K=1) limit and an independent K=2 reference. Each
generator is also fed to the scalable SA-AMG solver.
"""

using DiscStochSim, Catalyst, SparseArrays, LinearAlgebra, Printf

# ── parameters (Schlögl single-species kinetics + diffusion) ─────────────────────
c1 = 0.028; c2 = 3.2e-4; c3 = 19.5; c4 = 1.0; D = 2.0
graph = chain(3)                                   # 3-voxel chain

# ── reference: hand-built Schlögl graph generator over CartesianIndex{3} ──────────
e3(k) = CartesianIndex(ntuple(i -> i == k ? 1 : 0, 3))
function ref_schlogl()
    st = CartesianIndex{3}[]; pr = Function[]
    for k in 1:3
        let k = k, ek = e3(k)
            push!(st,  ek); push!(pr, (x,rv,t) -> c3)
            push!(st, -ek); push!(pr, (x,rv,t) -> c4 * max(0, x[k]))
            push!(st,  ek); push!(pr, (x,rv,t) -> begin nk=x[k]; c1*max(0,nk)*max(0,nk-1)/2 end)
            push!(st, -ek); push!(pr, (x,rv,t) -> begin nk=x[k]; c2*max(0,nk)*max(0,nk-1)*max(0,nk-2)/6 end)
        end
    end
    for (i, j) in graph.edges
        let i=i, j=j, ei=e3(i), ej=e3(j)
            push!(st, -ei+ej); push!(pr, (x,rv,t) -> D * max(0, x[i]))
            push!(st,  ei-ej); push!(pr, (x,rv,t) -> D * max(0, x[j]))
        end
    end
    DiscreteStochasticSystem{CartesianIndex{3}}(st, pr)
end
sysA = ref_schlogl()

# ── unified: same kinetics expressed as a Catalyst network ───────────────────────
rn = @reaction_network begin
    c3, ∅ --> X
    c4, X --> ∅
    c1, 2X --> 3X
    c2, 3X --> 2X
end
pmap = Dict(:c1 => c1, :c2 => c2, :c3 => c3, :c4 => c4)
P    = Float64[pmap[Symbol(string(p))] for p in Catalyst.get_ps(rn)]
sysB = build_rdme_joint(rn, graph; D = D)

# ── shared state space: full box [0,nmax]^3 ──────────────────────────────────────
nmax = 8
sp = StateSpace{CartesianIndex{3}, Float64}()
for n1 in 0:nmax, n2 in 0:nmax, n3 in 0:nmax
    add_state!(sp, CartesianIndex(n1, n2, n3))
end
@printf("state space: %d states (nmax=%d, K=3)\n", length(sp), nmax)

A1 = build_generator(sp, sysA, Float64[], 0.0)[1]
A2 = build_generator(sp, sysB, P,        0.0)[1]

# ── bit-for-bit agreement ────────────────────────────────────────────────────────
err = maximum(abs, A1 - A2)
@printf("max|A_schlogl - A_unified| = %.3e\n", err)
@printf("column-sum (mass) max|1'A_unified| = %.3e\n", maximum(abs, vec(sum(A2, dims = 1))))
@assert err < 1e-9 "unified generator does not match the Schlögl reference"
println("PASS: unified generator matches the hand-built Schlögl reference")

# ── scalable solver consumes it: one backward-Euler SA-AMG step ───────────────────
p0 = zeros(length(sp)); p0[sp.index[CartesianIndex(4, 4, 4)]] = 1.0
dt = 0.5
p1 = amg_be_step(A2, p0, dt; n_cycles = 30, rtol = 1e-10)
@printf("SA-AMG BE step: mass %.10f -> %.10f  (drift %.2e)\n",
        sum(p0), sum(p1), abs(sum(p1) - sum(p0)))
println("PASS: SA-AMG solved (I - dt*A_unified) p = p0")

# ═════════════════════════════════════════════════════════════════════════════════
#  Multi-species (S=2): A, B with annihilation A+B→∅ (couples species).
#  No reference RDME builder exists for S>1, so we check two independent ways.
# ═════════════════════════════════════════════════════════════════════════════════
println("\n── multi-species: A,B annihilation ──")
kA = 5.0; dA = 1.0; kB = 4.0; dB = 1.0; kab = 0.5; D2 = 1.5
rn2 = @reaction_network begin
    kA, ∅     --> A
    dA, A     --> ∅
    kB, ∅     --> B
    dB, B     --> ∅
    k,  A + B --> ∅
end
pm2 = Dict(:kA => kA, :dA => dA, :kB => kB, :dB => dB, :k => kab)
P2  = Float64[pm2[Symbol(string(p))] for p in Catalyst.get_ps(rn2)]

# Test 1 — K=1 must reproduce the well-mixed base generator bit-for-bit.
sysM1 = build_rdme_joint(rn2, chain(1); D = D2)     # no edges ⇒ pure base kinetics
base2 = DiscreteStochasticSystem(rn2)
m1 = 10
sp1 = StateSpace{CartesianIndex{2}, Float64}()
for a in 0:m1, b in 0:m1; add_state!(sp1, CartesianIndex(a, b)); end
Aw = build_generator(sp1, base2, P2, 0.0)[1]
Aj = build_generator(sp1, sysM1, P2, 0.0)[1]
@printf("K=1: max|A_wellmixed - A_joint| = %.3e\n", maximum(abs, Aw - Aj))
@assert maximum(abs, Aw - Aj) < 1e-9 "K=1 joint != well-mixed base"

# Test 2 — K=2, independent hand-built reference over CartesianIndex{4}.
# layout: comp (1,2)=voxel1 (A,B); comp (3,4)=voxel2 (A,B).
e4(k) = CartesianIndex(ntuple(i -> i == k ? 1 : 0, 4))
function ref_dss()
    st = CartesianIndex{4}[]; pr = Function[]
    for v in 1:2
        a = 2v - 1; b = 2v
        let a=a, b=b
            push!(st,  e4(a));        push!(pr, (x,rv,t) -> kA)
            push!(st, -e4(a));        push!(pr, (x,rv,t) -> dA * max(0, x[a]))
            push!(st,  e4(b));        push!(pr, (x,rv,t) -> kB)
            push!(st, -e4(b));        push!(pr, (x,rv,t) -> dB * max(0, x[b]))
            push!(st, -e4(a)-e4(b));  push!(pr, (x,rv,t) -> kab * max(0, x[a]) * max(0, x[b]))
        end
    end
    for (i, j) in ((1, 3), (2, 4))                  # species A: 1↔3, species B: 2↔4
        let i=i, j=j
            push!(st, -e4(i)+e4(j)); push!(pr, (x,rv,t) -> D2 * max(0, x[i]))
            push!(st,  e4(i)-e4(j)); push!(pr, (x,rv,t) -> D2 * max(0, x[j]))
        end
    end
    DiscreteStochasticSystem{CartesianIndex{4}}(st, pr)
end
sysM2  = build_rdme_joint(rn2, chain(2); D = D2)
sysRef = ref_dss()
m2 = 6
sp2 = StateSpace{CartesianIndex{4}, Float64}()
for a1 in 0:m2, b1 in 0:m2, a2 in 0:m2, b2 in 0:m2
    add_state!(sp2, CartesianIndex(a1, b1, a2, b2))
end
Ar = build_generator(sp2, sysRef, Float64[], 0.0)[1]
Au = build_generator(sp2, sysM2, P2,        0.0)[1]
@printf("K=2: %d states  max|A_ref - A_joint| = %.3e  mass max|1'A| = %.3e\n",
        length(sp2), maximum(abs, Ar - Au), maximum(abs, vec(sum(Au, dims = 1))))
@assert maximum(abs, Ar - Au) < 1e-9 "K=2 joint != independent reference"
println("PASS: multi-species joint matches well-mixed (K=1) and independent K=2 reference")

# scalable solver on the multi-species generator
q0 = zeros(length(sp2)); q0[sp2.index[CartesianIndex(3, 3, 3, 3)]] = 1.0
q1 = amg_be_step(Au, q0, 0.2; n_cycles = 30, rtol = 1e-10)
@printf("SA-AMG BE step (S=2): mass %.10f -> %.10f  (drift %.2e)\n",
        sum(q0), sum(q1), abs(sum(q1) - sum(q0)))
println("PASS: SA-AMG solved multi-species (I - dt*A) q = q0")
