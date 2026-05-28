"""
QTT-RDME: Quantized Tensor Train for the Reaction-Diffusion Master Equation
===========================================================================
Based on: Kazeev et al. "Quantized Tensor-Train representation of CME" (2106.15188)

Key ideas borrowed from Kazeev:
1. QTT encoding: molecule count n ∈ {0,...,2^q-1} → q BINARY sites (dim=2)
   - Physical dim fixed at 2 regardless of n_max
   - Cost scales as O(q) = O(log n_max) instead of O(n_max)
   - For Schlögl n_max=220: q=8 bits, local dim 2 vs 221 (110× cheaper per site)

2. TDVP time evolution instead of expv:
   - expv fails for CME: large negative diagonals (-K×k_b) make exp(H dt)≈0
   - TDVP solves local Krylov problems at each BINARY SITE (dim 2r², tiny)
   - These local problems always converge — no large-diagonal issue
   - This is the "DMRG to solve linear equations" from Kazeev

3. N_k as sum of binary projectors:
   n_k = Σ_j 2^{j-1} × b_{k,j}   →   ⟨n_k⟩ = Σ_j 2^{j-1} × ⟨P_{k,j}⟩
   No need for n_max×n_max matrix for the number operator!

QTT site ordering: voxel-major, bit-minor
   Sites: [b₁(v1), b₂(v1), ..., b_q(v1), b₁(v2), ..., b_q(vK)]
   Total sites M = K × q   (vs M = K in standard TT)

Run: julia --project=examples/itensor_env examples/qtt_rdme.jl
"""

using ITensors, ITensorMPS, Printf, LinearAlgebra

# ── Qbit site type: binary digit with CME-compatible operators ─────────────────
ITensors.space(::SiteType"Qbit") = 2   # always dim=2

function ITensors.op(::OpName"Id",    ::SiteType"Qbit", s::Index)
    O = ITensor(s',dag(s)); O[s'=>1,s=>1]=1.0; O[s'=>2,s=>2]=1.0; O
end
function ITensors.op(::OpName"Plus",  ::SiteType"Qbit", s::Index) # |1⟩⟨0|
    O = ITensor(s',dag(s)); O[s'=>2,s=>1]=1.0; O
end
function ITensors.op(::OpName"Minus", ::SiteType"Qbit", s::Index) # |0⟩⟨1|
    O = ITensor(s',dag(s)); O[s'=>1,s=>2]=1.0; O
end
function ITensors.op(::OpName"P",     ::SiteType"Qbit", s::Index) # |1⟩⟨1|
    O = ITensor(s',dag(s)); O[s'=>2,s=>2]=1.0; O
end
function ITensors.op(::OpName"Q",     ::SiteType"Qbit", s::Index) # |0⟩⟨0|
    O = ITensor(s',dag(s)); O[s'=>1,s=>1]=1.0; O
end

# ── Binary encoding helpers ────────────────────────────────────────────────────

to_bits(n::Int, q::Int) = [Bool((n >> (j-1)) & 1) for j in 1:q]

# Operator for a single bit transition b_in → b_out at a given site
function bit_op_name(b_in::Bool, b_out::Bool)
    !b_in && !b_out && return "Q"     # stay at 0
     b_in &&  b_out && return "P"     # stay at 1
    !b_in &&  b_out && return "Plus"  # 0 → 1
    return "Minus"                     # 1 → 0
end

# Add an n→m transition at voxel sites sk with given coefficient
function add_transition!(os, coeff, n::Int, m::Int, sk, q)
    bn = to_bits(n, q); bm = to_bits(m, q)
    ops = [(bit_op_name(bn[j], bm[j]), sk[j]) for j in 1:q]
    add!(os, coeff, Iterators.flatten(ops)...)
end

# ── Build H_eff = -L as QTT MPO via AutoMPO ───────────────────────────────────
"""
    build_qtt_heff(sites, Kx, Ky, q, k_b, k_d, D) → MPO

Build the QTT operator H_eff = -L (the CME generator negated) as an MPO
on M = Kx×Ky×q binary sites.

H_eff = -L ensures eigenvalues ≥ 0, enabling stable TDVP imaginary-time
evolution: e^{-H_eff t} P = e^{L t} P (converges toward stationary state).

Key insight from Kazeev: operators in QTT basis have small MPO bond width:
- Birth gain: n_max terms × q-site operators  → W = O(n_max × q) at transitions
- Number operator N_k: q terms × 1-site operators → W = O(1) per bit!
- Death gain: n_max terms × q-site operators

Column diffusion MPO bond dim: W = O(Ky × q) due to long-range interactions
→ this is the remaining bottleneck for 2D systems.
"""
function build_qtt_heff(sites, Kx::Int, Ky::Int, q::Int,
                         k_b::Float64, k_d::Float64, D::Float64)
    K    = Kx * Ky
    M    = K * q   # total sites
    os   = OpSum()
    nmax = (1 << q) - 1   # 2^q - 1 = max molecule count

    for k in 1:K
        # Binary sites for voxel k
        sk = [(k-1)*q + j for j in 1:q]

        # ── Birth: H_eff_birth = k_b(I - Plus_k) ─────────────────────────────
        # -k_b × Plus_k: gain transitions n→n+1 (absorbing FSP)
        for n in 0:nmax-1
            add_transition!(os, -k_b, n, n+1, sk, q)
        end
        # +k_b × I (loss diagonal, all states)
        add!(os, k_b, "Id", sk[1])

        # ── Death: H_eff_death = k_d(N_k - Minus_k) ──────────────────────────
        # +k_d × N_k = Σ_j 2^{j-1} × P_{k,j}  [KEY QTT SIMPLIFICATION!]
        # Number operator decomposes into q single-site projectors — no matrix!
        for j in 1:q
            add!(os, k_d * Float64(1 << (j-1)), "P", sk[j])
        end
        # -k_d × Minus_k: propensity-weighted transitions n→n-1
        for n in 1:nmax
            add_transition!(os, -k_d * Float64(n), n, n-1, sk, q)
        end
    end

    # ── Row diffusion (nearest-neighbor in MPS ordering) ────────────────────
    for i in 1:Kx, j in 1:Ky-1
        k1 = (i-1)*Ky + j;  k2 = k1+1
        _add_diff!(os, k1, k2, q, D, nmax)
    end

    # ── Column diffusion (distance q×Ky sites in MPS ordering) ──────────────
    # These are LONG-RANGE operators — AutoMPO threads through identity sites.
    # MPO bond dim W = O(n_max² × something) at crossing sites.
    # For small n_max this is manageable; for large n_max use carry-adder MPO.
    for i in 1:Kx-1, j in 1:Ky
        k1 = (i-1)*Ky + j;  k2 = i*Ky + j
        _add_diff!(os, k1, k2, q, D, nmax)
    end

    os   # return the OpSum
end

function _add_diff!(os, k1::Int, k2::Int, q::Int, D::Float64, nmax::Int)
    sk1 = [(k1-1)*q + j for j in 1:q]
    sk2 = [(k2-1)*q + j for j in 1:q]
    combined_sk = [sk1; sk2]

    # Hop k1→k2 with propensity D×n1:
    for n1 in 1:nmax, n2 in 0:nmax-1
        bn1_in  = to_bits(n1, q);   bn1_out = to_bits(n1-1, q)
        bn2_in  = to_bits(n2, q);   bn2_out = to_bits(n2+1, q)
        ops     = [(bit_op_name(bn1_in[j], bn1_out[j]), sk1[j]) for j in 1:q]
        ops2    = [(bit_op_name(bn2_in[j], bn2_out[j]), sk2[j]) for j in 1:q]
        add!(os, -D * Float64(n1), Iterators.flatten([ops; ops2])...)
    end
    # Diagonal loss from k1: +D × N_k1
    for j in 1:q
        add!(os, D * Float64(1 << (j-1)), "P", sk1[j])
    end

    # Hop k2→k1 with propensity D×n2 (symmetric):
    for n1 in 0:nmax-1, n2 in 1:nmax
        bn1_in  = to_bits(n1, q);   bn1_out = to_bits(n1+1, q)
        bn2_in  = to_bits(n2, q);   bn2_out = to_bits(n2-1, q)
        ops     = [(bit_op_name(bn1_in[j], bn1_out[j]), sk1[j]) for j in 1:q]
        ops2    = [(bit_op_name(bn2_in[j], bn2_out[j]), sk2[j]) for j in 1:q]
        add!(os, -D * Float64(n2), Iterators.flatten([ops; ops2])...)
    end
    for j in 1:q
        add!(os, D * Float64(1 << (j-1)), "P", sk2[j])
    end
end

# ── QTT state encoding / decoding ─────────────────────────────────────────────

"""
    qtt_product_state(sites, K, q, counts) → MPS

Build a product-state QTT MPS where voxel k has exactly counts[k] molecules.
Each voxel contributes q binary sites encoding counts[k] in binary.
"""
function qtt_product_state(sites, K::Int, q::Int, counts::Vector{Int})
    M = K * q
    psi = MPS(sites; linkdims=1)
    for k in 1:K
        n = counts[k]
        bits = to_bits(n, q)
        for j in 1:q
            s  = (k-1)*q + j
            v  = zeros(2)
            v[bits[j] ? 2 : 1] = 1.0   # index 1=bit0, index 2=bit1
            psi[s] = ITensor(v, sites[s])
        end
    end
    psi
end

"""
    qtt_means(psi, sites, K, q) → Vector{Float64}

Compute ⟨n_k⟩ for each voxel using L1 inner products (probability-correct).

For a PROBABILITY distribution P (not a quantum state):
  ⟨n_k⟩ = Σ_n P(n) × n_k = ⟨all_ones | N_k | P⟩   [L1 weighted mean]
         = Σ_j 2^{j-1} × ⟨all_ones | P_{k,j} | P⟩

NOT ⟨ψ|N_k|ψ⟩/⟨ψ|ψ⟩ which is L2-weighted (wrong for probability distributions).

The all-ones bra ⟨1| has components ones[n]=1 for all n, so:
  ⟨all_ones|P⟩ = Σ_n P(n) = L1 norm (probability conservation check)
"""
function build_ones_mps(sites)
    M = length(sites)
    psi_1 = MPS(sites; linkdims=1)
    for s in 1:M
        psi_1[s] = ITensor([1.0, 1.0], sites[s])
    end
    psi_1
end

function qtt_l1_norm(psi::MPS, ones_mps::MPS)
    real(inner(ones_mps, psi))
end

function qtt_means(psi::MPS, sites, K::Int, q::Int, ones_mps::MPS)
    L1_norm = qtt_l1_norm(psi, ones_mps)
    μ = zeros(K)
    for k in 1:K, j in 1:q
        s    = (k-1)*q + j
        # Apply P_s = |1⟩⟨1| at site s: use apply without explicit site index
        # ITensorMPS detects the site from the operator's site indices
        P_op  = op("P", sites[s])   # ITensor with primed+unprimed site[s] indices
        P_psi = apply(P_op, psi; cutoff=1e-14)
        μ[k] += Float64(1 << (j-1)) * real(inner(ones_mps, P_psi)) / L1_norm
    end
    μ
end

# ── Kazeev backward-Euler time step (Algorithm 2, Box 2) ─────────────────────
"""
    qtt_step!(psi, H_eff, dt; max_r, cutoff) → MPS

Kazeev Algorithm 2 (Box 2) backward-Euler step:
    solve (I - dt·L) p_{n+1} = p_n
where L = -H_eff is the CME generator and H_eff = -L ≥ 0.

Rewritten as: (I + dt·H_eff) p_{n+1} = p_n
since H_eff = -L → I - dt·L = I + dt·H_eff.

Solved via ITensorMPS `linsolve` (DMRG-based ALS for linear systems):
- (I + dt·H_eff) x = rhs  where rhs = p_n
- The system matrix has eigenvalues ≥ 1 for any dt > 0 (always well-conditioned!)
- DMRG sweeps on the QTT structure: local problems are 2r² × 2r² (tiny)
- This is the exact stability fix vs TDVP (which has a diverging reverse sweep)
"""
function qtt_step!(psi::MPS, H_eff::MPO, dt::Float64;
                    max_r::Int      = 8,
                    cutoff::Float64 = 1e-12,
                    nsweeps::Int    = 15,
                    ones_mps::Union{MPS,Nothing} = nothing)

    # Kazeev Algorithm 2: backward Euler via DMRG linsolve
    # Solve (I + dt × H_eff) p_{n+1} = p_n
    # = (I - dt × L) p_{n+1} = p_n  (since H_eff = -L)
    # Matrix I - dt×L has all eigenvalues ≥ 1 → always well-conditioned
    # Verified correct for K=1 (single voxel): exact QTT mean matches FSP ✓
    sites  = siteinds(psi)
    I_mpo  = MPO(sites, "Id")
    BE_mpo = add(I_mpo, dt * H_eff; maxdim=maxlinkdim(H_eff)+4, cutoff=1e-14)

    psi_new = linsolve(BE_mpo, psi, psi;
                        nsweeps     = nsweeps,
                        maxdim      = max_r,
                        cutoff      = cutoff,
                        outputlevel = 0)
    truncate!(psi_new; maxdim=max_r, cutoff=cutoff)

    # L1 normalize to correct for linsolve approximation error
    if ones_mps !== nothing
        l1 = real(inner(ones_mps, psi_new))
        l1 > 1e-10 && (psi_new = psi_new / l1)
    end
    psi_new
end

# ── Main ───────────────────────────────────────────────────────────────────────
const Kx = 5;  const Ky = 5   # 5×5 = 25 voxels (verify on small case first)
const k_b = 1.0; const k_d = 0.5; const D = 1.0
const q    = 4     # bits per voxel → n_max = 2^4-1 = 15
const max_r = 4
const dt   = 0.5;  const T_END = 8.0   # backward Euler: stable for any dt

K = Kx * Ky; M = K * q
nmax = (1 << q) - 1

@printf("=== QTT-RDME (%d×%d grid, q=%d bits, n_max=%d) ===\n", Kx, Ky, q, nmax)
@printf("k_b=%.1f  k_d=%.1f  D=%.1f  μ_ss=%.1f\n", k_b, k_d, D, k_b/k_d)
@printf("M = K×q = %d×%d = %d binary sites (vs %d standard TT sites)\n\n", K, q, M, K)

# Build sites and H_eff
sites = siteinds("Qbit", M)
@printf("Building QTT MPO (H_eff = -L)...\n"); flush(stdout)
t0 = time()
os_heff = build_qtt_heff(sites, Kx, Ky, q, k_b, k_d, D)
H_eff   = MPO(os_heff, sites)   # construct MPO from OpSum
@printf("  Done in %.1fs  |  max MPO bond W=%d\n\n", time()-t0, maxlinkdim(H_eff))

# IC: 1 molecule at centre
centre = (Kx÷2)*Ky + Ky÷2 + 1
counts = zeros(Int, K); counts[centre] = 1
psi    = qtt_product_state(sites, K, q, counts)
ones_mps = build_ones_mps(sites)   # for L1 inner products
@printf("IC: 1 molecule at voxel %d (row=%d col=%d)\n", centre, Kx÷2+1, Ky÷2+1)
@printf("    L1 norm (prob sum) = %.6f  (should be 1.0)\n\n",
        qtt_l1_norm(psi, ones_mps))

# Time loop
nsteps = round(Int, T_END/dt)
@printf("%-5s  %-7s  %-12s  %-12s  %-8s\n",
        "step", "t", "total_μ", "analytical", "r_max")
@printf("%s\n", "─"^55)

for step in 1:nsteps
    global psi = qtt_step!(psi, H_eff, dt; max_r=8, cutoff=1e-12, nsweeps=20, ones_mps)
    t  = step * dt
    μ   = qtt_means(psi, sites, K, q, ones_mps)
    l1n = qtt_l1_norm(psi, ones_mps)
    tgt = K * (k_b/k_d) * (1-exp(-k_d*t)) + exp(-k_d*t)
    @printf("%-5d  %-7.1f  %-12.4f  %-12.4f  %-8d  L1=%.4f\n",
            step, t, sum(μ), tgt, maxlinkdim(psi), l1n)
    flush(stdout)
end

println("\nFinal means (5×5 grid):")
μ_final = reshape(qtt_means(psi, sites, K, q, ones_mps), Kx, Ky)
for i in 1:Kx
    @printf("  [%s]\n", join([@sprintf("%5.3f", μ_final[i,j]) for j in 1:Ky], " "))
end
@printf("\nTotal: %.4f  target: %.4f  L1_norm=%.4f\n",
        sum(μ_final), K*(k_b/k_d)*(1-exp(-k_d*T_END))+exp(-k_d*T_END),
        qtt_l1_norm(psi, ones_mps))
