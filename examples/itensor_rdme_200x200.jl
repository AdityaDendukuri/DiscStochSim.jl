"""
Full joint 200×200 RDME via ITensors.jl — TEBD gate approach
=============================================================

Time evolution via Strang splitting + pre-computed matrix exponential gates:

  P(t+dt) = exp(dt/2 L_local) exp(dt L_row) exp(dt L_col) exp(dt/2 L_local) P(t)

Each factor is applied as a set of 2-site ITensor gates.
- L_local: per-voxel birth/death CME — single-site gates (no bond growth)
- L_row:   row diffusion (nearest-neighbor in MPS ordering) — 2-site gates
- L_col:   column diffusion (distance Ky in MPS ordering) — ITensors handles
           non-adjacent gates via SWAP operations internally

This avoids the Krylov-expv issue (huge negative diagonal of L makes
exp(H_krylov dt) ≈ 0), and instead uses exact 2×2 block matrix exponentials.

Bond dimension r controls the accuracy of the joint distribution approximation.
Run: julia --project=examples/itensor_env examples/itensor_rdme_200x200.jl
"""

using ITensors, ITensorMPS, Printf, LinearAlgebra

# ── CME site type (Doi-Peliti operators) ───────────────────────────────────────
ITensors.space(::SiteType"CME"; dim::Int=6) = dim

function ITensors.op(::OpName"N", ::SiteType"CME", s::Index)
    d = ITensors.dim(s); O = ITensor(s', dag(s))
    for n in 1:d; O[s'=>n, s=>n] = Float64(n-1); end; O   # n_count = index-1
end
function ITensors.op(::OpName"Id", ::SiteType"CME", s::Index)
    d = ITensors.dim(s); O = ITensor(s', dag(s))
    for n in 1:d; O[s'=>n, s=>n] = 1.0; end; O
end

# ── Local (single-site) CME generator: birth + death ──────────────────────────
function local_cme_generator(k_b::Float64, k_d::Float64, n_max::Int)
    nv = n_max + 1
    L  = zeros(Float64, nv, nv)
    for n in 0:n_max      # n = molecule count (0-indexed)
        ni = n + 1        # Julia index (1-indexed)
        # Birth ∅→A: gain at n+1 from n, loss at n
        n < n_max && (L[ni+1, ni] += k_b)    # gain: n → n+1
        L[ni, ni] -= k_b                       # loss: outflow at n
        # Death A→∅: gain at n-1 from n, loss at n
        n > 0 && (L[ni-1, ni] += k_d * n)    # gain: n → n-1
        L[ni, ni] -= k_d * n                   # loss
    end
    L
end

# ── Pair (2-site) diffusion CME generator ─────────────────────────────────────
function pair_diffusion_generator(D::Float64, n_max::Int)
    nv = n_max + 1; N = nv * nv
    L  = zeros(Float64, N, N)
    for n1 in 0:n_max, n2 in 0:n_max
        col = n1*nv + n2 + 1
        L[col, col] -= D * (n1 + n2)          # total outflow
        n1 > 0 && n2 < n_max &&
            (L[(n1-1)*nv + (n2+1) + 1, col] += D * n1)   # hop n1→n2
        n2 > 0 && n1 < n_max &&
            (L[(n1+1)*nv + (n2-1) + 1, col] += D * n2)   # hop n2→n1
    end
    L
end

# ── Gate constructors ──────────────────────────────────────────────────────────
function single_site_gate(mat::Matrix{Float64}, s::Index)
    d = ITensors.dim(s); nv = d
    G = ITensor(s', dag(s))
    for i in 1:nv, j in 1:nv
        G[s'=>i, s=>j] = mat[i, j]
    end
    G
end

function two_site_gate(mat::Matrix{Float64}, s1::Index, s2::Index)
    d = ITensors.dim(s1); nv = d
    G = ITensor(s1', s2', dag(s1), dag(s2))
    for n1 in 1:nv, n2 in 1:nv, m1 in 1:nv, m2 in 1:nv
        row = (n1-1)*nv + n2
        col = (m1-1)*nv + m2
        G[s1'=>n1, s2'=>n2, s1=>m1, s2=>m2] = mat[row, col]
    end
    G
end

# ── Main simulation ────────────────────────────────────────────────────────────
const Kx    = 50;  const Ky = 50   # 50×50 = 2500 voxels
const k_b   = 1.0; const k_d = 0.5
const D     = 1.0
const n_max = 6
const max_r = 4
const dt    = 0.5; const T_END = 10.0

M = Kx * Ky
@printf("=== ITensors TEBD joint RDME (%d×%d = %d voxels) ===\n", Kx, Ky, M)
@printf("k_b=%.1f  k_d=%.1f  D=%.1f  μ_ss=%.1f  r=%d  n_max=%d\n\n",
        k_b, k_d, D, k_b/k_d, max_r, n_max)

# Precompute gate matrices (done once)
@printf("Precomputing gate matrices...\n")
L_local = local_cme_generator(k_b, k_d, n_max)
L_pair  = pair_diffusion_generator(D, n_max)
e_local_half = exp(L_local * (dt/2))   # exact matrix exp for local reactions
e_pair_full  = exp(L_pair  *  dt)      # exact matrix exp for pair diffusion
@printf("  Local gate: %.4f  (sum per col should be 1)\n", maximum(sum(e_local_half; dims=1)))
@printf("  Pair gate:  %.4f  (same)\n", maximum(sum(e_pair_full;  dims=1)))

# Site indices
sites = siteinds("CME", M; dim=n_max+1)

# Build gate lists (once, not every step)
@printf("Building gate lists...\n")
t0 = time()

# Local reaction gates (dt/2 each)
local_gates = [single_site_gate(e_local_half, sites[j]) for j in 1:M]

# Row diffusion gates (nearest-neighbor in 1D ordering)
row_gates = ITensor[]
for i in 1:Kx, j in 1:Ky-1
    s = (i-1)*Ky + j   # left site (row-major)
    push!(row_gates, two_site_gate(e_pair_full, sites[s], sites[s+1]))
end

# Column diffusion gates (non-adjacent: distance Ky)
# ITensors' apply() handles these via SWAP gates automatically
col_gates = ITensor[]
for j in 1:Ky, i in 1:Kx-1
    s_top = (i-1)*Ky + j   # site (i,j)
    s_bot = i*Ky + j        # site (i+1,j)
    push!(col_gates, two_site_gate(e_pair_full, sites[s_top], sites[s_bot]))
end

@printf("  local=%d  row=%d  col=%d gates  (%.1fs)\n\n",
        length(local_gates), length(row_gates), length(col_gates), time()-t0)

# Initial MPS: product state, 1 molecule at centre
centre = (Kx÷2)*Ky + (Ky÷2) + 1
psi    = MPS(sites; linkdims=1)
for s in 1:M
    d = ITensors.dim(sites[s]); v = zeros(d)
    v[s == centre ? 2 : 1] = 1.0   # index 1=n=0, index 2=n=1
    psi[s] = ITensor(v, sites[s])
end
@printf("IC: 1 molecule at site %d (row=%d col=%d)  norm=%.6f\n\n",
        centre, Kx÷2+1, Ky÷2+1, norm(psi))

# Time evolution
nsteps = round(Int, T_END/dt)
@printf("%-5s  %-7s  %-12s  %-12s  %-8s  %-10s\n",
        "step", "t", "total_μ", "analytical", "r_max", "s/step")
@printf("%s\n", "─"^65)

t_total = Ref(0.0)
for step in 1:nsteps
    te = @elapsed begin
        # Strang splitting: local/2 → row → col → local/2
        global psi = apply(local_gates, psi; maxdim=max_r, cutoff=1e-12)
        global psi = apply(row_gates,   psi; maxdim=max_r, cutoff=1e-12)
        global psi = apply(col_gates,   psi; maxdim=max_r, cutoff=1e-12)
        global psi = apply(local_gates, psi; maxdim=max_r, cutoff=1e-12)
        normalize!(psi)
    end
    t_total[] += te

    t       = step * dt
    mu_vec  = real(expect(psi, "N"))
    total_μ = sum(mu_vec)
    r_now   = maxlinkdim(psi)
    target  = M * (k_b/k_d) * (1-exp(-k_d*t)) + exp(-k_d*t)

    @printf("%-5d  %-7.1f  %-12.4f  %-12.4f  %-8d  %-10.2f\n",
            step, t, total_μ, target, r_now, te)
    flush(stdout)
end

@printf("\nTotal: %.1fs  avg: %.2fs/step\n", t_total[], t_total[]/nsteps)

# Final means
mu_final = reshape(real(expect(psi, "N")), Kx, Ky)
@printf("\nFinal mean grid (centre 5×5 subgrid, rows %d-%d):\n",
        Kx÷2-1, Kx÷2+3)
for i in Kx÷2-1:Kx÷2+3
    @printf("  [%s]\n", join([@sprintf("%5.3f", mu_final[i,j]) for j in Ky÷2-2:Ky÷2+2], " "))
end
@printf("\nTotal mean: %.4f  |  max bond r=%d\n",
        sum(mu_final), maxlinkdim(psi))
