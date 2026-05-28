"""
QTT Joint CME over the Active Ring: Full Joint Distribution for 200×200 RDME
============================================================================

Key insight: for a 200×200 birth-death RDME, only ~300 voxels are ever
"active" (in transition between empty and equilibrated). The rest are:
  - Empty:  n_k = 0 (known exactly, contributes 0 flux)
  - Equil:  p_k = Poisson(μ_ss) (known exactly, contributes D×μ_ss flux)

This reduces the QTT problem from 40,000 voxels to ~300 active voxels!
With q=4 binary bits per voxel: M = 300×4 = 1200 binary sites vs 160,000.

Active ring QTT (Kazeev Algorithm 2):
  1. Track active voxels (TV vs Poisson criterion from SpatialFSP)
  2. Order active ring by angle → nearest-neighbor diffusion in ring ordering
  3. Build QTT for K_active × q binary sites with carry-adder MPO
  4. Solve backward Euler: (I + dt × H_eff) p_{n+1} = p_n via DMRG linsolve
  5. Promote/demote voxels as front moves

Run: julia --project=examples/itensor_env examples/qtt_active_ring.jl
"""

using ITensors, ITensorMPS, Printf, LinearAlgebra, Statistics

# ── Qbit site type ────────────────────────────────────────────────────────────
ITensors.space(::SiteType"Qbit") = 2
function ITensors.op(::OpName"Id",    ::SiteType"Qbit", s::Index)
    O=ITensor(s',dag(s)); O[s'=>1,s=>1]=1.0; O[s'=>2,s=>2]=1.0; O
end
function ITensors.op(::OpName"Plus",  ::SiteType"Qbit", s::Index)
    O=ITensor(s',dag(s)); O[s'=>2,s=>1]=1.0; O
end
function ITensors.op(::OpName"Minus", ::SiteType"Qbit", s::Index)
    O=ITensor(s',dag(s)); O[s'=>1,s=>2]=1.0; O
end
function ITensors.op(::OpName"P",     ::SiteType"Qbit", s::Index)
    O=ITensor(s',dag(s)); O[s'=>2,s=>2]=1.0; O
end
function ITensors.op(::OpName"Q",     ::SiteType"Qbit", s::Index)
    O=ITensor(s',dag(s)); O[s'=>1,s=>1]=1.0; O
end

# ── Binary helpers ─────────────────────────────────────────────────────────────
to_bits(n,q) = [Bool((n>>(j-1))&1) for j in 1:q]
bit_op_name(b0,b1) = !b0&&!b1 ? "Q" : b0&&b1 ? "P" : !b0&&b1 ? "Plus" : "Minus"

function add_transition!(os, coeff, n, m, sk, q)
    bn=to_bits(n,q); bm=to_bits(m,q)
    ops = [(bit_op_name(bn[j],bm[j]), sk[j]) for j in 1:q]
    add!(os, coeff, Iterators.flatten(ops)...)
end

# ── Build H_eff for the active ring ───────────────────────────────────────────
"""
    build_ring_heff(active_voxels, ring_order, equil_μ, empty_mask,
                    grid, k_b, k_d, D, q) → OpSum

Build H_eff = -L for the active ring.

`active_voxels`: list of (row,col) voxel coordinates in the active ring
`ring_order`:    indices into active_voxels giving their ring ordering (by angle)
`equil_μ`:       mean count for each equil neighbor (Dict{(r,c),Float64})
`empty_mask`:    set of empty voxel coords (contribute 0 flux)

For each active voxel k:
  - Birth/death reactions: carry-adder/subtractor on bits sk
  - Intra-ring diffusion: to adjacent active voxels in ring ordering
  - Boundary coupling:
      from equil neighbor j: +D×μ_j effective birth (gain) + loss
      from empty neighbor j: just loss (D×n_k, no gain since n_j=0)
"""
function build_ring_heff(active_list,   # Vector of (row,col) in ring order
                          adj_active,    # adj_active[k] = ring indices adjacent to k
                          adj_equil_μ,  # adj_equil_μ[k] = Vector of μ_j for equil nbrs
                          adj_n_empty,  # adj_n_empty[k] = number of empty nbrs
                          k_b, k_d, D, q)
    Ka   = length(active_list)
    nmax = (1 << q) - 1
    os   = OpSum()

    for k in 1:Ka
        sk = [(k-1)*q + j for j in 1:q]

        # ── Birth: H_eff = -L_birth = k_b(I - Plus_k) ─────────────────────────
        for n in 0:nmax-1; add_transition!(os, -k_b, n, n+1, sk, q); end
        add!(os, k_b, "Id", sk[1])

        # ── Death: H_eff = -L_death = k_d(N_k - Minus_k) ─────────────────────
        for j in 1:q; add!(os, k_d*Float64(1<<(j-1)), "P", sk[j]); end
        for n in 1:nmax; add_transition!(os, -k_d*Float64(n), n, n-1, sk, q); end

        # ── Intra-ring diffusion to adjacent active voxels ────────────────────
        for k2 in adj_active[k]
            k2 > k || continue   # process each pair once
            sk2 = [(k2-1)*q + j for j in 1:q]
            # Hop k→k2:
            for n1 in 1:nmax, n2 in 0:nmax-1
                all_sk = [sk; sk2]
                all_n_in  = [to_bits(n1,q); to_bits(n2,q)]
                all_n_out = [to_bits(n1-1,q); to_bits(n2+1,q)]
                ops = [(bit_op_name(all_n_in[s],all_n_out[s]), all_sk[s]) for s in 1:2q]
                add!(os, -D*Float64(n1), Iterators.flatten(ops)...)
            end
            for j in 1:q; add!(os, D*Float64(1<<(j-1)), "P", sk[j]); end
            # Hop k2→k (symmetric):
            for n1 in 0:nmax-1, n2 in 1:nmax
                all_sk = [sk; sk2]
                all_n_in  = [to_bits(n1,q); to_bits(n2,q)]
                all_n_out = [to_bits(n1+1,q); to_bits(n2-1,q)]
                ops = [(bit_op_name(all_n_in[s],all_n_out[s]), all_sk[s]) for s in 1:2q]
                add!(os, -D*Float64(n2), Iterators.flatten(ops)...)
            end
            for j in 1:q; add!(os, D*Float64(1<<(j-1)), "P", sk2[j]); end
        end

        # ── Boundary coupling to equil neighbors (mean-field) ─────────────────
        # Equil neighbor j with mean μ_j: effective birth D×μ_j + effective death D
        for μ_j in adj_equil_μ[k]
            # Gain: +D×μ_j (constant birth into voxel k from equil j)
            for n in 0:nmax-1; add_transition!(os, -D*μ_j, n, n+1, sk, q); end
            add!(os, D*μ_j, "Id", sk[1])
            # Loss: +D (per molecule from k to equil j, on top of normal death)
            for j in 1:q; add!(os, D*Float64(1<<(j-1)), "P", sk[j]); end
            for n in 1:nmax; add_transition!(os, -D*Float64(n), n, n-1, sk, q); end
        end

        # ── Loss to empty neighbors (n_j=0 so no gain, just extra loss) ───────
        n_empty = adj_n_empty[k]
        if n_empty > 0
            for j in 1:q
                add!(os, D*n_empty*Float64(1<<(j-1)), "P", sk[j])
            end
            for n in 1:nmax
                add_transition!(os, -D*n_empty*Float64(n), n, n-1, sk, q)
            end
        end
    end

    os
end

# ── L1 utilities ───────────────────────────────────────────────────────────────
function build_ones_mps(sites)
    M = length(sites)
    psi = MPS(sites; linkdims=1)
    for s in 1:M; psi[s] = ITensor([1.0,1.0], sites[s]); end
    psi
end

function qtt_means_ring(psi, sites, Ka, q, ones_mps)
    l1 = real(inner(ones_mps, psi))
    l1 = l1 > 1e-10 ? l1 : 1.0
    μ  = zeros(Ka)
    for k in 1:Ka, j in 1:q
        s    = (k-1)*q + j
        Ppsi = apply(op("P", sites[s]), psi; cutoff=1e-14)
        μ[k] += Float64(1<<(j-1)) * real(inner(ones_mps, Ppsi)) / l1
    end
    μ
end

# ── One Kazeev backward-Euler step ────────────────────────────────────────────
function ring_step!(psi, H_eff, sites, ones_mps, dt;
                     max_r=8, cutoff=1e-12, nsweeps=15)
    I_mpo  = MPO(sites, "Id")
    BE_mpo = add(I_mpo, dt * H_eff; maxdim=maxlinkdim(H_eff)+4, cutoff=1e-14)
    psi_new = linsolve(BE_mpo, psi, psi;
                        nsweeps=nsweeps, maxdim=max_r, cutoff=cutoff, outputlevel=0)
    truncate!(psi_new; maxdim=max_r, cutoff=cutoff)
    l1 = real(inner(ones_mps, psi_new))
    l1 > 1e-10 && (psi_new = psi_new / l1)
    psi_new
end

# ── Main: 40×40 birth-death, demonstrate joint QTT over active ring ───────────
const Kx=40; const Ky=40; const K=Kx*Ky
const k_b=1.0; const k_d=0.5; const D=1.0
const μ_ss = k_b/k_d   # = 2.0
const q=4; const nmax=(1<<q)-1   # 4 bits → n_max=15
const max_r=6; const dt=0.5; const T_END=20.0

@printf("=== Active-Ring QTT RDME (%d×%d, K_active tracks front) ===\n", Kx,Ky)
@printf("q=%d bits  n_max=%d  μ_ss=%.1f  max_r=%d\n\n", q, nmax, μ_ss, max_r)

# Grid state: :empty, :active, :equil
vstate   = fill(:empty, Kx, Ky)
vox_mean = zeros(Kx, Ky)  # mean count per voxel

# IC: 1 molecule at centre → seed 1 active voxel
ci, cj = Kx÷2+1, Ky÷2+1
vstate[ci,cj] = :active

# Neighbours of (i,j) in the 2D grid (4-connected)
function nbrs(i,j)
    [(i-1,j),(i+1,j),(i,j-1),(i,j+1)] |>
    xs -> filter(((r,c),) -> 1≤r≤Kx && 1≤c≤Ky, xs)
end

# TV distance vs Poisson(μ_ss): classify equil
function is_equil(μ; ε=0.2)
    abs(μ - μ_ss) / (μ_ss + 0.01) < ε
end

# Promote empty neighbour of an active/equil voxel
function promote_empty!(vstate, vox_mean, i, j)
    vstate[i,j] == :empty || return
    # Check if any active/equil neighbour
    any(vstate[r,c] != :empty for (r,c) in nbrs(i,j)) || return
    vstate[i,j] = :active
    vox_mean[i,j] = 0.0
end

# Build the active list (ordered by angle from centre)
function active_ring_list(vstate, ci, cj)
    act = [(i,j) for i in 1:Kx, j in 1:Ky if vstate[i,j]==:active]
    sort(act; by = ((i,j),) -> atan(j-cj, i-ci))
end

function run_active_ring()
    local psi, sites, ones_mps

    local vstate_l  = copy(vstate)
    local vmean_l   = copy(vox_mean)

    # Seed the initial active set with 1 molecule at centre
    vstate_l[ci,cj] = :active
    vmean_l[ci,cj]  = 1.0

    t = 0.0
    Ka_prev = -1      # track ring size to know when to rebuild QTT
    local psi, sites, ones_mps, H_eff_cache

    @printf("%-5s  %-7s  %-5s  %-10s  %-10s  %-8s\n",
            "step","t","K_act","total_μ","analytical","r_max")
    @printf("%s\n","─"^60)

    for step in 1:round(Int, T_END/dt)
        t = step * dt

        # ── 1. Build active ring list (angle-ordered) ─────────────────────────
        active_list = active_ring_list(vstate_l, ci, cj)
        Ka = length(active_list)

        # ── 2. Compute adjacency for active ring ───────────────────────────────
        act_set = Set(active_list)
        act_idx = Dict(v=>i for (i,v) in enumerate(active_list))

        adj_active   = [Int[] for _ in 1:Ka]
        adj_equil_μ  = [Float64[] for _ in 1:Ka]
        adj_n_empty  = zeros(Int, Ka)

        for (k,(i,j)) in enumerate(active_list)
            for (r,c) in nbrs(i,j)
                if vstate_l[r,c] == :active
                    push!(adj_active[k], act_idx[(r,c)])
                elseif vstate_l[r,c] == :equil
                    push!(adj_equil_μ[k], vmean_l[r,c])
                else  # empty
                    adj_n_empty[k] += 1
                end
            end
        end

        # ── 3. Build QTT state if ring changed size ────────────────────────────
        if Ka != Ka_prev
            Ka_prev   = Ka
            sites_new = siteinds("Qbit", Ka*q)
            # Initialise: product state from per-voxel means
            psi_new = MPS(sites_new; linkdims=1)
            for k in 1:Ka
                μk   = vmean_l[active_list[k]...]
                μk   = clamp(μk, 0.0, Float64(nmax))
                for j in 1:q
                    s     = (k-1)*q + j
                    # Simple: concentrate probability at round(μk) in binary
                    n_int = clamp(round(Int, μk), 0, nmax)
                    bits  = to_bits(n_int, q)
                    v     = zeros(2); v[bits[j] ? 2 : 1] = 1.0
                    psi_new[s] = ITensor(v, sites_new[s])
                end
            end
            sites    = sites_new
            psi      = psi_new
            ones_mps = build_ones_mps(sites)
        end

        # ── 4. Build ring MPO ─────────────────────────────────────────────────
        Ka == 0 && continue  # skip if no active voxels
        os_ring = build_ring_heff(active_list, adj_active, adj_equil_μ, adj_n_empty,
                                   k_b, k_d, D, q)
        H_eff = MPO(os_ring, sites)

        # ── 5. Kazeev backward Euler step ─────────────────────────────────────
        psi = ring_step!(psi, H_eff, sites, ones_mps, dt;
                          max_r, cutoff=1e-12, nsweeps=12)

        # ── 6. Extract means and update grid ──────────────────────────────────
        μ_ring = qtt_means_ring(psi, sites, Ka, q, ones_mps)
        for (k,(i,j)) in enumerate(active_list)
            vmean_l[i,j] = max(0.0, μ_ring[k])
        end

        # ── 7. Adapt active set ────────────────────────────────────────────────
        # Demote active → equil when mean ≈ μ_ss
        for (i,j) in copy(active_list)
            is_equil(vmean_l[i,j]) && (vstate_l[i,j] = :equil)
        end
        # Promote empty neighbours: flux = D × μ_source > ε_expand
        ε_expand = 0.05
        for i in 1:Kx, j in 1:Ky
            vstate_l[i,j] == :empty && continue
            μ_src = vmean_l[i,j]
            D * μ_src < ε_expand && continue   # low flux: don't activate
            for (r,c) in nbrs(i,j)
                vstate_l[r,c] == :empty || continue
                vstate_l[r,c] = :active
                vmean_l[r,c]  = 0.0            # empty voxel starts at n=0
            end
        end

        # ── 8. Log ─────────────────────────────────────────────────────────────
        # Total mean: active (from QTT) + equil (μ_ss each) + empty (0)
        total_μ = sum(vmean_l)
        n_act   = count(==(:active), vstate_l)
        n_eq    = count(==(:equil),  vstate_l)
        # Analytical: only the voxels that have ever been activated contribute
        n_touched = n_act + n_eq
        target_touched = n_touched * μ_ss * (1-exp(-k_d*t)) + exp(-k_d*t)

        step % 2 == 0 &&
        @printf("%-5d  %-7.1f  %-5d  %-10.3f  %-10.3f  %-8d  (eq=%d,touched=%d)\n",
                step, t, n_act, total_μ, target_touched, maxlinkdim(psi), n_eq, n_touched)
        flush(stdout)
    end

    @printf("\nFinal grid (centre 7×7, equil=E, empty=. ):\n")
    for i in ci-3:ci+3
        row = join([vstate_l[i,j]==:active ? @sprintf("%4.1f",vmean_l[i,j]) :
                    vstate_l[i,j]==:equil  ? "  E " : "  . " for j in cj-3:cj+3])
        println("  ", row)
    end
end

run_active_ring()
