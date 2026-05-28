using ITensors, ITensorMPS, LinearAlgebra, Printf
# Test: QTT-CME single voxel K=1 vs exact FSP solution

# Qbit site type
ITensors.space(::SiteType"Qbit") = 2
for op_name in ("Id","Plus","Minus","P","Q")
    @eval function ITensors.op(::OpName{Symbol($op_name)}, ::SiteType"Qbit", s::Index)
        O = ITensor(s', dag(s))
        if $op_name == "Id";    O[s'=>1,s=>1]=1.0; O[s'=>2,s=>2]=1.0
        elseif $op_name=="Plus";  O[s'=>2,s=>1]=1.0
        elseif $op_name=="Minus"; O[s'=>1,s=>2]=1.0
        elseif $op_name=="P";     O[s'=>2,s=>2]=1.0
        elseif $op_name=="Q";     O[s'=>1,s=>1]=1.0
        end; O
    end
end

to_bits(n, q) = [Bool((n>>(j-1))&1) for j in 1:q]
bit_op_name(b0, b1) = !b0&&!b1 ? "Q" : b0&&b1 ? "P" : !b0&&b1 ? "Plus" : "Minus"

# K=1 voxel, q=4 bits
K=1; q=4; k_b=1.0; k_d=0.5; n_max=15
M = K*q
sites = siteinds("Qbit", M)

# Build H_eff = -L for K=1 voxel
os = OpSum()
sk = collect(1:q)
for n in 0:n_max-1
    bits_n = to_bits(n,q); bits_n1 = to_bits(n+1,q)
    ops = [(bit_op_name(bits_n[j], bits_n1[j]), sk[j]) for j in 1:q]
    add!(os, -k_b, Iterators.flatten(ops)...)  # -k_b × Plus_k
end
add!(os, k_b, "Id", sk[1])   # +k_b × I
for j in 1:q; add!(os, k_d * Float64(1<<(j-1)), "P", sk[j]); end  # +k_d × N
for n in 1:n_max
    bits_n = to_bits(n,q); bits_n_m1 = to_bits(n-1,q)
    ops = [(bit_op_name(bits_n[j], bits_n_m1[j]), sk[j]) for j in 1:q]
    add!(os, -k_d*Float64(n), Iterators.flatten(ops)...)  # -k_d×n × Minus
end
H_eff = MPO(os, sites)
@printf("H_eff MPO bond dim = %d\n", maxlinkdim(H_eff))

# Build all-ones MPS
ones_mps = MPS(sites; linkdims=1)
for s in 1:M; ones_mps[s] = ITensor([1.0,1.0], sites[s]); end

# IC: 1 molecule at voxel 1
psi = MPS(sites; linkdims=1)
bits_1 = to_bits(1, q)
for j in 1:q
    v = zeros(2); v[bits_1[j] ? 2 : 1] = 1.0
    psi[j] = ITensor(v, sites[j])
end
@printf("IC L1 norm: %.6f  (should be 1.0)\n", real(inner(ones_mps, psi)))

# Verify H_eff conservation: L1(H_eff × psi) should be 0
Hpsi = apply(H_eff, psi; maxdim=8, cutoff=1e-14)
@printf("L1(H_eff × psi) = %.2e  (should be ≈0 for CME conservation)\n",
        real(inner(ones_mps, Hpsi)))

# Build exact 16×16 CME generator L for comparison
L_exact = zeros(n_max+1, n_max+1)
for n in 0:n_max
    n < n_max && (L_exact[n+2, n+1] += k_b)  # birth gain
    L_exact[n+1, n+1] -= k_b                   # birth loss
    n > 0 && (L_exact[n, n+1] += k_d*n)        # death gain
    L_exact[n+1, n+1] -= k_d*n                  # death loss
end
@printf("L_exact eigenvalues (largest 3 magnitudes): %s\n",
        join([@sprintf("%.2f", v) for v in sort(abs.(eigvals(L_exact)), rev=true)[1:3]], ", "))

# Backward Euler: solve (I + dt × H_eff) p_new = p_old
dt = 0.5
I_mpo  = MPO(sites, "Id")
BE_mpo = add(I_mpo, dt * H_eff; maxdim=maxlinkdim(H_eff)+2, cutoff=1e-14)
@printf("BE_mpo bond dim = %d\n", maxlinkdim(BE_mpo))

# Check L1(BE × psi) = 1 (backward Euler preserves probability)
BEpsi = apply(BE_mpo, psi; maxdim=16, cutoff=1e-14)
@printf("L1(BE × psi) = %.6f  (should be 1.0 for prob. conservation)\n",
        real(inner(ones_mps, BEpsi)))

# linsolve
@printf("\nSolving (I + dt×H_eff) p_new = p_old via linsolve...\n")
psi_new = linsolve(BE_mpo, psi, psi; nsweeps=20, maxdim=8, cutoff=1e-12, outputlevel=0)
@printf("L1(p_new) = %.6f  (should be 1.0)\n", real(inner(ones_mps, psi_new)))

# QTT means for p_new
l1 = real(inner(ones_mps, psi_new))
μ_qtt = sum(Float64(1<<(j-1)) *
            real(inner(ones_mps, apply(op("P", sites[j]), psi_new; cutoff=1e-14))) / l1
            for j in 1:q)
@printf("QTT mean <n> = %.4f\n", μ_qtt)

# Exact solution: (I + dt×L_exact)^{-1} × p_0
p0_exact = zeros(n_max+1); p0_exact[2] = 1.0  # n=1
BE_exact = I + dt * (-L_exact)  # = I - dt × L
p_new_exact = BE_exact \ p0_exact
μ_exact = sum(n * p_new_exact[n+1] for n in 0:n_max)
@printf("Exact mean <n> = %.4f\n", μ_exact)
@printf("L1 of exact solution = %.6f\n", sum(p_new_exact))
