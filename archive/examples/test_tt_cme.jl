using DiscStochSim, ExponentialUtilities, Printf, LinearAlgebra

# Parameters
k_b = 1.0; k_d = 0.5   # μ_ss = 2
D   = 2.0               # diffusion rate D/dx²
n_max = 10
dt  = 0.3; T = 3.0
nsteps = round(Int, T/dt)

L_local = birth_death_generator(k_b, k_d, n_max)

# Analytical mean for K voxels, uniform birth-death, IC = 1 molecule at voxel k0
# d⟨N_total⟩/dt = K*k_b - k_d*⟨N_total⟩  →  ⟨N_total⟩(T) = K*μ_ss*(1-e^{-k_d T}) + N0*e^{-k_d T}
analytical_total(K, T) = K * (k_b/k_d) * (1 - exp(-k_d*T)) + 1 * exp(-k_d*T)

println("=== TT-CME for birth-death RDME ===")
println("k_b=$k_b  k_d=$k_d  D=$D  n_max=$n_max  T=$T  dt=$dt")
println("μ_ss = $(k_b/k_d)")
println()

# ── 1D test: K=6, vary bond dimension ──────────────────────────────────────────
K = 6
ana = analytical_total(K, T)
@printf("1D chain K=%d  (analytical total mean = %.3f)\n\n", K, ana)
@printf("%-6s  %-8s  %-10s  %-10s  %-40s\n",
        "max_r", "r_actual", "total_μ", "err_vs_ODE", "per-voxel means")
@printf("%-6s  %-8s  %-10s  %-10s  %-40s\n", "─"^5, "─"^7, "─"^9, "─"^9, "─"^38)

for max_r in [1, 2, 4, 8, 16]
    tt = CMETensorTrain(K, n_max, max_r;
                        ic = CartesianIndex(0,0,0,1,0,0))
    for _ in 1:nsteps
        tt_step_1d!(tt, L_local, D, dt; krylov_m=25, tol=1e-14)
    end
    μ = tt_means(tt)
    r_act = tt_max_bond(tt)
    @printf("%-6d  %-8d  %-10.4f  %-10.4f  [%s]\n",
            max_r, r_act, sum(μ), abs(sum(μ) - ana),
            join([@sprintf("%.3f", m) for m in μ], " "))
end
println()

# ── 1D test: show norm evolution (should stay ≈1 throughout) ──────────────────
println("Norm evolution (r=4, should stay ≈1.0 if truncation error is small):")
tt_check = CMETensorTrain(K, n_max, 4; ic = CartesianIndex(0,0,0,1,0,0))
for step in 1:nsteps
    tt_step_1d!(tt_check, L_local, D, dt; krylov_m=25, tol=1e-14)
    step % 3 == 0 && @printf("  t=%.1f  norm=%.6f  r_max=%d\n",
                               step*dt, tt_norm(tt_check), tt_max_bond(tt_check))
end
println()

# ── 2D test: Kx=6, Ky=6 — compare r=1 (product) vs r>1 ────────────────────────
Kx = 6; Ky = 6; M = Kx*Ky
ana2d = analytical_total(M, T)
println("2D grid $(Kx)×$(Ky)  (analytical total mean = $(round(ana2d, digits=3)))")
println()

# IC: 1 molecule at centre voxel (row 3, col 3) in 0-indexed = index (3-1)*6+3=15
ic_2d = [zeros(n_max+1) for _ in 1:M]
for k in 1:M; ic_2d[k][1] = 1.0; end
centre = (Kx÷2)*Ky + (Ky÷2) + 1   # centre in row-major indexing
ic_2d[centre][1] = 0.0; ic_2d[centre][2] = 1.0

@printf("%-6s  %-8s  %-12s  %-12s\n", "max_r", "r_actual", "total_μ", "err_vs_ODE")
@printf("%-6s  %-8s  %-12s  %-12s\n", "─"^5, "─"^7, "─"^11, "─"^11)

for max_r in [1, 2, 4, 8]
    tt2d = CMETensorTrain(M, n_max, max_r; ic = ic_2d)
    t0 = time()
    for _ in 1:nsteps
        tt_step_2d!(tt2d, L_local, D, dt, Kx, Ky; krylov_m=20, tol=1e-14)
    end
    elapsed = time() - t0
    μ2d   = tt_means(tt2d)
    r_act = tt_max_bond(tt2d)
    @printf("%-6d  %-8d  %-12.4f  %-12.4f  (%.1fs)\n",
            max_r, r_act, sum(μ2d), abs(sum(μ2d) - ana2d), elapsed)
end

println()
println("Key result: cost = O(M × r² × n_max²) — LINEAR in M for fixed r")
println("For M=$(M), r=4, n_max=$n_max: ≈$(M * 4^2 * (n_max+1)^2 ÷ 1000)K operations/sweep")
println("For M=40000 (200×200), r=4, n_max=10: ≈$(40000 * 16 * 121 ÷ 1_000_000)M ops/sweep")
