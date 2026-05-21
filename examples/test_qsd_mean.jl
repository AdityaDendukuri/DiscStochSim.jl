using DiscStochSim, LinearAlgebra, Printf

k_b=0.5; k_d=0.1; D=1.0; n_max=14
model = BranchingDeathRDME(k_b, k_d, D, 1.0)

# Isolated voxel QSD
L_iso = DiscStochSim._single_voxel_generator(model, 0.0, 0.0, n_max)
p_iso = abs.(nullspace(Matrix(L_iso)')[:,1]); p_iso ./= sum(p_iso)
mu_iso = sum((0:n_max) .* p_iso)
@printf("Isolated voxel QSD mean: %.2f (%.1f%% of n_max=%d)\n", mu_iso, 100*mu_iso/n_max, n_max)

# Interior voxel with 4 equil neighbours — self-consistent mean
mu_sat = mu_iso
for _ in 1:100
    L2 = DiscStochSim._single_voxel_generator(model, 4*D*mu_sat, 4*D, n_max)
    p2 = abs.(nullspace(Matrix(L2)')[:,1]); p2 ./= sum(p2)
    mu_new = sum((0:n_max) .* p2)
    abs(mu_new - mu_sat) < 1e-8 && break
    mu_sat = 0.5*mu_sat + 0.5*mu_new   # damped iteration
end
@printf("Interior QSD mean (self-consistent): %.2f (%.1f%% of n_max)\n",
        mu_sat, 100*mu_sat/n_max)
@printf("Good threshold: >= %.1f%% × n_max = %.1f\n", 70.0, 0.70*n_max)
