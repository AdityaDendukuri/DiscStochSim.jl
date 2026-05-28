using DiscStochSim, Printf, Random
Random.seed!(1)
model = SchloglModel1D(2.0, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("Fixed pts: n_lo=%d n_un=%d n_hi=%d\n", n_lo, n_un, n_hi)

D = 2.0
in_rate   = D * (n_hi + 3*n_lo)
out_coeff = D * 4
@printf("in_rate=%.1f  out_coeff=%.1f  MF fixed pt=%.1f\n",
        in_rate, out_coeff, in_rate/out_coeff)

p = zeros(221); p[n_lo+1] = 1.0
krylov_m = 30; dt = 0.5
voxel_mean_f(q) = sum((n-1)*q[n] for n in 1:length(q))
p_hi_f(q) = sum(q[n_un+2:end])

for step in 1:10
    max_rate = in_rate + out_coeff * n_hi
    n_sub = max(1, ceil(Int, max_rate * dt / krylov_m))
    dt_sub = dt / n_sub
    for _ in 1:n_sub
        global p = DiscStochSim._expv_voxel(model, in_rate, out_coeff, p, dt_sub, krylov_m; n_include=n_hi)
        p .= max.(p, 0.0); p ./= sum(p)
    end
    mu = voxel_mean_f(p); phi = p_hi_f(p)
    @printf("  step %2d (t=%.1f): mean=%6.1f  P(hi)=%.3f  p_nun=%.3f  snap4=%-5s snap40=%-5s\n",
            step, step*dt, mu, phi, p[n_un+1],
            abs(mu-n_un)>4 ? "YES" : "no",
            abs(mu-n_un)>40 ? "YES" : "no")
end
