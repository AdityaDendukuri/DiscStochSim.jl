using DiscStochSim, ExponentialUtilities, Printf
# Test: product convolution prolongation in V-cycle vs direct expv for Ka=4

model  = SchloglModel1D(1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
@printf("n_lo=%d  n_un=%d  n_hi=%d\n\n", n_lo, n_un, n_hi)

K = 20; dx = 1.0; n_max = 200
Ka = 4  # fixed active strip — V-cycle will trigger

# IC: SYMMETRIC interface — equal prob for (hi,hi,lo,lo) and (lo,lo,hi,hi)
sp0 = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp0, CartesianIndex(n_hi, n_hi, n_lo, n_lo), 0.5)
add_state!(sp0, CartesianIndex(n_lo, n_lo, n_hi, n_hi), 0.5)

# Mean-field for the full grid
μ = [i <= K÷2 ? Float64(n_hi) : Float64(n_lo) for i in 1:K]
μ[K÷2-1:K÷2+2] .= (n_hi + n_lo) / 2.0  # transition zone

k_lo = K÷2 - 1; k_hi = k_lo + Ka - 1
lbc = k_lo - 1; rbc = k_hi + 1
ga  = VoxelGrid(Ka, dx, 0)
sys = build_active_schlogl_1d(model, ga, lbc, rbc, μ, Val(Ka))
bc  = schlogl_mixed_bc(n_max)
expand!(sp0, sys, bc; depth=1)

@printf("Initial state: %d states  p_sum=%.4f\n\n", length(sp0), sum(sp0.probs))

# ── Direct expv (reference) ────────────────────────────────────────────────────
sp_dir = deepcopy(sp0)
A, = build_generator(sp_dir, sys, Float64[], 0.0)
sp_dir.probs .= expv(0.5, A, sp_dir.probs; m=30)
renormalize!(sp_dir); prune_threshold!(sp_dir, 1e-8)

# ── V-cycle with product convolution ──────────────────────────────────────────
sp_vc  = deepcopy(sp0)
hier   = build_hierarchy(ga, 1)  # Ka=4 → K2=2 at coarse level
sp_vc  = multi_level_vcycle_schlogl(sp_vc, model, hier, lbc, rbc, μ, Float64[],
                                     0.0, 0.5;
                                     τ_pre = 0.05, τ_post = 0.05,
                                     krylov_m = 30,
                                     weight_tol = 1e-8, expand_coarse = true,
                                     coarse_r_depth = 2, coarse_d_depth = 1,
                                     coarse_n_max = 2*n_max, prune_tol = 1e-8)
renormalize!(sp_vc)

@printf("After dt=0.5:\n")
@printf("  Direct:   %d states  p_sum=%.4f\n", length(sp_dir), sum(sp_dir.probs))
@printf("  V-cycle:  %d states  p_sum=%.4f\n\n", length(sp_vc), sum(sp_vc.probs))

# Compare means per voxel
μ_dir = zeros(Ka); μ_vc = zeros(Ka)
for (ci,p) in zip(sp_dir.states, sp_dir.probs); t=Tuple(ci); for k in 1:Ka; μ_dir[k]+=p*t[k]; end; end
for (ci,p) in zip(sp_vc.states,  sp_vc.probs);  t=Tuple(ci); for k in 1:Ka; μ_vc[k] +=p*t[k]; end; end

@printf("Voxel means:\n")
@printf("  %-8s  %-10s  %-10s  %-10s\n", "voxel", "direct", "vcycle", "|Δμ|")
for k in 1:Ka
    @printf("  %-8d  %-10.2f  %-10.2f  %-10.4f\n", k, μ_dir[k], μ_vc[k], abs(μ_dir[k]-μ_vc[k]))
end

# Compare config probabilities P(m voxels > n_un)
p_m_dir = zeros(Ka+1); p_m_vc = zeros(Ka+1)
for (ci,p) in zip(sp_dir.states, sp_dir.probs); m=count(n->n>n_un, Tuple(ci)); p_m_dir[m+1]+=p; end
for (ci,p) in zip(sp_vc.states,  sp_vc.probs);  m=count(n->n>n_un, Tuple(ci)); p_m_vc[m+1] +=p; end

@printf("\nConfig probabilities P(m voxels in hi state):\n")
@printf("  %-4s  %-10s  %-10s  %-10s\n", "m", "direct", "vcycle", "|Δ|")
for m in 0:Ka
    (p_m_dir[m+1] > 1e-4 || p_m_vc[m+1] > 1e-4) || continue
    @printf("  %-4d  %-10.4f  %-10.4f  %-10.4f\n", m, p_m_dir[m+1], p_m_vc[m+1],
            abs(p_m_dir[m+1]-p_m_vc[m+1]))
end

@printf("\nTV distance in config probs: %.4f\n",
        0.5 * sum(abs(p_m_dir[m+1] - p_m_vc[m+1]) for m in 0:Ka))
