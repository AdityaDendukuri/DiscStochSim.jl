"""
Policy validation for K=12 adaptive level assignment.

Runs the K=12 adaptive mixed-FSP (same parameters as compare_k12_levels Run B)
and after every time step reports per-pair diagnostics:
  μ_L, μ_R   : voxel means (left and right voxel of each pair)
  α_intra    : |μ_L − μ_R| / (μ_L + μ_R)  within-pair asymmetry
  β_left/right : inter-pair boundary gradient
  α_eff      : max(α_intra, β_left, β_right)
  level      : 0=fine, 1=pair-coarse, 2=quad-coarse

The check: front-adjacent pairs must have α_eff above the α_lo=0.10 threshold
(preventing level-2 promotion), while bulk pairs must have α_eff < α_lo2=0.03
(permitting level-2 promotion).

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/audit_k12_policy.jl
"""

using DiscStochSim
using Printf

D = 0.005; K = 12; h = 1.0/K
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_buf = 30; n̄_low = 2*n_low; n̄_high = 2*n_high
n_max_coarse = 2*(n_high + n_buf); n_max_fine = n_high + n_buf

levels_init = [1,1,0,1,1,1]
off0 = DiscStochSim._mixed_offsets(levels_init)
op3  = off0[3]   # fine pair 3 offset
KM0  = DiscStochSim._mixed_dim(levels_init)

u0a = CartesianIndex(ntuple(i -> i==1 ? n̄_low  :
                                 i==2 ? n̄_low  :
                                 i==op3   ? n_low  :
                                 i==op3+1 ? n_high :
                                 i <= KM0 ? n̄_high : 0, KM0)...)
u0b = CartesianIndex(ntuple(i -> i==1 ? n̄_low  :
                                 i==2 ? n̄_low  :
                                 i==op3   ? n_high :
                                 i==op3+1 ? n_low  :
                                 i <= KM0 ? n̄_high : 0, KM0)...)

sp = StateSpace{CartesianIndex{KM0}, Float64}()
add_state!(sp, u0a, 0.5); add_state!(sp, u0b, 0.5)
levels = copy(levels_init)

mxsys = build_schlogl_mixed_system(model, g, levels)
function bc_mix(s)
    t = Tuple(s)
    t[1] ≤ n_max_coarse && t[2] ≤ n_max_coarse &&
    t[op3]   ≤ n_max_fine  &&
    t[op3+1] ≤ n_max_fine  &&
    all(t[k] ≤ n_max_coarse for k in 5:length(t))
end
expand!(sp, mxsys, bc_mix; depth=1)

cache = MixedSolverCache()
cache.mixed_systems[levels] = mxsys
sc = Ref(0)

# ── Diagnostic helper ─────────────────────────────────────────────────────────
"""Compute per-pair diagnostics from the current mixed state."""
function pair_diagnostics(sp, levels, model, grid)
    n_pairs = length(levels)
    K2      = 2 * n_pairs
    off     = DiscStochSim._mixed_offsets(levels)
    d_fine  = diffusion_rate(model.D, grid)

    # Voxel means (same reconstruction as select_coarsening_mask_mixed)
    μ = zeros(K2)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        j = 1
        while j <= n_pairs
            op = off[j]; l = levels[j]
            if l == 0
                μ[2j-1] += p * t[op]; μ[2j] += p * t[op+1]; j += 1
            elseif l == 1
                h = p * t[op] / 2; μ[2j-1] += h; μ[2j] += h; j += 1
            else
                q = p * t[op] / 4
                μ[2j-1] += q; μ[2j] += q; μ[2(j+1)-1] += q; μ[2(j+1)] += q
                j += 2
            end
        end
    end

    # α_intra per pair
    α = zeros(n_pairs)
    for j in 1:n_pairs
        tot = μ[2j-1] + μ[2j]
        α[j] = tot > 0.0 ? abs(μ[2j-1] - μ[2j]) / tot : 0.0
    end

    # β inter-pair boundary gradient
    β = zeros(n_pairs - 1)
    for j in 1:n_pairs-1
        tot = μ[2j] + μ[2j+1]
        β[j] = tot > 0.0 ? abs(μ[2j] - μ[2j+1]) / tot : 0.0
    end

    α_eff = [max(α[j], j > 1 ? β[j-1] : 0.0, j < n_pairs ? β[j] : 0.0)
             for j in 1:n_pairs]

    μ, α, β, α_eff
end

# Thresholds (must match two_level_step_mixed defaults)
α_lo  = 0.10   # below this → eligible for level-1 (hysteresis lower)
α_hi  = 0.25   # above this → forced fine (level-0)
α_lo2 = 0.03   # below this AND already level-1 → eligible for level-2

println("α thresholds:  α_lo=$(α_lo)  α_hi=$(α_hi)  α_lo2=$(α_lo2)  (for level-2 promotion)")
println()

# ── Column header ─────────────────────────────────────────────────────────────
println("─"^88)
@printf("%-5s  %-7s  %-7s  %-7s  %-8s  %-8s  %-8s  %-8s  %s\n",
        "step", "pair", "μ_L", "μ_R", "α_intra", "β_left", "β_right", "α_eff", "level")
println("─"^88)

# ── Time-stepping loop ────────────────────────────────────────────────────────
n_steps = 5
for step in 1:n_steps
    global sp, levels, cache

    sp, levels = two_level_step_mixed(sp, levels, model, g, 0.1;
                     cache               = cache,
                     step_count          = sc,
                     mask_check_interval = 5,
                     windowed            = true,
                     halo                = 0,
                     expand_mixed        = true,
                     mixed_expand_depth  = 1,
                     mixed_n_max         = n_max_coarse,
                     prune_tol           = 0.0,
                     reexpand_depth      = 0,
                     krylov_m            = 40)
    sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)

    μ, α, β, α_eff = pair_diagnostics(sp, levels, model, g)
    n_pairs = length(levels)

    for j in 1:n_pairs
        β_l = j > 1       ? β[j-1] : NaN
        β_r = j < n_pairs ? β[j]   : NaN
        lv  = levels[j]
        # flag: '**' if level-2, '!!' if α_eff > α_hi (should be fine)
        flag = lv == 2 ? "← L2" :
               lv == 0 ? "← fine" :
               α_eff[j] < α_lo2 ? "(L2 eligible)" : ""
        @printf("%-5d  pair%-2d   %-7.1f  %-7.1f  %-8.4f  %-8s  %-8s  %-8.4f  L%d %s\n",
                step, j,
                μ[2j-1], μ[2j],
                α[j],
                isnan(β_l) ? "  —  " : @sprintf("%.4f", β_l),
                isnan(β_r) ? "  —  " : @sprintf("%.4f", β_r),
                α_eff[j], lv, flag)
    end
    @printf("  |S|=%d\n\n", length(sp))
end
println("─"^88)

# ── Summary ───────────────────────────────────────────────────────────────────
println("\nFinal levels = $(levels)")
μ_f, α_f, β_f, α_eff_f = pair_diagnostics(sp, levels, model, g)
n_pairs = length(levels)

level2_pairs = [j for j in 1:n_pairs if levels[j] == 2]
fine_pairs   = [j for j in 1:n_pairs if levels[j] == 0]
l1_pairs     = [j for j in 1:n_pairs if levels[j] == 1]

println("\nLevel summary:")
@printf("  fine   (L0): pairs %s  α_eff = [%s]\n",
        string(fine_pairs), join([@sprintf("%.4f",α_eff_f[j]) for j in fine_pairs], ", "))
@printf("  coarse (L1): pairs %s  α_eff = [%s]\n",
        string(l1_pairs), join([@sprintf("%.4f",α_eff_f[j]) for j in l1_pairs], ", "))
@printf("  quad   (L2): pairs %s  α_eff = [%s]\n",
        string(level2_pairs), join([@sprintf("%.4f",α_eff_f[j]) for j in level2_pairs], ", "))

println()
for j in level2_pairs
    @printf("  ✓ Pair %d (L2): α_eff=%.4f < α_lo2=%.2f\n", j, α_eff_f[j], α_lo2)
end
if !isempty(fine_pairs)
    println("  Fine-pair neighbours:")
    for j in fine_pairs
        for nb in [j-1, j+1]
            1 ≤ nb ≤ n_pairs || continue
            @printf("    Pair %d (L%d, adjacent to fine pair %d): α_eff=%.4f  [α_lo=%.2f → level-2 blocked]\n",
                    nb, levels[nb], j, α_eff_f[nb], α_lo)
        end
    end
end
