"""
Anisotropic expansion benchmark: K=12 adaptive run.

Demonstrates two effects of anisotropic per-coordinate bounds:
  1. CORRECTNESS: the old uniform `mixed_n_max` caps L2 quad coords at
     n_max_coarse = 384 even though L2 IC states sit at 4*n_high ≈ 648.
     This silently prevents any upward expansion of the L2 state space.
  2. TIGHTNESS: fine (L0) coords get a smaller floor (n_max_fine ≈ 192)
     instead of the inflated uniform cap, reducing spurious high-count states.

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/compare_anisotropic.jl
"""

using DiscStochSim
using Printf

D = 0.005; K = 12; h = 1.0/K
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_buf = 30
n_max_fine   = n_high + n_buf       # 192: floor per fine voxel
n_max_coarse = 2*(n_high + n_buf)   # 384: old uniform cap

@printf("n_high=%d  n_max_fine=%d  n_max_coarse (old uniform cap)=%d\n",
        n_high, n_max_fine, n_max_coarse)
@printf("L2 IC = 4*n_high = %d  →  old cap %d is BELOW the IC (truncation bug)\n\n",
        4*n_high, n_max_coarse)

n̄_low = 2*n_low; n̄_high = 2*n_high
levels_init = [1,1,0,1,1,1]
off0 = DiscStochSim._mixed_offsets(levels_init)
op3  = off0[3]; KM0 = DiscStochSim._mixed_dim(levels_init)

function make_sp()
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
    sp
end

function run_mode(label, anisotropic; n_steps=5)
    sp     = make_sp()
    levels = copy(levels_init)
    mxsys  = build_schlogl_mixed_system(model, g, levels)

    off_i = DiscStochSim._mixed_offsets(levels); op = off_i[3]
    function bc_init(s)
        t = Tuple(s)
        t[1] ≤ n_max_coarse && t[2] ≤ n_max_coarse &&
        t[op] ≤ n_max_fine && t[op+1] ≤ n_max_fine &&
        all(t[k] ≤ n_max_coarse for k in 5:length(t))
    end
    expand!(sp, mxsys, bc_init; depth=1)

    cache = MixedSolverCache(); cache.mixed_systems[levels] = mxsys
    sc = Ref(0); sizes = Int[]

    for _ in 1:n_steps
        sp, levels = two_level_step_mixed(sp, levels, model, g, 0.1;
                         cache               = cache,
                         step_count          = sc,
                         mask_check_interval = 5,
                         windowed            = true,
                         halo                = 0,
                         expand_mixed        = true,
                         mixed_expand_depth  = 1,
                         mixed_n_max         = n_max_coarse,
                         anisotropic_expand  = anisotropic,
                         n_max_per_voxel     = n_max_fine,
                         c_sigma             = 6.0,
                         prune_tol           = 0.0,
                         reexpand_depth      = 0,
                         krylov_m            = 40)
        sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
        push!(sizes, length(sp))
    end
    sp, levels, sizes
end

# ── Run both modes ────────────────────────────────────────────────────────────
println("── |S| per step ─────────────────────────────────────────────────────────")
@printf("%-20s  %s\n", "mode", join([@sprintf("step %d", s) for s in 1:5], "  "))
println("─"^60)

sp_uni,   lev_uni,   sz_uni   = run_mode("uniform (old)",  false)
sp_aniso, lev_aniso, sz_aniso = run_mode("anisotropic",    true)

@printf("%-20s  %s\n", "uniform (old)",
        join([@sprintf("%6d", s) for s in sz_uni], "  "))
@printf("%-20s  %s\n", "anisotropic",
        join([@sprintf("%6d", s) for s in sz_aniso], "  "))

# ── Per-coordinate bounds at step 5 ──────────────────────────────────────────
println("\n── Per-coordinate bounds at t=0.5 ──────────────────────────────────────")
@printf("%-8s  %-8s  %-8s  %-12s  %-12s  %s\n",
        "pair", "level", "coord", "μ", "uniform cap", "aniso bound")
println("─"^72)

bnds = mixed_coord_bounds(sp_aniso, lev_aniso; n_max_per_voxel=n_max_fine, c_sigma=6.0)
off5 = DiscStochSim._mixed_offsets(lev_aniso)
n_pairs = length(lev_aniso)

# compute coord means from aniso run
μ_coord = zeros(length(bnds))
for (i, s) in enumerate(sp_aniso.states)
    t = Tuple(s); p = sp_aniso.probs[i]
    for k in eachindex(μ_coord); μ_coord[k] += p * t[k]; end
end

let j2 = 1
    while j2 <= n_pairs
        l  = lev_aniso[j2]; op = off5[j2]
        if l == 0
            @printf("%-8d  %-8s  %-8d  %-12.1f  %-12d  %d\n",
                    j2, "L0", op,   μ_coord[op],   n_max_coarse, bnds[op])
            @printf("%-8d  %-8s  %-8d  %-12.1f  %-12d  %d\n",
                    j2, "L0", op+1, μ_coord[op+1], n_max_coarse, bnds[op+1])
            j2 += 1
        elseif l == 1
            @printf("%-8d  %-8s  %-8d  %-12.1f  %-12d  %d\n",
                    j2, "L1", op, μ_coord[op], n_max_coarse, bnds[op])
            j2 += 1
        else
            @printf("%-8d  %-8s  %-8d  %-12.1f  %-12d  %d  ← L2: old cap BELOW μ!\n",
                    j2, "L2", op, μ_coord[op], n_max_coarse, bnds[op])
            j2 += 2
        end
    end
end

println("\nKey: L2 uniform cap ($(n_max_coarse)) < L2 mean (≈$(round(Int,μ_coord[end]))) → expansion was blocked upward")
println("     L0 uniform cap ($(n_max_coarse)) >> L0 mean (≈$(round(Int,μ_coord[off5[3]]))) → unnecessary slack removed by aniso")
