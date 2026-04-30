"""
Validate level-2 promotion: compare K=12 mixed solver runs at t=0.5.
  Run A: fixed levels=[1,1,0,1,1,1]  (level-1 coarse everywhere)
  Run B: adaptive, starts at [1,1,0,1,1,1], pairs 5-6 promoted to level-2

Metric: TV(P_A(n₅,n₆), P_B(n₅,n₆)) at t=0.5

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/compare_k12_levels.jl
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

# ── Shared IC ─────────────────────────────────────────────────────────────────
function make_ic(levels)
    off = DiscStochSim._mixed_offsets(levels)
    op  = off[3]   # fine pair 3 offset
    KM  = _mixed_dim(levels)
    u0a = CartesianIndex(ntuple(i -> i==1 ? n̄_low  :
                                     i==2 ? n̄_low  :
                                     i==op ? n_low  :
                                     i==op+1 ? n_high :
                                     i <= KM ? n̄_high : 0, KM)...)
    u0b = CartesianIndex(ntuple(i -> i==1 ? n̄_low  :
                                     i==2 ? n̄_low  :
                                     i==op ? n_high :
                                     i==op+1 ? n_low  :
                                     i <= KM ? n̄_high : 0, KM)...)
    sp = StateSpace{CartesianIndex{KM}, Float64}()
    add_state!(sp, u0a, 0.5); add_state!(sp, u0b, 0.5)
    sp, off, op
end

# ── Run solver for n_steps, return P(n₅,n₆) at step t05 ─────────────────────
function run_mixed(levels_init; mask_check_interval, windowed, n_steps=5, t05=5)
    sp, off, op3 = make_ic(levels_init)
    levels = copy(levels_init)
    mxsys  = build_schlogl_mixed_system(model, g, levels)

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

    mg_fine  = Dict{Tuple{Int,Int}, Float64}()   # P(n₅,n₆)
    mg_coarse = Dict{Int, Float64}()              # P(n̄₃) = P(n₅+n₆)
    for step in 1:n_steps
        sp, levels = two_level_step_mixed(sp, levels, model, g, 0.1;
                         cache               = cache,
                         step_count          = sc,
                         mask_check_interval = mask_check_interval,
                         windowed            = windowed,
                         halo                = 0,
                         expand_mixed        = true,
                         mixed_expand_depth  = 1,
                         mixed_n_max         = n_max_coarse,
                         prune_tol           = 0.0,   # no pruning: preserve symmetry
                         reexpand_depth      = 0,
                         krylov_m            = 40)
        sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
        if step == t05
            off2 = DiscStochSim._mixed_offsets(levels)
            op   = off2[3]
            for (i,s) in enumerate(sp.states)
                tt = Tuple(s); p = sp.probs[i]
                mg_fine[(tt[op], tt[op+1])] = get(mg_fine, (tt[op], tt[op+1]), 0.0) + p
                nb = tt[op] + tt[op+1]
                mg_coarse[nb] = get(mg_coarse, nb, 0.0) + p
            end
        end
    end
    mg_fine, mg_coarse, levels, length(sp.states)
end

# ── Run A: fixed level-1 ─────────────────────────────────────────────────────
@printf("Run A: fixed levels=[1,1,0,1,1,1]...\n"); flush(stdout)
@time mgf_A, mgc_A, lev_A, sA = run_mixed([1,1,0,1,1,1]; mask_check_interval=100_000, windowed=false)
@printf("  final levels=%s  |S|=%d  |P(n₅,n₆)|=%d\n\n", string(lev_A), sA, length(mgf_A)); flush(stdout)

# ── Run B: adaptive, level-2 ─────────────────────────────────────────────────
@printf("Run B: adaptive (windowed, halo=0)...\n"); flush(stdout)
@time mgf_B, mgc_B, lev_B, sB = run_mixed([1,1,0,1,1,1]; mask_check_interval=5, windowed=true)
@printf("  final levels=%s  |S|=%d  |P(n₅,n₆)|=%d\n\n", string(lev_B), sB, length(mgf_B)); flush(stdout)

# ── TV on fine P(n₅,n₆) ──────────────────────────────────────────────────────
fine_keys = union(keys(mgf_A), keys(mgf_B))
tv_fine = 0.5 * sum(abs(get(mgf_A, k, 0.0) - get(mgf_B, k, 0.0)) for k in fine_keys)
@printf("TV(P_A(n₅,n₆), P_B(n₅,n₆)) = %.6f\n", tv_fine)

# ── TV on coarse P(n̄₃) — symmetric, unaffected by mode-orientation ───────────
coarse_keys = union(keys(mgc_A), keys(mgc_B))
tv_coarse = 0.5 * sum(abs(get(mgc_A, k, 0.0) - get(mgc_B, k, 0.0)) for k in coarse_keys)
@printf("TV(P_A(n̄₃),  P_B(n̄₃))     = %.6f\n", tv_coarse)
@printf("  (TV coarse expected ≪ 0.01 if level-2 propensities are correct)\n\n")

# ── Symmetry check: P(a,b) vs P(b,a) for each run ────────────────────────────
function symm_error(mg)
    0.5 * sum(abs(get(mg,(a,b),0.0) - get(mg,(b,a),0.0)) for (a,b) in keys(mg))
end
@printf("Symmetry error A: %.6f  (0 = perfectly symmetric)\n", symm_error(mgf_A))
@printf("Symmetry error B: %.6f\n\n", symm_error(mgf_B))

# ── Top modes ────────────────────────────────────────────────────────────────
for (lbl, mg) in [("A", mgf_A), ("B", mgf_B)]
    top = sort(collect(mg), by=x->-x[2])[1:min(5,end)]
    @printf("Top modes %s:", lbl)
    for ((a,b),p) in top; @printf("  (%d,%d)=%.4f", a, b, p); end
    @printf("\n")
end
