"""
Speedup figure for talk: K=4 Schlögl direct FSP vs K=2 coarse FSP,
plus K=8 state-space probe showing intractability.

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_speedup.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0; const MED = 9f0; const LW = 1.4f0
const C_DIRECT = RGBf(0.000, 0.447, 0.698)
const C_COARSE = RGBf(0.835, 0.369, 0.000)
const C_K8     = RGBf(0.8,   0.2,   0.2  )

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

# ── Model ─────────────────────────────────────────────────────────────────────
D = 0.005; K = 4; h = 1.0/K
model = SchloglModel1D(D)
rates = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_max = 180
@printf("n_low=%d  n_uns=%d  n_high=%d  d=%.4f\n", n_low, n_uns, n_high, D/h^2)

u0_fine   = CartesianIndex(n_low, n_low, n_high, n_high)
u0_coarse = CartesianIndex(2*n_low, 2*n_high)

T_END = 2.0; dt = 0.2; n_steps = round(Int, T_END / dt)

fine_sys   = build_schlogl_rdme_system(model, fine_grid)
coarse_sys = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
bc_fine    = s -> rdme_bc(s, n_max)
bc_coarse  = s -> rdme_bc(s, 2*n_max)

# ── K=4 direct adaptive FSP ───────────────────────────────────────────────────
function run_direct(u0, sys, bc, n_steps, dt, T_END)
    sp = StateSpace{CartesianIndex{4}, Float64}()
    add_state!(sp, u0, 1.0)
    ts = Float64[]; szs = Int[]; wall = 0.0
    t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt, T_END - t_cur)
        expand!(sp, sys, bc; depth=1)
        A, = build_generator(sp, sys, rates, t_cur)
        t0 = time()
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=50), 0.0)
        wall += time() - t0
        t_cur += dt_step
        renormalize!(sp); prune_threshold!(sp, 1e-12)
        push!(ts, t_cur); push!(szs, length(sp))
        step % 2 == 0 && (@printf("  step %d  |S|=%d\n", step, length(sp)); flush(stdout))
    end
    ts, szs, wall
end

# ── K=2 coarse adaptive FSP ───────────────────────────────────────────────────
function run_coarse(u0, sys, bc, n_steps, dt, T_END)
    sp = StateSpace{CartesianIndex{2}, Float64}()
    add_state!(sp, u0, 1.0)
    ts = Float64[]; szs = Int[]; wall = 0.0
    t_cur = 0.0
    for step in 1:n_steps
        dt_step = min(dt, T_END - t_cur)
        expand!(sp, sys, bc; depth=1)
        A, = build_generator(sp, sys, rates, t_cur)
        t0 = time()
        sp.probs .= max.(expv(dt_step, A, sp.probs; m=30), 0.0)
        wall += time() - t0
        t_cur += dt_step
        renormalize!(sp); prune_threshold!(sp, 1e-12)
        push!(ts, t_cur); push!(szs, length(sp))
        step % 2 == 0 && (@printf("  step %d  |S|=%d\n", step, length(sp)); flush(stdout))
    end
    ts, szs, wall
end

# ── K=8 direct: probe how fast state space grows ─────────────────────────────
function probe_k8(u0, sys, bc; max_states=500_000)
    sp = StateSpace{CartesianIndex{8}, Float64}()
    add_state!(sp, u0, 1.0)
    expand!(sp, sys, bc; depth=1)
    sz_per_depth = Int[length(sp)]
    t0 = time()
    for depth in 2:6
        expand!(sp, sys, bc; depth=1)
        push!(sz_per_depth, length(sp))
        @printf("  K=8 depth %d: |S|=%d  (%.1fs)\n", depth, length(sp), time()-t0)
        flush(stdout)
        length(sp) > max_states && break
    end
    sz_per_depth
end

# ── Run ───────────────────────────────────────────────────────────────────────
println("\nK=4 direct adaptive FSP..."); flush(stdout)
ts_d, szs_d, wall_d = run_direct(u0_fine, fine_sys, bc_fine, n_steps, dt, T_END)

println("\nK=2 coarse adaptive FSP..."); flush(stdout)
ts_c, szs_c, wall_c = run_coarse(u0_coarse, coarse_sys, bc_coarse, n_steps, dt, T_END)

println("\nK=8 direct probe..."); flush(stdout)
K8 = 8; h8 = 1.0/K8
fine_grid8 = VoxelGrid(K8, h8, 0)
fine_sys8  = build_schlogl_rdme_system(model, fine_grid8)
bc8        = s -> rdme_bc(s, n_max)
u0_8       = CartesianIndex(n_low, n_low, n_low, n_low, n_high, n_high, n_high, n_high)
szs_k8     = probe_k8(u0_8, fine_sys8, bc8; max_states=300_000)

sz_d = szs_d[end]; sz_c = szs_c[end]
println()
println("="^60)
@printf("  K=4 direct:  |S|=%6d  expv=%.2fs\n", sz_d, wall_d)
@printf("  K=2 coarse:  |S|=%6d  expv=%.2fs  (%d× fewer states)\n",
        sz_c, wall_c, round(Int, sz_d / max(sz_c, 1)))
@printf("  K=8 probe:   |S|=%6d after limited expansion\n", szs_k8[end])
println("="^60)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(780, 280))

# Panel (a): |S| vs time
ax_a = Axis(fig[1,1];
    xlabel="time", ylabel="|S|  (state space size)",
    title="(a)  State space growth  [K=4 Schlögl, front IC]",
    titlesize=MED)

lines!(ax_a, ts_d, Float64.(szs_d);
       color=C_DIRECT, linewidth=LW,
       label="K=4 direct FSP  (ground truth)")
lines!(ax_a, ts_c, Float64.(szs_c);
       color=C_COARSE, linewidth=LW, linestyle=:dash,
       label="K=2 coarse FSP  (=K=8 V-cycle coarse level)")
axislegend(ax_a; position=:lt, framevisible=false, labelsize=SMALL-1)

reduction = round(Int, sz_d / max(sz_c, 1))
text!(ax_a, T_END * 0.5, Float64(sz_d) * 0.4;
      text="$(reduction)× fewer states at t=$(T_END)",
      fontsize=SMALL, color=C_COARSE)

# Panel (b): bar chart — final |S| for K=4 direct and K=2 coarse;
# K=8 bar uses the depth-6 probe (states from initial expansion only, before solving)
k8_bar = szs_k8[end]   # 258k — just from expansion, not yet solved

ax_b = Axis(fig[1,2];
    ylabel="|S|  (log₁₀ scale)",
    title="(b)  State space at t=$(T_END)  [K=4] / depth-6 expansion [K=8]",
    titlesize=MED, yscale=log10,
    xticks=([1, 2, 3], ["K=4\ndirect", "K=2\ncoarse\n(K=8 level)", "K=8\ndirect\n(probe)"]))

barplot!(ax_b, [1, 2, 3], Float64.([sz_d, sz_c, k8_bar]);
         color=[C_DIRECT, C_COARSE, C_K8], gap=0.3)

text!(ax_b, 3.0, Float64(k8_bar) * 2.5;
      text="expansion\nonly\n(no solve)",
      fontsize=SMALL-1, align=(:center,:bottom), color=C_K8)

colgap!(fig.layout, 12)
CairoMakie.save(joinpath(OUTDIR, "fig_speedup.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_speedup.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_speedup.{pdf,png}")
