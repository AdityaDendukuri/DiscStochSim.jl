"""
Empirical timestep benchmark: K=4 direct vs K=2 coarse Schlögl FSP.

Builds frozen fine/coarse state spaces from a short warmup, then compares full-horizon
integration error and runtime over a dt sweep using the same Krylov dimension for both.
This isolates the practical timestep advantage of the less-stiff coarse generator.

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_stiffness.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Statistics
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0
const MED   = 9f0
const LW    = 1.4f0
const C_FINE   = RGBf(0.000, 0.447, 0.698)
const C_COARSE = RGBf(0.835, 0.369, 0.000)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

tv_distance(p, q) = 0.5 * sum(abs, p .- q)

function propagate_expv(A, p0, t_final, dt; m::Int)
    p = copy(p0)
    t = 0.0
    while t < t_final - 1e-12
        h = min(dt, t_final - t)
        p .= max.(expv(h, A, p; m=m), 0.0)
        s = sum(p)
        s > 0 && (p ./= s)
        t += h
    end
    p
end

function elapsed_propagate(A, p0, t_final, dt; m::Int)
    t0 = time()
    p = propagate_expv(A, p0, t_final, dt; m)
    time() - t0, p
end

# ── Model ─────────────────────────────────────────────────────────────────────
D = 0.005
K = 4
h = 1.0 / K
model = SchloglModel1D(D)
rates = Float64[]
fine_grid   = VoxelGrid(K, h, 0)
coarse_grid = coarsen(fine_grid)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_max = 180

u0_fine   = CartesianIndex(n_low, n_low, n_high, n_high)
u0_coarse = CartesianIndex(2 * n_low, 2 * n_high)

fine_sys   = build_schlogl_rdme_system(model, fine_grid)
coarse_sys = build_schlogl_coarse_system(model, coarse_grid, fine_grid)
bc_fine    = s -> rdme_bc(s, n_max)
bc_coarse  = s -> rdme_bc(s, 2 * n_max)

d_fine   = D / h^2
d_coarse = D / (2h)^2
@printf("d_fine=%.4f  d_coarse=%.4f  ratio=%.0f×\n", d_fine, d_coarse, d_fine / d_coarse)

# ── Build fixed state spaces by short warmup ─────────────────────────────────
println("\nBuilding fine state space (5 warmup steps at dt=0.2)..."); flush(stdout)
sp_fine = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp_fine, u0_fine, 1.0)
for step in 1:5
    expand!(sp_fine, fine_sys, bc_fine; depth=1)
    A_tmp, = build_generator(sp_fine, fine_sys, rates, 0.0)
    sp_fine.probs .= max.(expv(0.2, A_tmp, sp_fine.probs; m=60), 0.0)
    renormalize!(sp_fine)
    prune_threshold!(sp_fine, 1e-12)
    @printf("  fine step %d: |S|=%d\n", step, length(sp_fine)); flush(stdout)
end

println("\nBuilding coarse state space (5 warmup steps at dt=0.2)..."); flush(stdout)
sp_coarse = StateSpace{CartesianIndex{2}, Float64}()
add_state!(sp_coarse, u0_coarse, 1.0)
for step in 1:5
    expand!(sp_coarse, coarse_sys, bc_coarse; depth=1)
    A_tmp, = build_generator(sp_coarse, coarse_sys, rates, 0.0)
    sp_coarse.probs .= max.(expv(0.2, A_tmp, sp_coarse.probs; m=60), 0.0)
    renormalize!(sp_coarse)
    prune_threshold!(sp_coarse, 1e-12)
    @printf("  coarse step %d: |S|=%d\n", step, length(sp_coarse)); flush(stdout)
end

sz_fine = length(sp_fine)
sz_coarse = length(sp_coarse)
sz_ratio = round(Int, sz_fine / max(sz_coarse, 1))
@printf("\nFine   |S|=%d\n", sz_fine)
@printf("Coarse |S|=%d  (%d× fewer)\n", sz_coarse, sz_ratio)

# Build generators once on frozen state spaces
A_fine, = build_generator(sp_fine, fine_sys, rates, 0.0)
A_coarse, = build_generator(sp_coarse, coarse_sys, rates, 0.0)

p0_fine = zeros(Float64, length(sp_fine))
p0_coarse = zeros(Float64, length(sp_coarse))
p0_fine[sp_fine.index[u0_fine]] = 1.0
p0_coarse[sp_coarse.index[u0_coarse]] = 1.0

# ── Reference and dt sweep ────────────────────────────────────────────────────
t_final = 2.0
dt_ref  = 0.01
m_ref   = 100
m_test  = 30
tol_tv  = 1e-4
dt_vals = [0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.60]

println("\nComputing high-accuracy references..."); flush(stdout)
p_ref_fine   = propagate_expv(A_fine, p0_fine, t_final, dt_ref; m=m_ref)
p_ref_coarse = propagate_expv(A_coarse, p0_coarse, t_final, dt_ref; m=m_ref)

err_fine   = Float64[]
err_coarse = Float64[]
time_fine  = Float64[]
time_coarse = Float64[]

println("\ndt sweep: full solve to T=2.0 with fixed Krylov m=30"); flush(stdout)
@printf("%-6s  %-10s  %-10s  %-10s  %-10s\n",
        "dt", "TV fine", "TV coarse", "t fine (s)", "t coarse"); flush(stdout)
for dt in dt_vals
    tf, p_f = elapsed_propagate(A_fine, p0_fine, t_final, dt; m=m_test)
    tc, p_c = elapsed_propagate(A_coarse, p0_coarse, t_final, dt; m=m_test)

    ef = tv_distance(p_f, p_ref_fine)
    ec = tv_distance(p_c, p_ref_coarse)
    push!(err_fine, ef); push!(err_coarse, ec)
    push!(time_fine, tf); push!(time_coarse, tc)

    @printf("%-6.2f  %-10.3e  %-10.3e  %-10.4f  %-10.4f\n",
            dt, ef, ec, tf, tc); flush(stdout)
end

function max_accepted_dt(dt_vals, errs, tol)
    accepted = [dt for (dt, err) in zip(dt_vals, errs) if err <= tol]
    isempty(accepted) ? NaN : maximum(accepted)
end

dtacc_fine   = max_accepted_dt(dt_vals, err_fine, tol_tv)
dtacc_coarse = max_accepted_dt(dt_vals, err_coarse, tol_tv)
speedup_dt   = dtacc_coarse / dtacc_fine

@printf("\nTolerance TV ≤ %.1e\n", tol_tv)
@printf("  max accepted dt (fine)   = %.2f\n", dtacc_fine)
@printf("  max accepted dt (coarse) = %.2f\n", dtacc_coarse)
@printf("  empirical timestep gain  = %.1f×\n", speedup_dt)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(780, 280))

ax_a = Axis(fig[1, 1];
    xlabel="time step  Δt",
    ylabel="TV error at T=2",
    title="(a)  Accuracy vs Δt  [fixed Krylov m=30]",
    titlesize=MED, xscale=log10, yscale=log10)

lines!(ax_a, dt_vals, err_fine; color=C_FINE, linewidth=LW, label="K=4 direct")
scatter!(ax_a, dt_vals, err_fine; color=C_FINE, markersize=5)
lines!(ax_a, dt_vals, err_coarse; color=C_COARSE, linewidth=LW, linestyle=:dash, label="K=2 coarse")
scatter!(ax_a, dt_vals, err_coarse; color=C_COARSE, markersize=5)
hlines!(ax_a, [tol_tv]; color=(:gray, 0.6), linestyle=:dot, linewidth=0.8)
text!(ax_a, dt_vals[1] * 1.1, tol_tv * 1.2;
      text="TV tol = $(tol_tv)", fontsize=SMALL, color=:gray)
axislegend(ax_a; position=:lt, framevisible=false, labelsize=SMALL - 1)

ax_b = Axis(fig[1, 2];
    xlabel="time step  Δt",
    ylabel="wall time to T=2 (s)",
    title="(b)  Runtime vs Δt",
    titlesize=MED, xscale=log10, yscale=log10)

lines!(ax_b, dt_vals, time_fine; color=C_FINE, linewidth=LW, label="K=4 direct")
scatter!(ax_b, dt_vals, time_fine; color=C_FINE, markersize=5)
lines!(ax_b, dt_vals, time_coarse; color=C_COARSE, linewidth=LW, linestyle=:dash, label="K=2 coarse")
scatter!(ax_b, dt_vals, time_coarse; color=C_COARSE, markersize=5)

text!(ax_b, dt_vals[2], maximum(time_fine) * 0.35;
      text="|S| ratio: $(sz_ratio)×\naccepted Δt gain: $(round(speedup_dt, digits=1))×",
      fontsize=SMALL, color=:black)

colgap!(fig.layout, 12)
CairoMakie.save(joinpath(OUTDIR, "fig_stiffness.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_stiffness.png"), fig; px_per_unit=2)
println("\nSaved: paper/figures/fig_stiffness.{pdf,png}")
