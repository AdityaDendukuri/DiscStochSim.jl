"""
Nonlinear closure validation for Schlögl pairs.

Runs a K=4 fine FSP with a front IC and compares exact conditional factorial moments
inside two pairs:
  - pair 1 = interface pair (n₁, n₂)
  - pair 2 = homogeneous bulk pair (n₃, n₄)

For each pair total N, compares exact conditional moments against the Binomial
closure used in the coarse propensities:
  E[n(n-1) | N]        vs  N(N-1)/4
  E[n(n-1)(n-2) | N]   vs  N(N-1)(N-2)/8

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_closure.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0
const MED   = 9f0
const LW    = 1.2f0
const C_IFACE = RGBf(0.800, 0.180, 0.160)
const C_BULK  = RGBf(0.000, 0.447, 0.698)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

fac2(n) = n * (n - 1)
fac3(n) = n * (n - 1) * (n - 2)

function pair_stats(sp, idx1, idx2; pmin=5e-5)
    # Aggregate pair marginal
    mg = Dict{Tuple{Int, Int}, Float64}()
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        key = (t[idx1], t[idx2])
        mg[key] = get(mg, key, 0.0) + sp.probs[i]
    end

    P_total = Dict{Int, Float64}()
    M2 = Dict{Int, Float64}()
    M3 = Dict{Int, Float64}()
    for ((n1, n2), p) in mg
        N = n1 + n2
        P_total[N] = get(P_total, N, 0.0) + p
        M2[N] = get(M2, N, 0.0) + p * fac2(n1)
        M3[N] = get(M3, N, 0.0) + p * fac3(n1)
    end
    for N in keys(P_total)
        pN = P_total[N]
        if pN > 0
            M2[N] /= pN
            M3[N] /= pN
        end
    end

    N_vals = sort([N for (N, pN) in P_total if pN > pmin])
    probs = [P_total[N] for N in N_vals]
    m2_exact = [M2[N] for N in N_vals]
    m3_exact = [M3[N] for N in N_vals]
    m2_bin = [N * (N - 1) / 4 for N in N_vals]
    m3_bin = [N * (N - 1) * (N - 2) / 8 for N in N_vals]
    err2 = abs.(m2_exact .- m2_bin)
    err3 = abs.(m3_exact .- m3_bin)
    rel2 = err2 ./ max.(m2_bin, 1.0)
    rel3 = err3 ./ max.(m3_bin, 1.0)

    (; N_vals, probs, m2_exact, m3_exact, m2_bin, m3_bin, err2, err3, rel2, rel3, P_total)
end

# ── Model ─────────────────────────────────────────────────────────────────────
D = 0.005
K = 4
h = 1.0 / K
model = SchloglModel1D(D)
rates = Float64[]
fine_grid = VoxelGrid(K, h, 0)
fine_sys = build_schlogl_rdme_system(model, fine_grid)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_max = 200
bc = s -> rdme_bc(s, n_max)
@printf("n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)

# Front inside pair 1; pair 2 stays homogeneous high bulk
u0 = CartesianIndex(n_low, n_high, n_high, n_high)
T_END = 1.5
dt = 0.15
n_steps = round(Int, T_END / dt)

sp = StateSpace{CartesianIndex{4}, Float64}()
add_state!(sp, u0, 1.0)

println("Running K=4 fine FSP ($(n_steps) steps)..."); flush(stdout)
for step in 1:n_steps
    expand!(sp, fine_sys, bc; depth=1)
    A, = build_generator(sp, fine_sys, rates, 0.0)
    sp.probs .= max.(expv(dt, A, sp.probs; m=50), 0.0)
    renormalize!(sp)
    prune_threshold!(sp, 1e-12)
    step % 3 == 0 && @printf("  step %d/%d  |S|=%d\n", step, n_steps, length(sp))
    flush(stdout)
end
@printf("Final |S|=%d\n\n", length(sp))

iface = pair_stats(sp, 1, 2)
bulk  = pair_stats(sp, 3, 4)

println("Key closure errors:")
for (label, stats, N_key) in [("interface", iface, n_low + n_high),
                              ("bulk",      bulk,  2 * n_high)]
    if N_key in stats.N_vals
        idx = findfirst(==(N_key), stats.N_vals)
        @printf("  %-10s N=%d: rel2=%.3f  rel3=%.3f  P(N)=%.4f\n",
                label, N_key, stats.rel2[idx], stats.rel3[idx], stats.probs[idx])
    end
end

function marker_sizes(probs)
    pmax = maximum(probs)
    [4.0f0 + 9.0f0 * sqrt(p / pmax) for p in probs]
end

ms_iface = marker_sizes(iface.probs)
ms_bulk  = marker_sizes(bulk.probs)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(780, 280))

ax_a = Axis(fig[1, 1];
    xlabel="pair total  N",
    ylabel="relative error",
    title="(a)  Second factorial-moment closure error",
    titlesize=MED, yscale=log10)

scatter!(ax_a, Float64.(bulk.N_vals), bulk.rel2;
         color=C_BULK, markersize=ms_bulk, label="bulk pair (n₃,n₄)")
scatter!(ax_a, Float64.(iface.N_vals), iface.rel2;
         color=C_IFACE, markersize=ms_iface, label="interface pair (n₁,n₂)")

for (N_key, label) in [(2 * n_low, "2n_lo"), (n_low + n_high, "n_lo+n_hi"), (2 * n_high, "2n_hi")]
    vlines!(ax_a, [Float64(N_key)]; color=(:gray, 0.4), linewidth=0.8, linestyle=:dot)
    text!(ax_a, Float64(N_key) + 1, maximum(vcat(iface.rel2, bulk.rel2)) * 0.08;
          text=label, fontsize=SMALL - 1, color=:gray)
end
axislegend(ax_a; position=:rt, framevisible=false, labelsize=SMALL - 1)

ax_b = Axis(fig[1, 2];
    xlabel="pair total  N",
    ylabel="relative error",
    title="(b)  Third factorial-moment closure error",
    titlesize=MED, yscale=log10)

scatter!(ax_b, Float64.(bulk.N_vals), bulk.rel3;
         color=C_BULK, markersize=ms_bulk)
scatter!(ax_b, Float64.(iface.N_vals), iface.rel3;
         color=C_IFACE, markersize=ms_iface)

for N_key in [2 * n_low, n_low + n_high, 2 * n_high]
    vlines!(ax_b, [Float64(N_key)]; color=(:gray, 0.4), linewidth=0.8, linestyle=:dot)
end

text!(ax_b, Float64(2 * n_high) * 0.80, maximum(vcat(iface.rel3, bulk.rel3)) * 0.12;
      text="bulk pair:\nsmall error", fontsize=SMALL, color=C_BULK)
text!(ax_b, Float64(n_low + n_high) * 0.72, maximum(vcat(iface.rel3, bulk.rel3)) * 0.65;
      text="interface pair:\nlarge error", fontsize=SMALL, color=C_IFACE)

colgap!(fig.layout, 12)
CairoMakie.save(joinpath(OUTDIR, "fig_closure.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_closure.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_closure.{pdf,png}")
