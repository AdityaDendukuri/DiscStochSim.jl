using DiscStochSim
using CairoMakie

# ── model ─────────────────────────────────────────────────────────────────────
K     = 4
D     = 1.0
model = SchloglModel1D(D)
fine_grid = VoxelGrid(K, 1.0, 0)
rates = Float64[]

n_low, n_unstable, n_high = schlogl_fixed_points(model)
println("Schlögl fixed points: low=$n_low  unstable=$n_unstable  high=$n_high")

u0 = CartesianIndex(n_unstable, 0, 0, 0)   # v1 at unstable point, rest empty

alg = RDMEMultigridFSP(
    dt                = 0.1,
    n_levels          = 1,
    n_max             = 250,
    save_every        = 5,       # snapshot every 0.5 time units
    max_states        = 500_000,
    expand_depth      = 2,       # fine expand! drives wavefront propagation
    prune_tol         = 1e-9,
    coarse_expand_depth = 2,
    τ_pre             = 0.5,
)

tspan = (0.0, 5.0)
println("Running V-cycle: IC=$u0, K=$K, D=$D ...")
sol = solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)
println("Done.  $(length(sol.t)) snapshots")
println("State space sizes: ", sol.state_space_sizes)
println()

mu = mean_voxel_counts(sol)
for (ti, t) in enumerate(sol.t)
    println("  t=$(round(t,digits=1))  |S|=$(sol.state_space_sizes[ti])  μ=$(round.(mu[ti,:],digits=1))")
end

# ── marginal helper ───────────────────────────────────────────────────────────
function voxel_marginals(sol, ti, n_max)
    states, probs = sol.snapshots[ti]
    P = zeros(Float64, K, n_max + 1)
    for (s, p) in zip(states, probs)
        t = Tuple(s)
        for k in 1:K
            n = t[k]; n <= n_max && (P[k, n+1] += p)
        end
    end
    P
end

# ── figure ────────────────────────────────────────────────────────────────────
# Pick 5 evenly spaced snapshots
tidx = round.(Int, range(1, length(sol.t), length=5))
tidx = unique(clamp.(tidx, 1, length(sol.t)))

n_fig_max = 200
P_all = [voxel_marginals(sol, i, n_fig_max) for i in tidx]
p_max = maximum(maximum(P) for P in P_all)
clims = (0.0, p_max ^ (1/3))

fig = Figure(size = (155 * length(tidx) + 90, 490), fontsize = 11)

Label(fig[0, 1:length(tidx)],
    "Schlögl RDME: K=$K, D=$D — IC: v1 at unstable point n=$n_unstable, v2–v4 empty  [color: P^{1/3}]",
    fontsize = 10, tellwidth = false)

hms = Any[]
for (col, ti) in enumerate(tidx)
    t_val = sol.t[ti]
    P     = P_all[col] .^ (1/3)

    ax = Axis(fig[1, col],
              title  = "t = $(round(t_val, digits=1))  (|S|=$(sol.state_space_sizes[ti]))",
              xlabel = "voxel k",
              ylabel = col == 1 ? "molecules n" : "",
              xticks = (1:K, ["v$k" for k in 1:K]),
              yticklabelsvisible = col == 1,
              ylabelvisible      = col == 1,
              titlesize = 11)

    hm = heatmap!(ax, 1:K, 0:n_fig_max, P,
                  colormap = :inferno, colorrange = clims,
                  lowclip = :black, interpolate = false)
    push!(hms, hm)

    hlines!(ax, [n_low, n_high]; color = (:white, 0.55), linestyle = :dash,  linewidth = 1.0)
    hlines!(ax, [n_unstable];   color = (:white, 0.30), linestyle = :dot,   linewidth = 0.8)
    for x in 1.5:1.0:K-0.5
        vlines!(ax, [x]; color = (:white, 0.15), linewidth = 0.5)
    end

    xlims!(ax, 0.5, K + 0.5); ylims!(ax, 0, n_fig_max)
end

cb_true = [0.0, 0.001, 0.01, 0.05, 0.15, 0.4, 0.8]
Colorbar(fig[1, length(tidx)+1], hms[end],
         label = "P(nₖ = n)",
         ticks = (cb_true .^ (1/3), string.(cb_true)),
         width = 14, ticksize = 4)

Label(fig[2, 1:length(tidx)],
    "n-level V-cycle (K=$K → K=$(K÷2)), dynamic-π prolongation.  " *
    "White dashed: stable states (n=$n_low, n=$n_high).  White dotted: unstable point (n=$n_unstable).  " *
    "v1 escapes the unstable point via reactions; diffusion carries probability to neighboring voxels.",
    fontsize = 8, tellwidth = false, color = :gray40)

rowsize!(fig.layout, 0, Auto(0.07))
rowsize!(fig.layout, 2, Auto(0.07))
colgap!(fig.layout, 4)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_schlogl_unstable_spread.pdf")
save(out, fig)
println("\nSaved: $out")
fig
