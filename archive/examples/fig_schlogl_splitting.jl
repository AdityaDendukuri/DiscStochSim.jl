using DiscStochSim
using CairoMakie

# ── model & simulation ────────────────────────────────────────────────────────
K     = 4
D     = 5.0
model = SchloglModel1D(D)
fine_grid = VoxelGrid(K, 1.0, 0)
rates = Float64[]

n_low, n_unstable, n_high = schlogl_fixed_points(model)

u0    = CartesianIndex(ntuple(_ -> n_unstable, Val(K)))
tspan = (0.0, 12.0)

alg = RDMEMultigridFSP(
    dt                = 0.1,
    n_levels          = 1,
    n_max             = 300,
    save_every        = 10,
    max_states        = 500_000,
    expand_depth      = 2,
    prune_tol         = 1e-10,
    coarse_expand_depth = 2,
    τ_pre             = 0.5,
)

println("Running Schlögl splitting simulation...")
sol = solve_rdme_multigrid(model, fine_grid, u0, tspan, rates, alg)
println("Done. $(length(sol.t)) snapshots, final |S|=$(sol.state_space_sizes[end])")

# ── marginal distribution helper ──────────────────────────────────────────────
function voxel_marginals(sol, ti, n_max)
    states, probs = sol.snapshots[ti]
    P = zeros(Float64, K, n_max + 1)   # P[k, n+1] = P(voxel k has n molecules)
    for (s, p) in zip(states, probs)
        t = Tuple(s)
        for k in 1:K
            n = t[k]
            n <= n_max && (P[k, n+1] += p)
        end
    end
    P
end

# ── figure setup ──────────────────────────────────────────────────────────────
# Deduplicate times (final save can produce a duplicate last snapshot)
unique_times = unique(round.(sol.t, digits=6))
snap_indices = [findfirst(t -> round(t, digits=6) == ut, sol.t) for ut in unique_times]
snap_indices = snap_indices[round.(sol.t[snap_indices], digits=1) .∈ Ref([0.0, 3.0, 6.0, 9.0, 12.0])]
n_panels = length(snap_indices)
n_max    = 200

# cube-root transform: reveals low-probability tails while keeping structure
cbrtP(P) = P .^ (1/3)

n_max  = 175   # captures both stable states (31, 162) with headroom
P_all  = [cbrtP(voxel_marginals(sol, i, n_max)) for i in snap_indices]
p_max  = maximum(maximum(P) for P in P_all)
clims  = (0.0, p_max)

cmap = :inferno

fig = Figure(size = (155 * n_panels + 90, 500), fontsize = 11)

ax_label = Label(fig[0, 1:n_panels],
    "Schlögl RDME: K=$K voxels, D=$D — unstable IC (n=$n_unstable) splitting toward low (n=$n_low) and high (n=$n_high)  [color: P(nₖ=n)^{1/3}]",
    fontsize = 10, tellwidth = false)

axes = Axis[]
hms  = Any[]

for (col, si) in enumerate(snap_indices)
    t_val = sol.t[si]
    P     = P_all[col]

    ax = Axis(fig[1, col],
              title   = "t = $(round(t_val, digits=1))",
              xlabel  = "voxel k",
              ylabel  = col == 1 ? "molecules n" : "",
              xticks  = (1:K, ["v$k" for k in 1:K]),
              yticklabelsvisible = col == 1,
              ylabelvisible      = col == 1,
              titlesize = 12)
    push!(axes, ax)

    # heatmap: x = voxel (1..K), y = n (0..n_max), color = P(n_k = n)
    hm = heatmap!(ax, 1:K, 0:n_max, P,
                  lowclip = :black,
                  colormap = cmap, colorrange = clims,
                  interpolate = false)
    push!(hms, hm)

    # stable-state and unstable-point markers
    hlines!(ax, [n_low, n_high];   color = (:white, 0.55), linestyle = :dash,  linewidth = 1)
    hlines!(ax, [n_unstable];      color = (:white, 0.30), linestyle = :dot,   linewidth = 0.8)
    # pair-boundary separator (between v2 and v3)
    vlines!(ax, [2.5]; color = (:white, 0.25), linestyle = :solid, linewidth = 0.6)

    xlims!(ax, 0.5, K + 0.5)
    ylims!(ax, 0, n_max)
end

# Colorbar ticks labelled as true probability (cube-root scale inverted)
cb_true = [0.0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.6, 1.0]
Colorbar(fig[1, n_panels + 1], hms[end],
         label = "P(nₖ = n)",
         ticks = (cb_true .^ (1/3), string.(cb_true)),
         width = 14, ticksize = 4)

# annotation row
caption = "n-level V-cycle (K=$K → K=$(K÷2)), " *
          "dynamic-π prolongation, τ_pre=0.5.  " *
          "White dashed: stable states (n=$n_low, n=$n_high).  " *
          "White dotted: unstable point (n=$n_unstable)."
Label(fig[2, 1:n_panels], caption, fontsize = 8, tellwidth = false,
      color = :gray40)

rowsize!(fig.layout, 0, Auto(0.08))
rowsize!(fig.layout, 2, Auto(0.08))
colgap!(fig.layout, 4)

out = joinpath(@__DIR__, "..", "paper", "figures", "fig_schlogl_splitting.pdf")
save(out, fig)
println("Saved: $out")
fig
