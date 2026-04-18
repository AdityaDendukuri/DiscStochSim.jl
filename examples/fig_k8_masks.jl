"""
K=8 runtime mask history via the deterministic RDME rate equation.

Evolves voxel means μ_k under:
  dμ_k/dt = f(μ_k) + d*(μ_{k-1} + μ_{k+1} - 2μ_k)
where f is the Schlögl rate equation and d = D/h².

From the mean-field trajectory we compute α_j and β_j at each step and show
which pairs the adaptive selector would mark fine.  This is a valid proxy for
the mask logic: the selector only sees voxel means, so the mean-field ODE
gives the same mask decisions as the full CME would on average.
"""

using DiscStochSim
using CairoMakie
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 9f0; const MED = 10f0; const LW = 1.5f0
const C_FINE   = RGBf(0.835, 0.369, 0.000)
const C_COARSE = RGBf(0.000, 0.447, 0.698)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.7f0)))

D     = 0.005
K     = 8
h     = 1.0 / K
d     = D / h^2          # diffusion rate per molecule

model = SchloglModel1D(D)
fp    = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
@printf("n_low=%d  n_uns=%d  n_high=%d  d=%.4f\n", n_low, n_uns, n_high, d)

# ── Deterministic ODE step (RK4) ─────────────────────────────────────────────
function rhs!(du, u)
    for k in 1:K
        rxn = schlogl_rate_eq(model, u[k])
        left  = k > 1 ? d * (u[k-1] - u[k]) : 0.0
        right = k < K ? d * (u[k+1] - u[k]) : 0.0
        du[k] = rxn + left + right
    end
end

function rk4_step!(u, dt)
    k1 = zeros(K); rhs!(k1, u)
    u2 = u .+ (dt/2).*k1
    k2 = zeros(K); rhs!(k2, u2)
    u3 = u .+ (dt/2).*k2
    k3 = zeros(K); rhs!(k3, u3)
    u4 = u .+ dt.*k3
    k4 = zeros(K); rhs!(k4, u4)
    u .+= (dt/6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
end

# ── Mask from voxel means ─────────────────────────────────────────────────────
function mask_from_means(μ; α_lo=0.10, α_hi=0.25)
    n_pairs = K ÷ 2
    α = zeros(n_pairs)
    β = zeros(n_pairs - 1)
    for j in 1:n_pairs
        tot = μ[2j-1] + μ[2j]
        α[j] = tot > 0.0 ? abs(μ[2j-1] - μ[2j]) / tot : 0.0
    end
    for j in 1:n_pairs-1
        tot = μ[2j] + μ[2j+1]
        β[j] = tot > 0.0 ? abs(μ[2j] - μ[2j+1]) / tot : 0.0
    end
    mask = trues(n_pairs)   # true = coarse
    for j in 1:n_pairs
        βl = j > 1       ? β[j-1] : 0.0
        βr = j < n_pairs ? β[j]   : 0.0
        α_eff = max(α[j], βl, βr)
        if   α_eff > α_hi; mask[j] = false
        end
    end
    mask, α, β
end

# ── Time integration ──────────────────────────────────────────────────────────
t_end   = 8.0
dt_ode  = 0.002      # small ODE step for accuracy
dt_rec  = 0.1        # record interval
n_rec   = round(Int, t_end / dt_rec)
steps_per_rec = round(Int, dt_rec / dt_ode)

n_pairs = K ÷ 2
ts      = Float64[]
μ_hist  = Matrix{Float64}(undef, 0, K)
masks   = Matrix{Bool}(undef, 0, n_pairs)
α_hist  = Matrix{Float64}(undef, 0, n_pairs)
β_hist  = Matrix{Float64}(undef, 0, n_pairs-1)

# IC: front between voxels 3 and 4 (inside pair 2)
u = Float64[n_low, n_low, n_low, n_high, n_high, n_high, n_high, n_high]

for rec in 0:n_rec
    m, α, β = mask_from_means(u)
    push!(ts, rec * dt_rec)
    μ_hist = vcat(μ_hist, u')
    masks  = vcat(masks,  m')
    α_hist = vcat(α_hist, α')
    β_hist = vcat(β_hist, β')
    for _ in 1:steps_per_rec
        rk4_step!(u, dt_ode)
    end
end

@printf("Final means: %s\n", string(round.(u; digits=1)))

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(580, 320))

# top panel: voxel means
ax1 = Axis(fig[1,1];
    ylabel = "mean count  μₖ",
    title  = "K=8 Schlögl — voxel means (deterministic)",
    xticklabelsvisible = false,
)
colors_vox = range(colorant"royalblue", colorant"firebrick"; length=K)
for k in 1:K
    lines!(ax1, ts, μ_hist[:,k]; color=colors_vox[k], linewidth=LW,
           label="v$k")
end
hlines!(ax1, [Float64(n_low), Float64(n_uns), Float64(n_high)];
        color=:black, linestyle=:dot, linewidth=0.6)
axislegend(ax1; position=:rc, framevisible=false, labelsize=SMALL-1,
           nbanks=2, tellwidth=false)

# bottom panel: mask heatmap (fine=1, coarse=0)
ax2 = Axis(fig[2,1];
    xlabel = "time",
    ylabel = "pair",
    title  = "Adaptive mask  (orange = fine, blue = coarse)",
    yticks = (1:n_pairs, ["pair $j" for j in 1:n_pairs]),
    limits = (nothing, (0.5, n_pairs+0.5)),
)
for j in 1:n_pairs
    for i in eachindex(ts)
        is_coarse = masks[i, j]
        scatter!(ax2, [ts[i]], [j];
                 color     = is_coarse ? C_COARSE : C_FINE,
                 marker    = is_coarse ? :rect    : :circle,
                 markersize = 6)
    end
end

rowsize!(fig.layout, 1, Relative(0.62))

CairoMakie.save(joinpath(OUTDIR, "fig_k8_masks.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_k8_masks.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_k8_masks.{pdf,png}")
