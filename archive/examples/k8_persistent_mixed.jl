"""
Persistent mixed-state solver for the K=8 Schlögl RDME.

Story:
  (a) Single-voxel stationary CME P(n): establishes local bistability
  (b) Spatial front profile over time: pair 2 kept fine throughout (never coarsened)
  (c) P(n̄₂) vs P(n₃,n₄): what the coarse sum discards vs what fine resolution preserves

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/k8_persistent_mixed.jl
"""

using DiscStochSim
using CairoMakie
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0; const MED = 9f0; const LW = 1.4f0
const C_FINE   = RGBf(0.835, 0.369, 0.000)
const C_COARSE = RGBf(0.000, 0.447, 0.698)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

# ── Model parameters ──────────────────────────────────────────────────────────
D = 0.005; K = 8; h = 1.0/K; n_pairs = K ÷ 2
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
@printf("n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high); flush(stdout)

n_buf  = 30
n̄_low  = 2 * n_low
n̄_high = 2 * n_high

# ══════════════════════════════════════════════════════════════════════════════
# (a) Single-voxel stationary CME — analytical via detailed balance
# ══════════════════════════════════════════════════════════════════════════════
let c1=model.c1, c2=model.c2, c3=model.c3, c4=model.c4
    λ(n) = c3 + c1*n*(n-1)/2
    μ(n) = c4*n + c2*n*(n-1)*(n-2)/6
    n_max = n_high + 60
    log_P = zeros(n_max+1)
    for n in 0:n_max-1
        mn1 = μ(n+1)
        log_P[n+2] = log_P[n+1] + log(λ(n)) - (mn1 > 0 ? log(mn1) : -700.0)
    end
    log_P .-= maximum(log_P)
    P = exp.(log_P); P ./= sum(P)
    global ns_sv = collect(0:n_max)
    global ps_sv = P
end
@printf("P peaks: n_low=%.4f  n_high=%.4f\n", ps_sv[n_low+1], ps_sv[n_high+1]); flush(stdout)

# ══════════════════════════════════════════════════════════════════════════════
# Persistent mixed solver — pair 2 fine, pairs 1/3/4 coarse
# State tuple: (n̄₁, n₃, n₄, n̄₃, n̄₄)   [KM = 5]
# ══════════════════════════════════════════════════════════════════════════════
# Fixed levels: pair 2 = fine (interface pair), others level-1 coarse
levels = [1, 0, 1, 1]
off    = DiscStochSim._mixed_offsets(levels)
op2    = off[2]    # index of n₃ in tuple; n₄ is at op2+1
KM     = _mixed_dim(levels)   # = 5
@printf("\nlevels=%s  KM=%d  pair2 offset=%d\n", string(levels), KM, op2); flush(stdout)

# IC: pair 2 is the interface pair, 50/50 mixture of (n_low,n_high) and (n_high,n_low).
# Both states have the same coarse sum n̄₂ = n_low+n_high = 193, so the coarse
# marginal P(n̄₂) cannot distinguish them — pair 2 must stay fine to resolve this.
u0a = CartesianIndex(n̄_low, n_low,  n_high, n̄_high, n̄_high)  # n₃=n_low,  n₄=n_high
u0b = CartesianIndex(n̄_low, n_high, n_low,  n̄_high, n̄_high)  # n₃=n_high, n₄=n_low
@printf("u0a = %s  (50%%)\n", string(u0a)); flush(stdout)
@printf("u0b = %s  (50%%)\n", string(u0b)); flush(stdout)

sp = StateSpace{CartesianIndex{KM}, Float64}()
add_state!(sp, u0a, 0.5)
add_state!(sp, u0b, 0.5)

mxsys = build_schlogl_mixed_system(model, g, levels)

# BC: coarse dimensions ≤ 2*(n_high+n_buf), fine dimensions ≤ n_high+n_buf
n_max_coarse = 2*(n_high + n_buf)
n_max_fine   = n_high + n_buf
function bc_mix(s)
    t = Tuple(s)
    t[1] ≤ n_max_coarse &&       # n̄₁
    t[op2]   ≤ n_max_fine   &&   # n₃
    t[op2+1] ≤ n_max_fine   &&   # n₄
    t[4]  ≤ n_max_coarse &&      # n̄₃
    t[5]  ≤ n_max_coarse         # n̄₄
end
expand!(sp, mxsys, bc_mix; depth=1)
@printf("Initial |S_mix|=%d\n\n", length(sp.states)); flush(stdout)

# ── helper: voxel means and stds ──────────────────────────────────────────────
function voxel_stats(sp, levels, off)
    μ  = zeros(K)
    μ2 = zeros(K)
    for (i,s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for j in 1:n_pairs
            o = off[j]
            if levels[j] > 0   # coarse (level 1 or 2): pair total at t[o]
                h2 = p*t[o]/2
                μ[2j-1] += h2;             μ[2j] += h2
                μ2[2j-1] += p*(t[o]/2)^2; μ2[2j] += p*(t[o]/2)^2
            else               # fine: individual voxel coords t[o], t[o+1]
                μ[2j-1] += p*t[o];    μ[2j] += p*t[o+1]
                μ2[2j-1] += p*t[o]^2; μ2[2j] += p*t[o+1]^2
            end
        end
    end
    μ, sqrt.(max.(μ2 .- μ.^2, 0.0))
end

# ── time loop ─────────────────────────────────────────────────────────────────
dt = 0.1; n_steps = 60     # t_final = 6
snap_steps = [0, 5, 15, 30, 60]   # t = 0, 0.5, 1.5, 3, 6

snaps = Dict{Int, Tuple{Vector{Float64},Vector{Float64}}}()
snaps[0] = voxel_stats(sp, levels, off)

ts       = Float64[]
S_sizes  = Int[]
levelss  = [copy(levels)]   # track levels at each step for diagnostic

cache = MixedSolverCache(); sc = Ref(0)
prev_wc = Ref(-1.0)

@printf("%-6s  %-8s  %-8s  %-20s\n", "step", "|S_mix|", "Σp", "levels"); flush(stdout)
for step in 1:n_steps
    global sp, levels, off, op2
    sp, levels = two_level_step_mixed(sp, levels, model, g, dt;
                                      cache               = cache,
                                      step_count          = sc,
                                      mask_check_interval = 100_000,
                                      expand_mixed        = true,
                                      mixed_expand_depth  = 1,
                                      mixed_n_max         = n_max_coarse,
                                      prune_tol           = 1e-5,
                                      reexpand_depth      = 0,
                                      krylov_m            = 40)
    off = DiscStochSim._mixed_offsets(levels)
    op2 = off[2]
    sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
    push!(ts, step*dt); push!(S_sizes, length(sp.states)); push!(levelss, copy(levels))
    step in snap_steps && (snaps[step] = voxel_stats(sp, levels, off))
    step % 10 == 0 && @printf("%-6d  %-8d  %-8.5f  %s\n",
        step, S_sizes[end], sum(sp.probs), string(levels)); flush(stdout)
    # Save joint marginal at step 5 (t=0.5) — bimodal structure clearest here
    if step == 5
        global mg_n3n4_early = Dict{Tuple{Int,Int}, Float64}()
        for (i,s) in enumerate(sp.states)
            tt = Tuple(s); p = sp.probs[i]
            key = (tt[op2], tt[op2+1])
            mg_n3n4_early[key] = get(mg_n3n4_early, key, 0.0) + p
        end
        @printf("  → saved P(n₃,n₄) at t=0.5: %d unique states\n", length(mg_n3n4_early))
        flush(stdout)
    end
    # Save intermediate joint marginal at step 30 (t=3) for FSP-vs-SSA comparison
    if step == 30
        global mg_n3n4_t3 = Dict{Tuple{Int,Int}, Float64}()
        for (i,s) in enumerate(sp.states)
            tt = Tuple(s); p = sp.probs[i]
            key = (tt[op2], tt[op2+1])
            mg_n3n4_t3[key] = get(mg_n3n4_t3, key, 0.0) + p
        end
        @printf("  → saved P(n₃,n₄) at t=3: %d unique states\n", length(mg_n3n4_t3))
        flush(stdout)
    end
end
@printf("\nFinal |S_mix|=%d\n\n", length(sp.states)); flush(stdout)

# ── P(n₃,n₄) marginal at t_final ─────────────────────────────────────────────
mg_n3n4 = Dict{Tuple{Int,Int}, Float64}()
for (i,s) in enumerate(sp.states)
    t = Tuple(s); p = sp.probs[i]
    key = (t[op2], t[op2+1])
    mg_n3n4[key] = get(mg_n3n4, key, 0.0) + p
end

# Use early-time marginal (t=0.5) for panels (c)/(d): bimodal structure is clear
mg_n3n4_plot = mg_n3n4_early
t_plot = 0.5

# Coarse-equivalent marginal P(n̄₂) — what coarse would see
mg_nbar2 = Dict{Int, Float64}()
for ((n3,n4), p) in mg_n3n4_plot
    nb = n3 + n4
    mg_nbar2[nb] = get(mg_nbar2, nb, 0.0) + p
end
@printf("P(n₃,n₄) at t=%.1f: %d unique states\n", t_plot, length(mg_n3n4_plot)); flush(stdout)

# ══════════════════════════════════════════════════════════════════════════════
# Figure
# ══════════════════════════════════════════════════════════════════════════════
fig = Figure(size=(780, 460))

# ── (a) Single-voxel P(n) ─────────────────────────────────────────────────────
ax_a = Axis(fig[1:2, 1];
    xlabel="n  (molecules / voxel)", ylabel="P_stat(n)",
    title="(a)  Single-voxel stationary CME", titlesize=MED)
barplot!(ax_a, ns_sv, ps_sv; color=(C_COARSE,0.7), strokewidth=0, gap=0, width=1)
vlines!(ax_a, [Float64(n_low), Float64(n_high)]; color=:black, linestyle=:dash, linewidth=0.8)
vlines!(ax_a, [Float64(n_uns)];                  color=:red,   linestyle=:dot,  linewidth=0.8)
text!(ax_a, Float64(n_low)+1,  maximum(ps_sv)*0.9; text="n_low",  fontsize=SMALL-1, align=(:left,:top))
text!(ax_a, Float64(n_high)-1, maximum(ps_sv)*0.9; text="n_high", fontsize=SMALL-1, align=(:right,:top))
text!(ax_a, Float64(n_uns)+1,  maximum(ps_sv)*0.5; text="n_uns",  fontsize=SMALL-1, color=:red, align=(:left,:bottom))

# ── (b) Spatial front profile ─────────────────────────────────────────────────
ax_b = Axis(fig[1, 2:4];
    xlabel="voxel k", ylabel="E[nₖ]  (molecules)",
    title="(b)  K=8 persistent mixed solver  —  pair 2 (v₃,v₄) always fine",
    titlesize=MED,
    xticks=(1:K, ["v$k" for k in 1:K]))

xs = Float64.(1:K)
snap_colors = resample_cmap(:viridis, length(snap_steps))
for (ci, ss) in enumerate(snap_steps)
    haskey(snaps, ss) || continue
    μs, σs = snaps[ss]
    t_val  = ss * dt
    clr    = snap_colors[ci]
    errorbars!(ax_b, xs, μs, σs; color=(clr,0.45), linewidth=0.9, whiskerwidth=4)
    scatterlines!(ax_b, xs, μs; color=clr, linewidth=LW, markersize=7,
                  label="t=$(round(t_val, digits=1))")
end
# Shade pair backgrounds
for j in 1:n_pairs
    c   = levels[j] > 0 ? C_COARSE : C_FINE
    lbl = levels[j] > 0 ? "C$j" : "F$j"
    poly!(ax_b, Rect(2j-1.5, -8.0, 2.0, Float64(n_high)+25); color=(c,0.08), strokewidth=0)
    text!(ax_b, Float64(2j-0.5), Float64(n_high)+12;
          text=lbl, fontsize=SMALL-2, color=c, align=(:center,:bottom), font=:bold)
end
hlines!(ax_b, [Float64(n_low), Float64(n_high)]; color=:black, linestyle=:dash, linewidth=0.7)
text!(ax_b, K+0.1, Float64(n_low);  text=" n_low",  fontsize=SMALL-1, align=(:left,:center))
text!(ax_b, K+0.1, Float64(n_high); text=" n_high", fontsize=SMALL-1, align=(:left,:center))
axislegend(ax_b; position=:rc, framevisible=false, labelsize=SMALL-1, nbanks=length(snap_steps))
xlims!(ax_b, 0.3, K+1.5); ylims!(ax_b, -5, n_high+30)

# ── (c) Coarse P(n̄₂) ─────────────────────────────────────────────────────────
ax_c = Axis(fig[2, 2];
    xlabel="n̄₂ = n₃+n₄", ylabel="probability",
    title="(c)  Coarse: P(n̄₂)  at t=$(t_plot)", titlesize=MED)
nb_ks = sort(collect(keys(mg_nbar2)))
nb_vs = [mg_nbar2[k] for k in nb_ks]
barplot!(ax_c, Float64.(nb_ks), nb_vs; color=(C_COARSE,0.75), strokewidth=0, gap=0, width=1)
vlines!(ax_c, [Float64(n̄_low), Float64(n̄_high)]; color=:black, linestyle=:dash, linewidth=0.8)
text!(ax_c, Float64(n̄_low)+2,  maximum(nb_vs)*0.85; text="2n_low",  fontsize=SMALL-2, align=(:left,:center))
text!(ax_c, Float64(n̄_high)-2, maximum(nb_vs)*0.85; text="2n_high", fontsize=SMALL-2, align=(:right,:center))

# ── (d) Fine P(n₃,n₄) ────────────────────────────────────────────────────────
ax_d = Axis(fig[2, 3:4];
    xlabel="n₃  (molecules in voxel 3)", ylabel="n₄  (molecules in voxel 4)",
    title="(d)  Fine: P(n₃,n₄)  at t=$(t_plot)  —  two fine modes, same coarse total",
    titlesize=MED)

bin_w  = 4
n_max_hm = n_high + n_buf
n_bins = div(n_max_hm, bin_w) + 1
hm_data = zeros(n_bins, n_bins)
for ((n3,n4), p) in mg_n3n4_plot
    i3 = clamp(div(n3, bin_w) + 1, 1, n_bins)
    i4 = clamp(div(n4, bin_w) + 1, 1, n_bins)
    hm_data[i3, i4] += p
end
bin_ctrs = Float64.(0:bin_w:bin_w*(n_bins-1))
heatmap!(ax_d, bin_ctrs, bin_ctrs, hm_data;
         colormap=:hot, colorrange=(0, max(maximum(hm_data)*0.95, 1e-10)))

# Crosshairs at fixed points
vlines!(ax_d, [Float64(n_low), Float64(n_high)]; color=(:white,0.4), linestyle=:dash, linewidth=0.9)
hlines!(ax_d, [Float64(n_low), Float64(n_high)]; color=(:white,0.4), linestyle=:dash, linewidth=0.9)

# Anti-diagonal: restriction maps both modes to the same coarse sum n₃+n₄ = n_low+n_high
x_diag = Float64.([0, n_low+n_high])
lines!(ax_d, x_diag, reverse(x_diag); color=RGBf(0.9,0.2,0.2), linestyle=:dot, linewidth=2.0)
# Label the diagonal near the centre
mid = (n_low + n_high) / 2
text!(ax_d, mid - 12, mid + 8; text="n₃+n₄=$(n_low+n_high)", fontsize=SMALL-1,
      color=RGBf(1.0,0.5,0.5), align=(:right,:bottom), rotation=Float32(-atan(1)))

# Label only the two occupied modes
text!(ax_d, Float64(n_low)+3,  Float64(n_high)+3;  text="(nₗ,nₕ)", fontsize=SMALL-1,
      color=:white, align=(:left,:bottom))
text!(ax_d, Float64(n_high)+3, Float64(n_low)+3;   text="(nₕ,nₗ)", fontsize=SMALL-1,
      color=:white, align=(:left,:bottom))

rowsize!(fig.layout, 1, Relative(0.50))
rowsize!(fig.layout, 2, Relative(0.40))
colsize!(fig.layout, 1, Relative(0.24))

CairoMakie.save(joinpath(OUTDIR, "fig_k8_persistent.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_k8_persistent.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_k8_persistent.{pdf,png}")

# ── Serialize FSP data for FSP-vs-SSA comparison ─────────────────────────────
using Serialization
fsp_data = Dict(
    "mg_n3n4"    => mg_n3n4,
    "mg_n3n4_t3" => mg_n3n4_t3,   # richer intermediate distribution at t=3
    "mg_nbar2"   => mg_nbar2,
    "snaps"    => snaps,
    "ts"       => ts,
    "S_sizes"  => S_sizes,
    "n_low"    => n_low,
    "n_uns"    => n_uns,
    "n_high"   => n_high,
    "model_c1" => model.c1,
    "model_c2" => model.c2,
    "model_c3" => model.c3,
    "model_c4" => model.c4,
    "D"        => D,
    "K"        => K,
    "h"        => h,
    "t_final"  => n_steps * dt,
)
serialize(joinpath(OUTDIR, "fsp_k8_data.jls"), fsp_data)
println("Serialized FSP data → paper/figures/fsp_k8_data.jls")
