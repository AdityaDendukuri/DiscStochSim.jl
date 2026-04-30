"""
Figures for the adaptive coarsening section.

Fig A (adaptive_criterion):  Two-panel.
  Left:  α_intra[j] and β_right[j] for scenarios A, B, C on K=8 grid.
         Shows which pairs are forced fine by each criterion.
  Right: Resulting KM and conceptual mixed-grid layout for scenarios A and B.

Fig B (adaptive_accuracy):   K=4 front accuracy comparison.
  Direct FSP (reference) vs. adaptive V-cycle vs. non-adaptive V-cycle.
  Voxel means and TV(direct FSP) over time on the Schlögl front IC.

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_adaptive.jl
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using Printf

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)
save_fig(name, fig) = begin
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=2)
    println("  Saved: paper/figures/$name.{pdf,png}")
end

# ── theme ──────────────────────────────────────────────────────────────────────
const SMALL = 9f0; const MED = 10f0; const LW = 1.5f0
const C_FINE   = RGBf(0.835, 0.369, 0.000)   # red-orange: fine pair
const C_COARSE = RGBf(0.000, 0.447, 0.698)   # blue: coarse pair
const C_ALPHA  = RGBf(0.000, 0.620, 0.451)   # green: α_intra bars
const C_BETA   = RGBf(0.800, 0.475, 0.655)   # purple: β_right bars
const C_THRESH = RGBf(0.4,   0.4,   0.4)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.7f0)))

# ── model ──────────────────────────────────────────────────────────────────────
D   = 0.005
K8  = 8;  h8 = 1.0/K8;  g8 = VoxelGrid(K8, h8, 0)
model = SchloglModel1D(D)
fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_max8 = n_high + 20
r8  = build_schlogl_rdme_system(model, g8)
rm8 = RDMEModel1D(model.D, 0.0, 0.0)

@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)

function diagnostics(u0; depth=1)
    sp = StateSpace{CartesianIndex{K8}, Float64}()
    add_state!(sp, u0, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max8, Tuple(s))
    expand!(sp, r8, bc; depth)
    α, φ = pair_admissibility(sp, rm8, g8)
    μ = zeros(K8)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s); p = sp.probs[i]
        for k in 1:K8; μ[k] += p * t[k]; end
    end
    n_pairs = K8 ÷ 2
    β = [abs(μ[2j]-μ[2j+1])/max(μ[2j]+μ[2j+1],1e-14) for j in 1:n_pairs-1]
    mask = select_coarsening_mask(sp, rm8, g8)
    α, β, mask
end

α_A, β_A, mask_A = diagnostics(
    CartesianIndex(n_low, n_low, n_low, n_high, n_high, n_high, n_high, n_high))
α_B, β_B, mask_B = diagnostics(
    CartesianIndex(n_low, n_low, n_low, n_low, n_high, n_high, n_high, n_high))
α_C, β_C, mask_C = diagnostics(
    CartesianIndex(n_low, n_low, n_low, n_low, n_low, n_low, n_low, n_low))

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE A: adaptive criterion diagnostic bars
# ══════════════════════════════════════════════════════════════════════════════
println("\n── Fig A: adaptive criterion ──")

let
    n_pairs = K8 ÷ 2
    pair_labels = ["P$j" for j in 1:n_pairs]
    # β padded to length n_pairs (no right boundary for last pair)
    pad(β) = [β; 0.0]

    fig = Figure(size=(580, 420))

    # helper: draw one scenario panel
    function panel!(pos, label, α, β, mask; show_ylabel=false)
        ax = Axis(fig[pos...];
            title = label,
            xticks = (1:n_pairs, pair_labels),
            yticks = 0.0:0.2:1.0,
            ylabel = show_ylabel ? "asymmetry  α, β" : "",
            ylabelsize = MED,
            limits = (0.5, n_pairs+0.5, -0.02, 1.0),
        )

        # threshold lines
        hlines!(ax, 0.10; color = C_THRESH, linestyle = :dot, linewidth = 0.9)
        hlines!(ax, 0.25; color = C_THRESH, linestyle = :dash, linewidth = 0.9)
        text!(ax, n_pairs + 0.45, 0.11; text = "αlo", fontsize = SMALL-1,
              align = (:right, :bottom), color = C_THRESH)
        text!(ax, n_pairs + 0.45, 0.26; text = "αhi", fontsize = SMALL-1,
              align = (:right, :bottom), color = C_THRESH)

        # bars: α_intra and β_right side by side
        w = 0.35
        for j in 1:n_pairs
            # α_intra bar
            barplot!(ax, [j - w/2], [α[j]];
                     width = w, color = mask[j] ? C_COARSE : C_FINE,
                     strokewidth = 0.5, strokecolor = :black)
            # β_right bar (only for j < n_pairs)
            if j < n_pairs
                barplot!(ax, [j + w/2], [β[j]];
                         width = w, color = C_BETA,
                         strokewidth = 0.5, strokecolor = :black)
            end
        end

        # mask annotation below each pair
        for j in 1:n_pairs
            clr = mask[j] ? C_COARSE : C_FINE
            lbl = mask[j] ? "C" : "F"
            text!(ax, Float64(j), -0.015; text = lbl, fontsize = SMALL,
                  align = (:center, :top), color = clr, font = :bold)
        end
        ax
    end

    panel!([1,1], "Scenario A — front inside pair 2",  α_A, β_A, mask_A; show_ylabel=true)
    panel!([1,2], "Scenario B — front on pair 2|3 boundary", α_B, β_B, mask_B)
    panel!([1,3], "Scenario C — homogeneous (all-low)", α_C, β_C, mask_C)

    # legend
    elem_ai  = PolyElement(color=C_COARSE, strokecolor=:black, strokewidth=0.5)
    elem_ai2 = PolyElement(color=C_FINE,   strokecolor=:black, strokewidth=0.5)
    elem_b   = PolyElement(color=C_BETA,   strokecolor=:black, strokewidth=0.5)
    Legend(fig[2,:],
           [elem_ai, elem_ai2, elem_b],
           ["α_intra (coarsened pair)", "α_intra (fine pair)", "β_right (inter-pair)"];
           orientation=:horizontal, framevisible=false, labelsize=SMALL,
           tellwidth=false)

    save_fig("fig_adaptive_criterion", fig)
end

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE B: K=4 accuracy comparison  (direct FSP vs adaptive vs non-adaptive)
# ══════════════════════════════════════════════════════════════════════════════
println("\n── Fig B: K=4 accuracy comparison ──")

K4 = 4; h4 = 1.0/K4; g4 = VoxelGrid(K4, h4, 0)
coarse_grid4 = coarsen(g4)
# Use a small cap — the goal is a method demo, not precise fixed-point numerics.
# Synthetic counts lo4=3, hi4=15 within [0, n_max4] give clear asymmetry.
n_max4 = 20
lo4, hi4 = 3, 15
r4 = build_schlogl_rdme_system(model, g4)

# IC: front within pair 1  →  (lo4, hi4, hi4, hi4)
u0_4 = CartesianIndex(lo4, hi4, hi4, hi4)

t_plot = range(0.0, 2.0, length=41)
dt     = t_plot[2] - t_plot[1]
n_steps = length(t_plot) - 1

function voxel_means(sp)
    mu = zeros(K4)
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        for k in 1:K4; mu[k] += sp.probs[i] * t[k]; end
    end
    mu
end

function tv_distance(sp_a, sp_b)
    overlap = 0.0
    for (i, s) in enumerate(sp_a.states)
        j = get(sp_b.index, s, 0)
        if j > 0; overlap += min(sp_a.probs[i], sp_b.probs[j]); end
    end
    1.0 - overlap
end

# ── Direct FSP (reference) ────────────────────────────────────────────────────
println("  Running direct K=4 FSP...")
let
    sp = StateSpace{CartesianIndex{K4}, Float64}()
    add_state!(sp, u0_4, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max4, Tuple(s))
    expand!(sp, r4, bc; depth=2)
    global sp_ref_all = [sp]
    global mu_ref = zeros(n_steps+1, K4)
    mu_ref[1,:] = voxel_means(sp)
    t_cur = 0.0
    for step in 1:n_steps
        dt_step = dt
        expand!(sp, r4, bc; depth=1)
        A, = build_generator(sp, r4, Float64[], t_cur)
        sp.probs .= expv(dt_step, A, sp.probs; m=40)
        sp.probs  = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
        t_cur += dt_step
        mu_ref[step+1,:] = voxel_means(sp)
        push!(sp_ref_all, deepcopy(sp))
    end
end
println("  Direct FSP done.")

# dynamic-π table (K=2 at h=1/K4)
pair_grid4 = VoxelGrid(2, h4, 0)
pair_sys4  = build_schlogl_rdme_system(model, pair_grid4)
sp_lo4 = let sp = StateSpace{CartesianIndex{2},Float64}()
    add_state!(sp, CartesianIndex(n_low, n_low), 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max4, Tuple(s))
    expand!(sp, pair_sys4, bc; depth=3)
    sp
end
A_lo4, = build_generator(sp_lo4, pair_sys4, Float64[], 0.0)
pi_table4 = compute_dynamic_pi(sp_lo4, A_lo4; n_max=n_max4)

# ── Non-adaptive V-cycle (Mult. prolong.) ─────────────────────────────────────
println("  Running non-adaptive V-cycle...")
let
    sp = StateSpace{CartesianIndex{K4}, Float64}()
    add_state!(sp, u0_4, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max4, Tuple(s))
    expand!(sp, r4, bc; depth=1)
    global mu_vc = zeros(n_steps+1, K4)
    global tv_vc = zeros(n_steps+1)
    mu_vc[1,:] = voxel_means(sp)
    tv_vc[1]   = tv_distance(sp, sp_ref_all[1])
    t_cur = 0.0
    for step in 1:n_steps
        sp = two_level_vcycle_schlogl_injection(sp, model, g4, coarse_grid4,
                                                pi_table4, Float64[], t_cur, dt;
                                                τ_pre=0.0, krylov_m=30,
                                                use_dynamic_pi=true,
                                                weight_tol=1e-8,
                                                coarse_n_max=2*n_max4,
                                                expand_coarse=true,
                                                coarse_expand_depth=1)
        sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
        t_cur += dt
        mu_vc[step+1,:] = voxel_means(sp)
        tv_vc[step+1]   = tv_distance(sp, sp_ref_all[step+1])
    end
end
println("  Non-adaptive done.")

# ── Adaptive V-cycle ──────────────────────────────────────────────────────────
println("  Running adaptive V-cycle...")
let
    sp = StateSpace{CartesianIndex{K4}, Float64}()
    add_state!(sp, u0_4, 1.0)
    bc = s -> all(c -> 0 ≤ c ≤ n_max4, Tuple(s))
    expand!(sp, r4, bc; depth=1)
    global mu_ad = zeros(n_steps+1, K4)
    global tv_ad = zeros(n_steps+1)
    global mask_ad = Vector{BitVector}(undef, n_steps+1)
    mu_ad[1,:]  = voxel_means(sp)
    tv_ad[1]    = tv_distance(sp, sp_ref_all[1])
    mask_ad[1]  = BitVector([true, true])
    prev_mask   = nothing
    t_cur       = 0.0
    for step in 1:n_steps
        sp, prev_mask = two_level_vcycle_adaptive(sp, model, g4, Float64[], t_cur, dt;
                                                   prev_mask,
                                                   krylov_m=30,
                                                   expand_mixed=true,
                                                   mixed_expand_depth=1,
                                                   mixed_n_max=2*n_max4,
                                                   binom_tol=1e-4)
        sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
        t_cur += dt
        mu_ad[step+1,:]  = voxel_means(sp)
        tv_ad[step+1]    = tv_distance(sp, sp_ref_all[step+1])
        mask_ad[step+1]  = copy(prev_mask)
    end
end
println("  Adaptive done.")

# ── Plot: panels (a)+(b) — voxel means and TV ─────────────────────────────────
let
    ts = collect(t_plot)
    fig = Figure(size=(580, 270))

    ax1 = Axis(fig[1,1];
        xlabel = "time", ylabel = "mean molecule count",
        title  = "(a)  Voxel means  (K=4, front IC)",
    )
    colors_vox = [RGBf(0.1,0.1,0.9), RGBf(0.1,0.7,0.3),
                  RGBf(0.9,0.5,0.1), RGBf(0.7,0.1,0.1)]
    for k in 1:K4
        lines!(ax1, ts, mu_ref[:,k]; color=colors_vox[k], linewidth=LW+0.4)
        lines!(ax1, ts, mu_vc[:,k];  color=colors_vox[k], linewidth=LW, linestyle=:dash)
        lines!(ax1, ts, mu_ad[:,k];  color=colors_vox[k], linewidth=LW, linestyle=:dot)
    end
    elem_r = LineElement(color=:black, linewidth=LW+0.4)
    elem_v = LineElement(color=:black, linewidth=LW, linestyle=:dash)
    elem_a = LineElement(color=:black, linewidth=LW, linestyle=:dot)
    Legend(fig[1,1], [elem_r, elem_v, elem_a],
           ["Direct FSP", "Non-adaptive V-cycle", "Adaptive V-cycle"];
           position=:rt, framevisible=false, labelsize=SMALL, tellwidth=false)

    ax2 = Axis(fig[1,2];
        xlabel = "time", ylabel = "TV(direct FSP)",
        title  = "(b)  Accuracy",
        yscale = log10,
    )
    lines!(ax2, ts[2:end], tv_vc[2:end]; color=C_COARSE, linewidth=LW,
           linestyle=:dash, label="Non-adaptive")
    lines!(ax2, ts[2:end], tv_ad[2:end]; color=C_FINE, linewidth=LW,
           linestyle=:dot, label="Adaptive")
    axislegend(ax2; position=:rb, framevisible=false, labelsize=SMALL)

    save_fig("fig_adaptive_accuracy", fig)
end

# ── Plot: panel (c) — mask evolution ──────────────────────────────────────────
let
    ts = collect(t_plot)
    fig = Figure(size=(580, 180))

    ax3 = Axis(fig[1,1];
        xlabel = "time",
        ylabel = "pair status",
        title  = "(c)  Adaptive mask evolution",
        yticks = ([0, 1], ["fine", "coarse"]),
        limits = (nothing, (-0.15, 1.15)),
    )
    for j in 1:2
        ys = [Float64(mask_ad[i][j]) for i in eachindex(mask_ad)]
        scatter!(ax3, ts, ys; markersize = 6,
                 marker = j == 1 ? :circle : :rect,
                 color  = j == 1 ? C_FINE : C_COARSE,
                 label  = "pair $j")
    end
    axislegend(ax3; position=:rt, framevisible=false, labelsize=SMALL)

    save_fig("fig_adaptive_mask", fig)
end

println("\nDone.")
