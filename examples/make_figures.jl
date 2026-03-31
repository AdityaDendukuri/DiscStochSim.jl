"""
Publication figures for the RDME multigrid paper.

Generates CairoMakie figures for Sections 5-6 of rdme_multigrid.tex:
  Fig 1 (schlogl_bistability):   Schlögl rate equation and potential well
  Fig 2 (prolongation_failure):  p(n1|nc) — Binomial vs dynamic-π vs truth
  Fig 3 (heatmap_comparison):    K=2 fine distribution heatmaps (3-panel)
  Fig 4 (msa_tv_decay):          TV time series for K=2 MSA methods
  Fig 5 (vcycle_means):          Voxel means over time, K=4 V-cycle benchmark
  Fig 6 (vcycle_statespace):     Fine state-space size over time

All figures are saved to paper/figures/ as PDF + PNG (300 dpi).

Run (from repo root):
  JULIA_DEPOT_PATH=/tmp/jl_depot:~/.julia \\
    julia --project=paper/figures_env examples/make_figures.jl

  Set SKIP_VCYCLE=1 to skip the expensive K=4 V-cycle section (~60s).

Note: requires the separate figures_env (paper/figures_env/) which uses
CairoMakie 0.15 (compatible with Julia ≥1.12).  On first run, the
JULIA_DEPOT_PATH trick redirects precompile caches to /tmp to avoid
permission issues with root-owned cache dirs on some macOS installs.
"""

using DiscStochSim
using CairoMakie
using ExponentialUtilities
using LinearAlgebra
using SparseArrays
using Printf

# ── Output directory ──────────────────────────────────────────────────────────
const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)
save_fig(name, fig) = begin
    CairoMakie.save(joinpath(OUTDIR, name * ".pdf"), fig)
    CairoMakie.save(joinpath(OUTDIR, name * ".png"), fig; px_per_unit=2)
    println("  Saved: paper/figures/$name.{pdf,png}")
end

# ── CairoMakie theme: LaTeX-like, paper-quality ───────────────────────────────
const SMALL = 9f0       # tick/caption font size
const MED   = 10f0      # axis label font size
const LARGE = 11f0      # title font size
const LW    = 1.6f0     # default line width
const MS    = 6f0       # marker size

# Color palette (colorblind-friendly: Wong 2011)
const C_BLUE   = RGBf(0.000, 0.447, 0.698)   # direct FSP / ground truth
const C_ORANGE = RGBf(0.902, 0.624, 0.000)   # Binomial
const C_GREEN  = RGBf(0.000, 0.620, 0.451)   # dynamic-π (coarse)
const C_RED    = RGBf(0.835, 0.369, 0.000)   # V-cycle Binom
const C_PURPLE = RGBf(0.800, 0.475, 0.655)   # V-cycle Dyn-π
const C_GRAY   = RGBf(0.4,   0.4,   0.4)     # reference / annotation

paper_theme = Theme(
    fontsize      = MED,
    Axis          = (
        spinewidth    = 0.8,
        xgridvisible  = false,
        ygridvisible  = false,
        xticklabelsize = SMALL,
        yticklabelsize = SMALL,
        xlabelsize     = MED,
        ylabelsize     = MED,
        titlesize      = MED,
    ),
    Lines = (linewidth = LW,),
    Scatter = (markersize = MS,),
    Legend = (
        labelsize  = SMALL,
        framewidth = 0.6,
        padding    = (4f0, 4f0, 4f0, 4f0),
    ),
)
set_theme!(paper_theme)

# ── Model parameters ─────────────────────────────────────────────────────────
const D     = 0.005
const K_MSA = 2
const h_MSA = 1.0 / K_MSA
const d_MSA = D / h_MSA^2

model = SchloglModel1D(D)
rates = Float64[]

println("="^65)
println("Schlögl RDME Figure Generation")
@printf("Model: c1=%.4f c2=%.2e c3=%.1f c4=%.1f D=%.4f\n",
        model.c1, model.c2, model.c3, model.c4, model.D)
println("="^65)

fp = schlogl_fixed_points(model; n_max=250)
n_low, n_uns, n_high = fp
@printf("Fixed points: n_low=%d  n_uns=%d  n_high=%d\n", n_low, n_uns, n_high)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Schlögl rate equation and potential
# ══════════════════════════════════════════════════════════════════════════════

println("\n── Fig 1: Schlögl bistability ──")

let
    ns = 0:0.5:220
    f_vals = [schlogl_rate_eq(model, n) for n in ns]

    # Potential  U(n) = -∫₀ⁿ f(m)dm  (numerical integration)
    Δn = step(ns)
    U_vals = cumsum(-f_vals) .* Δn
    U_vals .-= U_vals[argmin(U_vals)]     # shift so global minimum = 0

    fig = Figure(size = (420, 310))

    # ─── panel (a): rate equation ───────────────────────────────────────────
    ax1 = Axis(fig[1,1];
        xlabel = "molecule count  n",
        ylabel = "f(n)  (net rate)",
        title  = "(a) Schlögl rate equation",
        yticks = WilkinsonTicks(5),
    )
    hlines!(ax1, 0; color = C_GRAY, linestyle = :dash, linewidth = 0.8)
    lines!(ax1, collect(ns), f_vals; color = C_BLUE, linewidth = LW + 0.4)

    # mark fixed points
    scatter!(ax1, [Float64(n_low), Float64(n_high)], [0.0, 0.0];
             color = C_BLUE, marker = :circle, markersize = MS + 2,
             label = "stable")
    scatter!(ax1, [Float64(n_uns)], [0.0];
             color = :white, marker = :circle, markersize = MS + 2,
             strokewidth = 1.4, strokecolor = C_BLUE)
    scatter!(ax1, [Float64(n_uns)], [0.0];
             color = :transparent, marker = :circle, markersize = MS + 2,
             label = "unstable")

    # annotations
    text!(ax1, Float64(n_low), 8.0; text = "n_low=$(n_low)", fontsize = SMALL,
          align = (:center, :bottom))
    text!(ax1, Float64(n_high), 8.0; text = "n_high=$(n_high)", fontsize = SMALL,
          align = (:center, :bottom))
    text!(ax1, Float64(n_uns), -12.0; text = "n_uns=$(n_uns)", fontsize = SMALL,
          align = (:center, :top))

    axislegend(ax1; position = :rb, framevisible = false, labelsize = SMALL)

    # ─── panel (b): potential well ──────────────────────────────────────────
    ax2 = Axis(fig[1,2];
        xlabel = "molecule count  n",
        ylabel = "U(n)  (–∫f dn, arb.)",
        title  = "(b) Effective potential",
    )
    lines!(ax2, collect(ns), U_vals; color = C_BLUE, linewidth = LW + 0.4)

    # shade stable wells
    n_arr = collect(ns)
    idx_lo_l = findfirst(n_arr .>= 0);     idx_lo_r = findfirst(n_arr .>= n_uns)
    idx_hi_l = idx_lo_r;                    idx_hi_r = length(n_arr)
    band!(ax2, n_arr[idx_lo_l:idx_lo_r], zeros(idx_lo_r - idx_lo_l + 1),
          U_vals[idx_lo_l:idx_lo_r]; color = (C_BLUE, 0.15))
    band!(ax2, n_arr[idx_hi_l:idx_hi_r], zeros(idx_hi_r - idx_hi_l + 1),
          U_vals[idx_hi_l:idx_hi_r]; color = (C_BLUE, 0.15))

    scatter!(ax2, [Float64(n_low), Float64(n_high)],
             [U_vals[argmin(abs.(n_arr .- n_low))],
              U_vals[argmin(abs.(n_arr .- n_high))]];
             color = C_BLUE, markersize = MS)
    scatter!(ax2, [Float64(n_uns)],
             [U_vals[argmin(abs.(n_arr .- n_uns))]];
             color = :white, markersize = MS, strokewidth = 1.4,
             strokecolor = C_BLUE)

    save_fig("fig1_schlogl_bistability", fig)
end

# ══════════════════════════════════════════════════════════════════════════════
# K=2 COMPUTATIONS (needed for Figs 2, 3, 4)
# ══════════════════════════════════════════════════════════════════════════════

println("\n── Running K=2 direct FSP ──")

const n_max = 200
fine_grid_2   = VoxelGrid(K_MSA, h_MSA, 0)
coarse_grid_2 = coarsen(fine_grid_2)

# Direct FSP state space
sp_fine2 = StateSpace{CartesianIndex{2}, Float64}()
for n1 in 0:n_max, n2 in 0:n_max
    add_state!(sp_fine2, CartesianIndex(n1, n2), 0.0)
end
u0_2 = CartesianIndex(n_low, n_high)
sp_fine2.probs[sp_fine2.index[u0_2]] = 1.0

fine_sys2   = build_schlogl_rdme_system(model, fine_grid_2)
A_fine2, _  = build_generator(sp_fine2, fine_sys2, rates, 0.0)
@printf("  Generator: |S|=%d  nnz=%d\n", length(sp_fine2), nnz(A_fine2))

t_solve_msa = [0.5, 1.0, 2.0, 5.0]
probs_fine2 = Dict{Float64, Vector{Float64}}()
t0 = time()
let p = copy(sp_fine2.probs), t_prev = 0.0
    for t in t_solve_msa
        p = max.(expv(t - t_prev, A_fine2, p; m = 50), 0.0)
        p ./= sum(p)
        probs_fine2[t] = copy(p)
        t_prev = t
    end
end
@printf("  Solve to t=%.0f: %.2fs\n", maximum(t_solve_msa), time() - t0)

# dynamic-π table
println("\n── Computing dynamic-π table ──")
t0 = time()
pi_table = compute_dynamic_pi(sp_fine2, A_fine2; n_max = n_max)
@printf("  Stationary solve: %.2fs\n", time() - t0)

nc_front = n_low + n_high
@printf("  Dynamic-π at nc=%d: peak n1=%d (p=%.4f)  Binom peak n1=%d\n",
        nc_front, argmax(pi_table[nc_front+1]) - 1,
        maximum(pi_table[nc_front+1]), nc_front ÷ 2)

# ── Coarse MSA runs ───────────────────────────────────────────────────────────

println("\n── Running K=2 coarse MSA (Binomial + dynamic-π) ──")

const dt_c_msa  = 0.05
const n_steps_msa = round(Int, maximum(t_solve_msa) / dt_c_msa)
coarse_bc_msa = s -> rdme_bc(s, 2 * n_max)
nc0_msa = n_low + n_high

# Helper: run one coarse MSA method
function run_coarse_msa(coarse_sys; label="")
    sp = StateSpace{CartesianIndex{1}, Float64}()
    add_state!(sp, CartesianIndex(nc0_msa), 1.0)
    expand!(sp, coarse_sys, coarse_bc_msa; depth = 1)
    snaps = Dict{Float64, Tuple{Vector{CartesianIndex{1}}, Vector{Float64}}}()
    t0_ = time()
    let t_cur = 0.0
        for step in 1:n_steps_msa
            dt_step = min(dt_c_msa, maximum(t_solve_msa) - t_cur)
            expand!(sp, coarse_sys, coarse_bc_msa; depth = 1)
            A, _ = build_generator(sp, coarse_sys, rates, t_cur)
            sp.probs .= expv(dt_step, A, sp.probs; m = 30)
            t_cur += dt_step
            renormalize!(sp); prune_threshold!(sp, 1e-12)
            for t in t_solve_msa
                if abs(t_cur - t) < dt_c_msa / 2 && !haskey(snaps, t)
                    snaps[t] = (copy(sp.states), copy(sp.probs))
                end
            end
        end
    end
    @printf("  %s: %.2fs\n", label, time() - t0_)
    snaps
end

coarse_sys_binom = build_schlogl_coarse_system(model, coarse_grid_2, fine_grid_2)
coarse_sys_dyn   = build_schlogl_coarse_system_dynamic(model, coarse_grid_2, fine_grid_2, pi_table)

snaps_binom = run_coarse_msa(coarse_sys_binom; label = "Coarse Binom")
snaps_dyn   = run_coarse_msa(coarse_sys_dyn;   label = "Coarse Dyn-π")

# TV vs direct FSP helper
function tv_vs_fine2(sp_approx)
    tv = 0.0
    for (i, s) in enumerate(sp_fine2.states)
        idx = get(sp_approx.index, s, 0)
        p_approx = idx != 0 ? sp_approx.probs[idx] : 0.0
        tv += abs(sp_fine2.probs[i] - p_approx)
    end
    for (i, s) in enumerate(sp_approx.states)
        get(sp_fine2.index, s, 0) != 0 && continue
        tv += sp_approx.probs[i]
    end
    tv / 2
end

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Conditional distribution p(n1 | nc = n_low + n_high)
# ══════════════════════════════════════════════════════════════════════════════

println("\n── Fig 2: Prolongation failure at the front ──")

let
    t_diag = 0.5
    p_true = probs_fine2[t_diag]
    nc = nc_front

    # True conditional p(n1 | nc) from direct FSP
    cond_true = zeros(nc + 1)
    for (i, s) in enumerate(sp_fine2.states)
        t2 = Tuple(s)
        if t2[1] + t2[2] == nc
            cond_true[t2[1] + 1] += p_true[i]
        end
    end
    mass_nc = sum(cond_true)
    if mass_nc > 1e-10
        cond_true ./= mass_nc
    end

    # Binomial(nc, 1/2)
    function binom_dist(nc)
        logph = -nc * log(2.0)
        b = zeros(nc + 1)
        for n1 in 0:nc
            lbc = sum(log(i) for i in (nc - n1 + 1):nc; init = 0.0) -
                  sum(log(i) for i in 1:n1; init = 0.0)
            b[n1 + 1] = exp(logph + lbc)
        end
        b ./ sum(b)
    end
    cond_binom = binom_dist(nc)

    # Dynamic-π
    cond_dyn = pi_table[nc + 1]

    n1_vals = 0:nc

    fig = Figure(size = (400, 270))
    ax = Axis(fig[1, 1];
        xlabel = "n₁   (left voxel count, n₂ = $(nc) − n₁)",
        ylabel = "p(n₁ | nc = $(nc))",
        title  = "Conditional distribution at the front  (nc = n_low + n_high = $(nc))",
        xticks = [0, n_low, nc ÷ 2, n_high, nc],
        xticklabelsvisible = true,
    )

    # Shade: Binom concentrates in the wrong place
    vspan!(ax, nc ÷ 2 - 8, nc ÷ 2 + 8; color = (C_ORANGE, 0.12))
    vspan!(ax, n_low - 8, n_low + 8; color = (C_GREEN, 0.12))

    lines!(ax, collect(n1_vals), cond_binom;
           color = C_ORANGE, linewidth = LW, label = "Binomial($(nc), ½)")
    lines!(ax, collect(n1_vals), cond_dyn;
           color = C_GREEN, linewidth = LW, linestyle = :dash,
           label = "Dynamic-π (stationary)")
    lines!(ax, collect(n1_vals), cond_true;
           color = C_BLUE, linewidth = LW, linestyle = :dashdot,
           label = "True (direct FSP, t=$(t_diag))")

    # Annotate the Binom mode
    binom_mode = nc ÷ 2
    binom_peak = cond_binom[binom_mode + 1]
    text!(ax, Float64(binom_mode), binom_peak + 0.003;
          text = "n₁ = $(binom_mode)\n(nc/2, wrong)",
          align = (:center, :bottom), fontsize = SMALL, color = C_ORANGE)

    # Annotate n_low (true initial position)
    text!(ax, Float64(n_low) - 3, maximum(cond_true) + 0.003;
          text = "n_low = $(n_low)",
          align = (:right, :bottom), fontsize = SMALL, color = C_BLUE)

    xlims!(ax, -5, nc + 5)

    axislegend(ax; position = :rt, framevisible = true)

    save_fig("fig2_prolongation_failure", fig)
end

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: K=2 fine distribution heatmaps (restrict–prolong accuracy)
# ══════════════════════════════════════════════════════════════════════════════

println("\n── Fig 3: Restrict-prolong heatmaps ──")

let
    t_hm = 0.5
    p_ref = probs_fine2[t_hm]

    # Restrict to coarse
    coarse_dict = Dict{Int, Float64}()
    for (i, s) in enumerate(sp_fine2.states)
        nc_s = Tuple(s)[1] + Tuple(s)[2]
        coarse_dict[nc_s] = get(coarse_dict, nc_s, 0.0) + p_ref[i]
    end
    sp_rp_c = StateSpace{CartesianIndex{1}, Float64}()
    for (nc_s, p) in coarse_dict
        add_state!(sp_rp_c, CartesianIndex(nc_s), p)
    end
    renormalize!(sp_rp_c)

    # Prolong
    sp_binom_rp = prolong(sp_rp_c, Val(2); binom_tol = 1e-6);  renormalize!(sp_binom_rp)
    sp_dyn_rp   = prolong_dynamic(sp_rp_c, pi_table, n_max);    renormalize!(sp_dyn_rp)

    # TV distances
    function tv2(sp_a, p_true)
        tv = 0.0
        for (i, s) in enumerate(sp_fine2.states)
            idx = get(sp_a.index, s, 0)
            pa  = idx != 0 ? sp_a.probs[idx] : 0.0
            tv += abs(p_true[i] - pa)
        end
        for (i, s) in enumerate(sp_a.states)
            get(sp_fine2.index, s, 0) != 0 && continue
            tv += sp_a.probs[i]
        end
        tv / 2
    end
    tv_b = tv2(sp_binom_rp, p_ref)
    tv_d = tv2(sp_dyn_rp,   p_ref)
    @printf("  TV(Binom, FSP)   = %.4f\n", tv_b)
    @printf("  TV(Dyn-π, FSP)   = %.4f\n", tv_d)

    # Build probability matrices for plotting
    nplot = min(n_max, 200)
    function sp_to_matrix(sp)
        Z = zeros(nplot + 1, nplot + 1)
        for (i, s) in enumerate(sp.states)
            t2 = Tuple(s)
            n1, n2 = t2[1], t2[2]
            (n1 <= nplot && n2 <= nplot) || continue
            Z[n1 + 1, n2 + 1] = sp.probs[i]
        end
        Z
    end

    # For the reference: build a temporary sp with current probs
    sp_ref_sp = StateSpace{CartesianIndex{2}, Float64}()
    for i in eachindex(sp_fine2.states)
        add_state!(sp_ref_sp, sp_fine2.states[i], p_ref[i])
    end

    Z_true  = sp_to_matrix(sp_ref_sp)
    Z_binom = sp_to_matrix(sp_binom_rp)
    Z_dyn   = sp_to_matrix(sp_dyn_rp)

    cmax = maximum(Z_true) * 0.7   # clip for visibility

    fig = Figure(size = (560, 220))

    ax_labels = [
        "(a) Direct FSP (truth)",
        "(b) Binom–R–P  (TV=$(round(tv_b, digits=3)))",
        "(c) Dyn-π–R–P  (TV=$(round(tv_d, digits=3)))",
    ]
    mats = [Z_true, Z_binom, Z_dyn]

    for (j, (lbl, Z)) in enumerate(zip(ax_labels, mats))
        ax = Axis(fig[1, j];
            title  = lbl,
            xlabel = "n₁",
            ylabel = j == 1 ? "n₂" : "",
            yticklabelsvisible = (j == 1),
            xticks = [0, n_low, n_high, nplot],
            yticks = [0, n_low, n_high, nplot],
        )
        hm = heatmap!(ax, 0:nplot, 0:nplot, Z;
                      colormap = :viridis, colorrange = (0, cmax),
                      lowclip = :black)
        # Mark the initial condition
        scatter!(ax, [Float64(n_low)], [Float64(n_high)]; color = :white,
                 marker = :cross, markersize = 8, strokewidth = 1.5, strokecolor = :white)
        if j == 3
            Colorbar(fig[1, 4], hm; label = "p(n₁,n₂)", labelsize = SMALL,
                     ticklabelsize = SMALL, width = 10)
        end
    end

    colgap!(fig.layout, 4)
    save_fig("fig3_heatmap_comparison", fig)
end

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: TV time series — K=2 coarse MSA
# ══════════════════════════════════════════════════════════════════════════════

println("\n── Fig 4: MSA TV decay ──")

let
    tv_times = Float64[]
    tv_binom_vals = Float64[]
    tv_dyn_vals   = Float64[]

    for t in sort(collect(keys(snaps_binom)))
        # Set sp_fine2.probs to the FSP snapshot at t
        sp_fine2.probs .= probs_fine2[t]

        # Binom prolongation
        states_b, probs_b = snaps_binom[t]
        sp_tmp_b = StateSpace{CartesianIndex{1}, Float64}()
        for (i, s) in enumerate(states_b); add_state!(sp_tmp_b, s, probs_b[i]); end
        sp_prol_b = prolong(sp_tmp_b, Val(2); binom_tol = 1e-6); renormalize!(sp_prol_b)

        # Dyn-π prolongation
        states_d, probs_d = snaps_dyn[t]
        sp_tmp_d = StateSpace{CartesianIndex{1}, Float64}()
        for (i, s) in enumerate(states_d); add_state!(sp_tmp_d, s, probs_d[i]); end
        sp_prol_d = prolong_dynamic(sp_tmp_d, pi_table, n_max); renormalize!(sp_prol_d)

        push!(tv_times,       t)
        push!(tv_binom_vals,  tv_vs_fine2(sp_prol_b))
        push!(tv_dyn_vals,    tv_vs_fine2(sp_prol_d))
    end

    fig = Figure(size = (340, 240))
    ax = Axis(fig[1, 1];
        xlabel = "time  t",
        ylabel = "Total variation distance",
        title  = "K=2 MSA accuracy  (IC: n_low–n_high front)",
    )
    ylims!(ax, 0.0, 1.05)

    hlines!(ax, 1.0; color = C_GRAY, linestyle = :dot, linewidth = 0.8,
            label = "TV = 1 (completely wrong)")
    hlines!(ax, 0.0; color = C_GRAY, linestyle = :dot, linewidth = 0.8)

    scatterlines!(ax, tv_times, tv_binom_vals;
                  color = C_ORANGE, marker = :circle, linewidth = LW,
                  label = "Coarse Binomial-π")
    scatterlines!(ax, tv_times, tv_dyn_vals;
                  color = C_GREEN, marker = :utriangle, linewidth = LW,
                  label = "Coarse Dynamic-π")

    axislegend(ax; position = :rb, framevisible = true)

    save_fig("fig4_msa_tv_decay", fig)
end

# ══════════════════════════════════════════════════════════════════════════════
# K=4 V-CYCLE COMPUTATION  (Figs 5, 6)
# ══════════════════════════════════════════════════════════════════════════════

skip_vcycle = haskey(ENV, "SKIP_VCYCLE")
if skip_vcycle
    println("\n── Skipping K=4 V-cycle (SKIP_VCYCLE set) ──")
else
    println("\n── Running K=4 V-cycle benchmark ──")

    const K_vc    = 4
    const h_vc    = 1.0 / K_vc
    const n_max_vc = 180
    const dt_vc   = 0.2
    const t_solve_vc = [0.5, 1.0, 2.0]
    const t_max_vc   = maximum(t_solve_vc)
    const n_steps_vc = round(Int, t_max_vc / dt_vc)

    fine_grid_vc   = VoxelGrid(K_vc, h_vc, 0)
    coarse_grid_vc = coarsen(fine_grid_vc)

    u0_vc = CartesianIndex(n_low, n_low, n_high, n_high)

    # Ground truth: two independent K=2 pair FSPs
    pair_grid_vc = VoxelGrid(2, h_vc, 0)
    pair_sys_vc  = build_schlogl_rdme_system(model, pair_grid_vc)

    sp_p1 = StateSpace{CartesianIndex{2}, Float64}()
    sp_p2 = StateSpace{CartesianIndex{2}, Float64}()
    for n1 in 0:n_max_vc, n2 in 0:n_max_vc
        add_state!(sp_p1, CartesianIndex(n1, n2), 0.0)
        add_state!(sp_p2, CartesianIndex(n1, n2), 0.0)
    end
    sp_p1.probs[sp_p1.index[CartesianIndex(n_low,  n_low)]]  = 1.0
    sp_p2.probs[sp_p2.index[CartesianIndex(n_high, n_high)]] = 1.0
    A_p1, _ = build_generator(sp_p1, pair_sys_vc, rates, 0.0)
    A_p2, _ = build_generator(sp_p2, pair_sys_vc, rates, 0.0)

    probs_gt_vc = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    let p1 = copy(sp_p1.probs), p2 = copy(sp_p2.probs), tprev = 0.0
        for t in t_solve_vc
            p1 = max.(expv(t - tprev, A_p1, p1; m = 50), 0.0); p1 ./= sum(p1)
            p2 = max.(expv(t - tprev, A_p2, p2; m = 50), 0.0); p2 ./= sum(p2)
            probs_gt_vc[t] = (copy(p1), copy(p2))
            tprev = t
        end
    end
    @printf("  Ground truth (2×K=2): done\n")

    # dynamic-π table for K=4 pair size
    pi_table_vc = compute_dynamic_pi(sp_p1, A_p1; n_max = n_max_vc)

    # Init coarse (K=2) for coarse-only methods
    nc1_vc = n_low + n_low; nc2_vc = n_high + n_high
    coarse_bc_vc = s -> rdme_bc(s, 2 * n_max_vc)

    function init_coarse_vc()
        sp = StateSpace{CartesianIndex{2}, Float64}()
        add_state!(sp, CartesianIndex(nc1_vc, nc2_vc), 1.0)
        sp
    end

    # Voxel means helper
    function voxel_means_k4(states, probs)
        mu = zeros(4)
        for (i, s) in enumerate(states)
            t4 = Tuple(s)
            for k in 1:4; mu[k] += probs[i] * t4[k]; end
        end
        mu
    end

    # TV vs pair GT helper
    function tv_vs_gt_vc(fine_states, fine_probs, pp1, pp2)
        overlap = 0.0
        for (i, s) in enumerate(fine_states)
            t4 = Tuple(s)
            i1 = get(sp_p1.index, CartesianIndex(t4[1], t4[2]), 0)
            i2 = get(sp_p2.index, CartesianIndex(t4[3], t4[4]), 0)
            p_gt = (i1 != 0 ? pp1[i1] : 0.0) * (i2 != 0 ? pp2[i2] : 0.0)
            overlap += min(fine_probs[i], p_gt)
        end
        1.0 - overlap
    end

    τ_pre_vc = 1.0 / (D / h_vc^2 * n_low + 1.0)

    # Storage
    result_means = Dict{String, Dict{Float64, Vector{Float64}}}()
    result_tv    = Dict{String, Dict{Float64, Float64}}()
    result_sizes = Dict{String, Vector{Int}}()
    result_times = Dict{String, Float64}()
    for m in ["coarse_binom","coarse_dyn","vc_binom","vc_dyn"]
        result_means[m] = Dict{Float64, Vector{Float64}}()
        result_tv[m]    = Dict{Float64, Float64}()
        result_sizes[m] = Int[]
    end

    # ── Method 1: Coarse Binomial ─────────────────────────────────────────────
    @printf("  Running Coarse Binom... ")
    coarse_sys_vc_b = build_schlogl_coarse_system(model, coarse_grid_vc, fine_grid_vc)
    sp_cb = init_coarse_vc()
    expand!(sp_cb, coarse_sys_vc_b, coarse_bc_vc; depth = 1)
    t0_ = time()
    let t_cur = 0.0
        for step in 1:n_steps_vc
            dt_step = min(dt_vc, t_max_vc - t_cur)
            expand!(sp_cb, coarse_sys_vc_b, coarse_bc_vc; depth = 1)
            A, _ = build_generator(sp_cb, coarse_sys_vc_b, rates, t_cur)
            sp_cb.probs .= expv(dt_step, A, sp_cb.probs; m = 30)
            t_cur += dt_step
            renormalize!(sp_cb); prune_threshold!(sp_cb, 1e-12)
            push!(result_sizes["coarse_binom"], length(sp_cb))
            for t in t_solve_vc
                if abs(t_cur - t) < dt_vc / 2 && !haskey(result_means["coarse_binom"], t)
                    sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
                    for (i, s) in enumerate(sp_cb.states)
                        add_state!(sp_tmp, s, sp_cb.probs[i])
                    end
                    sp_prol = prolong(sp_tmp, Val(K_vc); binom_tol = 1e-6)
                    renormalize!(sp_prol)
                    result_means["coarse_binom"][t] = voxel_means_k4(sp_prol.states, sp_prol.probs)
                    pp1, pp2 = probs_gt_vc[t]
                    result_tv["coarse_binom"][t] = tv_vs_gt_vc(sp_prol.states, sp_prol.probs, pp1, pp2)
                end
            end
        end
    end
    result_times["coarse_binom"] = time() - t0_
    @printf("done (%.1fs)\n", result_times["coarse_binom"])

    # ── Method 2: Coarse Dynamic-π ────────────────────────────────────────────
    @printf("  Running Coarse Dyn-π... ")
    coarse_sys_vc_d = build_schlogl_coarse_system_dynamic(model, coarse_grid_vc, fine_grid_vc, pi_table_vc)
    sp_cd = init_coarse_vc()
    expand!(sp_cd, coarse_sys_vc_d, coarse_bc_vc; depth = 1)
    t0_ = time()
    let t_cur = 0.0
        for step in 1:n_steps_vc
            dt_step = min(dt_vc, t_max_vc - t_cur)
            expand!(sp_cd, coarse_sys_vc_d, coarse_bc_vc; depth = 1)
            A, _ = build_generator(sp_cd, coarse_sys_vc_d, rates, t_cur)
            sp_cd.probs .= expv(dt_step, A, sp_cd.probs; m = 30)
            t_cur += dt_step
            renormalize!(sp_cd); prune_threshold!(sp_cd, 1e-12)
            push!(result_sizes["coarse_dyn"], length(sp_cd))
            for t in t_solve_vc
                if abs(t_cur - t) < dt_vc / 2 && !haskey(result_means["coarse_dyn"], t)
                    sp_tmp = StateSpace{CartesianIndex{2}, Float64}()
                    for (i, s) in enumerate(sp_cd.states)
                        add_state!(sp_tmp, s, sp_cd.probs[i])
                    end
                    sp_prol = prolong_dynamic(sp_tmp, pi_table_vc, n_max_vc)
                    renormalize!(sp_prol)
                    result_means["coarse_dyn"][t] = voxel_means_k4(sp_prol.states, sp_prol.probs)
                    pp1, pp2 = probs_gt_vc[t]
                    result_tv["coarse_dyn"][t] = tv_vs_gt_vc(sp_prol.states, sp_prol.probs, pp1, pp2)
                end
            end
        end
    end
    result_times["coarse_dyn"] = time() - t0_
    @printf("done (%.1fs)\n", result_times["coarse_dyn"])

    # ── Method 3: V-cycle Binom ───────────────────────────────────────────────
    @printf("  Running V-cycle Binom... ")
    sp_vb = StateSpace{CartesianIndex{4}, Float64}()
    add_state!(sp_vb, u0_vc, 1.0)
    t0_ = time()
    let sp = sp_vb, t_cur = 0.0
        for step in 1:n_steps_vc
            dt_step = min(dt_vc, t_max_vc - t_cur)
            sp = two_level_vcycle_schlogl(sp, model, fine_grid_vc, coarse_grid_vc,
                                          pi_table_vc, rates, t_cur, dt_step;
                                          τ_pre = τ_pre_vc, krylov_m = 30,
                                          use_dynamic_pi = false, weight_tol = 1e-8,
                                          coarse_n_max = 2 * n_max_vc,
                                          expand_coarse = true, coarse_expand_depth = 1)
            t_cur += dt_step
            renormalize!(sp); prune_threshold!(sp, 1e-12)
            push!(result_sizes["vc_binom"], length(sp))
            for t in t_solve_vc
                if abs(t_cur - t) < dt_vc / 2 && !haskey(result_means["vc_binom"], t)
                    result_means["vc_binom"][t] = voxel_means_k4(sp.states, sp.probs)
                    pp1, pp2 = probs_gt_vc[t]
                    result_tv["vc_binom"][t] = tv_vs_gt_vc(sp.states, sp.probs, pp1, pp2)
                end
            end
        end
        global sp_vb = sp
    end
    result_times["vc_binom"] = time() - t0_
    @printf("done (%.1fs, %d fine states)\n", result_times["vc_binom"], length(sp_vb))

    # ── Method 4: V-cycle Dyn-π ───────────────────────────────────────────────
    @printf("  Running V-cycle Dyn-π... ")
    sp_vd = StateSpace{CartesianIndex{4}, Float64}()
    add_state!(sp_vd, u0_vc, 1.0)
    t0_ = time()
    let sp = sp_vd, t_cur = 0.0
        for step in 1:n_steps_vc
            dt_step = min(dt_vc, t_max_vc - t_cur)
            sp = two_level_vcycle_schlogl(sp, model, fine_grid_vc, coarse_grid_vc,
                                          pi_table_vc, rates, t_cur, dt_step;
                                          τ_pre = τ_pre_vc, krylov_m = 30,
                                          use_dynamic_pi = true, weight_tol = 1e-8,
                                          coarse_n_max = 2 * n_max_vc,
                                          expand_coarse = true, coarse_expand_depth = 1)
            t_cur += dt_step
            renormalize!(sp); prune_threshold!(sp, 1e-12)
            push!(result_sizes["vc_dyn"], length(sp))
            for t in t_solve_vc
                if abs(t_cur - t) < dt_vc / 2 && !haskey(result_means["vc_dyn"], t)
                    result_means["vc_dyn"][t] = voxel_means_k4(sp.states, sp.probs)
                    pp1, pp2 = probs_gt_vc[t]
                    result_tv["vc_dyn"][t] = tv_vs_gt_vc(sp.states, sp.probs, pp1, pp2)
                end
            end
        end
        global sp_vd = sp
    end
    result_times["vc_dyn"] = time() - t0_
    @printf("done (%.1fs, %d fine states)\n", result_times["vc_dyn"], length(sp_vd))

    # ── Print benchmark table ──────────────────────────────────────────────────
    println()
    println("  t     GT(μ1,μ2|μ3,μ4)         C-Binom TV   C-Dyn TV   VC-Binom TV   VC-Dyn TV")
    println("  ─────────────────────────────────────────────────────────────────────────────")
    for t in t_solve_vc
        pp1, pp2 = probs_gt_vc[t]
        μ1_gt = sum(Tuple(sp_p1.states[i])[1] * pp1[i] for i in eachindex(pp1))
        μ2_gt = sum(Tuple(sp_p1.states[i])[2] * pp1[i] for i in eachindex(pp1))
        μ3_gt = sum(Tuple(sp_p2.states[i])[1] * pp2[i] for i in eachindex(pp2))
        μ4_gt = sum(Tuple(sp_p2.states[i])[2] * pp2[i] for i in eachindex(pp2))
        tv_cb = get(result_tv["coarse_binom"], t, NaN)
        tv_cd = get(result_tv["coarse_dyn"],   t, NaN)
        tv_vb = get(result_tv["vc_binom"],     t, NaN)
        tv_vd = get(result_tv["vc_dyn"],       t, NaN)
        @printf("  %4.1f  (%4.1f,%4.1f|%4.1f,%4.1f)   %.4f       %.4f     %.4f        %.4f\n",
                t, μ1_gt, μ2_gt, μ3_gt, μ4_gt, tv_cb, tv_cd, tv_vb, tv_vd)
    end

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE 5: Voxel means over time
    # ══════════════════════════════════════════════════════════════════════════

    println("\n── Fig 5: K=4 voxel means ──")

    let
        times_vc = sort(collect(keys(probs_gt_vc)))

        fig = Figure(size = (520, 360))
        Label(fig[0, :], "K=4 V-cycle benchmark  (symmetric IC, D=$(D))",
              fontsize = MED)

        # Pair 1 (voxels 1,2) left panel; Pair 2 (voxels 3,4) right panel
        ax1 = Axis(fig[1, 1];
            xlabel = "time  t",
            ylabel = "mean molecule count",
            title  = "Pair 1  (voxels 1–2,  IC=$(n_low),$(n_low))",
        )
        ax2 = Axis(fig[1, 2];
            xlabel = "time  t",
            ylabel = "",
            yticklabelsvisible = false,
            title  = "Pair 2  (voxels 3–4,  IC=$(n_high),$(n_high))",
        )

        method_labels = ["Coarse Binom", "Coarse Dyn-π", "V-cyc Binom", "V-cyc Dyn-π"]
        method_keys   = ["coarse_binom", "coarse_dyn", "vc_binom", "vc_dyn"]
        colors        = [C_ORANGE, C_GREEN, C_RED, C_PURPLE]
        styles        = [:solid, :dash, :solid, :dash]

        # Ground truth
        gt_μ1 = [sum(Tuple(sp_p1.states[i])[1]*probs_gt_vc[t][1][i] for i in eachindex(sp_p1.states)) for t in times_vc]
        gt_μ2 = [sum(Tuple(sp_p1.states[i])[2]*probs_gt_vc[t][1][i] for i in eachindex(sp_p1.states)) for t in times_vc]
        gt_μ3 = [sum(Tuple(sp_p2.states[i])[1]*probs_gt_vc[t][2][i] for i in eachindex(sp_p2.states)) for t in times_vc]
        gt_μ4 = [sum(Tuple(sp_p2.states[i])[2]*probs_gt_vc[t][2][i] for i in eachindex(sp_p2.states)) for t in times_vc]

        # Prepend t=0 (delta on IC)
        t_plot = vcat([0.0], times_vc)
        for (ax, gt_a, gt_b, ic_a, ic_b) in [
                (ax1, gt_μ1, gt_μ2, Float64(n_low),  Float64(n_low)),
                (ax2, gt_μ3, gt_μ4, Float64(n_high), Float64(n_high))]
            lines!(ax, t_plot, vcat([ic_a], gt_a);
                   color = C_BLUE, linewidth = LW + 0.4, linestyle = :solid,
                   label = "GT (pair FSP)")
            lines!(ax, t_plot, vcat([ic_b], gt_b);
                   color = C_BLUE, linewidth = LW + 0.4, linestyle = :dashdot)
        end

        for (label, key, col, sty) in zip(method_labels, method_keys, colors, styles)
            ts = sort(collect(keys(result_means[key])))
            isempty(ts) && continue
            mu_mat = hcat([result_means[key][t] for t in ts]...)  # 4×nT

            lines!(ax1, vcat([0.0], ts), vcat([Float64(n_low)],  mu_mat[1, :]);
                   color = col, linewidth = LW, linestyle = sty, label = label)
            lines!(ax1, vcat([0.0], ts), vcat([Float64(n_low)],  mu_mat[2, :]);
                   color = col, linewidth = LW, linestyle = sty)
            lines!(ax2, vcat([0.0], ts), vcat([Float64(n_high)], mu_mat[3, :]);
                   color = col, linewidth = LW, linestyle = sty)
            lines!(ax2, vcat([0.0], ts), vcat([Float64(n_high)], mu_mat[4, :]);
                   color = col, linewidth = LW, linestyle = sty)
        end

        # Horizontal guides at fixed points
        for ax in [ax1, ax2]
            hlines!(ax, [Float64(n_low), Float64(n_high)];
                    color = (C_GRAY, 0.4), linestyle = :dot, linewidth = 0.8)
        end

        Legend(fig[2, :], ax1; orientation = :horizontal, framevisible = true,
               labelsize = SMALL, nbanks = 2)
        rowgap!(fig.layout, 6)
        colgap!(fig.layout, 8)

        save_fig("fig5_vcycle_means", fig)
    end

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE 6: Fine state-space size + TV over time (combined)
    # ══════════════════════════════════════════════════════════════════════════

    println("\n── Fig 6: State-space size + TV ──")

    let
        step_times = range(dt_vc, stop = t_max_vc, length = n_steps_vc)
        snap_times = sort(collect(keys(result_tv["coarse_binom"])))

        fig = Figure(size = (480, 340))

        # panel (a): state space size
        ax_sz = Axis(fig[1, 1];
            xlabel = "time  t",
            ylabel = "fine state-space  |S|",
            title  = "(a) Fine state-space growth",
        )

        for (label, key, col, sty) in [
                ("V-cyc Binom", "vc_binom", C_RED,    :solid),
                ("V-cyc Dyn-π", "vc_dyn",   C_PURPLE, :dash),
                ("Coarse Binom", "coarse_binom", C_ORANGE, :solid),
                ("Coarse Dyn-π", "coarse_dyn",   C_GREEN,  :dash)]
            sz = result_sizes[key]
            isempty(sz) && continue
            ts = collect(range(dt_vc, stop = t_max_vc, length = length(sz)))
            lines!(ax_sz, ts, Float64.(sz); color = col, linewidth = LW,
                   linestyle = sty, label = label)
        end

        axislegend(ax_sz; position = :lt, framevisible = true, labelsize = SMALL)

        # panel (b): TV over time
        ax_tv = Axis(fig[1, 2];
            xlabel = "time  t",
            ylabel = "TV vs pair-level ground truth",
            title  = "(b) Accuracy (TV distance)",
        )
        ylims!(ax_tv, 0.0, 1.05)
        hlines!(ax_tv, 1.0; color = (C_GRAY, 0.4), linestyle = :dot, linewidth = 0.8)
        hlines!(ax_tv, 0.0; color = (C_GRAY, 0.4), linestyle = :dot, linewidth = 0.8)

        for (label, key, col, sty, mk) in [
                ("Coarse Binom", "coarse_binom", C_ORANGE, :solid, :circle),
                ("Coarse Dyn-π", "coarse_dyn",   C_GREEN,  :dash,  :utriangle),
                ("V-cyc Binom",  "vc_binom",     C_RED,    :solid, :rect),
                ("V-cyc Dyn-π",  "vc_dyn",       C_PURPLE, :dash,  :diamond)]
            ts = sort(collect(keys(result_tv[key])))
            isempty(ts) && continue
            tv_vals = [result_tv[key][t] for t in ts]
            scatterlines!(ax_tv, ts, tv_vals; color = col, linewidth = LW,
                          linestyle = sty, marker = mk, label = label)
        end

        axislegend(ax_tv; position = :rb, framevisible = true, labelsize = SMALL)
        colgap!(fig.layout, 10)

        save_fig("fig6_vcycle_accuracy", fig)
    end

    println("\n── Timing summary ──")
    for (lbl, key) in [("Coarse Binom","coarse_binom"),("Coarse Dyn-π","coarse_dyn"),
                        ("V-cyc Binom", "vc_binom"),   ("V-cyc Dyn-π","vc_dyn")]
        @printf("  %-16s  %.1fs\n", lbl, result_times[key])
    end
end   # if !skip_vcycle

# ══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^65)
println("All figures saved to paper/figures/")
println("="^65)
