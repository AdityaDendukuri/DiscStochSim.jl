"""
fig_vcycle_demo.jl

Experimental demonstration of the V-cycle prolongation step.

IC: 70% mass at (n_hi, n_lo) = (164, 38)  and  30% at (n_hi−10, n_lo+10) = (154, 48).
Both modes lie on the anti-diagonal N = n_hi+n_lo = 202.
This joint is NOT a product of its marginals.

The marginals are: m₁(164)=0.70, m₁(154)=0.30; m₂(38)=0.70, m₂(48)=0.30.
Product-convolution at N=202 gives p(164|N)∝0.7²=0.49, p(154|N)∝0.3²=0.09 → normalises to
0.845 and 0.155 — the dominant mode is over-weighted and the minority under-weighted.

Exact multiplicative correction gives the true conditional (0.70 / 0.30) with zero error.
"""

using DiscStochSim, CairoMakie, ExponentialUtilities, Printf, LinearAlgebra

using DiscStochSim: StateSpace, add_state!, renormalize!, expand!, build_generator,
                    DiscreteStochasticSystem

_sg    = DiscStochSim._single_voxel_generator
_clean = DiscStochSim._clean!
_pml   = DiscStochSim._prolong_multilevel

# ── Model & parameters ─────────────────────────────────────────────────────
const D = 2.0
model   = SchloglModel1D(D, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)

const N_MAX = 220
const NMAX2 = 2 * N_MAX
const STOL  = 1e-10
const KM    = 30

@printf("Schlögl: n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)
@printf("Dominant N = %d+%d = %d\n", n_hi, n_lo, n_hi+n_lo)

# External boundary rates
in1 = D * Float64(n_hi);  in2 = D * Float64(n_lo)
out1 = D;                  out2 = D

# ── IC: correlated joint concentrated at N = n_hi+n_lo = 202 ──────────────
# Two modes on the same anti-diagonal → NOT a product of marginals.
# Δ controls the spread: Δ=0 → single mode; Δ=10 → two distinct modes
const Δ = 10        # spread (molecules shifted between voxels)
const W1 = 0.70     # weight of dominant mode
const W2 = 1 - W1   # weight of minority mode

N_dom = n_hi + n_lo   # = 202

sp = StateSpace{CartesianIndex{2}, Float64}()
add_state!(sp, CartesianIndex(n_hi,     n_lo),     W1)   # dominant
add_state!(sp, CartesianIndex(n_hi - Δ, n_lo + Δ), W2)   # minority
renormalize!(sp)

@printf("IC: %.0f%% at (%d,%d)  +  %.0f%% at (%d,%d)\n",
    100W1, n_hi, n_lo, 100W2, n_hi-Δ, n_lo+Δ)
@printf("Both on anti-diagonal N=%d — NOT a product of marginals\n", N_dom)

# ── Pre-smooth (minimal: just expand support) ──────────────────────────────
# Use very short τ_pre so the conditional structure is preserved.
# The failure of product-convolution is about the WITHIN-N conditional structure;
# fast diffusion would destroy it by equilibrating to Binomial(N,½).
pair_intra = DiscreteStochasticSystem{CartesianIndex{2}}(
    [CartesianIndex(-1,1), CartesianIndex(1,-1)],
    [(x,r,t) -> D*max(0.0, Float64(Tuple(x)[1])),
     (x,r,t) -> D*max(0.0, Float64(Tuple(x)[2]))])
bc2 = x -> all(c -> 0 ≤ c ≤ N_MAX, Tuple(x))

dt    = 0.5
τ_pre = 0.002 * dt   # very short — preserves conditional structure

A_intra, = build_generator(sp, pair_intra, nothing, 0.0)
sp.probs .= expv(τ_pre, A_intra, sp.probs; m=min(KM, size(A_intra,1)))
_clean(sp.probs, STOL)
renormalize!(sp)
@printf("Post-pre-smooth: %d states\n", length(sp))

# ── Restriction ────────────────────────────────────────────────────────────
p_pre = zeros(NMAX2+1)
for (ci, prob) in zip(sp.states, sp.probs)
    N = sum(Tuple(ci)); N ≤ NMAX2 && (p_pre[N+1] += prob)
end
p_pre ./= max(sum(p_pre), 1e-100)

@printf("p_pre: dominant N=%d has p=%.4f\n", N_dom, p_pre[N_dom+1])

# ── Coarse solve ───────────────────────────────────────────────────────────
in_rate   = in1 + in2
out_coeff = (out1 + out2) / 2
Q_c = _sg(model, in_rate, out_coeff, NMAX2; n_sites=2)
p_post = expv(Float64(dt), Q_c, p_pre; m=min(KM, size(Q_c,1)))
_clean(p_post, STOL)

# Which N values appeared/changed?
new_N = [N for N in 0:NMAX2 if p_post[N+1] > STOL && p_pre[N+1] ≤ STOL]
@printf("Coarse solve: %d new N values created\n", length(new_N))

# ── Old prolongation: product-convolution ──────────────────────────────────
function prolong_product(sp, p_post, n_max, stol)
    m1 = zeros(n_max+1); m2 = zeros(n_max+1)
    for (ci, prob) in zip(sp.states, sp.probs)
        n1, n2 = Tuple(ci); m1[n1+1] += prob; m2[n2+1] += prob
    end
    s = sum(m1); s > 0 && (m1 ./= s)
    s = sum(m2); s > 0 && (m2 ./= s)
    M = zeros(n_max+1, n_max+1)
    for N in 0:2*n_max
        p_post[N+1] > stol || continue
        Z = sum(m1[n1+1]*m2[N-n1+1] for n1 in max(0,N-n_max):min(N,n_max))
        Z < 1e-300 && continue
        coeff = p_post[N+1] / Z
        for n1 in max(0,N-n_max):min(N,n_max)
            M[n1+1, N-n1+1] += coeff * m1[n1+1] * m2[N-n1+1]
        end
    end
    M
end

joint_old = prolong_product(sp, p_post, N_MAX, STOL)

# ── New prolongation: multi-level stochastic ───────────────────────────────
sp_new = _pml(sp, p_pre, p_post, in1, in2, out1, out2,
              model, N_MAX, STOL, 200_000)
joint_new = zeros(N_MAX+1, N_MAX+1)
if sp_new !== nothing
    for (ci, prob) in zip(sp_new.states, sp_new.probs)
        n1, n2 = Tuple(ci); joint_new[n1+1, n2+1] += prob
    end
end

# ── Reference: exact multiplicative correction ─────────────────────────────
joint_correct = zeros(N_MAX+1, N_MAX+1)
for (ci, prob) in zip(sp.states, sp.probs)
    n1, n2 = Tuple(ci); N = n1+n2
    p_pre[N+1] > 1e-300 || continue
    joint_correct[n1+1, n2+1] = prob * p_post[N+1] / p_pre[N+1]
end

# ── Pre-smooth joint matrix ────────────────────────────────────────────────
joint_pre = zeros(N_MAX+1, N_MAX+1)
for (ci, prob) in zip(sp.states, sp.probs)
    n1, n2 = Tuple(ci); joint_pre[n1+1, n2+1] = prob
end

# ── Error metrics — compared only at EXISTING N values (where reference is defined) ──
# joint_correct is 0 at new N values (coarse solve created them); comparing there
# would penalise the new method for correctly populating new N via L/R.
# Instead, measure L1 only at N=202 (the existing dominant N).
function l1_at_N(M_test, M_ref, N, n_max)
    err = 0.0
    for n1 in max(0,N-n_max):min(N,n_max)
        n2 = N-n1; n2 ≤ n_max || continue
        err += abs(M_test[n1+1,n2+1] - M_ref[n1+1,n2+1])
    end
    err
end

err_old = l1_at_N(joint_old, joint_correct, N_dom, N_MAX)
err_new = l1_at_N(joint_new, joint_correct, N_dom, N_MAX)
@printf("\nL1 error at N=%d (existing, where reference is exact):\n", N_dom)
@printf("  Old (product-conv):  %.5f\n", err_old)
@printf("  New (multi-level):   %.5f\n", err_new)
@printf("  Improvement:         %.1f×\n", err_old / max(err_new, 1e-12))

# Conditional at dominant N
function cond_at_N(M, N, n_max)
    f = zeros(n_max+1)
    for n1 in max(0,N-n_max):min(N,n_max)
        n2 = N-n1; n2 ≤ n_max && (f[n1+1] = M[n1+1, n2+1])
    end
    s = sum(f); s > 1e-300 && (f ./= s); f
end

c_pre = cond_at_N(joint_pre,     N_dom, N_MAX)
c_old = cond_at_N(joint_old,     N_dom, N_MAX)
c_new = cond_at_N(joint_new,     N_dom, N_MAX)
c_ref = cond_at_N(joint_correct, N_dom, N_MAX)

@printf("\nConditional at N=%d (dominant):\n", N_dom)
@printf("  Pre-smooth true: p(n₁=%d)=%.3f  p(n₁=%d)=%.3f\n",
    n_hi, c_pre[n_hi+1], n_hi-Δ, c_pre[n_hi-Δ+1])
@printf("  Old prod-conv:   p(n₁=%d)=%.3f  p(n₁=%d)=%.3f\n",
    n_hi, c_old[n_hi+1], n_hi-Δ, c_old[n_hi-Δ+1])
@printf("  New multi-level: p(n₁=%d)=%.3f  p(n₁=%d)=%.3f\n",
    n_hi, c_new[n_hi+1], n_hi-Δ, c_new[n_hi-Δ+1])
@printf("  Reference:       p(n₁=%d)=%.3f  p(n₁=%d)=%.3f\n",
    n_hi, c_ref[n_hi+1], n_hi-Δ, c_ref[n_hi-Δ+1])

# ── Colour scale: max over all panels ─────────────────────────────────────
vmax = max(maximum(joint_pre), maximum(joint_old),
           maximum(joint_new), maximum(joint_correct))

# ── Zoom region around the two modes ──────────────────────────────────────
lo = n_hi - Δ - 8;  hi = n_hi + 8
lo2 = n_lo - 8;     hi2 = n_lo + Δ + 8
rng1 = lo:hi;       rng2 = lo2:hi2

function zoom(M) M[rng1 .+ 1, rng2 .+ 1] end

# ── Figure 1: 2D heatmaps ─────────────────────────────────────────────────
fig1 = Figure(size=(1400, 420))
Label(fig1[0,1:4],
    "V-cycle prolongation  ·  Schlögl pair  " *
    @sprintf("(%.0f%% at (%d,%d), %.0f%% at (%d,%d),  dt=%.2f)",
        100W1, n_hi, n_lo, 100W2, n_hi-Δ, n_lo+Δ, dt),
    fontsize=12, font=:bold)

function make_hm_ax(pos, title)
    Axis(fig1[pos...]; title, xlabel="n₁ (voxel 1)", ylabel="n₂ (voxel 2)",
         titlesize=11, xlabelsize=9, ylabelsize=9,
         xticks=[n_hi-Δ, n_hi], yticks=[n_lo, n_lo+Δ])
end

panels = [
    ("Pre-smooth joint\n(input to prolongation)", joint_pre, "IC (input)"),
    ("Old prolong\n(product-convolution)", joint_old,
     @sprintf("L1@N=%d: %.4f", N_dom, err_old)),
    ("New prolong\n(multi-level stochastic)", joint_new,
     @sprintf("L1@N=%d: %.2e", N_dom, err_new)),
    ("Reference\n(exact multiplicative corr.)", joint_correct, "reference"),
]

axs = [make_hm_ax([1,col], t) for (col,(t,_,_)) in enumerate(panels)]
hms = [heatmap!(axs[col], rng1, rng2, zoom(M);
                colormap=:hot, colorrange=(0, vmax), interpolate=false)
       for (col, (_,M,_)) in enumerate(panels)]

for (col, (_, _, lbl)) in enumerate(panels)
    # mark the two modes
    scatter!(axs[col], [n_hi, n_hi-Δ], [n_lo, n_lo+Δ];
             color=:cyan, marker=:circle, markersize=10, strokewidth=1.5,
             strokecolor=:white)
    text!(axs[col], n_hi+0.5, n_lo+0.5;
          text=@sprintf("%.2f", panels[col][2][n_hi+1, n_lo+1]),
          fontsize=9, color=:cyan, align=(:left,:bottom))
    text!(axs[col], n_hi-Δ+0.5, n_lo+Δ+0.5;
          text=@sprintf("%.2f", panels[col][2][n_hi-Δ+1, n_lo+Δ+1]),
          fontsize=9, color=:cyan, align=(:left,:bottom))
    text!(axs[col], lo+0.3, hi2-0.5; text=lbl, fontsize=8, color=:white,
          align=(:left,:top))
end

Colorbar(fig1[1,5], hms[end]; label="p(n₁,n₂)", labelsize=9,
         ticklabelsize=8, width=14)
colgap!(fig1.layout, 6)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
save(joinpath(outdir, "fig_vcycle_demo.png"), fig1; px_per_unit=3)
println("\nSaved → fig_vcycle_demo.png")

# ── Figure 2: conditional at N_dom ────────────────────────────────────────
fig2 = Figure(size=(760, 380))
ax = Axis(fig2[1,1];
    title=@sprintf("Conditional p(n₁ | N=%d): product-conv vs multi-level vs reference", N_dom),
    xlabel="n₁  (molecules in voxel 1)", ylabel=@sprintf("p(n₁ | N=%d)", N_dom),
    titlesize=11)

n_rng = (n_hi-Δ-3):(n_hi+3)

# True conditional: two spikes
barplot!(ax, n_rng, c_ref[n_rng .+ 1];
         color=(:black,0.25), label="Reference (exact)")
lines!(ax, n_rng, c_old[n_rng .+ 1]; color=:red, linewidth=2.5,
       linestyle=:dash, label="Old: product-convolution")
lines!(ax, n_rng, c_new[n_rng .+ 1]; color=:steelblue, linewidth=2.5,
       label="New: multi-level")
scatter!(ax, [n_hi, n_hi-Δ], [c_ref[n_hi+1], c_ref[n_hi-Δ+1]];
         color=:black, markersize=10, label=nothing)

# Annotations showing the error
for (n1, label_y) in [(n_hi, 0.05), (n_hi-Δ, 0.05)]
    p_ref = c_ref[n1+1]
    p_old = c_old[n1+1]
    p_new = c_new[n1+1]
    if p_ref > 0.01
        text!(ax, n1, label_y;
              text=@sprintf("truth=%.3f\nold=%.3f\nnew=%.3f", p_ref, p_old, p_new),
              fontsize=8, color=:black, align=(:center,:bottom))
    end
end

vlines!(ax, [n_hi, n_hi-Δ]; color=(:gray,0.4), linewidth=1, linestyle=:dot)
axislegend(ax; position=:lt, labelsize=10, framevisible=true)

save(joinpath(outdir, "fig_vcycle_conditional.png"), fig2; px_per_unit=3)
println("Saved → fig_vcycle_conditional.png")

# ── Figure 3: N-distribution ───────────────────────────────────────────────
fig3 = Figure(size=(600,320))
ax3 = Axis(fig3[1,1]; title="p_pre(N) vs p_post(N) after coarse CME solve",
    xlabel="N=n₁+n₂", ylabel="probability", titlesize=11)
N_show = findall(>(1e-5), max.(p_pre, p_post)) .- 1
barplot!(ax3, N_show, p_pre[N_show .+ 1];  color=(:steelblue,0.7), label="p_pre")
barplot!(ax3, N_show, p_post[N_show .+ 1]; color=(:red,0.5),       label="p_post")
vlines!(ax3, [N_dom]; color=:black, linewidth=1.5, linestyle=:dash)
text!(ax3, N_dom+0.5, maximum(p_pre)*0.8; text="N=$(N_dom)", fontsize=9)
axislegend(ax3; position=:rt, labelsize=9)
save(joinpath(outdir, "fig_vcycle_Ndist.png"), fig3; px_per_unit=3)
println("Saved → fig_vcycle_Ndist.png")
