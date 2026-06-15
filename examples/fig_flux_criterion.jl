#=
fig_flux_criterion.jl

The abstract's adaptive-refinement claim (gap 3):
  "combining total spatial flux with a per-molecule flux measure ... distinguish
   inactive regions from structural bottlenecks where few molecules mediate critical
   transitions."

Naive activation uses TOTAL flux Φ_tot = D·⟨n⟩ (mean-based; this is what
DynCME.activation_candidates currently does). Failure mode: a bistable source voxel
whose probability is mostly in the LOW basin but has a small HIGH-state tail has a
LOW mean — indistinguishable from a genuinely sub-critical (inactive) voxel with the
same mean — yet only the one with the high-state tail can nucleate its neighbour
across the separatrix. The mean averages the critical few molecules away.

Per-molecule (high-state) flux Φ_pm = D·Σ_{n>n_un} n·p(n) measures the flux carried
by the transition-competent molecules and separates the two cases.

GROUND TRUTH: a 2-voxel Schlögl joint solved EXACTLY (expv). Source = voxel 1 held
in a chosen distribution; target = voxel 2 starts pure-low. Nucleation = P(n_2 in
high basin) at time T. We scan two families of source distributions and plot them in
the (Φ_tot, Φ_pm) plane, coloured by whether they actually nucleate the neighbour.

Output: paper/figures/fig_flux_criterion.{png,pdf}
=#
using DiscStochSim, Printf, LinearAlgebra, SparseArrays, CairoMakie
using DiscStochSim: build_graph_rdme, chain_edges, LocalReaction,
                    StateSpace, add_state!, build_generator, build_basin_qsd
using ExponentialUtilities: expv

# ── single-voxel bistable Schlögl (roots ≈ 1 / 8 / 22) ──────────────────────────
c2=0.02; c1=0.2; c3=1.0667; c4=0.72; D=2.0
nmax = 28
schlogl_rxns = [
    LocalReaction{1}(( 1,), nv -> c3),
    LocalReaction{1}((-1,), nv -> c4 * nv[1]),
    LocalReaction{1}(( 1,), nv -> c1 * nv[1] * (nv[1]-1) / 2),
    LocalReaction{1}((-1,), nv -> c2 * nv[1] * (nv[1]-1) * (nv[1]-2) / 6),
]
birth(n) = c3 + c1*n*(n-1)/2
death(n) = c4*n + c2*n*(n-1)*(n-2)/6
Lsv = zeros(nmax+1, nmax+1)
for n in 0:nmax
    b=birth(n); d=death(n); Lsv[n+1,n+1]-=(b+d)
    n<nmax && (Lsv[n+2,n+1]+=b); n>0 && (Lsv[n,n+1]+=d)
end
bq = build_basin_qsd(Lsv, 8; n_max=nmax)
n_un = bq.n_un
qlo = zeros(nmax+1); qlo[(0:n_un).+1]       .= bq.qsd[1]    # low-basin QSD (full vector)
qhi = zeros(nmax+1); qhi[((n_un+1):nmax).+1] .= bq.qsd[2]   # high-basin QSD
@printf("separatrix n_un=%d  ⟨n⟩_lo=%.2f  ⟨n⟩_hi=%.2f\n", n_un, bq.mean[1], bq.mean[2])

# ── 2-voxel joint generator (exact) ─────────────────────────────────────────────
K2 = 2
sys = build_graph_rdme(K2, Val(1), chain_edges(K2), v -> schlogl_rxns, (D,))
bc  = x -> all(0 <= n <= nmax for n in Tuple(x))
sp = StateSpace{CartesianIndex{K2}, Float64}()
for n1 in 0:nmax, n2 in 0:nmax
    add_state!(sp, CartesianIndex(n1, n2), 0.0)
end
A2, = build_generator(sp, sys, nothing, 0.0; absorbing=false)
states = sp.states                                   # ordered CartesianIndex list
idx_of = Dict(s => i for (i,s) in enumerate(states))

# helper: true downstream nucleation prob given a source (voxel-1) distribution
T_nuc = 3.0
function nucleation(p1::Vector{Float64})
    p0 = zeros(length(states))
    for (i,s) in enumerate(states)
        p0[i] = p1[s[1]+1] * qlo[s[2]+1]              # source ⊗ target(pure low)
    end
    p0 ./= sum(p0)
    pT = expv(T_nuc, A2, p0)
    # P(voxel 2 in high basin)
    s_nuc = 0.0
    for (i,s) in enumerate(states); s[2] > n_un && (s_nuc += pT[i]); end
    s_nuc
end
flux_tot(p1) = D * sum((0:nmax) .* p1)
flux_pm(p1)  = D * sum(((n_un+1):nmax) .* p1[((n_un+1):nmax).+1])

# baseline spontaneous nucleation (source pure low) — subtract it off
base_nuc = nucleation(qlo)
@printf("baseline (pure-low source) nucleation at T=%.1f: %.3e\n", T_nuc, base_nuc)

# ── Family A: SUB-CRITICAL unimodal sources (peak below separatrix) ──────────────
#   many molecules, all below threshold → high Φ_tot possible, Φ_pm = 0, no nucleation
peaksA = 1.0:1.0:Float64(n_un)
famA = NamedTuple[]
for nc in peaksA
    p = exp.(-((0:nmax) .- nc).^2 ./ (2*1.2^2)); p ./= sum(p)
    push!(famA, (ftot=flux_tot(p), fpm=flux_pm(p), nuc=nucleation(p)-base_nuc, p=p, nc=nc))
end

# ── Family B: bistable MIXTURES (low basin + small high tail) ────────────────────
#   low mean but a few high-state molecules → low Φ_tot, Φ_pm > 0, DOES nucleate
αs = 0.04:0.04:0.6
famB = NamedTuple[]
for α in αs
    p = (1-α).*qlo .+ α.*qhi; p ./= sum(p)
    push!(famB, (ftot=flux_tot(p), fpm=flux_pm(p), nuc=nucleation(p)-base_nuc, α=α, p=p))
end

nuc_thresh = 0.05      # "nucleates" if it lifts neighbour high-basin prob > 5%
@printf("\nFamily A (sub-critical): Φ_tot∈[%.1f,%.1f]  Φ_pm≈0  max nuc=%.3f\n",
        minimum(a.ftot for a in famA), maximum(a.ftot for a in famA),
        maximum(a.nuc for a in famA))
@printf("Family B (bottleneck):  Φ_tot∈[%.1f,%.1f]  Φ_pm∈[%.1f,%.1f]  max nuc=%.3f\n",
        minimum(b.ftot for b in famB), maximum(b.ftot for b in famB),
        minimum(b.fpm for b in famB), maximum(b.fpm for b in famB),
        maximum(b.nuc for b in famB))

# matched-mean pair for panel (a): take the GENUINELY sub-critical voxel (high tail
# negligible, Φ_pm≈0) with the most molecules, then the Family-B bottleneck with the
# same total flux Φ_tot — the cleanest "same Φ_tot, opposite fate" contrast.
subcrit = [a for a in famA if a.fpm < 0.5]
mpk  = subcrit[argmax([a.ftot for a in subcrit])]
mref = famB[argmin([abs(b.ftot - mpk.ftot) for b in famB])]
@printf("matched pair: sub-critical(n=%.1f) Φtot=%.1f Φpm=%.1f nuc=%.3f | mixture(α=%.2f) Φtot=%.1f Φpm=%.1f nuc=%.3f\n",
        mpk.nc, mpk.ftot, mpk.fpm, mpk.nuc, mref.α, mref.ftot, mref.fpm, mref.nuc)

# ── Figure ──────────────────────────────────────────────────────────────────────
set_theme!(theme_minimal())
fig = Figure(size=(1080, 420), fontsize=15)

# (a) the matched-mean pair: same Φ_tot, different high tail → different fate
ax1 = Axis(fig[1,1], title="(a) two sources, same total flux Φ_tot",
           xlabel="count n", ylabel="P(n)")
vlines!(ax1, [n_un+0.5], color=:gray, linestyle=:dash, label="separatrix")
barplot!(ax1, 0:nmax, mpk.p, color=(:steelblue,0.6),
         label=@sprintf("sub-critical: nuc=%.2f", mpk.nuc))
barplot!(ax1, 0:nmax, mref.p, color=(:crimson,0.55),
         label=@sprintf("bottleneck: nuc=%.2f", mref.nuc))
axislegend(ax1, position=:rt, framevisible=false, labelsize=10)

# (b) decision plane: (Φ_tot, Φ_pm) coloured by TRUE nucleation
ax2 = Axis(fig[1,2], title="(b) activation decision plane (colour = true nucleation)",
           xlabel="total flux  Φ_tot = D·⟨n⟩",
           ylabel="per-molecule flux  Φ_pm = D·Σ_{n>n_un} n·p(n)")
allf = vcat(famA, famB)
sc = scatter!(ax2, [f.ftot for f in allf], [f.fpm for f in allf],
              color=[f.nuc for f in allf], colormap=:viridis,
              colorrange=(0, maximum(f.nuc for f in allf)), markersize=13,
              strokewidth=0.5, strokecolor=:black)
Colorbar(fig[1,3], sc, label="downstream nucleation P(n₂>n_un)")
# naive total-flux threshold: any vertical line misclassifies
ftot_thr = mref.ftot
vlines!(ax2, [ftot_thr], color=:red, linewidth=2, linestyle=:dash)
text!(ax2, ftot_thr, maximum(f.fpm for f in allf)*0.95,
      text="naive Φ_tot\nthreshold", color=:red, align=(:left,:top), fontsize=10)
# combined criterion: also activate if Φ_pm exceeds a small bound
fpm_thr = 0.5
hlines!(ax2, [fpm_thr], color=:black, linewidth=2)
text!(ax2, maximum(f.ftot for f in allf)*0.98, fpm_thr,
      text="per-molecule cut", color=:black, align=(:right,:bottom), fontsize=10)

save(joinpath(@__DIR__, "..", "paper", "figures", "fig_flux_criterion.png"), fig)
save(joinpath(@__DIR__, "..", "paper", "figures", "fig_flux_criterion.pdf"), fig)
println("saved paper/figures/fig_flux_criterion.{png,pdf}")
