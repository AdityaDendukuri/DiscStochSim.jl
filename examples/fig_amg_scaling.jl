"""
fig_amg_scaling.jl

The "multigrid actually scales" figure for the talk.

Data are the VERIFIED output of /tmp/ring_probe.jl (plain-aggregation Galerkin
AMG, V-cycle-preconditioned restarted GMRES) on the bistable Schlögl JOINT CME
over a K-voxel ring, NM=50 (51 states/voxel), D=0.15, M = I − 5A.

  K   |S|          levels  op-cplx  V-cycle ρ   GMRES iters (relres ~6e-9)
  2   2,601        2       1.25     0.7306      23
  3   132,651      5       1.39     0.9152      31
  4   6,765,201    8       1.48     0.8963      32

State space grows ×2600; GMRES iterations stay ~flat (23→32) — the canonical
multigrid size-independence. The standalone V-cycle stalls (ρ≈0.9, smoother-
limited), which is WHY it is wrapped as a GMRES preconditioner. GMRES drives the
true residual to ~1e-9, so accuracy is guaranteed regardless of coarse quality.

Output: paper/figures/fig_amg_scaling.{png,pdf}
"""

using CairoMakie

# ── verified data ─────────────────────────────────────────────────────────────
S     = [2601.0, 132651.0, 6765201.0]
Kvox  = [2, 3, 4]
iters = [23.0, 31.0, 32.0]
rho   = [0.7306, 0.9152, 0.8963]
opcx  = [1.25, 1.39, 1.48]

set_theme!(Theme(fontsize = 13))
fig = Figure(size = (920, 360))

# ── Panel (a): iterations vs problem size — the headline ──────────────────────
axa = Axis(fig[1, 1];
    title = "Iterations are size-independent",
    titlesize = 13,
    xlabel = "joint state-space size |S|", ylabel = "GMRES iterations",
    xscale = log10, xlabelsize = 12, ylabelsize = 12)
xlims!(axa, 1.5e3, 1.5e7)
ylims!(axa, 0, 45)
scatterlines!(axa, S, iters; color = :steelblue, linewidth = 2.5, markersize = 13,
    marker = :circle)
for i in eachindex(S)
    text!(axa, S[i], iters[i] + 3.2; text = "K=$(Kvox[i])", align = (:center, :bottom),
        fontsize = 11, color = :gray25)
end
# annotate the span
text!(axa, 3.0e3, 7;
    text = "|S| ×2600  ⇒  iters 23→32", align = (:left, :bottom),
    fontsize = 12, color = :steelblue, font = :bold)

# ── Panel (b): standalone V-cycle stalls → motivates the Krylov wrap ──────────
axb = Axis(fig[1, 2];
    title = "Standalone V-cycle stalls (smoother-limited)",
    titlesize = 13,
    xlabel = "joint state-space size |S|", ylabel = "convergence factor ρ / cycle",
    xscale = log10, xlabelsize = 12, ylabelsize = 12)
xlims!(axb, 1.5e3, 1.5e7)
ylims!(axb, 0, 1.08)
hlines!(axb, [1.0]; color = (:crimson, 0.6), linestyle = :dash, linewidth = 1.5)
text!(axb, 2.0e3, 1.005; text = "ρ = 1 (no progress)", align = (:left, :bottom),
    fontsize = 10, color = :crimson)
scatterlines!(axb, S, rho; color = :crimson, linewidth = 2.5, markersize = 13,
    marker = :rect)
text!(axb, 2.0e3, 0.10;
    text = "→ used as a GMRES preconditioner;\n   true residual driven to ~1e-9",
    align = (:left, :bottom), fontsize = 11, color = :gray25)

colgap!(fig.layout, 26)

for ext in ("png", "pdf")
    out = joinpath(@__DIR__, "..", "paper", "figures", "fig_amg_scaling.$ext")
    save(out, fig; px_per_unit = 2)
    println("wrote $out")
end
