"""
fig_vcycle_spatial.jl

Clean V-cycle figure:
  - Two voxel distributions (v1, v2) shown as step-function lines
  - Restrict to total count N (one merged distribution)
  - CME solve: dotted = p_pre, solid = p_post
  - Prolongate back: dotted = original voxel marginals, solid = prolongated
  - Second row shows the L1 intermediate level (coarse split)
"""

using DiscStochSim, CairoMakie, ExponentialUtilities, Printf, LinearAlgebra

using DiscStochSim: StateSpace, add_state!, renormalize!, expand!, build_generator,
                    DiscreteStochasticSystem

_sg  = DiscStochSim._single_voxel_generator
_cl  = DiscStochSim._clean!
_pml = DiscStochSim._prolong_multilevel
_bpf = DiscStochSim._build_pair_fibers

# ── Model ──────────────────────────────────────────────────────────────────
const D = 2.0
model = SchloglModel1D(D, 0.028, 3.2e-4, 21.0, 1.0)
n_lo, n_un, n_hi = schlogl_fixed_points(model)
const N_MAX=220; const NMAX2=2*N_MAX; const KM=30; const STOL=1e-8

pmf(μ,nmax) = (p=zeros(nmax+1); p[1]=exp(-μ);
               for n in 1:nmax; p[n+1]=p[n]*μ/n end; p./=sum(p); p)

# ── Per-voxel distributions at the bistable front ──────────────────────────
# v1: just behind the front (mostly hi, transitioning)
p1 = 0.65.*pmf(n_hi,N_MAX) .+ 0.25.*pmf(n_un+15,N_MAX) .+ 0.10.*pmf(n_lo,N_MAX)
p1 ./= sum(p1)

# v2: at the front (mostly lo, some hi)
p2 = 0.15.*pmf(n_hi,N_MAX) .+ 0.30.*pmf(n_un-5,N_MAX) .+ 0.55.*pmf(n_lo,N_MAX)
p2 ./= sum(p2)

in1=D*Float64(n_hi); in2=D*Float64(n_lo); out1=D; out2=D

# ── Build joint, pre-smooth, restrict, coarse-solve, prolong ───────────────
sp = StateSpace{CartesianIndex{2},Float64}()
for n1 in findall(>(STOL),p1).-1, n2 in findall(>(STOL),p2).-1
    w=p1[n1+1]*p2[n2+1]; w>STOL&&add_state!(sp,CartesianIndex(n1,n2),w)
end
renormalize!(sp)

intra = DiscreteStochasticSystem{CartesianIndex{2}}(
    [CartesianIndex(-1,1),CartesianIndex(1,-1)],
    [(x,r,t)->D*max(0.0,Float64(Tuple(x)[1])),
     (x,r,t)->D*max(0.0,Float64(Tuple(x)[2]))])
bc2 = x->all(c->0≤c≤N_MAX,Tuple(x))
expand!(sp,intra,bc2;depth=1)

A,=build_generator(sp,intra,nothing,0.0)
sp.probs .= expv(0.05,A,sp.probs; m=min(KM,size(A,1)))
_cl(sp.probs,STOL); renormalize!(sp)

p_pre=zeros(NMAX2+1)
for (ci,p) in zip(sp.states,sp.probs)
    N=sum(Tuple(ci)); N≤NMAX2&&(p_pre[N+1]+=p)
end
p_pre./=max(sum(p_pre),1e-100)

dt=0.5
Q_c=_sg(model,in1+in2,(out1+out2)/2,NMAX2;n_sites=2)
p_post=expv(dt,Q_c,p_pre; m=min(KM,size(Q_c,1)))
_cl(p_post,STOL)

sp_new=_pml(sp,p_pre,p_post,in1,in2,out1,out2,model,N_MAX,STOL,200_000)
p1_new=zeros(N_MAX+1); p2_new=zeros(N_MAX+1)
if sp_new!==nothing
    for (ci,p) in zip(sp_new.states,sp_new.probs)
        n1,n2=Tuple(ci); p1_new[n1+1]+=p; p2_new[n2+1]+=p
    end
    s=sum(p1_new); s>0&&(p1_new./=s)
    s=sum(p2_new); s>0&&(p2_new./=s)
end

# L1 fibers at dominant N
_f0,_f1=_bpf(sp,NMAX2,STOL)
N_dom=argmax(p_pre)-1
f1_pre=get(_f1,N_dom,Dict{Int,Float64}())
f1_new=Dict{Int,Float64}()
if sp_new!==nothing
    for (ci,p) in zip(sp_new.states,sp_new.probs)
        n1,n2=Tuple(ci); sum(Tuple(ci))==N_dom||continue
        m=n1÷2; f1_new[m]=get(f1_new,m,0.0)+p
    end
    s=sum(values(f1_new)); s>0&&map!(v->v/s,values(f1_new))
end

@printf("Schlögl: n_lo=%d n_un=%d n_hi=%d  dominant N=%d\n",n_lo,n_un,n_hi,N_dom)

# ── Plotting helpers ────────────────────────────────────────────────────────
# Draw a distribution as a step-function line (cleaner than bars for overlay)
function plot_dist!(ax, p, color; ls=:solid, lw=1.8, tol=5e-4, pad=6, label="")
    lo=max(1, findfirst(>(tol),p)-pad)
    hi=min(length(p), findlast(>(tol),p)+pad)
    lo≤hi || return
    ns=lo-1:hi-1
    stairs!(ax, ns, p[lo:hi]; color=color, linestyle=ls, linewidth=lw,
            step=:post, label=label)
end

# Axis style
function styled_ax(pos, title; kw...)
    Axis(pos; title, titlesize=10, titlefont=:bold,
         xlabel="n", xlabelsize=8, ylabelsize=8,
         xticklabelsize=7, yticklabelsize=7,
         xgridvisible=false, ygridvisible=false,
         spinewidth=0.8, kw...)
end

# Voxel background colour (P(hi) → RdBu-like)
phi1 = sum(p1[n_un+2:end]);  phi2 = sum(p2[n_un+2:end])

# ── Figure ─────────────────────────────────────────────────────────────────
fig = Figure(size=(1400, 480))

# ── Row labels ─────────────────────────────────────────────────────────────
Label(fig[0, 2:3], "L0 · per-voxel CME",   fontsize=10, font=:bold, color=:gray30)
Label(fig[0, 5],   "L2 · total count N",   fontsize=10, font=:bold, color=:gray30)
Label(fig[0, 7:8], "L0 · prolongated",     fontsize=10, font=:bold, color=:steelblue4)
Label(fig[1, 4],   "restrict\n▶▶▶\nΣ anti-diagonals", fontsize=9, color=:gray50,
      tellwidth=false, valign=:center)
Label(fig[1, 6],   "prolong\n◀◀◀\nhierarchical L/R", fontsize=9, color=:steelblue,
      tellwidth=false, valign=:center)

# ── L1 row labels ───────────────────────────────────────────────────────────
Label(fig[2, 2:3], "L1 · coarse split  p(⌊n₁/2⌋ | N=$N_dom)",
      fontsize=9, font=:bold, color=:mediumpurple4)
Label(fig[2, 7:8], "L1 · prolongated split",
      fontsize=9, font=:bold, color=:mediumpurple4)

# ═══════════ LEFT: voxel distributions (L0) ════════════════════════════════
ax_v1 = styled_ax(fig[1,2], "voxel  v₁";
    backgroundcolor=RGBAf(0.12,0.35,0.60,0.07))
plot_dist!(ax_v1, p1, :steelblue; lw=2.0, label="p₁(n)")
vlines!(ax_v1,[n_lo,n_hi]; color=(:black,0.18), linewidth=0.8, linestyle=:dot)
axislegend(ax_v1; position=:rt, labelsize=8, framevisible=false)
ylims!(ax_v1, 0, nothing)

ax_v2 = styled_ax(fig[1,3], "voxel  v₂";
    backgroundcolor=RGBAf(0.60,0.12,0.12,0.07))
plot_dist!(ax_v2, p2, :firebrick; lw=2.0, label="p₂(n)")
vlines!(ax_v2,[n_lo,n_hi]; color=(:black,0.18), linewidth=0.8, linestyle=:dot)
axislegend(ax_v2; position=:rt, labelsize=8, framevisible=false)
ylims!(ax_v2, 0, nothing)

# ═══════════ CENTRE: coarse total-count N (L2) ══════════════════════════════
ax_N = styled_ax(fig[1,5], "coarse variable  p(N)";
    xlabel="N = n₁+n₂",
    backgroundcolor=RGBAf(0.90,0.90,0.90,0.15))

N_show = findall(>(5e-4), max.(p_pre,p_post)) .- 1
stairs!(ax_N, N_show, p_pre[N_show.+1];
    color=:gray50, linestyle=:dash, linewidth=1.8, step=:post, label="p_pre  (before CME)")
stairs!(ax_N, N_show, p_post[N_show.+1];
    color=:darkorange, linestyle=:solid, linewidth=2.2, step=:post, label="p_post (after CME, dt=$dt)")
axislegend(ax_N; position=:rt, labelsize=8, framevisible=true)
ylims!(ax_N, 0, nothing)

# ═══════════ RIGHT: prolongated voxel distributions (L0) ════════════════════
ax_v1p = styled_ax(fig[1,7], "voxel  v₁  (prolongated)";
    backgroundcolor=RGBAf(0.12,0.35,0.60,0.07))
plot_dist!(ax_v1p, p1,     :gray60;  lw=1.4, ls=:dash,  label="original")
plot_dist!(ax_v1p, p1_new, :steelblue; lw=2.2, ls=:solid, label="prolongated")
vlines!(ax_v1p,[n_lo,n_hi]; color=(:black,0.18), linewidth=0.8, linestyle=:dot)
axislegend(ax_v1p; position=:rt, labelsize=8, framevisible=false)
ylims!(ax_v1p, 0, nothing)

ax_v2p = styled_ax(fig[1,8], "voxel  v₂  (prolongated)";
    backgroundcolor=RGBAf(0.60,0.12,0.12,0.07))
plot_dist!(ax_v2p, p2,     :gray60;  lw=1.4, ls=:dash,  label="original")
plot_dist!(ax_v2p, p2_new, :firebrick; lw=2.2, ls=:solid, label="prolongated")
vlines!(ax_v2p,[n_lo,n_hi]; color=(:black,0.18), linewidth=0.8, linestyle=:dot)
axislegend(ax_v2p; position=:rt, labelsize=8, framevisible=false)
ylims!(ax_v2p, 0, nothing)

# ═══════════ L1 ROW: coarse-split conditional at dominant N ═════════════════
# L1 left (from pre-smooth)
ax_L1l = Axis(fig[3, 2:3];
    title="pre-smooth fiber at N=$N_dom", titlesize=9,
    xlabel="bin  m = ⌊n₁/2⌋", xlabelsize=8,
    xticklabelsize=7, yticklabelsize=7,
    backgroundcolor=RGBAf(0.58,0.40,0.74,0.07),
    spinewidth=0.8, xgridvisible=false, ygridvisible=false)
if !isempty(f1_pre)
    ms=sort(collect(keys(f1_pre))); ps=[f1_pre[m] for m in ms]
    stairs!(ax_L1l, ms, ps; color=:mediumpurple, linewidth=2.0, step=:post,
            label="p(⌊n₁/2⌋|N=$N_dom)")
    vlines!(ax_L1l,[n_lo÷2,n_hi÷2]; color=(:black,0.18),linewidth=0.8,linestyle=:dot)
    axislegend(ax_L1l; position=:rt, labelsize=7, framevisible=false)
end
ylims!(ax_L1l, 0, nothing)

# L1 arrows (vertical connections)
Label(fig[3, 4], "L/R stochastic\noperators  ↑", fontsize=8, color=:mediumpurple4,
      tellwidth=false, valign=:center)
Label(fig[3, 6], "↑  within-bin\nrefinement",   fontsize=8, color=:mediumpurple4,
      tellwidth=false, valign=:center)

# L1 right (prolongated)
ax_L1r = Axis(fig[3, 7:8];
    title="prolongated fiber at N=$N_dom", titlesize=9,
    xlabel="bin  m = ⌊n₁/2⌋", xlabelsize=8,
    xticklabelsize=7, yticklabelsize=7,
    backgroundcolor=RGBAf(0.58,0.40,0.74,0.07),
    spinewidth=0.8, xgridvisible=false, ygridvisible=false)
if !isempty(f1_pre)
    ms=sort(collect(keys(f1_pre))); ps_pre=[f1_pre[m] for m in ms]
    stairs!(ax_L1r, ms, ps_pre; color=:gray60, linestyle=:dash, linewidth=1.4, step=:post,
            label="original")
end
if !isempty(f1_new)
    ms_new=sort(collect(keys(f1_new))); ps_new=[f1_new[m] for m in ms_new]
    stairs!(ax_L1r, ms_new, ps_new; color=:mediumpurple, linestyle=:solid, linewidth=2.0,
            step=:post, label="prolongated")
    vlines!(ax_L1r,[n_lo÷2,n_hi÷2]; color=(:black,0.18),linewidth=0.8,linestyle=:dot)
end
axislegend(ax_L1r; position=:rt, labelsize=7, framevisible=false)
ylims!(ax_L1r, 0, nothing)

# L1 centre: show which N this slice corresponds to
ax_L1c = Axis(fig[3, 5];
    title="N=$N_dom  slice", titlesize=9, xlabel="N",
    xticklabelsize=7, yticklabelsize=7,
    spinewidth=0.8, xgridvisible=false, ygridvisible=false)
scatter!(ax_L1c, [N_dom], [p_post[N_dom+1]]; color=:darkorange, markersize=12,
         label="p_post[$N_dom]")
stairs!(ax_L1c, N_show, p_post[N_show.+1];
    color=(:darkorange,0.3), linestyle=:solid, linewidth=1.0, step=:post)
vlines!(ax_L1c,[N_dom]; color=(:red,0.5), linewidth=1.2, linestyle=:dash)
axislegend(ax_L1c; position=:rt, labelsize=7, framevisible=false)
ylims!(ax_L1c, 0, nothing)

# ─── Vertical connectors between row 1 and row 3 ───────────────────────────
Label(fig[2, 4], "▼  restrict\n    bin n₁÷2", fontsize=8, color=:gray50, tellwidth=false)
Label(fig[2, 6], "▲  refine\n    within bin", fontsize=8, color=:steelblue, tellwidth=false)
Label(fig[2, 5], "▼  pick N-slice  ▲", fontsize=8, color=:gray50, tellwidth=false)

# ─── Layout ────────────────────────────────────────────────────────────────
colsize!(fig.layout, 2, Relative(0.17))  # v1
colsize!(fig.layout, 3, Relative(0.17))  # v2
colsize!(fig.layout, 4, Relative(0.06))  # arrow
colsize!(fig.layout, 5, Relative(0.20))  # coarse
colsize!(fig.layout, 6, Relative(0.06))  # arrow
colsize!(fig.layout, 7, Relative(0.17))  # v1'
colsize!(fig.layout, 8, Relative(0.17))  # v2'

rowsize!(fig.layout, 0, Auto())
rowsize!(fig.layout, 1, Relative(0.45))
rowsize!(fig.layout, 2, Auto())
rowsize!(fig.layout, 3, Relative(0.32))

rowgap!(fig.layout, 4); colgap!(fig.layout, 6)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out = joinpath(outdir, "fig_vcycle_spatial.png")
save(out, fig; px_per_unit=3)
println("Saved → $out")
