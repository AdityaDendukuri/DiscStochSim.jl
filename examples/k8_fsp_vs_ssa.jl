"""
FSP vs SSA comparison for the K=4 Schlögl RDME coarse-fine interface.

Minimal coarse-fine setting: K=4 voxels with h=1/8 (same voxel spacing as K=8),
mask=[1,0] — pair 1 coarse (nc₁), pair 2 fine (n₃,n₄).

The 3D mixed state (nc₁,n₃,n₄) is tractable for a numerically accurate FSP
reference (mixed_expand_depth=10, prune_tol=1e-8).  The SSA simulates the same
3D mixed CME exactly, so TV(FSP, SSA) → 0 as 1/√N.

Figure panels:
  (a) FSP exact P(n₃,n₄) at t=3   (b) SSA histogram N=10³   (c) TV vs N

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/k8_fsp_vs_ssa.jl
"""

using DiscStochSim
using CairoMakie
using Printf
using Random
using Statistics

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
const SMALL = 8f0; const MED = 9f0; const LW = 1.4f0
const C_FSP  = RGBf(0.000, 0.447, 0.698)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

# ── Model: same parameters as the K=8 study ───────────────────────────────────
D = 0.005; K = 4; h = 1.0/8    # h=1/8 → same d = D/h² = 0.32 as K=8
model = SchloglModel1D(D)
g     = VoxelGrid(K, h, 0)

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_buf = 30
@printf("n_low=%d  n_uns=%d  n_high=%d  d=%.4f\n",
        n_low, n_uns, n_high, D/h^2); flush(stdout)

# ── Mixed system: pair 1 coarse, pair 2 fine ─────────────────────────────────
mask = BitVector([true, false])     # [coarse, fine]
off  = DiscStochSim._mixed_offsets(mask)
op2  = off[2]                       # index of n₃ in the 3-tuple; n₄ at op2+1
KM   = K - sum(mask)               # = 3
@printf("mask=%s  KM=%d  op2=%d\n", string(mask), KM, op2); flush(stdout)

u0 = CartesianIndex(2*n_low, n_low, n_high)    # (62, 31, 162)
@printf("u0 = %s\n", string(u0)); flush(stdout)

sp    = StateSpace{CartesianIndex{KM}, Float64}()
add_state!(sp, u0, 1.0)
mxsys = build_schlogl_mixed_system(model, g, mask)

n_max_c = 2*(n_high + n_buf)    # coarse bound: 384
function bc_init(s); all(c -> 0 ≤ c ≤ n_max_c, Tuple(s)); end
expand!(sp, mxsys, bc_init; depth=1)
@printf("Initial |S|=%d\n", length(sp.states)); flush(stdout)

# ── Mixed FSP: 30 steps × dt=0.1 → t=3.0 ────────────────────────────────────
dt      = 0.1
n_steps = 30
t_comp  = n_steps * dt

cache = MixedSolverCache(); sc = Ref(0)

@printf("%-6s  %-8s  %-8s\n", "step", "|S|", "Σp"); flush(stdout)
for step in 1:n_steps
    global sp
    sp, _ = two_level_step_mixed(sp, mask, model, g, dt;
                                  cache               = cache,
                                  step_count          = sc,
                                  mask_check_interval = 100_000,
                                  expand_mixed        = true,
                                  mixed_expand_depth  = 10,
                                  mixed_n_max         = n_max_c,
                                  prune_tol           = 1e-8,
                                  reexpand_depth      = 1,
                                  krylov_m            = 60)
    sp.probs = max.(sp.probs, 0.0); sp.probs ./= sum(sp.probs)
    step % 5 == 0 && @printf("%-6d  %-8d  %-8.5f\n",
                              step, length(sp.states), sum(sp.probs)); flush(stdout)
end

# ── Marginal P(n₃,n₄) ────────────────────────────────────────────────────────
mg_fsp = Dict{Tuple{Int,Int}, Float64}()
for (i,s) in enumerate(sp.states)
    t_tup = Tuple(s); p = sp.probs[i]
    key = (t_tup[op2], t_tup[op2+1])
    mg_fsp[key] = get(mg_fsp, key, 0.0) + p
end
@printf("FSP P(n₃,n₄) at t=%.1f: %d unique states\n", t_comp, length(mg_fsp))
flush(stdout)

# ── SSA for the K=4 mixed CME ─────────────────────────────────────────────────
# State [nc₁, n₃, n₄] — 16 reactions matching build_schlogl_mixed_system(mask=[1,0])
function k4_mixed_ssa!(s, c1, c2, c3, c4, d, t_final, rng)
    rates = Vector{Float64}(undef, 16)

    function fill_rates!(r)
        nc1, n3, n4 = s[1], s[2], s[3]
        r[1]  = 2*c3;           r[2]  = c4*nc1
        r[3]  = c1*max(0, nc1*(nc1-1)) / 4
        r[4]  = c2*max(0, nc1*(nc1-1)*(nc1-2)) / 24
        r[5]  = c3;             r[6]  = c4*n3
        r[7]  = c1*max(0, n3*(n3-1)) / 2
        r[8]  = c2*max(0, n3*(n3-1)*(n3-2)) / 6
        r[9]  = c3;             r[10] = c4*n4
        r[11] = c1*max(0, n4*(n4-1)) / 2
        r[12] = c2*max(0, n4*(n4-1)*(n4-2)) / 6
        r[13] = d*n3;           r[14] = d*n4   # intra-pair2: v3↔v4
        r[15] = d*nc1/2;        r[16] = d*n3   # inter-pair:  pair1↔v3
    end

    function apply!(idx)
        if     idx==1;  s[1]+=1  elseif idx==2;  s[1]-=1
        elseif idx==3;  s[1]+=1  elseif idx==4;  s[1]-=1
        elseif idx==5;  s[2]+=1  elseif idx==6;  s[2]-=1
        elseif idx==7;  s[2]+=1  elseif idx==8;  s[2]-=1
        elseif idx==9;  s[3]+=1  elseif idx==10; s[3]-=1
        elseif idx==11; s[3]+=1  elseif idx==12; s[3]-=1
        elseif idx==13; s[2]-=1; s[3]+=1   # v3→v4
        elseif idx==14; s[3]-=1; s[2]+=1   # v4→v3
        elseif idx==15; s[1]-=1; s[2]+=1   # pair1→v3
        elseif idx==16; s[2]-=1; s[1]+=1   # v3→pair1
        end
    end

    t = 0.0
    fill_rates!(rates)
    while t < t_final
        a0 = sum(rates)
        a0 ≤ 0.0 && break
        t  += -log(rand(rng)) / a0
        t  > t_final && break
        rv = rand(rng)*a0; cs = 0.0; idx = 16
        for i in eachindex(rates)
            cs += rates[i]; if cs ≥ rv; idx = i; break; end
        end
        apply!(idx)
        fill_rates!(rates)
    end
end

c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4
d               = D / h^2
u0_ssa          = [2*n_low, n_low, n_high]
rng             = Random.default_rng()

N_max  = 20_000
N_vals = [50, 100, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000]

@printf("Running SSA N=%d ...", N_max); flush(stdout)
t0 = time()
all_samps = Vector{Tuple{Int,Int}}(undef, N_max)
for i in 1:N_max
    s = copy(u0_ssa)
    k4_mixed_ssa!(s, c1, c2, c3, c4, d, t_comp, rng)
    all_samps[i] = (s[2], s[3])    # (n₃, n₄)
end
@printf(" done (%.1fs)\n", time()-t0); flush(stdout)

# TV distance
function tv_distance(mg_fsp, samps)
    N  = length(samps)
    mg = Dict{Tuple{Int,Int},Float64}()
    for s in samps; mg[s] = get(mg, s, 0.0) + 1.0/N; end
    all_keys = union(keys(mg_fsp), keys(mg))
    0.5 * sum(abs(get(mg_fsp, k, 0.0) - get(mg, k, 0.0)) for k in all_keys)
end

n_boot = 20
tv_mean = Float64[]; tv_lo = Float64[]; tv_hi = Float64[]
for N in N_vals
    tvs = [tv_distance(mg_fsp, all_samps[randperm(rng, N_max)[1:N]]) for _ in 1:n_boot]
    push!(tv_mean, mean(tvs)); push!(tv_lo, quantile(tvs,0.1)); push!(tv_hi, quantile(tvs,0.9))
end

fsp_n3 = sum(n3*p for ((n3,_),p) in mg_fsp)
fsp_n4 = sum(n4*p for ((_,n4),p) in mg_fsp)
ssa_n3 = mean(s[1] for s in all_samps)
ssa_n4 = mean(s[2] for s in all_samps)
@printf("FSP mean (n₃,n₄) = (%.1f, %.1f)  range n₃=[%d,%d]  n₄=[%d,%d]\n",
    fsp_n3, fsp_n4,
    minimum(k[1] for k in keys(mg_fsp)), maximum(k[1] for k in keys(mg_fsp)),
    minimum(k[2] for k in keys(mg_fsp)), maximum(k[2] for k in keys(mg_fsp)))
@printf("SSA mean (n₃,n₄) = (%.1f, %.1f)  range n₃=[%d,%d]  n₄=[%d,%d]\n",
    ssa_n3, ssa_n4,
    minimum(s[1] for s in all_samps), maximum(s[1] for s in all_samps),
    minimum(s[2] for s in all_samps), maximum(s[2] for s in all_samps))
@printf("TV at N=%-6d: %.3f\n", N_vals[1], tv_mean[1])
@printf("TV at N=%-6d: %.3f\n", N_vals[end], tv_mean[end]); flush(stdout)

# ── Heatmap helper ────────────────────────────────────────────────────────────
function build_heatmap(mg, bin_w, n_max_hm)
    n_bins = div(n_max_hm, bin_w) + 1
    hm = zeros(n_bins, n_bins)
    for ((n3,n4), p) in mg
        i3 = clamp(div(n3, bin_w)+1, 1, n_bins)
        i4 = clamp(div(n4, bin_w)+1, 1, n_bins)
        hm[i3, i4] += p
    end
    Float64.(0:bin_w:bin_w*(n_bins-1)), hm
end

bin_w    = 5
n_max_hm = n_high + 30
bin_ctrs, hm_fsp_hm = build_heatmap(mg_fsp, bin_w, n_max_hm)

mg_ssa_1k = Dict{Tuple{Int,Int},Float64}()
for s in all_samps[1:1_000]; mg_ssa_1k[s] = get(mg_ssa_1k, s, 0.0) + 1e-3; end
_, hm_ssa_1k = build_heatmap(mg_ssa_1k, bin_w, n_max_hm)

hm_max = max(maximum(hm_fsp_hm) * 0.95, 1e-10)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(800, 320))

ax_a = Axis(fig[1,1];
    xlabel="n₃", ylabel="n₄",
    title="(a)  FSP exact P(n₃,n₄)  [K=4, t=$(Int(t_comp))]",
    titlesize=MED, aspect=1)
heatmap!(ax_a, bin_ctrs, bin_ctrs, hm_fsp_hm;
         colormap=:hot, colorrange=(0, hm_max))
scatter!(ax_a, [Float64(n_low)],  [Float64(n_high)];
         color=:cyan, marker=:cross, markersize=10, strokewidth=0)
scatter!(ax_a, [Float64(n_high)], [Float64(n_low)];
         color=:cyan, marker=:cross, markersize=10, strokewidth=0)
text!(ax_a, 5.0, Float64(n_max_hm)-5;
      text="$(length(mg_fsp)) support pts",
      fontsize=SMALL-1, color=:white, align=(:left,:top))

ax_b = Axis(fig[1,2];
    xlabel="n₃", ylabel="",
    title="(b)  SSA histogram  N=10³  t=$(Int(t_comp))",
    titlesize=MED, aspect=1)
heatmap!(ax_b, bin_ctrs, bin_ctrs, hm_ssa_1k;
         colormap=:hot, colorrange=(0, hm_max))
scatter!(ax_b, [Float64(n_low)],  [Float64(n_high)];
         color=:cyan, marker=:cross, markersize=10, strokewidth=0)
scatter!(ax_b, [Float64(n_high)], [Float64(n_low)];
         color=:cyan, marker=:cross, markersize=10, strokewidth=0)

ax_c = Axis(fig[1,3];
    xlabel="N  (SSA trajectories)",
    ylabel="TV distance to FSP",
    title="(c)  SSA convergence to FSP distribution",
    titlesize=MED,
    xscale=log10, yscale=log10)

ref_x   = [Float64(N_vals[1]), Float64(N_vals[end])]
ref_idx = findfirst(==(1_000), N_vals)
ref_y   = tv_mean[ref_idx] .* sqrt.(Float64(N_vals[ref_idx]) ./ ref_x)
lines!(ax_c, ref_x, ref_y;
       color=:gray, linewidth=0.9, linestyle=:dash, label="∝ N^{-1/2}")
band!(ax_c, Float64.(N_vals), tv_lo, tv_hi; color=(C_FSP, 0.25))
scatterlines!(ax_c, Float64.(N_vals), tv_mean;
              color=C_FSP, linewidth=LW, markersize=6, label="TV(SSA, FSP)")
axislegend(ax_c; position=:lb, framevisible=false, labelsize=SMALL-1)

colgap!(fig.layout, 10)
CairoMakie.save(joinpath(OUTDIR, "fig_fsp_vs_ssa.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_fsp_vs_ssa.png"), fig; px_per_unit=2)
println("Saved: paper/figures/fig_fsp_vs_ssa.{pdf,png}")
