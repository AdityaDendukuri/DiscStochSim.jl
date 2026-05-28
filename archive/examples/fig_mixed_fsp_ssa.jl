"""
Mixed-FSP vs SSA: K=4 Schlögl with pair 1 coarse (nc₁), pair 2 fine (n₃, n₄).

Uses direct FSP on the 3D mixed CME and compares the interface marginal P(n₃,n₄,t)
against SSA of the identical mixed generator. Independent SSA batches are used for
each ensemble size N, and stepwise FSP diagnostics report boundary/pruning losses.

Physical state at t=0: [31, 31, 31, 162] (front between voxels 3 and 4)
Mixed state: (nc₁, n₃, n₄) = (62, 31, 162)

Run: JULIA_DEPOT_PATH=/tmp/ds_depot julia --project examples/fig_mixed_fsp_ssa.jl
"""

using DiscStochSim
using ExponentialUtilities
using CairoMakie
using Printf, Random, Statistics

const OUTDIR = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(OUTDIR)

const SMALL = 8f0
const MED   = 9f0
const LW    = 1.4f0
const C_FSP = RGBf(0.000, 0.447, 0.698)

set_theme!(Theme(fontsize=MED,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3f0, tickwidth=0.6f0)))

function boundary_mass(sp, upper_bound; halo=2)
    sum((sp.probs[i] for i in eachindex(sp.states)
         if any(c >= upper_bound - halo for c in Tuple(sp.states[i]))); init=0.0)
end

function marginal_n3n4(sp, op2)
    mg = Dict{Tuple{Int, Int}, Float64}()
    for (i, s) in enumerate(sp.states)
        t = Tuple(s)
        key = (t[op2], t[op2 + 1])
        mg[key] = get(mg, key, 0.0) + sp.probs[i]
    end
    mg
end

function tv_distance(mg_fsp, mg_ssa)
    all_keys = union(keys(mg_fsp), keys(mg_ssa))
    0.5 * sum(abs(get(mg_fsp, k, 0.0) - get(mg_ssa, k, 0.0)) for k in all_keys)
end

function samples_to_marginal(samps)
    N = length(samps)
    mg = Dict{Tuple{Int, Int}, Float64}()
    invN = 1.0 / N
    for s in samps
        mg[s] = get(mg, s, 0.0) + invN
    end
    mg
end

function build_heatmap(mg, bin_w, n_max_hm)
    n_bins = div(n_max_hm, bin_w) + 1
    hm = zeros(n_bins, n_bins)
    for ((n3, n4), p) in mg
        i3 = clamp(div(n3, bin_w) + 1, 1, n_bins)
        i4 = clamp(div(n4, bin_w) + 1, 1, n_bins)
        hm[i3, i4] += p
    end
    Float64.(0:bin_w:bin_w * (n_bins - 1)), hm
end

# ── Model ─────────────────────────────────────────────────────────────────────
D = 0.005
K = 4
h = 1.0 / 8   # use K=8 voxel spacing so the local interface dynamics match the K=8 front
model = SchloglModel1D(D)
g = VoxelGrid(K, h, 0)
rates = Float64[]

fp = schlogl_fixed_points(model; n_max=300)
n_low, n_uns, n_high = fp
n_buf = 40
d = D / h^2
@printf("n_low=%d  n_uns=%d  n_high=%d  d=%.4f\n", n_low, n_uns, n_high, d)

# ── Mixed system: pair 1 coarse, pair 2 fine ──────────────────────────────────
mask = BitVector([true, false])
KM = K - sum(mask)   # = 3 (nc₁, n₃, n₄)
u0 = CartesianIndex(2 * n_low, n_low, n_high)   # (62, 31, 162)
@printf("u0 = (%d, %d, %d)\n", Tuple(u0)...)

mxsys = build_schlogl_mixed_system(model, g, mask)
upper_bound = 2 * (n_high + n_buf)
bc = s -> all(c -> 0 <= c <= upper_bound, Tuple(s))
off = DiscStochSim._mixed_offsets(mask)
op2 = off[2]   # index of n₃ in the mixed tuple

# ── Direct FSP for the 3D mixed CME ───────────────────────────────────────────
sp = StateSpace{CartesianIndex{KM}, Float64}()
add_state!(sp, u0, 1.0)

dt = 0.1
n_steps = 5
t_comp = n_steps * dt
expand_depth = 16
prune_tol = 1e-9

println("\nRunning direct FSP for mixed CME..."); flush(stdout)
@printf("%-6s  %-8s  %-10s  %-10s  %-10s\n",
        "step", "|S|", "mass pre", "pruned", "bdry mass"); flush(stdout)

max_pruned = 0.0
max_bdry = 0.0
for step in 1:n_steps
    expand!(sp, mxsys, bc; depth=expand_depth)
    A, = build_generator(sp, mxsys, rates, 0.0)
    sp.probs .= max.(expv(dt, A, sp.probs; m=60), 0.0)

    mass_pre = sum(sp.probs)
    bdry = boundary_mass(sp, upper_bound)

    pruned = 0.0
    for i in eachindex(sp.probs)
        sp.probs[i] < prune_tol && (pruned += sp.probs[i])
    end

    renormalize!(sp)
    prune_threshold!(sp, prune_tol)

    global max_pruned = max(max_pruned, pruned)
    global max_bdry = max(max_bdry, bdry)

    step % 5 == 0 && @printf("%-6d  %-8d  %-10.6f  %-10.2e  %-10.2e\n",
                              step, length(sp), mass_pre, pruned, bdry)
    flush(stdout)
end

@printf("\nFSP diagnostics: max pruned mass = %.2e, max boundary mass = %.2e\n",
        max_pruned, max_bdry)

mg_fsp = marginal_n3n4(sp, op2)
@printf("FSP P(n₃,n₄) at t=%.1f: %d support points\n", t_comp, length(mg_fsp))
fsp_n3 = sum(n3 * p for ((n3, _), p) in mg_fsp)
fsp_n4 = sum(n4 * p for ((_, n4), p) in mg_fsp)
fsp_cov = sum((n3 - fsp_n3) * (n4 - fsp_n4) * p for ((n3, n4), p) in mg_fsp)
@printf("FSP means/cov: ⟨n₃⟩=%.2f  ⟨n₄⟩=%.2f  Cov=%.2f\n", fsp_n3, fsp_n4, fsp_cov)

# ── SSA for the identical mixed CME ──────────────────────────────────────────
c1, c2, c3, c4 = model.c1, model.c2, model.c3, model.c4

function mixed_ssa!(s, t_final, rng, d, c1, c2, c3, c4)
    # State s = [nc₁, n₃, n₄]; rates match build_schlogl_mixed_system(mask=[true,false])
    rates_ssa = Vector{Float64}(undef, 16)

    function fill_rates!()
        nc1, n3, n4 = s[1], s[2], s[3]
        rates_ssa[1]  = 2 * c3
        rates_ssa[2]  = c4 * nc1
        rates_ssa[3]  = c1 * max(0, nc1 * (nc1 - 1)) / 4
        rates_ssa[4]  = c2 * max(0, nc1 * (nc1 - 1) * (nc1 - 2)) / 24

        rates_ssa[5]  = c3
        rates_ssa[6]  = c4 * n3
        rates_ssa[7]  = c1 * max(0, n3 * (n3 - 1)) / 2
        rates_ssa[8]  = c2 * max(0, n3 * (n3 - 1) * (n3 - 2)) / 6

        rates_ssa[9]  = c3
        rates_ssa[10] = c4 * n4
        rates_ssa[11] = c1 * max(0, n4 * (n4 - 1)) / 2
        rates_ssa[12] = c2 * max(0, n4 * (n4 - 1) * (n4 - 2)) / 6

        rates_ssa[13] = d * n3
        rates_ssa[14] = d * n4
        rates_ssa[15] = d * nc1 / 2
        rates_ssa[16] = d * n3
    end

    function apply!(idx)
        if idx == 1 || idx == 3
            s[1] += 1
        elseif idx == 2 || idx == 4
            s[1] -= 1
        elseif idx == 5 || idx == 7
            s[2] += 1
        elseif idx == 6 || idx == 8
            s[2] -= 1
        elseif idx == 9 || idx == 11
            s[3] += 1
        elseif idx == 10 || idx == 12
            s[3] -= 1
        elseif idx == 13
            s[2] -= 1; s[3] += 1
        elseif idx == 14
            s[3] -= 1; s[2] += 1
        elseif idx == 15
            s[1] -= 1; s[2] += 1
        elseif idx == 16
            s[2] -= 1; s[1] += 1
        end
    end

    t = 0.0
    fill_rates!()
    while t < t_final
        a0 = sum(rates_ssa)
        a0 <= 0.0 && break
        t += -log(rand(rng)) / a0
        t > t_final && break
        rv = rand(rng) * a0
        cs = 0.0
        idx = 16
        for i in eachindex(rates_ssa)
            cs += rates_ssa[i]
            if cs >= rv
                idx = i
                break
            end
        end
        apply!(idx)
        fill_rates!()
    end
    nothing
end

function run_ssa_ensemble(N, seed, u0_ssa, t_final, d, c1, c2, c3, c4)
    rng = MersenneTwister(seed)
    samps = Vector{Tuple{Int, Int}}(undef, N)
    for i in 1:N
        s = copy(u0_ssa)
        mixed_ssa!(s, t_final, rng, d, c1, c2, c3, c4)
        samps[i] = (s[2], s[3])
    end
    samps
end

u0_ssa = [2 * n_low, n_low, n_high]
N_vals = [100, 200, 500, 1_000, 2_000, 5_000, 10_000]
n_boot = 8

println("\nRunning independent SSA batches..."); flush(stdout)
tv_mean = Float64[]
tv_lo = Float64[]
tv_hi = Float64[]
mg_ssa_show = nothing
show_N = 1_000

for N in N_vals
    tvs = Float64[]
    for rep in 1:n_boot
        samps = run_ssa_ensemble(N, 10_000 * N + rep, u0_ssa, t_comp, d, c1, c2, c3, c4)
        mg_ssa = samples_to_marginal(samps)
        push!(tvs, tv_distance(mg_fsp, mg_ssa))
        if N == show_N && rep == 1
            global mg_ssa_show = mg_ssa
            ssa_n3 = mean(s[1] for s in samps)
            ssa_n4 = mean(s[2] for s in samps)
            ssa_cov = mean((s[1] - ssa_n3) * (s[2] - ssa_n4) for s in samps)
            @printf("SSA N=%d means/cov: ⟨n₃⟩=%.2f  ⟨n₄⟩=%.2f  Cov=%.2f\n",
                    N, ssa_n3, ssa_n4, ssa_cov)
        end
    end
    push!(tv_mean, mean(tvs))
    push!(tv_lo, quantile(tvs, 0.1))
    push!(tv_hi, quantile(tvs, 0.9))
    @printf("  N=%-6d  TV=%.4f  [%.4f, %.4f]\n", N, tv_mean[end], tv_lo[end], tv_hi[end])
    flush(stdout)
end

bin_w = 5
n_max_hm = n_high + 30
bin_ctrs, hm_fsp = build_heatmap(mg_fsp, bin_w, n_max_hm)
_, hm_ssa = build_heatmap(mg_ssa_show, bin_w, n_max_hm)
hm_max = max(maximum(hm_fsp) * 0.95, 1e-10)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(800, 280))

ax_a = Axis(fig[1, 1];
    xlabel="n₃", ylabel="n₄",
    title="(a)  Mixed FSP  P(n₃,n₄)  [t=$(round(t_comp, digits=1))]",
    titlesize=MED, aspect=1)
heatmap!(ax_a, bin_ctrs, bin_ctrs, hm_fsp; colormap=:hot, colorrange=(0, hm_max))
text!(ax_a, 5.0, Float64(n_max_hm) - 5;
      text="$(length(mg_fsp)) support pts",
      fontsize=SMALL - 1, color=:white, align=(:left, :top))

ax_b = Axis(fig[1, 2];
    xlabel="n₃", ylabel="",
    title="(b)  Mixed SSA histogram  N=$(show_N)",
    titlesize=MED, aspect=1)
heatmap!(ax_b, bin_ctrs, bin_ctrs, hm_ssa; colormap=:hot, colorrange=(0, hm_max))

ax_c = Axis(fig[1, 3];
    xlabel="N  (SSA trajectories)",
    ylabel="TV(FSP, SSA)",
    title="(c)  Mixed SSA convergence to mixed FSP",
    titlesize=MED, xscale=log10, yscale=log10)

ref_idx = findfirst(==(show_N), N_vals)
ref_y = tv_mean[ref_idx] .* sqrt.(Float64(N_vals[ref_idx]) ./ Float64.(N_vals))
lines!(ax_c, Float64.(N_vals), ref_y;
       color=:gray, linewidth=0.9, linestyle=:dash, label="∝ N^{-1/2}")
band!(ax_c, Float64.(N_vals), tv_lo, tv_hi; color=(C_FSP, 0.25))
scatterlines!(ax_c, Float64.(N_vals), tv_mean;
              color=C_FSP, linewidth=LW, markersize=6, label="TV(SSA, FSP)")
axislegend(ax_c; position=:lb, framevisible=false, labelsize=SMALL - 1)

colgap!(fig.layout, 10)
CairoMakie.save(joinpath(OUTDIR, "fig_mixed_fsp_ssa.pdf"), fig)
CairoMakie.save(joinpath(OUTDIR, "fig_mixed_fsp_ssa.png"), fig; px_per_unit=2)
println("\nSaved: paper/figures/fig_mixed_fsp_ssa.{pdf,png}")
