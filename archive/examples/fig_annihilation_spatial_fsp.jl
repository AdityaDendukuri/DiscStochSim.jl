"""
2D A+B→∅ Annihilation — Spatially Adaptive Per-Voxel FSP
==========================================================

A starts concentrated at corner (1,1), B at corner (K,K).
They diffuse toward each other and annihilate on contact: A + B → ∅.

Each active voxel stores the full joint distribution P(n_A, n_B, t),
evolved under the per-voxel 2-species CME with mean-field diffusion
coupling to neighbours. The nonlinear annihilation term k_r·n_A·n_B
is treated exactly within each voxel.

Run:  julia --project examples/fig_annihilation_spatial_fsp.jl
      SPATIAL_K=10 SPATIAL_KR=0.5 julia --project examples/fig_annihilation_spatial_fsp.jl
"""

using DiscStochSim
using ExponentialUtilities
using SparseArrays
using CairoMakie, Printf

# ── Parameters ────────────────────────────────────────────────────────────────
const K     = parse(Int,     get(ENV, "SPATIAL_K",    "6"))
const D     = parse(Float64, get(ENV, "SPATIAL_D",    "2.0"))
const k_r   = parse(Float64, get(ENV, "SPATIAL_KR",   "0.5"))
const N_A0  = parse(Int,     get(ENV, "SPATIAL_NA",   "10"))
const N_B0  = parse(Int,     get(ENV, "SPATIAL_NB",   "10"))
const n_max = parse(Int,     get(ENV, "SPATIAL_NMAX", "12"))
const dt    = parse(Float64, get(ENV, "SPATIAL_DT",   "0.1"))
const t_end = parse(Float64, get(ENV, "SPATIAL_TEND", "8.0"))
const ε_flux = 0.03   # activate empty neighbour when cumulative flux exceeds this
const krylov_m = 25

@printf("K=%d×%d  D=%.2f  k_r=%.2f  N_A0=%d  N_B0=%d  n_max=%d\n",
        K, K, D, k_r, N_A0, N_B0, n_max)

# ── Linear index for 2-species state (n_A, n_B) ───────────────────────────────
const NV = n_max + 1            # states per species
const NS = NV * NV              # states per voxel

lin(nA, nB) = nA * NV + nB + 1  # 1-based linear index

function voxel_means(p::Vector{Float64})
    μA = 0.0; μB = 0.0
    for nA in 0:n_max, nB in 0:n_max
        pr = p[lin(nA, nB)]
        μA += nA * pr; μB += nB * pr
    end
    μA, μB
end

# ── Per-voxel 2-species generator ────────────────────────────────────────────
"""
Build sparse generator for one voxel.
  Reactions:
    A inflow  (mean-field): ∅ → A  at rate μ_A_in
    A outflow (per mol):    A → ∅  at rate d·n_nb·n_A
    B inflow  (mean-field): ∅ → B  at rate μ_B_in
    B outflow (per mol):    B → ∅  at rate d·n_nb·n_B
    Annihilation:           A+B→∅  at rate k_r·n_A·n_B
"""
function build_voxel_gen(d, n_nb, μ_A_in, μ_B_in)
    Is = Int[]; Js = Int[]; Vs = Float64[]
    add!(row, col, val) = (push!(Is, row); push!(Js, col); push!(Vs, val))

    for nA in 0:n_max, nB in 0:n_max
        j = lin(nA, nB); diag = 0.0

        # A inflow
        if nA < n_max && μ_A_in > 0
            r = μ_A_in; add!(lin(nA+1, nB), j, r); diag -= r
        end
        # A outflow
        if nA > 0
            r = d * n_nb * nA; add!(lin(nA-1, nB), j, r); diag -= r
        end
        # B inflow
        if nB < n_max && μ_B_in > 0
            r = μ_B_in; add!(lin(nA, nB+1), j, r); diag -= r
        end
        # B outflow
        if nB > 0
            r = d * n_nb * nB; add!(lin(nA, nB-1), j, r); diag -= r
        end
        # Annihilation
        if nA > 0 && nB > 0
            r = k_r * nA * nB; add!(lin(nA-1, nB-1), j, r); diag -= r
        end

        add!(j, j, diag)
    end
    sparse(Is, Js, Vs, NS, NS)
end

# ── Grid helpers ──────────────────────────────────────────────────────────────
function grid_neighbors(i, j)
    nb = Tuple{Int,Int}[]
    i > 1 && push!(nb, (i-1, j))
    i < K && push!(nb, (i+1, j))
    j > 1 && push!(nb, (i, j-1))
    j < K && push!(nb, (i, j+1))
    nb
end

# ── Initial state ─────────────────────────────────────────────────────────────
# Initialise ALL voxels active from t=0 to avoid conservation loss:
# molecules leaving a voxel must always have an active destination.
active = Dict{Tuple{Int,Int}, Vector{Float64}}()
p_empty = zeros(NS); p_empty[lin(0, 0)] = 1.0
for i in 1:K, j in 1:K
    active[(i,j)] = copy(p_empty)
end
p0_A = zeros(NS); p0_A[lin(N_A0, 0)] = 1.0
p0_B = zeros(NS); p0_B[lin(0, N_B0)] = 1.0
active[(1, 1)] = p0_A
active[(K, K)] = p0_B

# ── Recording ─────────────────────────────────────────────────────────────────
snap_ts    = unique([0.0; filter(<(t_end), [3.0, 7.0, 12.0]); t_end])
snap_steps = Set(max(0, round(Int, ts/dt)) for ts in snap_ts)

snap_μA = Dict{Int, Matrix{Float64}}()
snap_μB = Dict{Int, Matrix{Float64}}()

t_hist   = Float64[0.0]
NA_hist  = Float64[Float64(N_A0)]
NB_hist  = Float64[Float64(N_B0)]
nact_hist = Int[length(active)]

function record_snap!(step, active)
    μA = zeros(K, K); μB = zeros(K, K)
    for ((i,j), p) in active
        a, b = voxel_means(p); μA[i,j] = a; μB[i,j] = b
    end
    snap_μA[step] = μA; snap_μB[step] = μB
end

record_snap!(0, active)

# ── Time loop ─────────────────────────────────────────────────────────────────
n_steps = round(Int, t_end / dt)
println("Running $n_steps steps on $(K)×$(K) grid...")

for step in 1:n_steps
    global active

    μA_map = Dict{Tuple{Int,Int}, Float64}()
    μB_map = Dict{Tuple{Int,Int}, Float64}()
    for ((i,j), p) in active
        a, b = voxel_means(p)
        μA_map[(i,j)] = a; μB_map[(i,j)] = b
    end

    new_active = Dict{Tuple{Int,Int}, Vector{Float64}}()
    for ((i,j), p) in active
        nbs    = grid_neighbors(i, j)
        n_nb   = length(nbs)
        μ_A_in = D * sum(get(μA_map, nb, 0.0) for nb in nbs; init=0.0)
        μ_B_in = D * sum(get(μB_map, nb, 0.0) for nb in nbs; init=0.0)
        Q      = build_voxel_gen(D, n_nb, μ_A_in, μ_B_in)
        p_new  = expv(Float64(dt), Q, p; m = min(krylov_m, NS))
        p_new .= max.(0.0, p_new)
        s = sum(p_new); s > 0.0 && (p_new ./= s)
        new_active[(i,j)] = p_new
    end
    active = new_active

    NA = sum(get(μA_map, v, 0.0) for v in keys(active); init=0.0)
    NB = sum(get(μB_map, v, 0.0) for v in keys(active); init=0.0)
    t  = step * dt
    push!(t_hist, t); push!(NA_hist, NA); push!(NB_hist, NB)
    push!(nact_hist, length(active))
    step in snap_steps && record_snap!(step, active)
    step % 20 == 0 && @printf("  t=%5.1f  active=%3d  E[NA]=%.2f  E[NB]=%.2f\n",
                               t, length(active), NA, NB)
end
println("Done.")

# ── Figure ────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.7, xgridvisible=false, ygridvisible=false,
          ticksize=3, tickwidth=0.6f0)))

n_snaps = length(snap_ts)
fig = Figure(size=(220*n_snaps + 120, 800))

Label(fig[0, 1:n_snaps],
      "A+B→∅ Annihilation  |  $(K)×$(K) RDME  D=$D  k_r=$k_r  " *
      "N_A0=$N_A0 (corner (1,1))  N_B0=$N_B0 (corner ($K,$K))  " *
      "|  per-voxel joint FSP: P(n_A,n_B,t)";
      fontsize=11, tellwidth=false)

μA_max = max(maximum(maximum(m) for m in values(snap_μA)), 0.1)
μB_max = max(maximum(maximum(m) for m in values(snap_μB)), 0.1)

sorted_snaps = sort(collect(keys(snap_μA)))

# ── Rows 1–2: E[nA] and E[nB] heatmaps ──────────────────────────────────────
for (col, step) in enumerate(sorted_snaps)
    ts = t_hist[step + 1]

    # Row 1: E[n_A]
    ax_a = Axis(fig[1, col]; aspect=DataAspect(),
                title="t = $(round(ts; digits=1))", titlesize=11,
                xticksvisible=false, yticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false)
    hm_a = heatmap!(ax_a, 1:K, 1:K, snap_μA[step]';
                    colormap=:Reds, colorrange=(0.0, μA_max))
    col == n_snaps && Colorbar(fig[1, col+1], hm_a; width=12, label="E[nA]",
                               ticks=([0, μA_max/2, μA_max],
                                      ["0", "$(round(μA_max/2,digits=1))", "$N_A0"]))

    # Row 2: E[n_B]
    ax_b = Axis(fig[2, col]; aspect=DataAspect(),
                xticksvisible=false, yticksvisible=false,
                xticklabelsvisible=false, yticklabelsvisible=false)
    hm_b = heatmap!(ax_b, 1:K, 1:K, snap_μB[step]';
                    colormap=:Blues, colorrange=(0.0, μB_max))
    col == n_snaps && Colorbar(fig[2, col+1], hm_b; width=12, label="E[nB]",
                               ticks=([0, μB_max/2, μB_max],
                                      ["0", "$(round(μB_max/2,digits=1))", "$N_B0"]))
end

# ── Row 3 left: total molecule counts over time ───────────────────────────────
ax_N = Axis(fig[3, 1:max(2, n_snaps÷2)];
            xlabel="time", ylabel="total molecules",
            title="E[N_A](t) and E[N_B](t) — FSP vs mean-field",
            titlesize=10)
lines!(ax_N, t_hist, NA_hist; color=:crimson,  linewidth=2.5, label="E[N_A] FSP")
lines!(ax_N, t_hist, NB_hist; color=:royalblue, linewidth=2.5, label="E[N_B] FSP")

# Mean-field ODE reference
function run_mf(K, D, k_r, N_A0, N_B0, t_end, dt_mf)
    μA = zeros(K, K); μB = zeros(K, K)
    μA[1,1] = Float64(N_A0); μB[K,K] = Float64(N_B0)
    ts = Float64[0.0]; NAs = Float64[N_A0]; NBs = Float64[N_B0]
    t = 0.0; s = 0
    while t < t_end
        dA = similar(μA); dB = similar(μB)
        fill!(dA, 0.0); fill!(dB, 0.0)
        for i in 1:K, j in 1:K
            rxn = k_r * μA[i,j] * μB[i,j]
            lap_A = (i>1 ? μA[i-1,j] : 0.0) + (i<K ? μA[i+1,j] : 0.0) +
                    (j>1 ? μA[i,j-1] : 0.0) + (j<K ? μA[i,j+1] : 0.0) -
                    ((i>1)+(i<K)+(j>1)+(j<K)) * μA[i,j]
            lap_B = (i>1 ? μB[i-1,j] : 0.0) + (i<K ? μB[i+1,j] : 0.0) +
                    (j>1 ? μB[i,j-1] : 0.0) + (j<K ? μB[i,j+1] : 0.0) -
                    ((i>1)+(i<K)+(j>1)+(j<K)) * μB[i,j]
            dA[i,j] = D * lap_A - rxn
            dB[i,j] = D * lap_B - rxn
        end
        μA .+= dt_mf .* dA; μB .+= dt_mf .* dB
        μA .= max.(0.0, μA); μB .= max.(0.0, μB)
        t += dt_mf; s += 1
        s % 5 == 0 && (push!(ts, t); push!(NAs, sum(μA)); push!(NBs, sum(μB)))
    end
    ts, NAs, NBs
end

print("Running mean-field ODE..."); flush(stdout)
ts_mf, NA_mf, NB_mf = run_mf(K, D, k_r, N_A0, N_B0, t_end, dt/10)
println(" done.")
lines!(ax_N, ts_mf, NA_mf; color=:crimson,  linewidth=1.5, linestyle=:dash, label="E[N_A] MF")
lines!(ax_N, ts_mf, NB_mf; color=:royalblue, linewidth=1.5, linestyle=:dash, label="E[N_B] MF")
axislegend(ax_N; position=:rt, labelsize=9, framevisible=false)

# ── Row 3 right: active voxel count ───────────────────────────────────────────
ax_act = Axis(fig[3, max(2, n_snaps÷2)+1:n_snaps];
              xlabel="time", ylabel="active voxels",
              title="Active voxel count (FSP tracks joint P(nA,nB) here)",
              titlesize=10)
lines!(ax_act, t_hist, Float64.(nact_hist); color=:darkorange, linewidth=2.5)
hlines!(ax_act, [Float64(K*K)]; color=(:black,0.3), linestyle=:dot, linewidth=1)
text!(ax_act, t_end*0.6, K^2*0.97; text="K²=$(K^2)", fontsize=9, color=:black)

rowgap!(fig.layout, 1, 6); rowgap!(fig.layout, 2, 10)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_annihilation_spatial_fsp.pdf")
save(out, fig)
println("Saved → $out")
