"""
anim_fsp_20x20.jl — FSP on 20×20 grid: adaptive three-phase state space

Uses SpatialFSP (BirthDeathRDME) to show R-step state space expansion.

IC: N0 molecules at center.  Wavefront spreads outward from center.

Three phases tracked per voxel:
  EMPTY      — not yet reached; zero probability assigned
  ACTIVE     — wavefront ring; full per-voxel CME distribution tracked
  EQUIL      — equilibrated; collapsed to Poisson(μ_ss)

Left panel : FSP mean ⟨n_k⟩ heatmap  (= mean-field for this linear system)
Right panel: phase map  (dark=empty, orange=active/fine, blue=equil/coarse)

For linear systems (∅→X, X→∅, diffusion), per-voxel Poisson is exact.
State space at peak ≈ (ring width) × K × n_max — far smaller than the full
joint space C(N_total + K² − 1, K² − 1).

Run: julia --project examples/anim_fsp_20x20.jl
"""

using DiscStochSim, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const K    = 100         # K×K grid
const k_b  = 1.0         # birth rate per voxel
const k_d  = 0.5         # death rate per molecule  →  μ_ss = 2 per voxel
const D    = 1.0
const dx   = 1.0
const N0   = 1           # IC molecules at center
const dt   = 0.3
const T_END = 80.0
const N_FRM = 100

const n_max      = 16    # truncation per voxel  (4 × μ_ss + margin)
const ε_expand   = 0.15  # activate when expected flux ≥ this
const ε_equil    = 0.07  # equilibrate when within 7% of μ_ss

model  = BirthDeathRDME(k_b, k_d, D, dx)
μ_ss   = ss_mean(model)   # = 2.0
center = CartesianIndex(K÷2 + 1, K÷2 + 1)

@printf("Grid %d×%d  μ_ss=%.1f  T=%.0f  N_frames=%d\n", K, K, μ_ss, T_END, N_FRM)
@printf("State space upper bound: active_voxels × n_max = %d × %d = %d\n",
        K^2, n_max, K^2 * n_max)

# ── Initialise SpatialFSP ──────────────────────────────────────────────────────
state = SpatialFSP(model, K, K;
                   n_max,
                   ε_expand,
                   ε_equil,
                   use_joint_active  = false,   # per-voxel block mode
                   use_block_vcycle  = true)
set_ic!(state, center, N0)

# ── Pre-compute frames ─────────────────────────────────────────────────────────
t_frames      = collect(range(0.0, T_END; length = N_FRM))
mean_frames   = Vector{Matrix{Float64}}(undef, N_FRM)
status_frames = Vector{Matrix{Int}}(undef, N_FRM)
kc_frames     = zeros(Int, N_FRM)   # "active" voxels = fine-grained FSP voxels

mean_frames[1]   = mean_grid(state)
status_frames[1] = status_grid(state)
kc_frames[1]     = n_active(state)

# ── RDME matrix renderer ───────────────────────────────────────────────────────
# Block size per voxel: active → n_max+1 states, equil → 1 state (collapsed),
# empty → not shown.  Variable block sizes show the coarsening directly.
const IMGSIZE = 280
const COL_DIAG  = RGBAf(0.27, 0.51, 0.71, 1.0)
const COL_COUPL = RGBAf(0.89, 0.45, 0.13, 1.0)
const COL_ZERO  = RGBAf(0.94, 0.94, 0.94, 1.0)
const COL_SEP   = RGBAf(0.30, 0.30, 0.30, 1.0)

function vox_neighbors(idx::CartesianIndex{2}, Kx, Ky)
    i, j = Tuple(idx)
    ns = CartesianIndex{2}[]
    i > 1  && push!(ns, CartesianIndex(i-1, j))
    i < Kx && push!(ns, CartesianIndex(i+1, j))
    j > 1  && push!(ns, CartesianIndex(i, j-1))
    j < Ky && push!(ns, CartesianIndex(i, j+1))
    ns
end

function render_rdme_matrix(s::SpatialFSP, nm::Int)
    active = sort!(collect(keys(s.dists)),       by = Tuple)
    equil  = sort!(collect(keys(s.equil_dists)), by = Tuple)

    isempty(active) && isempty(equil) && return fill(COL_ZERO, IMGSIZE, IMGSIZE)

    # State layout:
    #   Active voxels  → each gets (nm+1) states  (fine CME blocks)
    #   Equil region   → ALL equil voxels merged into ONE state
    #                    (the whole equil domain is one coarse super-voxel)
    # At full equilibrium: 0 active + 1 equil state → 1×1 matrix.
    n_act_states   = length(active) * (nm + 1)
    n_equil_states = isempty(equil) ? 0 : 1     # one merged super-voxel
    total          = n_act_states + n_equil_states

    # Map state index → active voxel (equil states map to nothing)
    stov = Vector{Union{CartesianIndex{2}, Nothing}}(undef, total)
    for (i, v) in enumerate(active)
        base = (i - 1) * (nm + 1)
        for si in 1:nm+1;  stov[base + si] = v;  end
    end
    if n_equil_states == 1;  stov[total] = nothing;  end  # equil super-voxel

    # Active-active adjacency
    active_set = Set(active)
    adj_aa = Set{NTuple{2, CartesianIndex{2}}}()
    for v in active, nb in vox_neighbors(v, s.K_x, s.K_y)
        nb ∈ active_set || continue
        push!(adj_aa, Tuple(v) < Tuple(nb) ? (v, nb) : (nb, v))
    end

    # Active voxels that border the equil region (coupled to equil super-voxel)
    equil_set = Set(equil)
    active_at_equil_border = Set{CartesianIndex{2}}()
    for v in active, nb in vox_neighbors(v, s.K_x, s.K_y)
        nb ∈ equil_set && push!(active_at_equil_border, v)
    end

    img = fill(COL_ZERO, IMGSIZE, IMGSIZE)
    for r in 1:IMGSIZE
        sr = clamp(floor(Int, (r - 0.5) / IMGSIZE * total) + 1, 1, total)
        vr = stov[sr]
        for c in 1:IMGSIZE
            sc = clamp(floor(Int, (c - 0.5) / IMGSIZE * total) + 1, 1, total)
            vc = stov[sc]

            if isnothing(vr) && isnothing(vc)
                img[r, c] = COL_DIAG    # equil super-voxel diagonal block
            elseif !isnothing(vr) && !isnothing(vc)
                if vr == vc
                    img[r, c] = COL_DIAG
                else
                    key = Tuple(vr) < Tuple(vc) ? (vr, vc) : (vc, vr)
                    key ∈ adj_aa && (img[r, c] = COL_COUPL)
                end
            else
                # One side is equil super-voxel, other is active
                av = isnothing(vr) ? vc : vr
                av ∈ active_at_equil_border && (img[r, c] = COL_COUPL)
            end
        end
    end

    # Separator between active and equil super-voxel
    if n_act_states > 0 && n_equil_states > 0
        sep = clamp(round(Int, n_act_states / total * IMGSIZE), 1, IMGSIZE)
        for q in 1:IMGSIZE
            img[sep, q] = COL_SEP
            img[q, sep] = COL_SEP
        end
    end

    img
end

mat_frames = Vector{Matrix{RGBAf}}(undef, N_FRM)
mat_frames[1] = render_rdme_matrix(state, n_max)

println("Stepping FSP …")
let frame_idx = 2
    for step in 1:round(Int, T_END / dt)
        step!(state, dt)
        t_now = step * dt
        if frame_idx <= N_FRM && t_now >= t_frames[frame_idx] - 1e-9
            mean_frames[frame_idx]   = mean_grid(state)
            status_frames[frame_idx] = status_grid(state)
            kc_frames[frame_idx]     = n_active(state)
            mat_frames[frame_idx]    = render_rdme_matrix(state, n_max)
            if frame_idx % 20 == 0
                @printf("  frame %3d  t=%5.1f  empty=%3d  active=%3d  equil=%3d\n",
                        frame_idx, t_now,
                        n_empty(state), n_active(state), n_equilibrated(state))
            end
            frame_idx += 1
        end
        frame_idx > N_FRM && break
    end
end
println("Done.  Peak active voxels: $(maximum(kc_frames))")

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1350, 530))

MU_MAX = μ_ss * 2.5

mu_obs  = Observable(mean_frames[1]')
sg_obs  = Observable(Float32.(status_frames[1]'))
mat_obs = Observable(mat_frames[1])

ax_kw = (aspect=DataAspect(),
    xticks=(1:5:K, string.(1:5:K)), yticks=(1:5:K, string.(1:5:K)),
    yreversed=true, xaxisposition=:top, titlesize=11,
    xticklabelsize=8, yticklabelsize=8)

# Panel A: FSP mean heatmap
ax_mu = Axis(fig[1,1]; title="FSP mean  ⟨n_k⟩  (exact CME, V-cycle solver)", ax_kw...)
hm = heatmap!(ax_mu, 1:K, 1:K, mu_obs;
    colormap=:inferno, colorrange=(0, MU_MAX), highclip=:white, lowclip=:black)
scatter!(ax_mu, [center[2]], [center[1]]; marker=:cross, markersize=10,
    color=:cyan, strokewidth=0)
Colorbar(fig[1,1][1,2], hm; label="⟨n⟩", width=10, labelsize=8, ticklabelsize=6)

# Panel B: phase map
PHASE_CMAP = cgrad([colorant"#1a1a2e", colorant"#e07b39", colorant"#4472c4"])
ax_ph = Axis(fig[1,2]; title="Adaptive state space  (orange=active, blue=equil)", ax_kw...)
heatmap!(ax_ph, 1:K, 1:K, sg_obs; colormap=PHASE_CMAP, colorrange=(0,2))
scatter!(ax_ph, [center[2]], [center[1]]; marker=:cross, markersize=10,
    color=:white, strokewidth=0)

# Panel C: RDME block matrix with variable block sizes
ax_mx = Axis(fig[1,3];
    title   = "RDME generator  (block size ∝ voxel state space)",
    aspect  = DataAspect(),
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    yreversed=true, xaxisposition=:top, titlesize=11)
image!(ax_mx, 0.5..K+0.5, 0.5..K+0.5, mat_obs)

colsize!(fig.layout, 1, Aspect(1,1.0))
colsize!(fig.layout, 2, Aspect(1,1.0))
colsize!(fig.layout, 3, Aspect(1,1.0))

t_label  = Label(fig[0,1:2], "t = 0.0"; fontsize=13, halign=:left)
kc_label = Label(fig[0,3],
    "active: $(kc_frames[1])  |  states: $(kc_frames[1]*(n_max+1))";
    fontsize=11, halign=:center)

Legend(fig[2,1:3],
    [PolyElement(color=colorant"#1a1a2e", strokecolor=:gray40, strokewidth=1),
     PolyElement(color=colorant"#e07b39"),
     PolyElement(color=colorant"#4472c4"),
     PolyElement(color=RGBf(0.27,0.51,0.71)),
     PolyElement(color=RGBf(0.89,0.45,0.13))],
    ["Empty  (0 states)",
     "Active  ($(n_max+1) states/voxel, V-cycle CME)",
     "Equil  (1 state/voxel, Poisson collapsed)",
     "Diagonal block  (local CME)",
     "Off-diagonal block  (diffusion coupling)"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(16,11))

rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,4)
colgap!(fig.layout,1,8); colgap!(fig.layout,2,8)

# ── Record ─────────────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "..", "paper", "figures", "anim_fsp_100x100.gif")
println("Rendering → $out")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]        = mean_frames[f]'
    sg_obs[]        = Float32.(status_frames[f]')
    mat_obs[]       = mat_frames[f]
    t_label.text[]  = "t = $(round(t_frames[f], digits=1))"
    na = kc_frames[f]
    kc_label.text[] = "active: $na  |  FSP states: $(na*(n_max+1)) + $(400-na) equil"
end
println("Saved → $out")
println("Kc (active): $(kc_frames[1]) → $(maximum(kc_frames)) → $(kc_frames[end])")
