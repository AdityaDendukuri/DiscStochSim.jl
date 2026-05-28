"""
schlogl_ring_anim.jl

Animated GIF: Schlogl bistable RDME wave on a 2D ring membrane.
Coloring encodes bimodality directly:
  R = P(n ≤ n_un)  — probability in lo basin  → red
  G = sin(π·P_hi)  — peaks when P_lo≈P_hi≈0.5 → green glow at wavefront
  B = P(n > n_un)  — probability in hi basin  → blue

The bimodal wavefront glows bright green between the red (lo) and blue (hi) regions.

Output: paper/figures/anim_schlogl_ring.gif
"""

using DiscStochSim, CairoMakie, Printf

# ── Grid ───────────────────────────────────────────────────────────────────────
const K_ang    = 800
const K_rad    = 50
const SEED_CTR = K_ang ÷ 2 + 1
const SEED_HALF = 1    # single-voxel nucleation site
const R_IN = 0.50
const R_OUT = 1.00

# ── Simulation ─────────────────────────────────────────────────────────────────
function run_wave(c3, T_snaps; dt=0.5)
    model = SchloglModel1D(2.0, 0.028, 3.2e-4, c3, 1.0)
    n_lo, n_un, n_hi = schlogl_fixed_points(model)
    @printf("  n_lo=%d  n_un=%d  n_hi=%d\n", n_lo, n_un, n_hi)

    s = SpatialFSP(model, K_ang, K_rad;
                   n_max=n_hi+30, n_max_eq=n_hi+30,
                   ε_expand=0.15, ε_equil=0.12,
                   use_block_vcycle=false)
    nm = s.n_max_eq

    seed_j = K_rad ÷ 2 + 1   # radial midpoint of ring
    for i in 1:K_ang, j in 1:K_rad
        hi = abs(i - SEED_CTR) ≤ SEED_HALF && abs(j - seed_j) ≤ SEED_HALF
        p  = zeros(nm+1); p[min((hi ? n_hi : n_lo)+1, nm+1)] = 1.0
        s.equil_dists[CartesianIndex(i,j)] = p
    end

    # Snapshot: P_lo and P_hi matrices (lo basin = n≤n_un, hi basin = n>n_un)
    p0_lo = zeros(nm+1); p0_lo[n_lo+1] = 1.0   # default = equil-lo
    function snap()
        PLo = zeros(Float32, K_ang, K_rad)
        PHi = zeros(Float32, K_ang, K_rad)
        for i in 1:K_ang, j in 1:K_rad
            idx = CartesianIndex(i,j)
            p   = get(s.dists, idx, get(s.equil_dists, idx, p0_lo))
            n   = min(n_un+1, length(p))
            PLo[i,j] = Float32(sum(@view p[1:n]))
            PHi[i,j] = Float32(sum(@view p[n+1:end]))
        end
        PLo, PHi
    end

    frames    = Dict{Float64, Tuple{Matrix{Float32}, Matrix{Float32}}}()
    remaining = Set(T_snaps)
    0.0 ∈ T_snaps && (frames[0.0] = snap(); delete!(remaining, 0.0))

    t0 = time()
    for step_i in 1:round(Int, maximum(T_snaps)/dt)
        step!(s, dt)
        t = round(step_i*dt; digits=6)
        for ts in collect(remaining)
            if t ≥ ts - 1e-9
                frames[ts] = snap()
                delete!(remaining, ts)
                step_i % 40 == 0 &&
                    @printf("  t=%5.0f  active=%3d  (%.0fs)\n",
                            ts, n_active(s), time()-t0)
                flush(stdout)
            end
        end
        isempty(remaining) && break
    end
    @printf("  Done %.0fs\n", time()-t0); flush(stdout)
    frames, n_lo, n_un, n_hi
end

# ── Bimodality colour: R=P_lo, G=sin(π·P_hi), B=P_hi ──────────────────────────
# Pure lo  (P_hi=0): red   (1, 0,   0)
# Bimodal  (P_hi=½): green (½, 1,   ½)   ← wavefront glows
# Pure hi  (P_hi=1): blue  (0, 0,   1)
@inline function bimodal_rgb(phi_smooth::Float32)
    G = Float32(sin(π * phi_smooth))
    RGBf(1f0 - phi_smooth, G, phi_smooth)
end

# Circular Gaussian blur in angular direction so wavefront band is visible.
# The true bimodal zone is 1-2 voxels; σ≈15 makes it ~45 voxels wide (≈5% of ring).
function smooth_phi!(PHi::Matrix{Float32}; σ::Float32=15f0)
    K, R = size(PHi)
    rad  = round(Int, 3σ)
    w    = Float32[exp(-0.5f0*(k/σ)^2) for k in -rad:rad]
    w  ./= sum(w)
    tmp  = similar(PHi, K)
    for j in 1:R
        fill!(tmp, 0f0)
        for (ki, wi) in enumerate(w)
            shift = ki - rad - 1
            for i in 1:K
                tmp[i] += wi * PHi[mod1(i + shift, K), j]
            end
        end
        PHi[:,j] .= tmp
    end
end

const _BG_OUT  = RGBf(0.90f0, 0.90f0, 0.92f0)
const _BG_IN   = RGBf(0.82f0, 0.82f0, 0.86f0)
const NIMG_ANIM = 800

function make_ring_image(PLo::Matrix{Float32}, PHi::Matrix{Float32}; Nimg=NIMG_ANIM)
    img   = fill(_BG_OUT, Nimg, Nimg)
    half  = (Nimg + 1.0) / 2.0
    scale = half / 1.35
    @inbounds for py in 1:Nimg, px in 1:Nimg
        x = (px - half) / scale
        y = (py - half) / scale
        r = sqrt(x*x + y*y)
        if r < R_IN
            img[px, py] = _BG_IN
        elseif r ≤ R_OUT
            θ = atan(y, x)
            i = mod1(round(Int, θ * K_ang / (2π)) + SEED_CTR, K_ang)
            j = clamp(ceil(Int, (r-R_IN)/(R_OUT-R_IN)*K_rad), 1, K_rad)
            img[px, py] = bimodal_rgb(PHi[i,j])
        end
    end
    img
end

# ── Run ────────────────────────────────────────────────────────────────────────
const N_FRAMES = 80
const T_MAX    = 550.0
const DT_FRAME = T_MAX / N_FRAMES

T_snaps = unique(sort([0.0; [round(DT_FRAME*k; digits=4) for k in 1:N_FRAMES]]))

@printf("Simulating %d frames to t=%.0f …\n", N_FRAMES, T_MAX); flush(stdout)
frames, n_lo, n_un, n_hi = run_wave(21.0, T_snaps)

# ── Build custom bimodality colorbar ───────────────────────────────────────────
# Gradient: red → green → blue  along P_hi = 0 → 0.5 → 1
_bimo_cmap = let p = range(0f0, 1f0; length=256)
    [RGBAf(1f0-v, sin(Float32(π)*v), v, 1f0) for v in p]
end

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(backgroundcolor=:white))
fig = Figure(size=(940, 940))

Label(fig[1, 1:2];
      text="Schlogl bistable RDME — ring membrane  (c3=21, D=2)",
      fontsize=14, color=RGBf(0.12,0.12,0.18), font=:bold,
      halign=:center, tellwidth=false)

ax = Axis(fig[2, 1];
          aspect=DataAspect(),
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false,
          leftspinevisible=false, rightspinevisible=false,
          topspinevisible=false, bottomspinevisible=false,
          backgroundcolor=RGBf(0.90f0, 0.90f0, 0.92f0))

PLo0, PHi0 = frames[0.0]
smooth_phi!(PHi0)
img_obs = Observable(make_ring_image(PLo0, PHi0))
image!(ax, -1.35..1.35, -1.35..1.35, img_obs)

# Seed dot marker at ring midpoint
R_mid  = (R_IN + R_OUT) / 2
scatter!(ax, [R_mid], [0.0]; color=(:white, 0.85), markersize=7, strokecolor=:black, strokewidth=1)

xlims!(ax, -1.35, 1.35); ylims!(ax, -1.35, 1.35)

# Colorbar
Colorbar(fig[2, 2];
         colormap=_bimo_cmap,
         limits=(0.0, 1.0),
         label="P(hi) = P(n > n_un)",
         ticks=([0.0, 0.5, 1.0],
                ["0\n(lo basin\nred)", "0.5\n(bimodal\ngreen)", "1\n(hi basin\nblue)"]),
         vertical=true, labelsize=10, ticklabelsize=9, width=18,
         tellheight=false)

# Annotation
Label(fig[3, 1:2];
      text="G = sin(π·P_hi): wavefront bimodality glows green between red (lo) and blue (hi) regions",
      fontsize=10, color=RGBf(0.30,0.30,0.40), halign=:center, tellwidth=false)

t_lbl = Label(fig[0, 1:2]; text="t = 0",
              fontsize=12, color=RGBf(0.25,0.25,0.40),
              halign=:center, tellwidth=false)

rowsize!(fig.layout, 2, Fixed(800))
colsize!(fig.layout, 2, Fixed(70))
colgap!(fig.layout, 6)

# ── Record ─────────────────────────────────────────────────────────────────────
outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "anim_schlogl_ring.gif")

@printf("Rendering %d frames → %s\n", length(T_snaps), out); flush(stdout)
t_render = time()

record(fig, out, T_snaps; framerate=15) do ts
    PLo, PHi  = frames[ts]
    smooth_phi!(PHi)          # blur wavefront band to make it visible
    img_obs[] = make_ring_image(PLo, PHi)
    t_lbl.text[] = @sprintf("t = %.0f", ts)
end

@printf("Saved in %.0fs → %s\n", time()-t_render, out)
