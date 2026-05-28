"""
schlogl_ring.jl

Schlogl bistable RDME on a 2D ring membrane — polarisome analog.
The ring is represented as a K_ang × K_rad rectangular grid; each voxel (i,j)
maps to angle θ_i and radius R_j in the annulus visualization.

  c3 = 21: bistable wave fills the whole ring
  c3 = 17: wave stalls → localized polarisome-like patch

Layout: 2 rows × 5 columns
  Row 1: ring diagrams at t1,t2,t3 (c3=21) + t_final (c3=17) + colorbar
  Row 2: ⟨n⟩ vs angle profiles (averaged over ring width)

Output: paper/figures/fig_schlogl_ring.png
"""

using DiscStochSim, CairoMakie, Printf, Statistics

# ── Grid ───────────────────────────────────────────────────────────────────────
const K_ang  = 800    # voxels along ring circumference
const K_rad  = 50     # voxels across ring thickness
const SEED_CTR   = K_ang ÷ 2 + 1      # angular centre of seed (index 401)
const SEED_HALF  = 40                  # seed half-width in voxels (= ±18°)

# Polar mapping: voxel (i,j) → (θ, R)
const R_IN   = 0.50   # inner ring radius (display units)
const R_OUT  = 1.00   # outer ring radius
θ_of(i)   = 2π * (i - SEED_CTR) / K_ang   # ∈ (-π, π]
R_of(j)   = R_IN + (j - 0.5) * (R_OUT - R_IN) / K_rad

# ── Run simulation ─────────────────────────────────────────────────────────────
function run_wave(c3, T_snaps; dt=0.5)
    model = SchloglModel1D(2.0, 0.028, 3.2e-4, c3, 1.0)
    n_lo, n_un, n_hi = schlogl_fixed_points(model)
    @printf("  c3=%.0f  n_lo=%d  n_un=%d  n_hi=%d\n", c3, n_lo, n_un, n_hi)

    s = SpatialFSP(model, K_ang, K_rad;
                   n_max        = n_hi + 30,
                   n_max_eq     = n_hi + 30,
                   ε_expand     = 0.15,
                   ε_equil      = 0.12,
                   use_block_vcycle = false)
    nm = s.n_max_eq

    # All voxels → equil_dists: seed arc = equil-hi delta, rest = equil-lo delta
    for i in 1:K_ang, j in 1:K_rad
        hi  = abs(i - SEED_CTR) ≤ SEED_HALF
        p   = zeros(nm + 1)
        p[min((hi ? n_hi : n_lo) + 1, nm + 1)] = 1.0
        s.equil_dists[CartesianIndex(i, j)] = p
    end

    # Profile snapshot: (K_ang × K_rad) matrix of ⟨n⟩
    p0 = zeros(nm + 1); p0[1] = 1.0
    function snap_profile()
        prof = zeros(K_ang, K_rad)
        for i in 1:K_ang, j in 1:K_rad
            idx = CartesianIndex(i, j)
            p   = get(s.dists, idx, get(s.equil_dists, idx, p0))
            prof[i, j] = sum((n-1)*p[n] for n in 1:length(p))
        end
        prof
    end

    profiles  = Dict{Float64, Matrix{Float64}}()
    remaining = Set(T_snaps)
    0.0 ∈ T_snaps && (profiles[0.0] = snap_profile(); delete!(remaining, 0.0))

    t0 = time()
    for step_i in 1:round(Int, maximum(T_snaps) / dt)
        step!(s, dt)
        t = round(step_i * dt; digits=6)
        for ts in collect(remaining)
            if t ≥ ts - 1e-9
                profiles[ts] = snap_profile()
                delete!(remaining, ts)
                @printf("  t=%5.0f  active=%4d  equil=%5d  (%.0fs)\n",
                        ts, n_active(s), n_equilibrated(s), time()-t0)
                flush(stdout)
            end
        end
        isempty(remaining) && break
    end
    profiles, n_lo, n_un, n_hi
end

# ── Colour helpers ─────────────────────────────────────────────────────────────
const _RDBU   = CairoMakie.to_colormap(:RdBu)
const _BG_OUT = RGBf(0.93f0, 0.93f0, 0.95f0)   # outside ring
const _BG_IN  = RGBf(0.86f0, 0.86f0, 0.89f0)   # ring interior disc
const NIMG    = 1400                              # pixel resolution per panel

rdbu_rgb(t) = let c = _RDBU[clamp(floor(Int, t*(length(_RDBU)-1))+1, 1, length(_RDBU))]
    RGBf(c.r, c.g, c.b)
end

# ── Ring image: pixel-remap (fast — one image! per panel, no poly! loops) ──────
# Each pixel (px,py) maps to Cartesian (x,y) → polar (r,θ) → voxel (i,j).
function make_ring_image(prof, n_lo, n_hi)
    img   = fill(_BG_OUT, NIMG, NIMG)
    half  = (NIMG + 1.0) / 2.0
    scale = half / 1.35    # display range ±1.35 units
    @inbounds for py in 1:NIMG, px in 1:NIMG
        x = (px - half) / scale
        y = (py - half) / scale
        r = sqrt(x*x + y*y)
        if r < R_IN
            img[px, py] = _BG_IN
        elseif r ≤ R_OUT
            θ = atan(y, x)   # ∈ (-π, π]
            i = mod1(round(Int, θ * K_ang / (2π)) + SEED_CTR, K_ang)
            j = clamp(ceil(Int, (r - R_IN) / (R_OUT - R_IN) * K_rad), 1, K_rad)
            v = clamp((prof[i, j] - n_lo) / (n_hi - n_lo), 0.0, 1.0)
            img[px, py] = rdbu_rgb(v)
        end
    end
    img
end

function draw_ring!(ax, prof, n_lo, n_hi)
    image!(ax, -1.35..1.35, -1.35..1.35, make_ring_image(prof, n_lo, n_hi))
    θ_s   = SEED_HALF * 2π / K_ang
    θ_arc = range(-θ_s, θ_s; length=120)
    lines!(ax, (R_OUT+0.08) .* cos.(θ_arc), (R_OUT+0.08) .* sin.(θ_arc);
           color=(:black, 0.5), linewidth=2.5)
end

# ── Angle profile (averaged over radial voxels) ────────────────────────────────
function angle_prof(prof)
    ang  = [θ_of(i) for i in 1:K_ang]
    vals = [mean(prof[i, :]) for i in 1:K_ang]
    perm = sortperm(ang)
    ang[perm], vals[perm]
end

# ── Run ────────────────────────────────────────────────────────────────────────
const D = 2.0

@printf("c3 = 21 (wave fills ring):\n"); flush(stdout)
T21  = [0.0, 90.0, 210.0, 360.0]
p21, n_lo21, n_un21, n_hi21 = run_wave(21.0, T21)

@printf("\nc3 = 14 (lo wins — collapsing polarisome):\n"); flush(stdout)
T17  = [0.0, 60.0]
p17, n_lo17, n_un17, n_hi17 = run_wave(14.0, T17)

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11,
    Axis=(spinewidth=0.6, xgridvisible=false, ygridvisible=false)))

fig = Figure(size=(1300, 700))

Label(fig[0, 1:5];
      text="Schlogl bistable wave on a ring membrane  (D=2, ring=$K_ang×$K_rad voxels, adaptive FSP)",
      fontsize=13, color=RGBf(0.15,0.15,0.20), font=:bold,
      halign=:center, tellwidth=false)

ring_axkw = (aspect=DataAspect(),
             xticksvisible=false, yticksvisible=false,
             xticklabelsvisible=false, yticklabelsvisible=false,
             leftspinevisible=false, rightspinevisible=false,
             topspinevisible=false, bottomspinevisible=false,
             backgroundcolor=RGBf(0.95, 0.95, 0.97))

t_labels21 = ["t = 0  (IC)", "t = 90", "t = 210", "t = 360"]

# Row 1: ring diagrams — c3=21 snapshots
for (c, ts) in enumerate(T21)
    ax = Axis(fig[1, c]; ring_axkw...,
              title=t_labels21[c], titlesize=11,
              titlecolor=RGBf(0.2,0.2,0.35))
    draw_ring!(ax, p21[ts], n_lo21, n_hi21)
    xlims!(ax, -1.35, 1.35); ylims!(ax, -1.35, 1.35)
end

# c3=17 in column 5
ax17 = Axis(fig[1, 5]; ring_axkw...,
            title="t = 60  (c3=14, lo wins)", titlesize=11,
            titlecolor=RGBf(0.6,0.15,0.15))
draw_ring!(ax17, p17[60.0], n_lo17, n_hi17)
xlims!(ax17, -1.35, 1.35); ylims!(ax17, -1.35, 1.35)

# Colorbar
Colorbar(fig[1:2, 6];
         colormap=:RdBu, limits=(0.0, 1.0),
         label="(⟨n⟩ − n_lo) / (n_hi − n_lo)\nred = lo basin   blue = hi basin",
         ticks=([0.0, Float64(n_un21-n_lo21)/(n_hi21-n_lo21), 1.0],
                ["n_lo\n=$n_lo21", "n_un\n=$n_un21", "n_hi\n=$n_hi21"]),
         vertical=true, labelsize=9, ticklabelsize=8, width=16)

# Row 2: angle profiles
prof_kw = (xticks=([-π,-π/2,0.0,π/2,π], ["-π","-π/2","0","π/2","π"]),
           xticklabelsvisible=true, xticksvisible=true)

for (c, ts) in enumerate(T21)
    ax = Axis(fig[2, c]; prof_kw...,
              xlabel = c == 2 ? "angle θ around ring" : "",
              ylabel = c == 1 ? "⟨n⟩" : "",
              yticklabelsvisible = c == 1, yticksvisible = c == 1,
              yticks = ([n_lo21, n_un21, n_hi21], ["$n_lo21", "$n_un21", "$n_hi21"]))
    ang, vals = angle_prof(p21[ts])
    band!(ax, ang, fill(Float64(n_lo21), length(ang)), vals; color=(:steelblue, 0.2))
    lines!(ax, ang, vals; color=:steelblue, linewidth=1.8)
    hlines!(ax, [Float64(n_lo21), Float64(n_un21), Float64(n_hi21)];
            color=(:gray40,0.4), linestyle=:dash, linewidth=0.9)
    xlims!(ax, -π, π); ylims!(ax, 0, n_hi21+20)
end

ax_p17 = Axis(fig[2, 5]; prof_kw...,
              yticklabelsvisible=false, yticksvisible=false,
              yticks=([n_lo17, n_un17, n_hi17], ["$n_lo17","$n_un17","$n_hi17"]))
ang17, vals17 = angle_prof(p17[60.0])

band!(ax_p17, ang17, fill(Float64(n_lo17), length(ang17)), vals17; color=(:crimson, 0.2))
lines!(ax_p17, ang17, vals17; color=:crimson, linewidth=1.8)
hlines!(ax_p17, [Float64(n_lo17), Float64(n_un17), Float64(n_hi17)];
        color=(:gray40,0.4), linestyle=:dash, linewidth=0.9)
xlims!(ax_p17, -π, π); ylims!(ax_p17, 0, n_hi17+20)

rowsize!(fig.layout, 1, Fixed(360))
rowsize!(fig.layout, 2, Fixed(160))
colgap!(fig.layout, 4); rowgap!(fig.layout, 6)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "fig_schlogl_ring.png")
save(out, fig; px_per_unit=3)
@printf("Saved → %s\n", out)
