using DiscStochSim, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const Kx   = 200; const Ky = 200
const k_b  = 1.0; const k_d = 0.5   # μ_ss = 2
const D    = 1.0
const n_max = 10
const max_r = 2        # TT bond dimension — r=2 is near-exact for birth-death
const dt   = 0.5
const T_END = 100.0
const N_FRM = 100

# ── Build solver ───────────────────────────────────────────────────────────────
println("Building RowTTSolver: $(Kx)×$(Ky) grid, r=$max_r, n_max=$n_max")
solver = RowTTSolver(Kx, Ky, k_b, k_d, D; n_max, max_r)

# IC: 1 molecule at centre
ci = Kx÷2 + 1; cj = Ky÷2 + 1
set_ic!(solver, ci, cj, 1)

@printf("Storage: %.1f MB (vs %.1e for full joint)\n",
        Kx * Ky * max_r^2 * (n_max+1) * 8 / 1e6,
        float(n_max+1)^(Kx*Ky))

# ── Colour helpers ─────────────────────────────────────────────────────────────
μ_ss = k_b / k_d
const _MAGMA = CairoMakie.to_colormap(:magma)
_magma(t)    = let c=_MAGMA[clamp(floor(Int,t*(length(_MAGMA)-1))+1,1,length(_MAGMA))];
               RGBf(c.r,c.g,c.b) end
const _mu_max = Ref(μ_ss * 2.5)

# ── Pre-compute frames ─────────────────────────────────────────────────────────
t_frames   = collect(range(0.0, T_END; length=N_FRM))
mu_frames  = Vector{Matrix{RGBf}}(undef, N_FRM)
bond_frames = Vector{Matrix{Float64}}(undef, N_FRM)

function make_frame(s)
    mg = mean_grid(s)
    bg = max_bond_grid(s)   # Kx × (Ky-1)
    # Mean heatmap
    mu_img = zeros(RGBf, Kx, Ky)
    for i in 1:Kx, j in 1:Ky
        mu_img[i, j] = _magma(clamp(mg[i,j] / _mu_max[], 0.0, 1.0))
    end
    mu_img, bg
end

mu_frames[1], bond_frames[1] = make_frame(solver)
println("Pre-computing $N_FRM frames…"); flush(stdout)

let fi = 2
    for step in 1:round(Int, T_END/dt)
        step!(solver, dt; krylov_m=20, tol=1e-12)
        t = step * dt
        if fi ≤ N_FRM && t ≥ t_frames[fi] - 1e-9
            mu_frames[fi], bond_frames[fi] = make_frame(solver)
            if fi % 10 == 0
                @printf("  frame %3d  t=%6.1f  total_μ=%.3f  max_r=%d\n",
                        fi, t, sum(solver.means), maximum(tt_max_bond(tt) for tt in solver.row_tts))
            end
            fi += 1
        end
        fi > N_FRM && break
    end
end

# ── Figure ─────────────────────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1100, 560))

mu_obs   = Observable(mu_frames[1])
bond_obs = Observable(bond_frames[1])

ax_kw = (aspect=DataAspect(), yreversed=true, xaxisposition=:top, titlesize=10,
         xticksvisible=false, yticksvisible=false,
         xticklabelsvisible=false, yticklabelsvisible=false)

ax_mu = Axis(fig[1,1];
    title = "Row-TT exact joint ⟨n_{k,j}⟩  ($(Kx)×$(Ky)  r=$max_r  μ_ss=$(μ_ss))",
    ax_kw...)
image!(ax_mu, 0.5..Kx+0.5, 0.5..Ky+0.5, mu_obs)

ax_bond = Axis(fig[1,2];
    title = "Within-row bond dimension  (correlations along each row)",
    ax_kw...)
hm = heatmap!(ax_bond, 1:Ky-1, 1:Kx,
              @lift($bond_obs');
              colormap=:viridis, colorrange=(1.0, Float64(max_r)))
Colorbar(fig[1,3], hm; label="bond dim r", width=14)

colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, Ky/(Ky-1)))

t_label  = Label(fig[0,1:2], "t = 0.0"; fontsize=13, halign=:left)
inf_label = Label(fig[2,1:2],
    "Row-TT: exact joint P(n^{k,1},...,n^{k,200}) per row | " *
    "mean-field between rows | storage=$(round(Kx*Ky*max_r^2*(n_max+1)*8/1e6,digits=1))MB";
    fontsize=9, halign=:center)

rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 3)
colgap!(fig.layout, 1, 8)

outdir = joinpath(@__DIR__, "..", "paper", "figures")
mkpath(outdir)
out = joinpath(outdir, "anim_tt_200x200.gif")
println("Rendering → $out"); flush(stdout)
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]   = mu_frames[f]
    bond_obs[] = bond_frames[f]
    t_label.text[] = "t = $(round(t_frames[f], digits=1))"
end
println("Saved → $out")
