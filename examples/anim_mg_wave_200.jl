"""
anim_mg_wave_200.jl

BranchingDeath RDME wave on 200×200 grid using multigrid L=3.
Single heatmap panel: plasma colormap of ⟨n_k⟩, white boundary around active ring.
"""

using DiscStochSim, CairoMakie, Printf

const k_b = 0.5; const k_d = 0.1; const D = 1.0
const K = 200; const N0 = 3; const dt = 0.4; const T_END = 80.0
const N_FRM = 80; const n_max = 14

const ci = K÷2 + 1; const cj = K÷2 + 1

model = BranchingDeathRDME(k_b, k_d, D, 1.0)
state = SpatialFSP(model, K, K;
    n_max, ε_expand=999.0, ε_equil=0.15,
    organic_activation=true, ε_organic=1e-6,
    use_multigrid=true, multigrid_levels=3, use_block_vcycle=false)

for di in -4:4, dj in -4:4
    di^2+dj^2 <= 16 || continue
    1<=ci+di<=K && 1<=cj+dj<=K || continue
    set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
end

@printf("Model: A→2A  k_b=%.1f  k_d=%.1f  D=%.1f  K=%d  MG L=3\n", k_b,k_d,D,K)
@printf("v_KPP ≈ %.3f  T_END=%.0f  frames=%d\n\n", 2*sqrt(D*(k_b-k_d)), T_END, N_FRM)

# ── Frame helpers ──────────────────────────────────────────────────────────────
const PLASMA = CairoMakie.to_colormap(:plasma)
plasma(t) = let c=PLASMA[clamp(floor(Int,t*(length(PLASMA)-1))+1,1,length(PLASMA))];
            RGBAf(c.r,c.g,c.b,1f0) end

function boundary_segs(sg)
    xs = Float32[]; ys = Float32[]
    for i in 1:K, j in 1:K
        sg[i,j] == 0 && continue
        # horizontal edge below
        if i < K && sg[i+1,j]==0
            push!(xs,j-0.5f0,j+0.5f0,NaN32); push!(ys,i+0.5f0,i+0.5f0,NaN32)
        end
        # horizontal edge above
        if i > 1 && sg[i-1,j]==0
            push!(xs,j-0.5f0,j+0.5f0,NaN32); push!(ys,i-0.5f0,i-0.5f0,NaN32)
        end
        # vertical edge right
        if j < K && sg[i,j+1]==0
            push!(xs,j+0.5f0,j+0.5f0,NaN32); push!(ys,i-0.5f0,i+0.5f0,NaN32)
        end
        # vertical edge left
        if j > 1 && sg[i,j-1]==0
            push!(xs,j-0.5f0,j-0.5f0,NaN32); push!(ys,i-0.5f0,i+0.5f0,NaN32)
        end
    end
    Point2f.(xs, ys)
end

const mu_max = Ref(1.0)

function make_frame(s)
    mg = mean_grid(s)
    sg = status_grid(s)
    cur = maximum(mg; init=1e-3)
    mu_max[] = max(mu_max[], cur * 1.02)
    img = zeros(RGBAf, K, K)
    sc  = max(mu_max[], 1e-3)
    for j in 1:K, i in 1:K
        v = clamp(sqrt(mg[i,j]/sc), 0f0, 1f0)
        img[j,i] = plasma(v)
    end
    img, boundary_segs(sg), length(s.dists), length(s.equil_dists)
end

# ── Pre-compute all frames ─────────────────────────────────────────────────────
t_frames   = collect(range(0.0, T_END; length=N_FRM))
img_frames = Vector{Matrix{RGBAf}}(undef, N_FRM)
seg_frames = Vector{Vector{Point2f}}(undef, N_FRM)
nact       = zeros(Int, N_FRM)
neq        = zeros(Int, N_FRM)

img_frames[1], seg_frames[1], nact[1], neq[1] = make_frame(state)
println("Pre-computing $N_FRM frames…"); flush(stdout)
t0 = time()

let fi = 2
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step_i * dt
        if fi <= N_FRM && t >= t_frames[fi] - 1e-9
            img_frames[fi], seg_frames[fi], nact[fi], neq[fi] = make_frame(state)
            if fi % 10 == 0
                @printf("  frame %3d / %d  t=%5.1f  active=%4d  equil=%5d  %.0fs\n",
                        fi, N_FRM, t, nact[fi], neq[fi], time()-t0)
                flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
end
@printf("Done in %.1fs\n\n", time()-t0)

# ── Render animation ───────────────────────────────────────────────────────────
fig = Figure(size=(640, 680), backgroundcolor=:black)
ax  = Axis(fig[1,1];
    aspect=DataAspect(), yreversed=true,
    backgroundcolor=:black,
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false,
    leftspinevisible=false, rightspinevisible=false,
    topspinevisible=false, bottomspinevisible=false)

img_obs = Observable(img_frames[1])
seg_obs = Observable(seg_frames[1])

image!(ax, 0.5..K+0.5, 0.5..K+0.5, img_obs)
linesegments!(ax, seg_obs; color=(:white,0.75), linewidth=0.8)

t_lbl  = Label(fig[2,1],
    @sprintf("t = 0.0   active: %d   equil: %d", nact[1], neq[1]);
    fontsize=11, color=:gray80, halign=:center)
Label(fig[0,1], "BranchingDeath RDME  200×200  Multigrid L=3";
    fontsize=12, color=:white, font=:bold)

rowgap!(fig.layout, 1, 2); rowgap!(fig.layout, 2, 4)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out    = joinpath(outdir, "anim_mg_wave_200.gif")
println("Rendering → $out"); flush(stdout)
@time record(fig, out, 1:N_FRM; framerate=15) do f
    img_obs[] = img_frames[f]
    seg_obs[] = seg_frames[f]
    t_lbl.text[] = @sprintf("t = %4.1f   active: %d   equil: %d",
                             t_frames[f], nact[f], neq[f])
end
println("Saved → $out")
