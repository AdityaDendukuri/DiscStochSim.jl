using DiscStochSim, CairoMakie, Printf

const k_b=0.5; const k_d=0.1; const D=1.0
const K=200; const N0=3; const dt=0.4; const T_END=200.0; const N_FRM=150
const n_max=14

model = BranchingDeathRDME(k_b, k_d, D, 1.0)
v_wave = 2*sqrt(D*(k_b-k_d))
@printf("Multigrid: A→2A k_b=%.1f k_d=%.1f D=%.1f  v_wave≈%.2f\n",k_b,k_d,D,v_wave)

state = SpatialFSP(model, K, K;
    n_max, ε_expand=0.02, ε_equil=0.15,
    use_multigrid=true, use_block_vcycle=false)

ci,cj = K÷2+1, K÷2+1
for di in -4:4, dj in -4:4
    di^2+dj^2<=16 || continue
    1<=ci+di<=K && 1<=cj+dj<=K || continue
    set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
end

const _PLASMA = CairoMakie.to_colormap(:plasma)
_plasma(t) = let c=_PLASMA[clamp(floor(Int,t*(length(_PLASMA)-1))+1,1,length(_PLASMA))];
             RGBf(c.r,c.g,c.b) end
const _mu_max = Ref(1.0)

function make_img(s)
    mg = mean_grid(s)
    cur_max = maximum(mg; init=1e-3)
    _mu_max[] = max(_mu_max[], cur_max*1.05)
    img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        v = clamp(sqrt(mg[i,j]/max(_mu_max[],1e-3)), 0.0, 1.0)
        img[j,i] = _plasma(v)
    end
    img
end

println("Pre-computing $N_FRM frames…")
t_frames = collect(range(0.0, T_END; length=N_FRM))
imgs = Vector{Matrix{RGBf}}(undef, N_FRM)
imgs[1] = make_img(state)

let fi = 2
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step_i*dt
        if fi <= N_FRM && t >= t_frames[fi]-1e-9
            imgs[fi] = make_img(state)
            act=length(state.dists); eq=length(state.equil_dists)
            fi%10==0 && @printf("  frame %3d  t=%5.1f  act=%5d  eq=%5d\n",fi,t,act,eq)
            fi+=1
        end
        fi > N_FRM && break
    end
end

outdir = joinpath(@__DIR__,"..","paper","figures"); mkpath(outdir)
out = joinpath(outdir,"anim_branching_mg.gif")
println("Rendering → $out")
fig = Figure(size=(520,520))
ax  = Axis(fig[1,1]; aspect=DataAspect(), yreversed=true,
           title="Multigrid RDME  (L0→L1→L2)  200×200")
img_obs = Observable(imgs[1])
image!(ax, 0.5..K+0.5, 0.5..K+0.5, img_obs)
@time record(fig, out, 1:length(imgs); framerate=12) do f
    img_obs[] = imgs[f]
end
println("Saved → $out")
