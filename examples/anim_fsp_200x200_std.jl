using DiscStochSim, Catalyst, CairoMakie, Printf

# ── 1. Declare the experiment ──────────────────────────────────────────────────
rn = @reaction_network begin
    k_b, ∅ → A
    k_d, A → ∅
end

exp = RDMEExperiment(rn;
    grid      = rect_grid(200),
    diffusion = (A = 1.0,),
    params    = (k_b = 1.0, k_d = 0.5),
    ic        = CenterIC(:A => 1))

println(exp)

# ── 2. Build state ─────────────────────────────────────────────────────────────
state = build_state(exp; n_max = 12, ε_expand = 0.15, ε_equil = 0.07,
                         use_block_vcycle = true)

const K        = 200
const μ_ss     = 2.0
const dt       = 0.4
const T_END    = 120.0
const N_FRM    = 120
const n_max    = 12
const IMGSIZE  = 300

# ── 3. Colour helpers ──────────────────────────────────────────────────────────
const _MAGMA  = CairoMakie.to_colormap(:magma)
_magma(t)     = let c = _MAGMA[clamp(floor(Int, t*(length(_MAGMA)-1))+1, 1, length(_MAGMA))];
                    RGBf(c.r, c.g, c.b) end
const _mu_max = Ref(μ_ss * 2.5)

const C_EMPTY  = RGBf(0.04, 0.02, 0.06)
const C_ACT    = RGBf(0.88, 0.48, 0.22)
const C_EQ     = RGBf(0.27, 0.45, 0.72)
const COL_DIAG = RGBAf(0.27, 0.51, 0.71, 1.0)
const COL_COUPL= RGBAf(0.89, 0.45, 0.13, 1.0)
const COL_ZERO = RGBAf(0.94, 0.94, 0.94, 1.0)
const COL_SEP  = RGBAf(0.20, 0.20, 0.20, 1.0)

# ── 4. Frame builder ───────────────────────────────────────────────────────────
function make_frames(s)
    mg = mean_grid(s)
    sg = status_grid(s)

    mu_img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        mu_img[j,i] = _magma(clamp(mg[i,j]/_mu_max[], 0.0, 1.0))
    end

    ph_img = zeros(RGBf, K, K)
    for j in 1:K, i in 1:K
        ph_img[j,i] = sg[i,j]==1 ? C_ACT : sg[i,j]==2 ? C_EQ : C_EMPTY
    end

    act   = [(i,j) for i in 1:K, j in 1:K if sg[i,j]==1]
    n_a   = length(act)
    has_eq= any(==(2), sg)
    total = n_a*n_max + (has_eq ? 1 : 0)

    mat_img = fill(COL_ZERO, IMGSIZE, IMGSIZE)
    if total > 0
        aidx = Dict(v=>i for (i,v) in enumerate(act))
        stov = Vector{Int}(undef, total)
        for (vi,_) in enumerate(act); for s in (vi-1)*n_max+1:vi*n_max; stov[s]=vi; end; end
        has_eq && (stov[total] = 0)
        adj = Set{NTuple{2,Int}}(); eq_b = Set{Int}()
        for (vi,(i,j)) in enumerate(act)
            for (ni,nj) in ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
                1≤ni≤K && 1≤nj≤K || continue
                nb = get(aidx,(ni,nj),nothing)
                !isnothing(nb) && push!(adj, vi<nb ? (vi,nb) : (nb,vi))
                sg[ni,nj]==2 && push!(eq_b, vi)
            end
        end
        for r in 1:IMGSIZE
            sr = clamp(floor(Int,(r-.5)/IMGSIZE*total)+1,1,total); vr = stov[sr]
            for c in 1:IMGSIZE
                sc = clamp(floor(Int,(c-.5)/IMGSIZE*total)+1,1,total); vc = stov[sc]
                if vr==vc;             mat_img[r,c]=COL_DIAG
                elseif vr==0||vc==0;   av=vr==0 ? vc : vr; av>0&&av∈eq_b&&(mat_img[r,c]=COL_COUPL)
                else (min(vr,vc),max(vr,vc))∈adj&&(mat_img[r,c]=COL_COUPL)
                end
            end
        end
        if n_a>0 && has_eq
            sep = clamp(round(Int, n_a*n_max/total*IMGSIZE), 1, IMGSIZE)
            for q in 1:IMGSIZE; mat_img[sep,q]=COL_SEP; mat_img[q,sep]=COL_SEP; end
        end
    end

    mu_img, ph_img, mat_img, n_active(s)
end

# ── 5. Pre-compute frames ──────────────────────────────────────────────────────
t_frames    = collect(range(0.0, T_END; length=N_FRM))
mu_frames   = Vector{Matrix{RGBf}}(undef, N_FRM)
ph_frames   = Vector{Matrix{RGBf}}(undef, N_FRM)
mat_frames  = Vector{Matrix{RGBAf}}(undef, N_FRM)
nact_frames = zeros(Int, N_FRM)

mu_frames[1], ph_frames[1], mat_frames[1], nact_frames[1] = make_frames(state)
println("Pre-computing $N_FRM frames…"); flush(stdout)

let fi = 2
    for step in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step*dt
        if fi ≤ N_FRM && t ≥ t_frames[fi]-1e-9
            mu_frames[fi], ph_frames[fi], mat_frames[fi], nact_frames[fi] = make_frames(state)
            fi%10==0 && @printf("  frame %3d  t=%6.1f  empty=%6d  active=%4d  equil=%6d\n",
                                fi, t, n_empty(state), n_active(state), n_equilibrated(state))
            fi += 1
        end
        fi > N_FRM && break
    end
end
@printf("Peak active: %d\n", maximum(nact_frames))

# ── 6. Figure & animation ──────────────────────────────────────────────────────
set_theme!(Theme(fontsize=11, Axis=(spinewidth=0.8, xgridvisible=false, ygridvisible=false)))
fig = Figure(size=(1400, 540))

mu_obs  = Observable(mu_frames[1])
ph_obs  = Observable(ph_frames[1])
mat_obs = Observable(mat_frames[1])

ax_kw = (aspect=DataAspect(), yreversed=true, xaxisposition=:top, titlesize=10,
         xticksvisible=false, yticksvisible=false,
         xticklabelsvisible=false, yticklabelsvisible=false)

ax_mu = Axis(fig[1,1]; title="Exact FSP ⟨n_k⟩  ($(K)×$(K)  μ_ss=$(μ_ss))", ax_kw...)
image!(ax_mu, 0.5..K+0.5, 0.5..K+0.5, mu_obs)

ax_ph = Axis(fig[1,2]; title="Adaptive state space  (orange=active, blue=equil)", ax_kw...)
image!(ax_ph, 0.5..K+0.5, 0.5..K+0.5, ph_obs)

ax_mx = Axis(fig[1,3]; title="RDME block matrix  (active=$(n_max) states, equil=1)",
    aspect=DataAspect(), yreversed=true, xaxisposition=:top, titlesize=10,
    xticksvisible=false, yticksvisible=false,
    xticklabelsvisible=false, yticklabelsvisible=false)
image!(ax_mx, 0.5..IMGSIZE+0.5, 0.5..IMGSIZE+0.5, mat_obs)

for col in 1:3; colsize!(fig.layout, col, Aspect(1,1.0)); end
t_label  = Label(fig[0,1],   "t = 0.0"; fontsize=13, halign=:left)
kc_label = Label(fig[0,2:3], "";         fontsize=11, halign=:center)
Legend(fig[2,1:3],
    [PolyElement(color=c) for c in
     [C_EMPTY, C_ACT, C_EQ,
      RGBf(COL_DIAG.r,COL_DIAG.g,COL_DIAG.b),
      RGBf(COL_COUPL.r,COL_COUPL.g,COL_COUPL.b)]],
    ["Empty","Active (per-voxel CME + V-cycle)","Equil (Poisson mean)",
     "Diagonal block","Diffusion coupling"],
    orientation=:horizontal, framevisible=false, labelsize=9, patchsize=(14,10))
rowgap!(fig.layout,1,5); rowgap!(fig.layout,2,4)
colgap!(fig.layout,1,8); colgap!(fig.layout,2,8)

outdir = joinpath(@__DIR__,"..","paper","figures")
mkpath(outdir)
out = joinpath(outdir, "anim_fsp_200x200.gif")
println("Rendering → $out")
@time record(fig, out, 1:N_FRM; framerate=15) do f
    mu_obs[]  = mu_frames[f]
    ph_obs[]  = ph_frames[f]
    mat_obs[] = mat_frames[f]
    t_label.text[]  = "t = $(round(t_frames[f], digits=1))"
    na = nact_frames[f]
    kc_label.text[] = "active: $na ($(na*n_max) FSP states) | equil: $(n_equilibrated(state))"
end
println("Saved → $out | Peak active: $(maximum(nact_frames))")
