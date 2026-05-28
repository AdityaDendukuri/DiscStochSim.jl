"""
anim_mg_structure_200.jl

Three-panel 200×200 BranchingDeath RDME wave with Multigrid L=3.

Panel 1 — plasma heatmap of ⟨n_k⟩ with white wavefront boundary
Panel 2 — multigrid block-structure map (voxel coloured by deepest MG level)
Panel 3 — RDME block matrix for horizontal center-row slice near wavefront:
          diagonal blocks = CME for each voxel (colour = MG level),
          off-diagonal blocks = diffusion coupling between adjacent voxels.
"""

using DiscStochSim, CairoMakie, Printf

# ── Parameters ─────────────────────────────────────────────────────────────────
const k_b  = 0.5;  const k_d = 0.1;  const D = 1.0
const K    = 200;  const N0  = 3;    const dt = 0.4
const T_END = 80.0; const N_FRM = 80; const n_max = 14
const N_LEV = 3        # multigrid levels
const N_SHOW = 13      # voxels in matrix panel slice  (13×15 = 195 ≈ K=200)
const NS = n_max + 1   # states per voxel

const ci = K÷2 + 1; const cj = K÷2 + 1

# ── Colours ────────────────────────────────────────────────────────────────────
const COL_EMPTY  = RGBf(0.04, 0.02, 0.06)
const COL_EQUIL  = RGBf(0.28, 0.28, 0.32)
const COL_EDGE   = RGBf(0.95, 0.18, 0.25)
const COL_UNPAIR = RGBf(0.75, 0.72, 0.08)

const LEV_COLS = [
    [RGBf(0.98,0.58,0.10), RGBf(0.72,0.36,0.04)],
    [RGBf(0.10,0.82,0.78), RGBf(0.04,0.52,0.54)],
    [RGBf(0.20,0.92,0.28), RGBf(0.08,0.58,0.18)],
    [RGBf(0.40,0.50,0.98), RGBf(0.18,0.28,0.72)],
]
lev_col(lev, grp_hue) = LEV_COLS[min(lev, length(LEV_COLS))][(grp_hue % 2) + 1]

const MAT_DIFF = RGBf(0.70, 0.70, 0.88)   # diffusion off-diagonal
const MAT_SEP  = RGBf(0.00, 0.00, 0.00)   # block separator

const PLASMA = CairoMakie.to_colormap(:plasma)
plasma(t) = let c=PLASMA[clamp(floor(Int,t*(length(PLASMA)-1))+1,1,length(PLASMA))];
            RGBAf(c.r,c.g,c.b,1f0) end
const mu_peak = Ref(1.0)

# ── Inline grid_neighbors ──────────────────────────────────────────────────────
@inline function gnb(idx::CartesianIndex{2}, Kx, Ky)
    i,j = Tuple(idx)
    ns = CartesianIndex{2}[]
    i>1  && push!(ns, CartesianIndex(i-1,j))
    i<Kx && push!(ns, CartesianIndex(i+1,j))
    j>1  && push!(ns, CartesianIndex(i,j-1))
    j<Ky && push!(ns, CartesianIndex(i,j+1))
    ns
end

# ── Multigrid hierarchy builder ────────────────────────────────────────────────
function mg_structure(s, n_lev, row_pass)
    active_keys = collect(keys(s.dists))
    active_set  = Set(active_keys)
    equil_set   = Set(keys(s.equil_dists))
    interior = CartesianIndex{2}[]
    edge     = CartesianIndex{2}[]
    for idx in active_keys
        if all(nb -> nb in active_set || nb in equil_set, gnb(idx,s.K_x,s.K_y))
            push!(interior, idx)
        else
            push!(edge, idx)
        end
    end
    Δ_par  = row_pass ? CartesianIndex(0,1) : CartesianIndex(1,0)
    Δ_perp = row_pass ? CartesianIndex(1,0) : CartesianIndex(0,1)
    Δdir(ℓ) = isodd(ℓ) ? Δ_par : Δ_perp
    pair_delta(ℓ) = CartesianIndex(Tuple(Δdir(ℓ+1)) .* 2^((ℓ-1)÷2))
    level_reps   = [interior]
    level_voxels = [[[v] for v in interior]]
    level_pairs  = Vector{Tuple{Int,Int}}[]
    for ℓ in 1:n_lev
        prev = level_reps[end]
        isempty(prev) && break
        ridx = Dict(r=>i for (i,r) in enumerate(prev))
        δ    = pair_delta(ℓ)
        used = Set{Int}(); pairs_ℓ = Tuple{Int,Int}[]
        new_reps = CartesianIndex{2}[]
        new_vox  = Vector{Vector{CartesianIndex{2}}}()
        for (i,rep) in enumerate(prev)
            i in used && continue
            j = get(ridx, rep+δ, 0)
            (j==0 || j in used) && continue
            push!(pairs_ℓ,(i,j)); push!(used,i); push!(used,j)
            push!(new_reps, min(rep,rep+δ))
            push!(new_vox, vcat(level_voxels[end][i], level_voxels[end][j]))
        end
        isempty(pairs_ℓ) && break
        push!(level_pairs, pairs_ℓ)
        push!(level_reps,  new_reps)
        push!(level_voxels, new_vox)
    end
    edge, interior, level_voxels, level_pairs
end

# ── Build voxel-depth dict (deepest MG level each voxel belongs to) ────────────
function build_vox_depth(level_voxels, level_pairs)
    n_built = length(level_pairs)
    vox_depth = Dict{CartesianIndex{2}, Int}()
    vox_grp   = Dict{CartesianIndex{2}, Int}()
    for ℓ in 1:n_built
        covered = Set{Int}()
        if ℓ < n_built
            for (i2,j2) in level_pairs[ℓ+1]
                push!(covered,i2); push!(covered,j2)
            end
        end
        for g in 1:length(level_voxels[ℓ+1])
            g in covered && continue
            vxs = level_voxels[ℓ+1][g]
            ri, rj = Tuple(vxs[1])
            hue = ri * 3 + rj * 7
            for v in vxs
                vox_depth[v] = ℓ
                vox_grp[v]   = hue
            end
        end
    end
    vox_depth, vox_grp
end

# ── Build RDME block matrix panel ─────────────────────────────────────────────
function build_matrix_panel(s, edge_set, vox_depth)
    active_set = Set(keys(s.dists))
    equil_set  = Set(keys(s.equil_dists))

    # Find wavefront in center row (rightward scan)
    wf_j = K
    for j in (K÷2+1):K
        idx = CartesianIndex(ci, j)
        if !(idx in active_set) && !(idx in equil_set)
            wf_j = j; break
        end
    end
    j_end   = min(K, wf_j + 3)
    j_start = max(1, j_end - N_SHOW + 1)
    voxels  = [CartesianIndex(ci,j) for j in j_start:j_end]
    N = length(voxels)

    # mat[x_col, y_row] — first dim = x (column direction), second = y (row direction)
    mat = fill(COL_EMPTY, N_SHOW*NS, N_SHOW*NS)

    function blk_color(idx)
        is_act = idx in active_set
        is_eq  = idx in equil_set
        (!is_act && !is_eq) && return COL_EMPTY
        is_eq  && return RGBf(0.35, 0.35, 0.40)
        idx in edge_set && return COL_EDGE
        lv = get(vox_depth, idx, 0)
        lv == 1 && return LEV_COLS[1][1]
        lv == 2 && return LEV_COLS[2][1]
        lv == 3 && return LEV_COLS[3][1]
        COL_UNPAIR
    end

    for (k, idx) in enumerate(voxels)
        k > N_SHOW && break
        Rk = (k-1)*NS+1 : k*NS
        col = blk_color(idx)
        # Diagonal CME block
        for c in Rk, r in Rk
            mat[c, r] = col
        end
        # Off-diagonal diffusion to right neighbor
        if k < N && k < N_SHOW
            nb = voxels[k+1]
            if nb in active_set || nb in equil_set
                Rkp1 = k*NS+1 : (k+1)*NS
                for c in Rkp1, r in Rk;  mat[c,r] = MAT_DIFF; end
                for c in Rk,   r in Rkp1; mat[c,r] = MAT_DIFF; end
            end
        end
    end

    # 1-px separator between each voxel block
    for k in 1:N_SHOW-1
        sp = k * NS
        for x in 1:N_SHOW*NS
            mat[sp, x] = MAT_SEP
            mat[x, sp] = MAT_SEP
        end
    end
    mat
end

# ── FrameData ─────────────────────────────────────────────────────────────────
struct FrameData
    heat    :: Matrix{RGBAf}
    bnd     :: Vector{Point2f}
    struct_ :: Matrix{RGBf}
    grp     :: Vector{Point2f}
    mat     :: Matrix{RGBf}
    n_act   :: Int
    n_eq    :: Int
end

# ── Build one frame ────────────────────────────────────────────────────────────
function make_frame(s, row_pass)
    K_x, K_y = s.K_x, s.K_y

    # Panel 1: plasma heatmap
    mg  = mean_grid(s)
    cur = maximum(mg; init=1e-3)
    mu_peak[] = max(mu_peak[], cur * 1.02)
    sc = max(mu_peak[], 1e-3)
    heat = [plasma(clamp(sqrt(mg[i,j]/sc), 0f0, 1f0)) for j in 1:K_y, i in 1:K_x]

    sg = status_grid(s)
    bxs = Float32[]; bys = Float32[]
    for i in 1:K_x, j in 1:K_y
        sg[i,j]==0 && continue
        i>1   && sg[i-1,j]==0 && (push!(bxs,j-.5f0,j+.5f0,NaN32); push!(bys,i-.5f0,i-.5f0,NaN32))
        i<K_x && sg[i+1,j]==0 && (push!(bxs,j-.5f0,j+.5f0,NaN32); push!(bys,i+.5f0,i+.5f0,NaN32))
        j>1   && sg[i,j-1]==0 && (push!(bxs,j-.5f0,j-.5f0,NaN32); push!(bys,i-.5f0,i+.5f0,NaN32))
        j<K_y && sg[i,j+1]==0 && (push!(bxs,j+.5f0,j+.5f0,NaN32); push!(bys,i-.5f0,i+.5f0,NaN32))
    end
    bnd_pts = Point2f.(bxs, bys)

    # Panel 2: multigrid structure map
    struct_img = fill(COL_EMPTY, K_y, K_x)
    for idx in keys(s.equil_dists)
        i,j = Tuple(idx); struct_img[j,i] = COL_EQUIL
    end

    edge, interior, level_voxels, level_pairs = mg_structure(s, N_LEV, row_pass)
    edge_set = Set(edge)
    for idx in edge
        i,j = Tuple(idx); struct_img[j,i] = COL_EDGE
    end

    vox_depth, vox_grp = build_vox_depth(level_voxels, level_pairs)

    for idx in interior
        vi, vj = Tuple(idx)
        if haskey(vox_depth, idx)
            struct_img[vj,vi] = lev_col(vox_depth[idx], vox_grp[idx])
        else
            struct_img[vj,vi] = COL_UNPAIR
        end
    end

    # White lines between same-group voxels
    n_built = length(level_pairs)
    gx = Float32[]; gy = Float32[]
    for ℓ in 1:n_built
        covered = Set{Int}()
        if ℓ < n_built
            for (i2,j2) in level_pairs[ℓ+1]; push!(covered,i2); push!(covered,j2); end
        end
        for g in 1:length(level_voxels[ℓ+1])
            g in covered && continue
            vxs = level_voxels[ℓ+1][g]
            vset = Set(vxs)
            for v in vxs
                vi,vj = Tuple(v)
                for nb in gnb(v, K_x, K_y)
                    nb in vset || continue
                    ni,nj = Tuple(nb)
                    push!(gx, Float32(vj), Float32(nj), NaN32)
                    push!(gy, Float32(vi), Float32(ni), NaN32)
                end
            end
        end
    end
    grp_pts = Point2f.(gx, gy)

    # Panel 3: RDME block matrix
    mat_img = build_matrix_panel(s, edge_set, vox_depth)

    n_act = length(s.dists); n_eq = length(s.equil_dists)
    heat, bnd_pts, struct_img, grp_pts, mat_img, n_act, n_eq
end

# ── Build SpatialFSP state ─────────────────────────────────────────────────────
model = BranchingDeathRDME(k_b, k_d, D, 1.0)
state = SpatialFSP(model, K, K;
    n_max, ε_expand=999.0, ε_equil=0.15,
    organic_activation=true, ε_organic=1e-6,
    use_multigrid=true, multigrid_levels=N_LEV, use_block_vcycle=false)

for di in -4:4, dj in -4:4
    di^2+dj^2 <= 16 || continue
    1<=ci+di<=K && 1<=cj+dj<=K || continue
    set_ic!(state, CartesianIndex(ci+di,cj+dj), N0)
end

@printf("BranchingDeath RDME 200×200  MG L=%d  T_END=%.0f  frames=%d\n",
        N_LEV, T_END, N_FRM)

# ── Pre-compute all frames ─────────────────────────────────────────────────────
t_frames = collect(range(0.0, T_END; length=N_FRM))
frames   = Vector{FrameData}(undef, N_FRM)
frames[1] = FrameData(make_frame(state, true)...)

println("Pre-computing $N_FRM frames…"); flush(stdout)
t0 = time()
let fi = 2
    for step_i in 1:round(Int, T_END/dt)
        step!(state, dt)
        t = step_i * dt
        row_pass = isodd(step_i)
        if fi <= N_FRM && t >= t_frames[fi] - 1e-9
            frames[fi] = FrameData(make_frame(state, row_pass)...)
            if fi % 10 == 0
                @printf("  frame %3d/%d  t=%5.1f  active=%5d  equil=%5d  %.0fs\n",
                        fi, N_FRM, t, frames[fi].n_act, frames[fi].n_eq, time()-t0)
                flush(stdout)
            end
            fi += 1
        end
        fi > N_FRM && break
    end
end
@printf("Done in %.1fs\n\n", time()-t0)

# ── Figure layout ──────────────────────────────────────────────────────────────
ax_kw = (aspect=DataAspect(), yreversed=true,
          xticksvisible=false, yticksvisible=false,
          xticklabelsvisible=false, yticklabelsvisible=false,
          leftspinevisible=false, rightspinevisible=false,
          topspinevisible=false, bottomspinevisible=false)

fig = Figure(size=(1900, 730), backgroundcolor=:black)

Label(fig[0,1], "⟨n_k⟩  heatmap";
      fontsize=11, color=:gray90, halign=:center)
Label(fig[0,2], "multigrid block structure  (■red=edge · ■orange=L1 · ■teal=L2 · ■green=L3 · ■gray=equil)";
      fontsize=9,  color=:gray70, halign=:center)
Label(fig[0,3], "RDME block matrix — center-row slice near wavefront  (diagonal=CME · off-diag=diffusion)";
      fontsize=9,  color=:gray70, halign=:center)

ax_heat   = Axis(fig[1,1]; backgroundcolor=:black,              ax_kw...)
ax_struct = Axis(fig[1,2]; backgroundcolor=RGBf(0.04,0.02,0.06), ax_kw...)
ax_mat    = Axis(fig[1,3]; backgroundcolor=RGBf(0.04,0.02,0.06), ax_kw...)

heat_obs   = Observable(frames[1].heat)
bnd_obs    = Observable(frames[1].bnd)
struct_obs = Observable(frames[1].struct_)
grp_obs    = Observable(frames[1].grp)
mat_obs    = Observable(frames[1].mat)

image!(ax_heat,   0.5..K+0.5, 0.5..K+0.5, heat_obs)
linesegments!(ax_heat, bnd_obs; color=(:white,0.75), linewidth=0.8)
image!(ax_struct, 0.5..K+0.5, 0.5..K+0.5, struct_obs)
linesegments!(ax_struct, grp_obs; color=(:white,0.55), linewidth=0.7)
image!(ax_mat,    0.5..K+0.5, 0.5..K+0.5, mat_obs)

t_lbl = Label(fig[2,1:3],
    @sprintf("t = 0.0   active: %d   equil: %d", frames[1].n_act, frames[1].n_eq);
    fontsize=11, color=:gray70, halign=:center)

rowgap!(fig.layout, 1, 2); rowgap!(fig.layout, 2, 4)
colgap!(fig.layout, 1, 6); colgap!(fig.layout, 2, 6)

outdir = joinpath(@__DIR__, "..", "paper", "figures"); mkpath(outdir)
out    = joinpath(outdir, "anim_mg_structure_200.gif")
println("Rendering → $out"); flush(stdout)

@time record(fig, out, 1:N_FRM; framerate=15) do f
    fr = frames[f]
    heat_obs[]   = fr.heat
    bnd_obs[]    = fr.bnd
    struct_obs[] = fr.struct_
    grp_obs[]    = fr.grp
    mat_obs[]    = fr.mat
    t_lbl.text[] = @sprintf("t = %4.1f   active: %d   equil: %d",
                             t_frames[f], fr.n_act, fr.n_eq)
end
println("Saved → $out")
