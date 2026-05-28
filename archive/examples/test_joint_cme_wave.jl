using DiscStochSim, ExponentialUtilities, Printf, Statistics, Random

const k_b=0.5; const k_d=0.1; const D=1.0; const K=40; const N0=3; const dt=0.4
const ci=21; const cj=21; const T_END=14.0
const t_int = 0:round(Int,T_END)

function pure_joint_step!(jf::JointRDMEFSP, dt::Float64)
    isempty(jf.fine_voxels) && return
    DiscStochSim._diffuse_expand!(jf)
    DiscStochSim._rstep_expand!(jf)    # no kemeny_restrict
    sys = DiscStochSim._build_joint_rdme_system(jf)
    A, = build_generator(jf.ss, sys, nothing, jf.t; absorbing=false)
    km = min(jf.krylov_m, length(jf.ss)); km = max(km, 4)
    jf.ss.probs .= expv(dt, A, jf.ss.probs; m=km)
    jf.ss.probs  = max.(0.0, jf.ss.probs)
    prune_threshold!(jf.ss, jf.ε_flux)
    jf.t += dt
end

function emaxr_independent(s::SpatialFSP, ci, cj)
    isempty(s.dists) && return 0.0
    data = [(sqrt(Float64(idx[1]-ci)^2+Float64(idx[2]-cj)^2), 1.0-p[1])
            for (idx,p) in s.dists]
    isempty(data) && return 0.0
    sort!(data; by=first, rev=true)
    R_max = data[1][1]; dr = 0.5; E = 0.0
    for r in 0:dr:R_max+dr
        log_pall = sum(log(max(1e-300,1-pk)) for (rk,pk) in data if rk>=r; init=0.0)
        E += (1.0 - exp(log_pall)) * dr
    end
    E
end

function run_ssa(seed::Int)
    Random.seed!(seed)
    n = zeros(Int32,K,K); n[ci,cj] = N0
    vrate(i,j) = begin
        nv=n[i,j]; nv==0 && return 0.0
        nb=Int(i>1)+Int(i<K)+Int(j>1)+Int(j<K)
        (nv<10 ? k_b*nv : 0.0) + (k_d+D*nb)*nv
    end
    rates=[vrate(i,j) for i in 1:K, j in 1:K]; R=sum(rates)
    front=Float64[]; t=0.0
    for ti in t_int
        while t<Float64(ti) && R>1e-15
            τ=-log(rand())/R; t+=τ
            u=rand()*R; cum=0.0; ei=ej=ci
            for i in 1:K, j in 1:K; cum+=rates[i,j]; if cum>=u; ei=i; ej=j; break end end
            nv=n[ei,ej]; nv==0 && continue
            nb=Int(ei>1)+Int(ei<K)+Int(ej>1)+Int(ej<K)
            birth=nv<10 ? k_b*nv : 0.0
            u2=rand()*(birth+(k_d+D*nb)*nv)
            if u2<birth; n[ei,ej]+=1
            elseif u2<birth+k_d*nv; n[ei,ej]-=1
            else
                dirs=Tuple{Int,Int}[]; ei>1&&push!(dirs,(-1,0)); ei<K&&push!(dirs,(1,0))
                ej>1&&push!(dirs,(0,-1)); ej<K&&push!(dirs,(0,1))
                di,dj=rand(dirs); ni,nj=ei+di,ej+dj
                n[ei,ej]-=1; n[ni,nj]+=1
                R-=rates[ni,nj]; rates[ni,nj]=vrate(ni,nj); R+=rates[ni,nj]
            end
            R-=rates[ei,ej]; rates[ei,ej]=vrate(ei,ej); R+=rates[ei,ej]
        end
        mx=0.0
        for i in 1:K, j in 1:K; n[i,j]>0 && (mx=max(mx,sqrt(Float64(i-ci)^2+Float64(j-cj)^2))); end
        push!(front,mx)
    end
    front
end

# ── Run joint CME ─────────────────────────────────────────────────────────────
println("Running pure joint CME...")
model = BranchingDeathRDME(k_b,k_d,D,1.0)
jf = JointRDMEFSP(model, K, K; n_max=10, ε_flux=1e-5, krylov_m=30)
set_ic!(jf, CartesianIndex(ci,cj), N0)
fronts_joint = Float64[wavefront_radius_mean(jf,ci,cj)]
t0=time()
for si in 1:round(Int,T_END/dt)
    pure_joint_step!(jf, dt)
    t2=si*dt
    if isapprox(t2, round(t2); atol=dt/2)
        push!(fronts_joint, wavefront_radius_mean(jf,ci,cj))
        @printf("  t=%.0f  K_f=%2d  |S|=%4d  E[max_r]=%.2f  sum_p=%.5f  (%.0fs)\n",
            t2, n_fine(jf), n_joint_states(jf), fronts_joint[end], sum(jf.ss.probs), time()-t0)
    end
end

# ── Run per-voxel with organic activation ─────────────────────────────────────
println("\nRunning per-voxel + organic activation...")
state = SpatialFSP(model, K, K; n_max=10, ε_equil=0.15,
                   ε_expand=999.0, organic_activation=true, ε_organic=1e-6)
set_ic!(state, CartesianIndex(ci,cj), N0)
fronts_pv = Float64[emaxr_independent(state,ci,cj)]
t0=time()
for si in 1:round(Int,T_END/dt)
    step!(state, dt); t2=si*dt
    if isapprox(t2, round(t2); atol=dt/2)
        push!(fronts_pv, emaxr_independent(state,ci,cj))
        @printf("  t=%.0f  active=%3d  E[max_r]=%.2f  (%.0fs)\n",
            t2, n_active(state), fronts_pv[end], time()-t0)
    end
end

# ── SSA ───────────────────────────────────────────────────────────────────────
println("\nRunning 20 SSA trajectories...")
ssa_fronts=[run_ssa(s) for s in 1:20]
ssa_mean=mean(ssa_fronts)

println("\n  t    SSA    JointCME  PerVoxelOrganic")
for (i,ti) in enumerate(t_int)
    jv = i<=length(fronts_joint) ? fronts_joint[i] : NaN
    pv = i<=length(fronts_pv)    ? fronts_pv[i]    : NaN
    @printf("  %2d   %5.2f    %5.2f       %5.2f\n", ti, ssa_mean[i], jv, pv)
end
