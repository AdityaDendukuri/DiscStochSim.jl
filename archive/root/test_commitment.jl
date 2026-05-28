using DiscStochSim, Random
model = SchloglModel1D(0.1); K=4; d=0.1
n_low,n_uns,n_high = schlogl_fixed_points(model)
c1=model.c1; c2=model.c2; c3=model.c3; c4=model.c4

function ssa_final(seed, t_end)
    rng = Random.MersenneTwister(seed)
    n = fill(n_uns, K); t = 0.0
    while t < t_end
        props = Float64[]
        for k in 1:K
            nk = n[k]
            push!(props, c3, c4*nk, c1*nk*(nk-1)/2, c2*nk*(nk-1)*(nk-2)/6)
            k>1 && push!(props, d*nk)
            k<K && push!(props, d*nk)
        end
        a0 = sum(props); a0 == 0 && break
        t += -log(rand(rng))/a0
        r = rand(rng)*a0; cum = 0.0; fired = 0
        for i in eachindex(props); cum += props[i]; fired = i; cum >= r && break; end
        idx = 0
        for k in 1:K
            base = idx; idx += 4 + (k > 1 ? 1 : 0) + (k < K ? 1 : 0)
            fired <= idx || continue
            lf = fired - base
            lf==1 && (n[k]+=1); lf==2 && (n[k]=max(0,n[k]-1))
            lf==3 && (n[k]+=1); lf==4 && (n[k]=max(0,n[k]-1))
            if lf==5; nb = k > 1 ? k-1 : k+1; n[k]=max(0,n[k]-1); n[nb]+=1; end
            lf==6 && (n[k]=max(0,n[k]-1); n[k+1]+=1)
            break
        end
    end
    n
end

for tend in [40, 100, 200, 400]
    t0 = time()
    finals = [ssa_final(i, Float64(tend)) for i in 1:200]
    dt = time() - t0
    committed = [all(n -> n < 50 || n > 140, f) for f in finals]
    hi_frac = sum(count(x->x>n_uns, f) for f in finals) / (200*K)
    println("t=$(tend): $(sum(committed))/200 committed  hi_frac=$(round(hi_frac,digits=2))  (200 traj in $(round(dt,digits=1))s → 3000 traj ≈ $(round(dt*15,digits=0))s)")
end
